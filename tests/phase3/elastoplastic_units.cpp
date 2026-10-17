#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <memory>

#include "mechanics/elastic2d.hpp"
#include "mechanics/elastoplastic2d.hpp"
#include "mesh/gmsh_reader.hpp"
#include "phase3/phase3_common.hpp"

// Supporting tests of the global elastoplastic solver (docs/learning/phase3.md section 9). No reference row.

using namespace tmfem;

namespace {

double factorial(int n) { return n <= 1 ? 1.0 : n * factorial(n - 1); }

// Bar of P2-b in plane strain: u_x = 0 on the left, u_y = 0 at the corner, traction 100 MPa on the right.
Elastic2DProblem bar_problem(int t_nodes) {
  const auto zero = [](const Eigen::Vector2d&) { return 0.0; };
  const Eigen::Vector2d t(100e6, 0.0);
  return {read_gmsh41(test::source_dir() / "cases/bar" / ("bar_t" + std::to_string(t_nodes) + ".msh")),
          test::phase2_material(),
          PlaneCondition::PlaneStrain,
          "domain",
          {{"left", 0, zero}, {"corner", 1, zero}},
          {{"right", [t](const Eigen::Vector2d&) { return t; }}}};
}

}  // namespace

TEST_CASE("Six-point triangle rule is exact to degree 4", "[phase3]") {
  const QuadRuleTri rule = triangle_rule(6);
  for (int i = 0; i <= 4; ++i) {
    for (int j = 0; i + j <= 4; ++j) {
      double sum = 0.0;
      for (std::size_t q = 0; q < rule.w.size(); ++q) sum += rule.w[q] * std::pow(rule.xi[q](0), i) * std::pow(rule.xi[q](1), j);
      const double exact = factorial(i) * factorial(j) / factorial(i + j + 2);
      INFO("xi^" << i << " eta^" << j);
      CHECK(std::abs(sum - exact) <= 1e-14);
    }
  }
}

TEST_CASE("Elastoplastic solver in the elastic limit equals the phase 2 solver", "[phase3]") {
  for (const int t : {3, 6}) {
    INFO("T" << t);
    const Elastic2DProblem e = bar_problem(t);
    const Elastic2DSolution ref = solve_elastic2d(e);
    // Yield stress 1e30 Pa: the J2 model never leaves its elastic branch.
    const auto model = test::j2(std::make_shared<LinearHardening>(1e30, 0.0));
    ElastoplasticSolver solver({e.mesh, model, e.plane, e.region, e.displacements, e.tractions});
    const StepReport rep = solver.step_to(1.0);
    CHECK(rep.corrections() == 1);
    const double diff = (solver.displacement_m() - ref.u_m).cwiseAbs().maxCoeff();
    CHECK(diff <= 1e-10 * ref.u_m.cwiseAbs().maxCoeff());
    CHECK((solver.reactions_N() - ref.reactions_N).cwiseAbs().maxCoeff() <= 1e-8 * ref.reactions_N.cwiseAbs().maxCoeff());
  }
}

TEST_CASE("Global tangent matches finite differences of the internal force", "[phase3]") {
  // A plastic state of the Kirsch plate (plane strain, level 0, source Ludwik), then the linearisation at
  // 1.2 times the converged displacement, where many points are well inside the plastic regime.
  const auto model = test::j2(test::ludwik_hardening(test::source_ludwik()));
  ElastoplasticSolver solver(test::kirsch_plastic_problem(6, 0, PlaneCondition::PlaneStrain, model));
  solver.step_to(0.8);
  solver.step_to(1.0);
  const Eigen::VectorXd u = 1.2 * solver.displacement_m();
  const ElastoplasticSolver::Linearization lin = solver.linearize(u);
  int n_plastic = 0;
  for (const PlaneResponse& r : lin.response) n_plastic += r.plastic ? 1 : 0;
  REQUIRE(n_plastic > 0);
  Eigen::VectorXd v2(u.size());
  for (Eigen::Index i = 0; i < u.size(); ++i) v2(i) = std::sin(0.37 * static_cast<double>(i) + 1.0);
  for (const Eigen::VectorXd& dir : {Eigen::VectorXd(u), Eigen::VectorXd(v2 * u.cwiseAbs().maxCoeff())}) {
    const double h = 1e-7;
    const Eigen::VectorXd fd = (solver.linearize(u + h * dir).f_int_N - solver.linearize(u - h * dir).f_int_N) / (2.0 * h);
    const Eigen::VectorXd Kv = lin.K_T * dir;
    CHECK((fd - Kv).cwiseAbs().maxCoeff() <= 1e-6 * Kv.cwiseAbs().maxCoeff());
  }
}
