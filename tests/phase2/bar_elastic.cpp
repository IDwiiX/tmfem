#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>

#include "mechanics/elastic2d.hpp"
#include "mesh/gmsh_reader.hpp"
#include "phase2/kirsch_common.hpp"
#include "support/paths.hpp"
#include "support/validation_record.hpp"

// P2-b uniform bar in tension (V2.2): end displacement sigma L / E.
// Reference: docs/REFERENCE_VALUES.md row P2-b, tolerance 1e-8 relative.
// Bar 100 mm x 20 mm (cases/bar/bar.geo), u_x = 0 on the left edge, u_y = 0 at the corner (0, 0),
// traction (sigma, 0) on the right edge. The row is plane stress; the plane-strain value
// (1 - nu^2) sigma L / E, the lateral contraction and the reaction are extra checks.

using namespace tmfem;

namespace {

constexpr double kTolBar = 1e-8;      // P2-b: 1e-8 relative
constexpr double kSigma_Pa = 100e6;   // applied tension
constexpr double kL_m = 0.100, kH_m = 0.020;

void run_bar(int t_nodes, PlaneCondition plane) {
  const ElasticMaterial m = test::phase2_material();
  const auto zero = [](const Eigen::Vector2d&) { return 0.0; };
  const Elastic2DProblem p{read_gmsh41(test::source_dir() / ("cases/bar/bar_t" + std::to_string(t_nodes) + ".msh")),
                           m,
                           plane,
                           "domain",
                           {{"left", 0, zero}, {"corner", 1, zero}},
                           {{"right", [](const Eigen::Vector2d&) { return Eigen::Vector2d(kSigma_Pa, 0.0); }}}};
  const Elastic2DSolution sol = solve_elastic2d(p);

  const bool ps = plane == PlaneCondition::PlaneStress;
  const double eps_xx = ps ? kSigma_Pa / m.E_Pa : (1 - m.nu * m.nu) * kSigma_Pa / m.E_Pa;
  const double eps_yy = ps ? -m.nu * kSigma_Pa / m.E_Pa : -m.nu * (1 + m.nu) * kSigma_Pa / m.E_Pa;
  const double u_end = eps_xx * kL_m;

  double worst_end = 0, worst_any = 0;
  for (const int node : p.mesh.nodes_of("right")) {
    worst_end = std::max(worst_end, std::abs(sol.u_m(2 * node) - u_end) / u_end);
  }
  for (int i = 0; i < p.mesh.n_nodes(); ++i) {  // the whole field is linear, so every node is exact
    const double ux = eps_xx * p.mesh.x_m(i, 0), uy = eps_yy * p.mesh.x_m(i, 1);
    worst_any = std::max(worst_any, std::hypot(sol.u_m(2 * i) - ux, sol.u_m(2 * i + 1) - uy) / u_end);
  }
  double reaction_x = 0;  // N per m of thickness on the x-dofs of the left edge
  for (std::size_t k = 0; k < sol.constrained_dofs.size(); ++k) {
    if (sol.constrained_dofs[k] % 2 == 0) reaction_x += sol.reactions_N(static_cast<Eigen::Index>(k));
  }
  const double reaction_ref = -kSigma_Pa * kH_m;
  const char* pname = ps ? "plane stress" : "plane strain";
  std::printf("P2-b T%d %s: u_end ref %.10e m, max rel error on the end %.2e, whole field %.2e, reaction %.6e N/m\n",
              t_nodes, pname, u_end, worst_end, worst_any, reaction_x);
  INFO("T" << t_nodes << " " << pname);
  CHECK(worst_end < kTolBar);
  CHECK(worst_any < kTolBar);
  CHECK(std::abs(reaction_x - reaction_ref) < kTolBar * std::abs(reaction_ref));
  if (ps) {
    test::record("P2-b", "T" + std::to_string(t_nodes) + " max rel error of end displacement", 0.0, worst_end,
                 "1e-8 relative", "P2-b uniform bar end displacement");
  }
}

}  // namespace

TEST_CASE("P2-b uniform bar end displacement", "[phase2][P2-b]") {
  for (const int t : {3, 6}) {
    run_bar(t, PlaneCondition::PlaneStress);
    run_bar(t, PlaneCondition::PlaneStrain);
  }
}

TEST_CASE("Missing rigid-body constraint is detected", "[phase2][mechanics]") {
  // Bar with u_x = 0 on the left edge but no u_y support: rigid translation in y is free, K_ff singular.
  const auto zero = [](const Eigen::Vector2d&) { return 0.0; };
  for (const int t : {3, 6}) {
    const Elastic2DProblem p{read_gmsh41(test::source_dir() / ("cases/bar/bar_t" + std::to_string(t) + ".msh")),
                             test::phase2_material(),
                             PlaneCondition::PlaneStress,
                             "domain",
                             {{"left", 0, zero}},
                             {{"right", [](const Eigen::Vector2d&) { return Eigen::Vector2d(kSigma_Pa, 0.0); }}}};
    CHECK_THROWS_AS(solve_elastic2d(p), std::runtime_error);
  }
}
