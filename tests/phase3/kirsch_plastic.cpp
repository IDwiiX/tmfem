#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstdio>
#include <limits>
#include <string>
#include <vector>

#include "io/csv_writer.hpp"
#include "mechanics/elastoplastic2d.hpp"
#include "mechanics/neuber.hpp"
#include "phase3/phase3_common.hpp"
#include "support/validation_record.hpp"

// P3-f (V3.4): elastoplastic hole plate at a remote 140 MPa, report only. Reference: docs/REFERENCE_VALUES.md
// row P3-f ("report"; the source has only an Abaqus von Mises map). Quarter Kirsch model with the elastic
// Kirsch tractions scaled by lambda (kickoff decision Q13): once a plastic zone forms these tractions only
// approximate an infinite plate, but the zone is small next to the 5a boundary. Tabulated hardening, S ramped
// 0 -> 140 MPa in 28 steps of 5 MPa, T6 levels 1-3 (level 4 in a [slow] test).
// Asserted: every step converges, global equilibrium, and every plastic point lies on its yield surface.
// Reported: hole-edge strain and stress, plastic zone, and Neuber's estimate on the same curve
// (docs/learning/phase3.md section 12). Neuber assumes a uniaxial notch root, which holds only in plane stress.

using namespace tmfem;

namespace {

constexpr int kSteps = 28;
constexpr double kKtS_Pa = 3.0 * test::kKirschS_Pa;  // elastic peak 420 MPa
constexpr double kTolEquilibrium = 1e-8;
constexpr double kTolYield = 1e-8;

struct PlateResult {
  int level = 0, nodes = 0, max_corrections = 0;
  double eps_yy_edge = 0, syy_gp = 0, szz_gp = 0, vm_gp = 0, alpha_gp = 0, gp_distance = 0;
  double syy_max = 0, syy_max_r = 0, alpha_max = 0, zone_extent = 0, zone_area = 0;
  double residual_syy_gp = std::numeric_limits<double>::quiet_NaN();
};

std::string plane_name(PlaneCondition plane) { return plane == PlaneCondition::PlaneStress ? "plane stress" : "plane strain"; }

// Global equilibrium in y: reactions on sym_y plus the applied y forces sum to zero.
void check_equilibrium(const ElastoplasticSolver& s) {
  const Eigen::VectorXd f = s.load_factor() * s.traction_load_N();
  const Eigen::VectorXd R = s.reactions_N();
  double sum = 0.0, scale = 0.0;
  for (Eigen::Index i = 1; i < f.size(); i += 2) {
    sum += f(i);
    scale += std::abs(f(i));
  }
  for (std::size_t k = 0; k < s.constrained_dofs().size(); ++k) {
    if (s.constrained_dofs()[k] % 2 == 1) sum += R(static_cast<Eigen::Index>(k));
  }
  CHECK(std::abs(sum) <= kTolEquilibrium * scale);
}

std::size_t nearest_point(const ElastoplasticSolver& s, const Eigen::Vector2d& x) {
  std::size_t best = 0;
  for (std::size_t q = 1; q < s.quad_points().size(); ++q) {
    if ((s.quad_points()[q].x_m - x).norm() < (s.quad_points()[best].x_m - x).norm()) best = q;
  }
  return best;
}

PlateResult run_plate(PlaneCondition plane, int level, int steps, bool unload, bool write_gauss) {
  const auto h = test::tabulated_hardening();
  ElastoplasticSolver solver(test::kirsch_plastic_problem(6, level, plane, test::j2(h)));
  PlateResult r;
  r.level = level;
  r.nodes = solver.problem().mesh.n_nodes();
  for (int k = 1; k <= steps; ++k) {
    r.max_corrections = std::max(r.max_corrections, solver.step_to(static_cast<double>(k) / steps).corrections());
  }
  check_equilibrium(solver);

  const Eigen::Vector2d edge(test::kKirschA_m, 0.0);  // theta = 90 deg from the load axis
  r.eps_yy_edge = (*strain_at_point(solver.problem().mesh, "domain", solver.displacement_m(), edge))(1);
  const std::size_t g = nearest_point(solver, edge);
  const QuadPoint& qg = solver.quad_points()[g];
  r.gp_distance = (qg.x_m - edge).norm();
  r.syy_gp = qg.stress_Pa(1);
  r.szz_gp = qg.sigma_zz_Pa;
  r.vm_gp = test::equivalent(qg.stress_Pa, qg.sigma_zz_Pa);
  r.alpha_gp = qg.state.alpha;
  std::vector<double> x, y, al, vm, syy, szz;
  for (const QuadPoint& q : solver.quad_points()) {
    const double q_eq = test::equivalent(q.stress_Pa, q.sigma_zz_Pa);
    if (q.plastic) {
      const double sy = h->yield_stress_Pa(q.state.alpha);
      CHECK(std::abs(q_eq - sy) <= kTolYield * sy);
    }
    if (q.stress_Pa(1) > r.syy_max) {  // the field maximum lies at the plastic-zone boundary, not at the edge
      r.syy_max = q.stress_Pa(1);
      r.syy_max_r = q.x_m.norm();
    }
    r.alpha_max = std::max(r.alpha_max, q.state.alpha);
    if (q.state.alpha > 0.0) {
      r.zone_extent = std::max(r.zone_extent, q.x_m.norm() - test::kKirschA_m);
      r.zone_area += q.weight_m2;
    }
    x.push_back(q.x_m(0));
    y.push_back(q.x_m(1));
    al.push_back(q.state.alpha);
    vm.push_back(q_eq / test::kPaPerMPa);
    syy.push_back(q.stress_Pa(1) / test::kPaPerMPa);
    szz.push_back(q.sigma_zz_Pa / test::kPaPerMPa);
  }
  if (write_gauss) {
    const std::string tag = plane == PlaneCondition::PlaneStress ? "plane_stress" : "plane_strain";
    write_csv(test::results_dir() / ("phase3/kirsch_plastic_" + tag + "_l" + std::to_string(level) + "_gauss.csv"),
              {"x_m", "y_m", "alpha", "von_mises_MPa", "syy_MPa", "szz_MPa"}, {x, y, al, vm, syy, szz},
              {"P3-f quadrature points at S = 140 MPa, " + plane_name(plane) + ", tabulated hardening, T6 level " +
               std::to_string(level)});
  }
  if (unload) {
    for (int k = steps - 1; k >= 0; --k) solver.step_to(static_cast<double>(k) / steps);
    r.residual_syy_gp = solver.quad_points()[g].stress_Pa(1);
  }
  return r;
}

void print(const PlateResult& r) {
  const double M = test::kPaPerMPa;
  std::printf("  level %d (%5d nodes): eps_yy(a,0) = %.6f | nearest Gauss point (%.2e m from the edge): "
              "syy = %.2f, szz = %.2f, vm = %.2f MPa, alpha = %.6f | max syy = %.2f MPa at r = %.3f a, "
              "max alpha = %.6f | zone r - a = %.3f a, area %.4f a^2 | max %d corrections\n",
              r.level, r.nodes, r.eps_yy_edge, r.gp_distance, r.syy_gp / M, r.szz_gp / M, r.vm_gp / M, r.alpha_gp,
              r.syy_max / M, r.syy_max_r / test::kKirschA_m, r.alpha_max, r.zone_extent / test::kKirschA_m,
              r.zone_area / (test::kKirschA_m * test::kKirschA_m), r.max_corrections);
}

void report(PlaneCondition plane, const std::vector<int>& levels, const std::string& test_name, bool extras) {
  const NeuberEstimate n = neuber_uniaxial(*test::tabulated_hardening(), test::phase2_material().E_Pa, kKtS_Pa);
  const bool neuber_applies = plane == PlaneCondition::PlaneStress;
  std::printf("P3-f elastoplastic hole plate, %s, S = 140 MPa, tabulated hardening, %d steps\n",
              plane_name(plane).c_str(), kSteps);
  std::printf("  Neuber on the same curve: sigma = %.2f MPa, eps = %.6f, alpha = %.6f%s\n", n.sigma_Pa / 1e6, n.eps,
              n.alpha, neuber_applies ? "" : "  (not applicable here: biaxial notch root in plane strain)");
  std::vector<PlateResult> results;
  for (const int level : levels) {
    results.push_back(run_plate(plane, level, kSteps, extras && level == levels.back(), level == levels.back()));
    print(results.back());
  }
  const PlateResult& f = results.back();
  if (results.size() > 1) {
    const PlateResult& c = results[results.size() - 2];
    std::printf("  change from level %d to %d: eps_yy(a,0) %+.2f%%, max alpha %+.2f%%, zone extent %+.2f%%\n", c.level,
                f.level, 100 * (f.eps_yy_edge / c.eps_yy_edge - 1), 100 * (f.alpha_max / c.alpha_max - 1),
                100 * (f.zone_extent / c.zone_extent - 1));
  }
  if (extras) {
    const PlateResult h = run_plate(plane, 2, 2 * kSteps, false, false);
    const PlateResult& base = results[0].level == 2 ? results[0] : results[std::min<std::size_t>(1, results.size() - 1)];
    std::printf("  load-step sensitivity, level 2, %d vs %d steps: eps_yy(a,0) %+.3f%%, max alpha %+.3f%%\n",
                2 * kSteps, kSteps, 100 * (h.eps_yy_edge / base.eps_yy_edge - 1), 100 * (h.alpha_max / base.alpha_max - 1));
    std::printf("  after unloading to S = 0: residual syy at the Gauss point nearest the edge = %.2f MPa\n",
                f.residual_syy_gp / 1e6);
  }
  // The Neuber comparison is at the notch root: edge strain, stress at the quadrature point nearest the edge, and
  // alpha there (the maximum alpha is at the edge). In plane strain Neuber does not apply: no reference (NaN).
  const std::string p = plane_name(plane) + ", T6 level " + std::to_string(f.level) + ": ";
  const std::string ref = neuber_applies ? " (reference = Neuber)" : "";
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const auto neuber = [&](double v) { return neuber_applies ? v : nan; };
  test::record("P3-f", p + "hole-edge eps_yy" + ref, neuber(n.eps), f.eps_yy_edge, "report only", test_name);
  test::record("P3-f", p + "sigma_yy at the point nearest the edge [MPa]" + ref, neuber(n.sigma_Pa / 1e6),
               f.syy_gp / 1e6, "report only", test_name);
  test::record("P3-f", p + "max alpha (at the edge)" + ref, neuber(n.alpha), f.alpha_max, "report only", test_name);
  char where[96];
  std::snprintf(where, sizeof where, "max sigma_yy in the field [MPa], at r = %.3f a", f.syy_max_r / test::kKirschA_m);
  test::record("P3-f", p + where, nan, f.syy_max / 1e6, "report only", test_name);
  test::record("P3-f", p + "plastic zone extent (r - a) / a", nan, f.zone_extent / test::kKirschA_m, "report only",
               test_name);
}

}  // namespace

TEST_CASE("P3-f elastoplastic hole plate, plane strain (report)", "[phase3][P3-f]") {
  report(PlaneCondition::PlaneStrain, {1, 2, 3}, "P3-f elastoplastic hole plate, plane strain (report)", true);
}

TEST_CASE("P3-f elastoplastic hole plate, plane strain, level 4 (report)", "[phase3][P3-f][slow]") {
  report(PlaneCondition::PlaneStrain, {4}, "P3-f elastoplastic hole plate, plane strain, level 4 (report)", false);
}

// Stage B: the thin plate of the source's Abaqus model. Here the notch root is uniaxial (sigma_rr = 0 on the
// free hole, sigma_zz = 0), so Neuber's rule applies and its difference from FE is reported.
TEST_CASE("P3-f elastoplastic hole plate, plane stress (report)", "[phase3][P3-f]") {
  report(PlaneCondition::PlaneStress, {1, 2, 3}, "P3-f elastoplastic hole plate, plane stress (report)", true);
}

TEST_CASE("P3-f elastoplastic hole plate, plane stress, level 4 (report)", "[phase3][P3-f][slow]") {
  report(PlaneCondition::PlaneStress, {4}, "P3-f elastoplastic hole plate, plane stress, level 4 (report)", false);
}
