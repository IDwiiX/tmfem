#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "io/csv_writer.hpp"
#include "mechanics/elastic2d.hpp"
#include "mechanics/error_norms.hpp"
#include "mechanics/kirsch_analytic.hpp"
#include "phase2/kirsch_common.hpp"
#include "support/convergence.hpp"
#include "support/paths.hpp"
#include "support/validation_record.hpp"

// P2-f Kirsch convergence (V2.4): displacement L2 order 2 (T3) and 3 (T6); energy-norm order 1 (T3)
// and 2 (T6). Reference: docs/REFERENCE_VALUES.md row P2-f, tolerance +/- 0.15 on the observed order.
// Five Gmsh levels (-clscale 1 to 1/16, not nested), h = sqrt(area / n_cells) measured on each mesh,
// errors against the exact Kirsch field (docs/learning/phase2.md section 10), both plane conditions.
// The order is the least-squares slope over levels 1-4. Disclosed change: the first run used levels
// 0-3 only and gave a T6 L2 order of 3.25-3.27, with pairwise orders 3.53, 3.17, 3.12 falling toward 3,
// so level 0 is pre-asymptotic (error-norm quadrature was checked; T6 energy and all T3 orders were
// within tolerance). After review, level 4 was added and level 0 dropped from the fit, for
// T3 and T6 alike. Level 0 is still computed and printed.

using namespace tmfem;

namespace {

constexpr double kTolOrder = 0.15;  // P2-f: +/- 0.15
constexpr int kFitFirstLevel = 1;   // levels 1..4 enter the fit (see the note above)

void run_convergence(int t_nodes, double order_L2, double order_E, const std::string& test_name) {
  const KirschPlate k = test::kirsch_plate();
  const ElasticMaterial m = test::phase2_material();
  std::vector<double> col_plane, col_level, col_fit, col_h, col_L2, col_E;
  for (const auto plane : {PlaneCondition::PlaneStress, PlaneCondition::PlaneStrain}) {
    const std::string pname = plane == PlaneCondition::PlaneStress ? "plane stress" : "plane strain";
    std::vector<double> h, eL2, eE;
    for (int level = 0; level < test::kKirschLevels; ++level) {
      const Elastic2DProblem p = test::kirsch_problem(t_nodes, level, plane);
      const Elastic2DSolution sol = solve_elastic2d(p);
      const ElasticErrors e = elastic_errors(
          p, sol.u_m, [&](const Eigen::Vector2d& x) { return kirsch_displacement(k, m, plane, x); },
          [&](const Eigen::Vector2d& x) { return kirsch_stress_cartesian(k, x); });
      h.push_back(e.h_m);
      eL2.push_back(e.L2 / e.L2_exact);
      eE.push_back(e.energy / e.energy_exact);
      col_plane.push_back(plane == PlaneCondition::PlaneStress ? 0.0 : 1.0);
      col_level.push_back(level);
      col_fit.push_back(level >= kFitFirstLevel ? 1.0 : 0.0);
      col_h.push_back(e.h_m);
      col_L2.push_back(eL2.back());
      col_E.push_back(eE.back());
    }
    const std::string label = "T" + std::to_string(t_nodes) + " " + pname;
    test::print_orders(("P2-f " + label + ", relative L2 displacement error").c_str(), h, eL2);
    test::print_orders(("P2-f " + label + ", relative energy-norm error").c_str(), h, eE);
    const auto fit = [](const std::vector<double>& v) { return std::vector<double>(v.begin() + kFitFirstLevel, v.end()); };
    const double pL2 = test::observed_order(fit(h), fit(eL2)), pE = test::observed_order(fit(h), fit(eE));
    std::printf("  asserted orders (fit over levels %d-%d): L2 %.4f, energy %.4f\n", kFitFirstLevel,
                test::kKirschLevels - 1, pL2, pE);
    INFO(label);
    CHECK(std::abs(pL2 - order_L2) <= kTolOrder);
    CHECK(std::abs(pE - order_E) <= kTolOrder);
    test::record("P2-f", label + " L2 order (levels 1-4)", order_L2, pL2, "+/- 0.15", test_name);
    test::record("P2-f", label + " energy order (levels 1-4)", order_E, pE, "+/- 0.15", test_name);
  }
  write_csv(test::results_dir() / "phase2" / ("convergence_t" + std::to_string(t_nodes) + ".csv"),
            {"plane_strain", "level", "in_fit", "h_m", "L2_rel", "energy_rel"},
            {col_plane, col_level, col_fit, col_h, col_L2, col_E},
            {"P2-f Kirsch convergence, T" + std::to_string(t_nodes) + "; plane_strain 0 = plane stress"});
}

}  // namespace

TEST_CASE("P2-f Kirsch convergence T3", "[phase2][P2-f]") {
  run_convergence(3, 2.0, 1.0, "P2-f Kirsch convergence T3");
}

TEST_CASE("P2-f Kirsch convergence T6", "[phase2][P2-f]") {
  run_convergence(6, 3.0, 2.0, "P2-f Kirsch convergence T6");
}
