#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

#include "mechanics/elastoplastic2d.hpp"
#include "mechanics/material_point.hpp"
#include "mesh/gmsh_reader.hpp"
#include "phase3/phase3_common.hpp"
#include "support/validation_record.hpp"

// P3-e (V3.3): a uniform plane-strain bar in FE equals the material-point driver with eps_zz = 0, to 1e-6
// relative at every load step. Reference: docs/REFERENCE_VALUES.md row P3-e.
// Bar 100 mm x 20 mm (cases/bar, T3 and T6): u_x = 0 on the left edge, u_y = 0 at the corner, u_x = lambda
// 0.05 L on the right edge in 20 equal steps. The state is uniform, so every quadrature point must equal the
// driver run with eps_xx = lambda 0.05, eps_zz = 0, sigma_yy = sigma_xy = 0 and the same increments (the
// path is not proportional, so the increments must match; docs/learning/phase3.md section 10).

using namespace tmfem;

namespace {

constexpr double kTol = 1e-6;  // P3-e: 1e-6 relative
constexpr double kL_m = 0.100;
constexpr double kEpsEnd = 0.05;
constexpr int kSteps = 20;

double run_bar(int t_nodes, std::shared_ptr<const J2Plasticity> model, const char* label) {
  std::vector<double> eps_xx;
  for (int k = 1; k <= kSteps; ++k) eps_xx.push_back(kEpsEnd * k / kSteps);
  const DriverResult drv = MaterialPointDriver(*model).run(plane_strain_bar_steps(eps_xx));

  const auto zero = [](const Eigen::Vector2d&) { return 0.0; };
  const auto pull = [](const Eigen::Vector2d&) { return kEpsEnd * kL_m; };
  Elastoplastic2DProblem p{read_gmsh41(test::source_dir() / "cases/bar" / ("bar_t" + std::to_string(t_nodes) + ".msh")),
                           model,
                           PlaneCondition::PlaneStrain,
                           "domain",
                           {{"left", 0, zero}, {"corner", 1, zero}, {"right", 0, pull}},
                           {}};
  ElastoplasticSolver solver(std::move(p));
  double worst = 0.0;
  for (int k = 0; k < kSteps; ++k) {
    const StepReport rep = solver.step_to(static_cast<double>(k + 1) / kSteps);
    const Voigt4& s_ref = drv.stress_Pa[static_cast<std::size_t>(k)];
    const Voigt4& e_ref = drv.strain[static_cast<std::size_t>(k)];
    const double a_ref = drv.state[static_cast<std::size_t>(k)].alpha;
    const double sxx = std::abs(s_ref(0));
    for (std::size_t q = 0; q < solver.quad_points().size(); ++q) {
      const QuadPoint& pt = solver.quad_points()[q];
      const Eigen::Vector3d e = solver.quad_strain(q);
      const double errs[] = {std::abs(pt.stress_Pa(0) - s_ref(0)) / sxx, std::abs(pt.sigma_zz_Pa - s_ref(2)) / sxx,
                             std::abs(pt.stress_Pa(1)) / sxx,           std::abs(pt.stress_Pa(2)) / sxx,
                             std::abs(e(1) - e_ref(1)) / std::abs(e_ref(1)),
                             a_ref > 0.0 ? std::abs(pt.state.alpha - a_ref) / a_ref : pt.state.alpha};
      for (const double err : errs) {
        worst = std::max(worst, err);
        INFO(label << " T" << t_nodes << " step " << k + 1 << " point " << q);
        CHECK(err <= kTol);
      }
    }
    if (k == kSteps - 1) {
      std::printf("  %s T%d: final sigma_xx = %.3f MPa, sigma_zz = %.3f MPa, eps_yy = %.6f, alpha = %.6f, "
                  "last step %d corrections\n",
                  label, t_nodes, s_ref(0) / test::kPaPerMPa, s_ref(2) / test::kPaPerMPa, e_ref(1), a_ref,
                  rep.corrections());
    }
  }
  return worst;
}

}  // namespace

TEST_CASE("P3-e plane-strain bar FE equals material-point driver", "[phase3][P3-e]") {
  std::printf("P3-e plane-strain bar, u_x = lambda 0.05 L, %d steps\n", kSteps);
  double worst = 0.0;
  for (const int t : {3, 6}) {
    worst = std::max(worst, run_bar(t, test::j2(test::tabulated_hardening()), "tabulated"));
    worst = std::max(worst, run_bar(t, test::j2(test::ludwik_hardening(test::source_ludwik())), "Ludwik (source)"));
  }
  std::printf("  max relative difference FE vs driver over all steps and points: %.3e\n", worst);
  test::record("P3-e", "max relative difference FE vs driver (T3, T6, tabulated, Ludwik)", 0.0, worst, "1e-6 relative",
               "P3-e plane-strain bar FE equals material-point driver");
}
