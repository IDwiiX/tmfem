#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "io/csv_writer.hpp"
#include "mechanics/tensile_fit.hpp"
#include "phase3/phase3_common.hpp"
#include "support/validation_record.hpp"

// P3-a: Ludwik parameters from the true tensile curve (docs/learning/phase3.md section 11).
// Reference: docs/REFERENCE_VALUES.md row P3-a. The log-log fit with sigma_y = 310.56 MPa fixed must give
// m = 0.385 and K = 548 MPa (0.5%, kickoff decision Q9); the source's 0.4077 / 559.27 MPa came from a
// different choice of points and are reported next to ours. The deliverable is the 3-parameter least-squares
// fit with its uncertainty: it is reported, and only its implementation is asserted.

using namespace tmfem;

namespace {

constexpr double kTolLogLog = 0.005;      // P3-a log-log reproduction: 0.5% relative
constexpr double kRefM = 0.385;           // P3-a target m
constexpr double kRefK_Pa = 548e6;        // P3-a target K
constexpr double kTolRecovery = 1e-8;     // synthetic-data parameter recovery (implementation)
constexpr double kTolGradient = 1e-10;    // stationarity of the least-squares solution

const char* const kNames[3] = {"sigma_y0 [MPa]", "K [MPa]", "m"};

Eigen::Vector3d as_vector(const LudwikParams& p) { return {p.sigma_y0_Pa, p.K_Pa, p.m}; }
Eigen::Vector3d display_units() { return {1.0 / test::kPaPerMPa, 1.0 / test::kPaPerMPa, 1.0}; }

}  // namespace

TEST_CASE("P3-a Ludwik log-log fit with fixed sigma_y", "[phase3][P3-a]") {
  const test::PlasticPoints p = test::plastic_points();
  const LudwikParams fit = fit_ludwik_loglog(p.alpha, p.sigma_Pa, test::kSigmaStart_Pa);
  const LudwikParams src = test::source_ludwik();
  std::printf("P3-a log-log fit, sigma_y = 310.5611 MPa fixed, 7 points with alpha > 1e-6\n");
  std::printf("  m = %.5f (target 0.385, source 0.4077)   K = %.2f MPa (target 548, source 559.27)\n", fit.m,
              fit.K_Pa / test::kPaPerMPa);
  CHECK(std::abs(fit.m - kRefM) / kRefM <= kTolLogLog);
  CHECK(std::abs(fit.K_Pa - kRefK_Pa) / kRefK_Pa <= kTolLogLog);
  const std::string name = "P3-a Ludwik log-log fit with fixed sigma_y";
  test::record("P3-a", "log-log m (sigma_y fixed)", kRefM, fit.m, "0.5%", name);
  test::record("P3-a", "log-log K [MPa] (sigma_y fixed)", kRefK_Pa / test::kPaPerMPa, fit.K_Pa / test::kPaPerMPa, "0.5%",
               name);
  test::record("P3-a", "source m vs our log-log m (report only)", src.m, fit.m, "report only", name);
}

TEST_CASE("P3-a Ludwik least-squares fit with uncertainty", "[phase3][P3-a]") {
  const test::PlasticPoints p = test::plastic_points();
  const LudwikParams src = test::source_ludwik();

  SECTION("implementation: exact parameters recovered from synthetic Ludwik data") {
    const LudwikParams truth{310e6, 520e6, 0.36};
    std::vector<double> s;
    for (const double a : p.alpha) s.push_back(ludwik_stress_Pa(truth, a));
    const LudwikFit fit = fit_ludwik_least_squares(p.alpha, s, {300e6, 600e6, 0.45});
    CHECK(std::abs(fit.params.sigma_y0_Pa - truth.sigma_y0_Pa) <= kTolRecovery * truth.sigma_y0_Pa);
    CHECK(std::abs(fit.params.K_Pa - truth.K_Pa) <= kTolRecovery * truth.K_Pa);
    CHECK(std::abs(fit.params.m - truth.m) <= kTolRecovery * truth.m);
  }

  SECTION("real data: fit, uncertainty and the source values side by side") {
    const LudwikParams start = fit_ludwik_loglog(p.alpha, p.sigma_Pa, test::kSigmaStart_Pa);
    const LudwikFit fit = fit_ludwik_least_squares(p.alpha, p.sigma_Pa, start);
    const LudwikFit fit2 = fit_ludwik_least_squares(p.alpha, p.sigma_Pa, start, {false, true, true});
    CHECK(fit.gradient_ratio <= kTolGradient);
    CHECK(fit2.gradient_ratio <= kTolGradient);
    CHECK(fit.covariance.llt().info() == Eigen::Success);

    const Eigen::Vector3d u = display_units();
    const Eigen::Vector3d v = as_vector(fit.params).cwiseProduct(u), se = fit.standard_error.cwiseProduct(u),
                          ci = fit.ci95_halfwidth.cwiseProduct(u), jk = fit.jackknife_se.cwiseProduct(u),
                          v2 = as_vector(fit2.params).cwiseProduct(u), ci2 = fit2.ci95_halfwidth.cwiseProduct(u),
                          vs = as_vector(src).cwiseProduct(u), vl = as_vector(start).cwiseProduct(u);
    std::printf("P3-a least squares on stress, n = %d plastic points, 3 parameters free (dof %d, t = %.4f)\n", fit.n,
                fit.dof, fit.t95);
    std::printf("  %-15s %10s %9s %10s %12s | %10s %9s | %10s %10s\n", "parameter", "LSQ", "SE", "95% CI +/-",
                "jackknife SE", "LSQ sy0 fix", "95% CI", "source", "log-log");
    for (int k = 0; k < 3; ++k) {
      std::printf("  %-15s %10.4f %9.4f %10.4f %12.4f | %10.4f %9.4f | %10.4f %10.4f\n", kNames[k], v(k), se(k), ci(k),
                  jk(k), v2(k), ci2(k), vs(k), vl(k));
    }
    std::printf("  residual s = %.3f MPa (3 free), %.3f MPa (sigma_y0 fixed); gradient ratio %.2e\n",
                fit.s_Pa / test::kPaPerMPa, fit2.s_Pa / test::kPaPerMPa, fit.gradient_ratio);
    std::printf("  residuals (MPa):");
    for (const double r : fit.residuals_Pa) std::printf(" %+.2f", r / test::kPaPerMPa);
    std::printf("\n  The jackknife SE of sigma_y0 is large because the yield point alone fixes the knee: there are\n"
                "  no data between alpha = 0 and 0.0086 (docs/learning/phase3.md section 11).\n");

    const std::string name = "P3-a Ludwik least-squares fit with uncertainty";
    for (int k = 0; k < 3; ++k) {
      char q[160];
      std::snprintf(q, sizeof q, "LSQ %s, 95%% CI +/- %.4g, jackknife SE %.4g (reference = source)", kNames[k], ci(k), jk(k));
      test::record("P3-a", q, vs(k), v(k), "report only", name);
    }
    test::record("P3-a", "LSQ residual standard deviation s [MPa]", std::nan(""), fit.s_Pa / test::kPaPerMPa,
                 "report only", name);

    // Curves for the figure: data, LSQ fit with its 95% band, source fit, log-log fit.
    std::vector<double> a, s_fit, band, s_src, s_ll;
    for (int i = 0; i <= 360; ++i) {
      const double ai = 0.18 * i / 360.0;
      a.push_back(ai);
      s_fit.push_back(ludwik_stress_Pa(fit.params, ai) / test::kPaPerMPa);
      band.push_back(ludwik_band_halfwidth_Pa(fit, ai) / test::kPaPerMPa);
      s_src.push_back(ludwik_stress_Pa(src, ai) / test::kPaPerMPa);
      s_ll.push_back(ludwik_stress_Pa(start, ai) / test::kPaPerMPa);
    }
    write_csv(test::results_dir() / "phase3/tensile_fit_curves.csv",
              {"alpha", "sigma_lsq_MPa", "band95_halfwidth_MPa", "sigma_source_MPa", "sigma_loglog_MPa"},
              {a, s_fit, band, s_src, s_ll}, {"Ludwik fits to tensile_true.csv (P3-a)"});
    std::vector<double> s_MPa;
    for (const double s : p.sigma_Pa) s_MPa.push_back(s / test::kPaPerMPa);
    write_csv(test::results_dir() / "phase3/tensile_fit_data.csv", {"alpha", "sigma_MPa", "residual_MPa"},
              {p.alpha, s_MPa, [&] {
                 std::vector<double> r;
                 for (const double x : fit.residuals_Pa) r.push_back(x / test::kPaPerMPa);
                 return r;
               }()},
              {"plastic points of tensile_true.csv; residual = data - LSQ fit"});
  }
}
