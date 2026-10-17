#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "io/csv_writer.hpp"
#include "mechanics/material_point.hpp"
#include "phase3/phase3_common.hpp"
#include "support/validation_record.hpp"

// P3-b and P3-c: the material-point driver in uniaxial stress (docs/learning/phase3.md sections 3 and 10).

using namespace tmfem;

namespace {

constexpr double kTolP3b = 0.005;          // P3-b: 0.5% relative
constexpr double kRefEpsYY = 0.03158;      // P3-b target eps_yy (Abaqus)
constexpr double kRefEpsLateral = -0.01535;  // P3-b target eps_xx = eps_zz (Abaqus)
constexpr double kSigmaA_Pa = 457.2405e6;  // the data point of P3-b
constexpr double kTolRoundTrip = 1e-10;    // tabulated law reproduces its own data
constexpr double kTolClosedForm = 1e-8;    // P3-c (1): driver vs scalar closed-form Ludwik
constexpr double kBandFactor = 2.0;        // P3-c (2): within +/- 2 s of the fit (kickoff decision Q10)

// Uniaxial Ludwik in closed form: sigma = sigma_y0 + K (eps - sigma / E)^m solved for sigma by bisection.
double uniaxial_ludwik_Pa(const LudwikParams& p, double E, double eps) {
  if (E * eps <= p.sigma_y0_Pa) return E * eps;
  double lo = p.sigma_y0_Pa, hi = E * eps;
  for (int it = 0; it < 200 && hi - lo > 1e-15 * hi; ++it) {
    const double mid = 0.5 * (lo + hi);
    const double g = mid - ludwik_stress_Pa(p, std::max(eps - mid / E, 0.0));
    (g < 0.0 ? lo : hi) = mid;
  }
  return 0.5 * (lo + hi);
}

}  // namespace

TEST_CASE("P3-b uniaxial strains at 457.24 MPa, tabulated hardening", "[phase3][P3-b]") {
  // Reference: docs/REFERENCE_VALUES.md row P3-b (eps_yy 0.03158, eps_xx = eps_zz -0.01535, 0.5%).
  // Nearly tautological: 457.24 MPa is a data point and the tabulated law goes through it, so eps_yy is a
  // round trip. The informative part is the lateral strain: -nu sigma / E elastic plus -alpha / 2 plastic.
  const auto model = test::j2(test::tabulated_hardening());
  const MaterialPointDriver driver(*model);
  const DriverResult r = driver.run(uniaxial_stress_steps(kSigmaA_Pa, 20));
  const Voigt4 eps = r.strain.back();
  std::printf("P3-b at sigma_yy = 457.2405 MPa: eps_yy = %.6f (target 0.03158), eps_xx = %.7f, eps_zz = %.7f "
              "(target -0.01535), alpha = %.7f\n",
              eps(1), eps(0), eps(2), r.state.back().alpha);
  CHECK(std::abs(eps(1) - kRefEpsYY) / kRefEpsYY <= kTolP3b);
  CHECK(std::abs(eps(0) - kRefEpsLateral) / std::abs(kRefEpsLateral) <= kTolP3b);
  CHECK(std::abs(eps(2) - kRefEpsLateral) / std::abs(kRefEpsLateral) <= kTolP3b);
  CHECK(std::abs(eps(0) - eps(2)) <= 1e-14);
  const std::string name = "P3-b uniaxial strains at 457.24 MPa, tabulated hardening";
  test::record("P3-b", "eps_yy at 457.24 MPa", kRefEpsYY, eps(1), "0.5%", name);
  test::record("P3-b", "eps_xx at 457.24 MPa", kRefEpsLateral, eps(0), "0.5%", name);
  test::record("P3-b", "eps_zz at 457.24 MPa", kRefEpsLateral, eps(2), "0.5%", name);

  // V3.1: the driver reproduces the true curve at every plastic data point.
  const test::PlasticPoints p = test::plastic_points();
  std::vector<double> e_out, s_out;
  for (std::size_t i = 0; i < p.eps.size(); ++i) {
    const DriverResult ri = driver.run(uniaxial_stress_steps(p.sigma_Pa[i], 20));
    const double e = ri.strain.back()(1);
    e_out.push_back(e);
    s_out.push_back(p.sigma_Pa[i] / test::kPaPerMPa);
    INFO("data point " << i << ": sigma = " << p.sigma_Pa[i] / test::kPaPerMPa << " MPa");
    if (i == 0) {
      // alpha was clipped from -3.2e-6 to 0 here, so the strain differs from the data by that amount.
      std::printf("  yield point: eps = %.7f vs data %.7f (clipped alpha, difference %.2e)\n", e, p.eps[i], e - p.eps[i]);
      CHECK(std::abs(e - p.eps[i]) <= 4e-6);
    } else {
      CHECK(std::abs(e - p.eps[i]) <= kTolRoundTrip * p.eps[i]);
    }
  }
  write_csv(test::results_dir() / "phase3/material_point_tabulated.csv", {"eps_yy", "sigma_yy_MPa"}, {e_out, s_out},
            {"uniaxial material point, tabulated hardening, at each plastic data stress (P3-b)"});
}

TEST_CASE("P3-c Ludwik material point along the true curve", "[phase3][P3-c]") {
  // Reference: docs/REFERENCE_VALUES.md row P3-c ("within the fit's residual band"; kickoff decision Q10):
  // (1) the driver equals the closed-form uniaxial Ludwik solution, (2) every data point with eps_p > 0 lies
  // within +/- 2 s of the least-squares Ludwik curve, s being the fit's residual standard deviation.
  const test::PlasticPoints p = test::plastic_points();
  const LudwikParams start = fit_ludwik_loglog(p.alpha, p.sigma_Pa, test::kSigmaStart_Pa);
  const LudwikFit fit = fit_ludwik_least_squares(p.alpha, p.sigma_Pa, start);
  const double E = test::phase2_material().E_Pa;
  std::vector<double> eps;
  for (std::size_t i = 1; i < p.eps.size(); ++i) eps.push_back(p.eps[i]);  // alpha > 0

  double worst_band = 0.0;
  for (const auto& [label, params] : {std::make_pair("LSQ fit", fit.params), std::make_pair("source", test::source_ludwik())}) {
    const auto model = test::j2(test::ludwik_hardening(params));
    const DriverResult r = MaterialPointDriver(*model).run(uniaxial_strain_steps(eps));
    std::printf("P3-c Ludwik (%s) at the data strains; s = %.3f MPa\n  eps_yy     data MPa   driver MPa  (driver-data)/s\n",
                label, fit.s_Pa / test::kPaPerMPa);
    for (std::size_t k = 0; k < eps.size(); ++k) {
      const double s = r.stress_Pa[k](1), data = p.sigma_Pa[k + 1];
      const double z = (s - data) / fit.s_Pa;
      std::printf("  %.6f %10.3f %11.3f %+8.2f\n", eps[k], data / test::kPaPerMPa, s / test::kPaPerMPa, z);
      INFO(label << ", eps_yy = " << eps[k]);
      CHECK(std::abs(s - uniaxial_ludwik_Pa(params, E, eps[k])) <= kTolClosedForm * s);
      CHECK(std::abs(r.stress_Pa[k](0)) <= 1e-4);
      if (std::string(label) == "LSQ fit") {
        CHECK(std::abs(z) <= kBandFactor);
        worst_band = std::max(worst_band, std::abs(z));
      }
    }
  }
  test::record("P3-c", "max |driver - data| / s, LSQ Ludwik, 7 points", 0.0, worst_band, "+/- 2 s",
               "P3-c Ludwik material point along the true curve");
}
