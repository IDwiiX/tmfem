#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>

#include "io/csv_writer.hpp"
#include "io/json_file.hpp"
#include "support/convergence.hpp"
#include "support/paths.hpp"
#include "support/validation_record.hpp"
#include "thermal/diffusion_analytic.hpp"
#include "thermal/heat1d.hpp"
#include "thermal/theta_integrator.hpp"

// P1-h thin film on a zero-flux half-space (V1.2): dose conserved to 1e-10 relative and the peak
// decaying as t^-1/2 (slope -0.5 +/- 0.01). Reference: docs/REFERENCE_VALUES.md row P1-h.
// Parameters: ex. 5 drive-in at 1200 C (references/data/diffusion_cases.json), dose from the
// predeposition formula (P1-e).

using namespace tmfem;

namespace {

struct ThinFilmRun {
  std::vector<double> t_s, peak, peak_exact, dose;
};

constexpr double kT0_s = 1.0;      // start from the analytic profile at t0 > 0 (not a delta)
constexpr double kTend_s = 3600.0;
constexpr double kDomain_m = 2e-5; // >= 10 sqrt(D t_end) = 1.04e-5 m

double drive_in_D_m2_s() {
  const auto ex = read_json(test::source_dir() / "references/data/diffusion_cases.json").at("ex5_boron_in_silicon");
  return ex.at("D_1200C_cm2_per_s").get<double>() * 1e-4;
}

double film_dose_per_m2() {
  const auto ex = read_json(test::source_dir() / "references/data/diffusion_cases.json").at("ex5_boron_in_silicon");
  return predeposition_dose(ex.at("Cs_at_per_cm3").get<double>() * 1e6, ex.at("D_1100C_cm2_per_s").get<double>() * 1e-4,
                            ex.at("t_s").get<double>());
}

// Geometric mesh (self-similar, matching the self-similar solution) and BE with dt = c t.
ThinFilmRun run(int n_elements, double ratio, double dt_over_t) {
  const double D = drive_in_D_m2_s();
  const double M = film_dose_per_m2();
  Heat1DProblem p{graded_interval(kDomain_m, n_elements, ratio), D, 1.0, FluxBC{[](double) { return 0.0; }},
                  FluxBC{[](double) { return 0.0; }}, {}};
  ThetaIntegrator be(build_operators(p, MassType::Consistent), 1.0, dt_over_t * kT0_s);
  const auto& Mmat = be.operators().M;
  const Eigen::VectorXd ones = Eigen::VectorXd::Ones(p.mesh.n_nodes());
  Eigen::VectorXd C = nodal_values(p.mesh, [&](double x) { return thin_film_profile(x, kT0_s, D, M); });
  ThinFilmRun r;
  double t = kT0_s;
  auto store = [&] {
    r.t_s.push_back(t);
    r.peak.push_back(C(0));
    r.peak_exact.push_back(thin_film_profile(0.0, t, D, M));
    r.dose.push_back(ones.dot(Mmat * C));
  };
  store();
  while (t < kTend_s * (1.0 - 1e-12)) {
    be.set_dt(std::min(dt_over_t * t, kTend_s - t));
    C = be.step(C, t);
    t += be.dt();
    store();
  }
  return r;
}

double peak_slope(const ThinFilmRun& r, double t_from) {
  std::vector<double> t, c;
  for (std::size_t i = 0; i < r.t_s.size(); ++i) {
    if (r.t_s[i] >= t_from) {
      t.push_back(r.t_s[i]);
      c.push_back(r.peak[i]);
    }
  }
  return test::observed_order(t, c);
}

}  // namespace

TEST_CASE("P1-h thin film dose conservation", "[phase1][P1-h]") {
  constexpr double kTol = 1e-10;  // REFERENCE_VALUES P1-h, dose conservation
  const ThinFilmRun r = run(400, 1.02, 0.01);
  double worst = 0.0;
  for (double d : r.dose) worst = std::max(worst, std::abs(d - r.dose.front()) / r.dose.front());
  std::printf("P1-h dose: initial FE %.10e at/m2 (analytic M %.10e), max relative drift %.3e over %zu steps\n",
              r.dose.front(), film_dose_per_m2(), worst, r.t_s.size() - 1);
  test::record("P1-h", "max relative dose drift", 0.0, worst, "1e-10", "P1-h thin film dose conservation");
  CHECK(worst <= kTol);
}

TEST_CASE("P1-h thin film peak decays as t^-1/2", "[phase1][P1-h]") {
  constexpr double kTol = 0.01;      // REFERENCE_VALUES P1-h, slope -0.5 +/- 0.01
  constexpr double kFitFrom_s = 10.0;
  const ThinFilmRun r = run(400, 1.02, 0.01);
  const ThinFilmRun fine = run(800, std::sqrt(1.02), 0.005);
  const double slope = peak_slope(r, kFitFrom_s);
  const double slope_fine = peak_slope(fine, kFitFrom_s);
  std::printf("P1-h peak slope over [%g, %g] s: %.5f (refined %.5f); peak at t_end FE %.6e, exact %.6e\n",
              kFitFrom_s, kTend_s, slope, slope_fine, r.peak.back(), r.peak_exact.back());
  test::record("P1-h", "peak decay slope", -0.5, slope, "+/- 0.01", "P1-h thin film peak decays as t^-1/2");
  CHECK(std::abs(slope + 0.5) <= kTol);
  CHECK(std::abs(slope - slope_fine) < kTol / 5.0);

  std::vector<double> drift;
  for (double d : r.dose) drift.push_back((d - r.dose.front()) / r.dose.front());
  write_csv(test::results_dir() / "phase1/thin_film_decay.csv", {"t_s", "C0_fe_per_m3", "C0_exact_per_m3", "dose_rel_drift"},
            {r.t_s, r.peak, r.peak_exact, drift}, {"P1-h thin film, BE dt = 0.01 t, graded mesh 400 elements"});
}
