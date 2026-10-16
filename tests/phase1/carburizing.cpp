#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstdio>
#include <string>

#include "io/csv_writer.hpp"
#include "io/json_file.hpp"
#include "support/paths.hpp"
#include "support/validation_record.hpp"
#include "thermal/diffusion_analytic.hpp"
#include "thermal/heat1d.hpp"
#include "thermal/theta_integrator.hpp"

// P1-b carburizing (V1.1): time for C to reach 0.8 wt% at 0.5 mm, against the erf solution.
// Reference: docs/REFERENCE_VALUES.md row P1-b, target 25 370 s, tolerance 0.5%.

using namespace tmfem;

namespace {

struct Carburizing {
  double C0, Cs, D_m2_s, x_m, C_target;
};

Carburizing carburizing() {
  const auto ex = read_json(test::source_dir() / "references/data/diffusion_cases.json").at("ex2_carburizing");
  return {ex.at("C0_wt").get<double>(), ex.at("Cs_wt").get<double>(), ex.at("D_m2_per_s").get<double>(),
          ex.at("x_m").get<double>(), ex.at("C_target_wt").get<double>()};
}

constexpr double kDomain_m = 10e-3;
constexpr double kTmax_s = 12.0 * 3600.0;  // 10 sqrt(D t_max) = 8.3 mm < kDomain_m
constexpr double kProfileTimes_s[2] = {3600.0, 25370.0};

// Marches until C(x_target) crosses C_target and t >= the last profile time; the crossing time
// is linearly interpolated between steps. Profiles at kProfileTimes_s go to `profiles` if given.
double fe_crossing_time(const Carburizing& c, int n_elements, double dt_s, std::vector<std::vector<double>>* profiles) {
  Heat1DProblem p{uniform_interval(kDomain_m, n_elements), c.D_m2_s, 1.0,
                  DirichletBC{[Cs = c.Cs](double) { return Cs; }}, FluxBC{[](double) { return 0.0; }}, {}};
  const int probe = nearest_node(p.mesh, c.x_m);
  REQUIRE(std::abs(p.mesh.x_m[static_cast<std::size_t>(probe)] - c.x_m) < 1e-12);
  // Lumped mass: monotone for the step at the surface (docs/learning/phase1.md section 5).
  const ThetaIntegrator be(build_operators(p, MassType::Lumped), 1.0, dt_s);
  Eigen::VectorXd C = Eigen::VectorXd::Constant(p.mesh.n_nodes(), c.C0);
  C(0) = c.Cs;
  double t_cross = -1.0;
  const int n_max = static_cast<int>(std::lround(kTmax_s / dt_s));
  for (int n = 0; n < n_max; ++n) {
    const double t = n * dt_s;
    const Eigen::VectorXd C_next = be.step(C, t);
    if (t_cross < 0.0 && C(probe) < c.C_target && C_next(probe) >= c.C_target) {
      t_cross = t + dt_s * (c.C_target - C(probe)) / (C_next(probe) - C(probe));
    }
    C = C_next;
    const double t_next = (n + 1) * dt_s;
    if (profiles) {
      for (int k = 0; k < 2; ++k) {
        if (std::abs(t_next - kProfileTimes_s[k]) < 1e-6 * dt_s) {
          (*profiles)[static_cast<std::size_t>(2 * k + 1)] = std::vector<double>(C.data(), C.data() + C.size());
        }
      }
    }
    if (t_cross > 0.0 && t_next >= kProfileTimes_s[1]) break;
  }
  if (profiles) (*profiles)[0] = p.mesh.x_m;
  return t_cross;
}

}  // namespace

TEST_CASE("P1-b carburizing time to 0.8 wt% at 0.5 mm", "[phase1][P1-b]") {
  constexpr double kTol = 5e-3;  // REFERENCE_VALUES P1-b
  constexpr double kTarget_s = 25370.0;
  const Carburizing c = carburizing();
  REQUIRE(kDomain_m >= 10.0 * std::sqrt(c.D_m2_s * kTmax_s));  // SPEC V1.1 semi-infinite condition

  const double t_exact = time_to_reach(c.C_target, c.x_m, c.D_m2_s, c.Cs, c.C0);
  test::record("P1-b", "time analytic [s]", kTarget_s, t_exact, "0.5%", "P1-b carburizing time to 0.8 wt% at 0.5 mm");
  CHECK(std::abs(t_exact - kTarget_s) / kTarget_s <= kTol);

  // Resolution fixed a priori (h = 5 um, dt = 5 s), then halved in both to show convergence.
  std::vector<std::vector<double>> profiles(5);
  const double t_fe = fe_crossing_time(c, 2000, 5.0, &profiles);
  const double t_fine = fe_crossing_time(c, 4000, 2.5, nullptr);
  std::printf("P1-b crossing time: analytic %.2f s, FE %.2f s (%.3e rel), refined FE %.2f s\n", t_exact, t_fe,
              (t_fe - kTarget_s) / kTarget_s, t_fine);
  test::record("P1-b", "time finite element [s]", kTarget_s, t_fe, "0.5%", "P1-b carburizing time to 0.8 wt% at 0.5 mm");
  REQUIRE(t_fe > 0.0);
  CHECK(std::abs(t_fe - kTarget_s) / kTarget_s <= kTol);
  CHECK(std::abs(t_fe - t_fine) / t_fine < kTol / 5.0);

  // Profile CSV for docs/figures (erf profile vs numeric).
  for (int k = 0; k < 2; ++k) {
    std::vector<double> exact;
    for (double x : profiles[0]) exact.push_back(erf_profile(x, kProfileTimes_s[k], c.D_m2_s, c.Cs, c.C0));
    profiles[static_cast<std::size_t>(2 * k + 2)] = exact;
  }
  write_csv(test::results_dir() / "phase1/carburizing_profile.csv",
            {"x_m", "C_fe_t3600s", "C_exact_t3600s", "C_fe_t25370s", "C_exact_t25370s"}, profiles,
            {"P1-b carburizing, C in wt%, lumped BE, h = 5 um, dt = 5 s"});
}
