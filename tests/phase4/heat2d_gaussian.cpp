#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "io/csv_writer.hpp"
#include "mesh/gmsh_reader.hpp"
#include "phase4/phase4_common.hpp"
#include "support/convergence.hpp"
#include "support/validation_record.hpp"
#include "thermal/heat2d.hpp"
#include "thermal/theta_integrator.hpp"

// P4-e 2D Gaussian heat spreading (V4.4): T0 exp(-r^2/(2 s0^2)) in the plane spreads as
// T0 (s0^2/s^2) exp(-r^2/(2 s^2)), s^2 = s0^2 + 2 a t, with conserved energy (docs/learning/phase4.md section 10).
// Reference: docs/REFERENCE_VALUES.md row P4-e: variance within 0.5%, energy 1e-10 relative.
// P4-h: L2 error of T at t_end, observed order T3: 2, T6: 3, +/- 0.15 over levels 1-4, with
// dt_l = dt_0 (h_l/h_0)^((p+1)/2) (section 10, "Separating time and space error"). The asserted slope is the
// combined space-time order along that path. The time part is about 2% (T3) and half (T6) of the error on
// every level (0.45 to 0.53 for T6, numerics review), so the slope is p + 1 only if both parts have their
// theoretical orders: a T6 spatial order of 2.5 would show as about 2.6.
// Quarter square [0, L]^2, L = 12 mm = 6 s_end, insulated (symmetry on x = 0, y = 0; the outer edges see
// T/T0 < 1e-7). Copper-like k = 400 W/(m K), rho c = 3.44e6 J/(m^3 K); s0 = 1 mm; s_end = 2 s0; Crank-Nicolson.
// Variance = int r^2 T dA / (2 int T dA), per coordinate. P4-e runs on the finest level (4), like P4-c/d.
// Disclosed change, approved after review: the first run used level 2 and gave a T3 variance error of
// 0.55% at the first check. The cause is the nodal interpolation of the initial Gaussian, which shifts the
// discrete variance by a constant O(h^2) amount (16.9, 4.38, 1.11, 0.28, 0.07% of s0^2 on levels 0-4); the
// growth 2 a t itself was exact. T6 was 3e-8 at level 2.

using namespace tmfem;

namespace {

constexpr double kTolVar = 0.005;    // P4-e: 0.5%
constexpr double kTolEnergy = 1e-10;  // P4-e
constexpr double kTolOrder = 0.15;   // P4-h
constexpr double kK = 400.0, kRhoC = 3.44e6, kA = kK / kRhoC, kS0 = 1e-3, kT0 = 100.0;
const double kTend = 1.5 * kS0 * kS0 / kA;  // s_end^2 = 4 s0^2

double exact(const Eigen::Vector2d& x, double t) {
  const double s2 = kS0 * kS0 + 2 * kA * t;
  return kT0 * kS0 * kS0 / s2 * std::exp(-x.squaredNorm() / (2 * s2));
}

Heat2DProblem gaussian_problem(int t_nodes, int level) {
  Heat2DProblem p;
  p.mesh = read_gmsh41(test::square_mesh(t_nodes, level));
  p.k_W_mK = kK;
  p.rho_c_J_m3K = kRhoC;
  return p;
}

// T_h at t_end with n steps of CN from the nodal interpolant; observer gets every step.
Eigen::VectorXd run(const Heat2DProblem& p, int n, const StepObserver& obs = nullptr) {
  const Eigen::VectorXd T0 = nodal_values(p.mesh, [](const Eigen::Vector2d& x) { return exact(x, 0.0); });
  return march(ThetaIntegrator(build_operators(p), 0.5, kTend / n), T0, 0.0, n, obs);
}

}  // namespace

TEST_CASE("P4-e Gaussian spreading: variance and energy", "[phase4][P4-e]") {
  for (const int t_nodes : {3, 6}) {
    const Heat2DProblem p = gaussian_problem(t_nodes, test::kP4Levels - 1);
    const ThermalOperators ops = build_operators(p);
    const Eigen::VectorXd one = Eigen::VectorXd::Ones(p.mesh.n_nodes());
    const int n = 90;
    double E0 = 0.0, worst_E = 0.0, worst_var = 0.0;
    std::vector<double> col_t, col_var, col_ref;
    const Eigen::VectorXd T_end = run(p, n, [&](double t, const Eigen::VectorXd& T) {
      const double E = one.dot(ops.M * T);  // M holds rho c: J per m of depth (T_ref = 0)
      if (t == 0.0) E0 = E;
      worst_E = std::max(worst_E, std::abs(E - E0) / E0);
      const int step = static_cast<int>(std::lround(t / (kTend / n)));
      if (step > 0 && step % 30 == 0) {
        const double var = integrate_field(p.mesh, "domain", T, [](const Eigen::Vector2d& x) { return x.squaredNorm(); }) /
                           (2 * integrate_field(p.mesh, "domain", T, [](const Eigen::Vector2d&) { return 1.0; }));
        const double ref = kS0 * kS0 + 2 * kA * t;
        const double rel = std::abs(var - ref) / ref;
        worst_var = std::max(worst_var, rel);
        col_t.push_back(t), col_var.push_back(var), col_ref.push_back(ref);
        std::printf("P4-e T%d t = %.4e s: variance %.6e m^2, exact %.6e, rel error %.2e\n", t_nodes, t, var, ref, rel);
      }
    });
    const std::string tag = "t" + std::to_string(t_nodes);
    write_csv(test::results_dir() / "phase4" / ("gaussian_variance_" + tag + ".csv"), {"t_s", "variance_m2", "exact_m2"},
              {col_t, col_var, col_ref}, {"P4-e Gaussian spreading, finest level"});
    std::vector<double> xs, Ts, Tx;  // nodes on y = 0 at t_end, for the figure
    for (int i = 0; i < p.mesh.n_nodes(); ++i) {
      if (p.mesh.x_m(i, 1) != 0.0) continue;
      xs.push_back(p.mesh.x_m(i, 0)), Ts.push_back(T_end(i)), Tx.push_back(exact(p.mesh.x_m.row(i).transpose(), kTend));
    }
    write_csv(test::results_dir() / "phase4" / ("gaussian_profile_" + tag + ".csv"), {"x_m", "T_fe_K", "T_exact_K"},
              {xs, Ts, Tx}, {"P4-e Gaussian spreading at t_end, nodes on y = 0 (unsorted)"});
    const double E_exact = kRhoC * kT0 * 2 * M_PI * kS0 * kS0 / 4;  // quarter of the plane
    std::printf("P4-e T%d: energy drift %.2e (discrete), initial energy vs closed form %.2e (interpolation, report)\n",
                t_nodes, worst_E, (E0 - E_exact) / E_exact);
    INFO("T" << t_nodes);
    CHECK(worst_var < kTolVar);
    CHECK(worst_E < kTolEnergy);
    test::record("P4-e", "T" + std::to_string(t_nodes) + " max rel variance error", 0.0, worst_var, "0.5%",
                 "P4-e Gaussian spreading: variance and energy");
    test::record("P4-e", "T" + std::to_string(t_nodes) + " energy drift", 0.0, worst_E, "1e-10 relative",
                 "P4-e Gaussian spreading: variance and energy");
  }
}

TEST_CASE("P4-h Gaussian spreading convergence", "[phase4][P4-h]") {
  constexpr int kN0 = 10;  // steps on level 0 (section 10)
  for (const auto& [t_nodes, order] : {std::pair{3, 2.0}, std::pair{6, 3.0}}) {
    const double expo = (t_nodes == 3 ? 2.0 : 3.0) / 2.0;  // (p + 1)/2
    std::vector<double> h, err, lv, steps;
    for (int level = 0; level < test::kP4Levels; ++level) {
      const Heat2DProblem p = gaussian_problem(t_nodes, level);
      const int n = static_cast<int>(std::lround(kN0 * std::pow(2.0, level * expo)));
      const Eigen::VectorXd T = run(p, n);
      const FieldError e = field_l2_error(p.mesh, "domain", T, [](const Eigen::Vector2d& x) { return exact(x, kTend); });
      h.push_back(e.h_m);
      err.push_back(e.L2 / e.L2_exact);
      lv.push_back(level);
      steps.push_back(n);
      if (level == 2) {  // time part of the error: (4/3) ||T(dt) - T(dt/2)|| for CN (section 10)
        const Eigen::VectorXd T2 = run(p, 2 * n);
        const FieldError d = field_l2_error(p.mesh, "domain", T - T2, [](const Eigen::Vector2d&) { return 0.0; });
        const double e_time = 4.0 / 3.0 * d.L2 / (e.L2 / err.back());  // relative to ||T_exact||
        std::printf("P4-h T%d level 2: error %.4e with %d steps; time part (4/3)||T(dt) - T(dt/2)|| = %.4e, "
                    "%.0f%% of the error\n", t_nodes, err.back(), n, e_time, 100.0 * e_time / err.back());
      }
    }
    const std::string label = "T" + std::to_string(t_nodes);
    test::print_orders(("P4-h " + label + ", relative L2 error of T at t_end").c_str(), h, err);
    const std::vector<double> hf(h.begin() + test::kP4FitFirstLevel, h.end()), ef(err.begin() + test::kP4FitFirstLevel, err.end());
    const double p_obs = test::observed_order(hf, ef);
    std::printf("  asserted order (levels 1-4): %.4f\n", p_obs);
    INFO(label);
    CHECK(std::abs(p_obs - order) <= kTolOrder);
    test::record("P4-h", label + " L2 order (levels 1-4)", order, p_obs, "+/- 0.15", "P4-h Gaussian spreading convergence");
    write_csv(test::results_dir() / "phase4" / ("gaussian_convergence_t" + std::to_string(t_nodes) + ".csv"),
              {"level", "steps", "h_m", "L2_rel"}, {lv, steps, h, err}, {"P4-h Gaussian spreading, Crank-Nicolson"});
  }
}
