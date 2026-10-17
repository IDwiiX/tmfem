#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>

#include "mesh/gmsh_reader.hpp"
#include "phase4/phase4_common.hpp"
#include "support/validation_record.hpp"
#include "thermal/heat2d.hpp"
#include "thermal/theta_integrator.hpp"

// P4-f 2D heat with source and convection (docs/learning/phase4.md sections 4 and 10).
// Reference: docs/REFERENCE_VALUES.md row P4-f, tolerance 1e-10 relative for all three parts:
//  (i)   steady slab on [0, L], uniform q, Robin h at x = 0 and x = L, insulated in y:
//        T = T_inf + qL/(2h) + q x (L - x)/(2k), in the T6 space on straight edges, so T6 is exact;
//  (ii)  discrete energy balance 1^T M (T^n - T^0) = sum dt [theta P^{k+1} + (1 - theta) P^k],
//        P = 1^T F - 1^T A T, over a transient with time-varying q and T_inf, a curved Robin rim and a face sink.
//        This is an identity of the theta step (multiply it by 1^T), true for any M, A, F: it checks the
//        integrator and the bookkeeping, not the integrals (numerics review). The same test therefore also
//        asserts load/matrix consistency, which can fail: with q = 0 and T_inf uniform in space, T = T_inf 1
//        solves the semi-discrete system exactly iff F(t) = A T_inf(t) 1 (curved Robin rim and face sink);
//  (iii) uniform plate cooled through its faces only: T^n - T_inf = G^n (T^0 - T_inf),
//        G = (1 - (1 - theta) lambda dt)/(1 + theta lambda dt), lambda = 2 h_f/(rho c t_z).
//        Posed in excess temperature (T_inf = 0, T0 = 100 K) with k = 0.4 W/(m K). G^n depends on neither.
//        Disclosed change, approved by the user: the first run used T_inf = 300 K, k = 400 W/(m K) and gave
//        1e-8 to 6e-6. Both causes are floating point, not discretisation: subtracting 300 K from T to get a
//        0.03 K excess loses about 4 digits, and k dt/(rho c h^2) = 5.6e4 makes conduction dominate a solve in
//        which it plays no role (K 1 = 0); the error grew in proportion to k. See docs/learning/phase4.md section 4.

using namespace tmfem;

namespace {

constexpr double kTol = 1e-10;  // P4-f
constexpr double kK = 1.0, kRhoC = 3.44e6;  // (i) uses a low k so the parabola is 18 K, not 0.4 K

}  // namespace

TEST_CASE("P4-f (i) steady slab with source and Robin faces", "[phase4][P4-f]") {
  const double L = test::kSquare_L_m, q = 1e6, h = 100.0, T_inf = 300.0;
  for (const int t : {3, 6}) {
    Heat2DProblem p;
    p.mesh = read_gmsh41(test::square_mesh(t, 1));
    p.k_W_mK = kK;
    p.rho_c_J_m3K = kRhoC;
    p.source_W_m3 = [q](const Eigen::Vector2d&, double) { return q; };
    p.convection = {{"left", h, [T_inf](double) { return T_inf; }}, {"right", h, [T_inf](double) { return T_inf; }}};
    const Eigen::VectorXd T = solve_steady(build_operators(p), 0.0);
    double err = 0.0;
    for (int i = 0; i < p.mesh.n_nodes(); ++i) {
      const double x = p.mesh.x_m(i, 0);
      err = std::max(err, std::abs(T(i) - (T_inf + q * L / (2 * h) + q * x * (L - x) / (2 * kK))));
    }
    const double rel = err / (q * L / (2 * h) + q * L * L / (8 * kK));
    std::printf("P4-f (i) T%d: max |T - T_exact| / max(T - T_inf) = %.2e%s\n", t, rel, t == 3 ? " (T3: not exact, report)" : "");
    if (t == 6) {
      CHECK(rel < kTol);
      test::record("P4-f", "(i) steady slab T6 max rel error", 0.0, rel, "1e-10 relative", "P4-f (i) steady slab with source and Robin faces");
    }
  }
}

TEST_CASE("P4-f (ii) discrete energy balance with source, curved Robin rim and face sink", "[phase4][P4-f]") {
  const double b = test::kDisk_b_m, k = 400.0, tau = 1.0;
  const auto T_inf = [tau](double t) { return 293.15 + 20.0 * std::sin(2 * M_PI * t / tau); };
  for (const int t_nodes : {3, 6}) {
    Heat2DProblem p;
    p.mesh = read_gmsh41(test::disk_mesh(t_nodes, 2));
    p.k_W_mK = k;
    p.rho_c_J_m3K = kRhoC;
    p.source_W_m3 = [b, tau](const Eigen::Vector2d& x, double t) {
      return 1e8 * std::exp(-x.squaredNorm() / (2 * (b / 4) * (b / 4))) * (1 + 0.5 * std::cos(2 * M_PI * t / tau));
    };
    p.convection = {{"rim", 50.0, T_inf}};
    p.face_cooling = FaceCooling{10.0, 0.005, T_inf};
    const ThermalOperators ops = build_operators(p);
    const Eigen::VectorXd one = Eigen::VectorXd::Ones(p.mesh.n_nodes());
    {  // load/matrix consistency: F_no_source(t) = A T_inf(t) 1, relative to sum |A_ij| T_inf
      Heat2DProblem p0 = p;
      p0.source_W_m3 = nullptr;
      const ThermalOperators ops0 = build_operators(p0);
      const double sum_abs = ops0.A.cwiseAbs().sum();
      double worst = 0.0;
      for (const double t : {0.0, 0.13, 0.4, 0.77}) {
        worst = std::max(worst, (ops0.F(t) - ops0.A * (T_inf(t) * one)).lpNorm<1>() / (sum_abs * T_inf(t)));
      }
      std::printf("P4-f (ii) T%d: max |F - A T_inf 1|_1 / (sum|A_ij| T_inf) = %.2e (Robin rim + face sink loads)\n",
                  t_nodes, worst);
      CHECK(worst < 1e-13);
    }
    for (const double theta : {1.0, 0.5}) {
      const int n = 50;
      const double dt = tau / n;
      const ThetaIntegrator integ(ops, theta, dt);
      const Eigen::VectorXd T0 = Eigen::VectorXd::Constant(p.mesh.n_nodes(), 293.15);
      const auto power = [&](double t, const Eigen::VectorXd& T) { return one.dot(ops.F(t)) - one.dot(ops.A * T); };
      double sum = 0.0, scale = 0.0, worst = 0.0;
      Eigen::VectorXd T = T0;
      for (int s = 0; s < n; ++s) {
        const double t0 = s * dt;
        const Eigen::VectorXd T1 = integ.step(T, t0);
        const double P0 = power(t0, T), P1 = power(t0 + dt, T1);
        sum += dt * (theta * P1 + (1 - theta) * P0);
        scale += dt * (std::abs(one.dot(ops.F(t0 + dt))) + std::abs(one.dot(ops.A * T1)));
        T = T1;
        worst = std::max(worst, std::abs(one.dot(ops.M * (T - T0)) - sum));
      }
      const double rel = worst / scale;
      std::printf("P4-f (ii) T%d theta %.1f: stored energy %.6e J/m, max |balance residual| / sum |power dt| = %.2e\n",
                  t_nodes, theta, one.dot(ops.M * (T - T0)), rel);
      INFO("T" << t_nodes << " theta " << theta);
      CHECK(rel < kTol);
      test::record("P4-f", "(ii) energy balance T" + std::to_string(t_nodes) + " theta " + std::to_string(theta).substr(0, 3),
                   0.0, rel, "1e-10 relative", "P4-f (ii) discrete energy balance with source, curved Robin rim and face sink");
    }
  }
}

TEST_CASE("P4-f (iii) uniform plate cooled through its faces: amplification G^n", "[phase4][P4-f]") {
  const double h_f = 10.0, t_z = 0.002, T_inf = 0.0, T0 = 100.0;
  const double lambda = 2 * h_f / (kRhoC * t_z);
  for (const int t_nodes : {3, 6}) {
    Heat2DProblem p;
    p.mesh = read_gmsh41(test::square_mesh(t_nodes, 1));
    p.k_W_mK = 0.4;
    p.rho_c_J_m3K = kRhoC;
    p.face_cooling = FaceCooling{h_f, t_z, [T_inf](double) { return T_inf; }};
    const ThermalOperators ops = build_operators(p);
    for (const double theta : {1.0, 0.5}) {
      const double dt = 0.5 / lambda;
      const int n = 20;
      const double G = (1 - (1 - theta) * lambda * dt) / (1 + theta * lambda * dt);
      const Eigen::VectorXd T = march(ThetaIntegrator(ops, theta, dt),
                                      Eigen::VectorXd::Constant(p.mesh.n_nodes(), T0), 0.0, n, nullptr);
      const double ref = std::pow(G, n) * (T0 - T_inf);
      const double rel = (T.array() - T_inf - ref).abs().maxCoeff() / std::abs(ref);
      std::printf("P4-f (iii) T%d theta %.1f: G = %.10f, max |T - T_inf - G^n (T0 - T_inf)| / |G^n (T0 - T_inf)| = %.2e\n",
                  t_nodes, theta, G, rel);
      INFO("T" << t_nodes << " theta " << theta);
      CHECK(rel < kTol);
      test::record("P4-f", "(iii) amplification T" + std::to_string(t_nodes) + " theta " + std::to_string(theta).substr(0, 3),
                   0.0, rel, "1e-10 relative", "P4-f (iii) uniform plate cooled through its faces: amplification G^n");
    }
  }
}
