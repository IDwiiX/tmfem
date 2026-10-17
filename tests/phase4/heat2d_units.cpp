#include <catch2/catch_test_macros.hpp>

#include <cmath>

#include "mesh/gmsh_reader.hpp"
#include "phase4/phase4_common.hpp"
#include "thermal/heat2d.hpp"
#include "thermal/theta_integrator.hpp"

// Unit checks of the 2D heat operators (docs/learning/phase4.md section 3): partition of unity gives
// 1^T M 1 = rho c area, 1^T H 1 = h length, 1^T F_q = int q dA; K annihilates constants; M, A symmetric.
// Quarter square [0, L]^2 (straight edges), T3 and T6.

using namespace tmfem;

TEST_CASE("2D heat operators: sums, constants and symmetry", "[phase4][units]") {
  const double L = test::kSquare_L_m, rho_c = 3.5e6, k = 50.0, h = 25.0;
  for (const int t : {3, 6}) {
    Heat2DProblem p;
    p.mesh = read_gmsh41(test::square_mesh(t, 1));
    p.k_W_mK = k;
    p.rho_c_J_m3K = rho_c;
    p.source_W_m3 = [](const Eigen::Vector2d& x, double) { return 1e6 * x(0) * x(0) * x(1); };  // degree 3
    const ThermalOperators no_bc = build_operators(p);
    p.convection.push_back({"right", h, [](double) { return 300.0; }});
    const ThermalOperators ops = build_operators(p);

    const Eigen::VectorXd one = Eigen::VectorXd::Ones(p.mesh.n_nodes());
    INFO("T" << t);
    CHECK(std::abs(one.dot(ops.M * one) - rho_c * L * L) <= 1e-12 * rho_c * L * L);
    // Round-off in K 1 and 1^T A 1 comes from cancelling entries of size k, so it scales with sum |A_ij|.
    const double sum_abs = ops.A.cwiseAbs().sum();
    CHECK((no_bc.A * one).lpNorm<1>() <= 1e-14 * sum_abs);                        // K 1 = 0
    CHECK(std::abs(one.dot(ops.A * one) - h * L) <= 1e-14 * sum_abs);             // only H survives
    const double int_q = 1e6 * (L * L * L / 3) * (L * L / 2);
    CHECK(std::abs(one.dot(no_bc.F(0.0)) - int_q) <= 1e-12 * int_q);              // 6-point rule, degree 3
    CHECK(std::abs(one.dot(ops.F(0.0)) - int_q - h * 300.0 * L) <= 1e-12 * (int_q + h * 300.0 * L));
    const Eigen::SparseMatrix<double> Mt = ops.M.transpose(), At = ops.A.transpose();
    CHECK((ops.M - Mt).norm() <= 1e-14 * ops.M.norm());
    CHECK((ops.A - At).norm() <= 1e-14 * ops.A.norm());
  }
}

TEST_CASE("2D heat: imposed temperature and imposed flux reproduce exact 1D profiles", "[phase4][units]") {
  // Quarter square, insulated top and bottom, so the solution depends on x only.
  // (a) T = T_D on x = 0 and x = L, uniform q: T = T_D + q x (L - x)/(2k), quadratic, exact in T6.
  // (b) flux q_in into x = 0, T = T_D on x = L, no source: T = T_D + q_in (L - x)/k, linear, exact in T3 and T6.
  const double L = test::kSquare_L_m, k = 2.0, q = 1e6, q_in = 5e3, T_D = 300.0;
  for (const int t : {3, 6}) {
    Heat2DProblem a;
    a.mesh = read_gmsh41(test::square_mesh(t, 1));
    a.k_W_mK = k;
    a.rho_c_J_m3K = 1e6;
    a.source_W_m3 = [q](const Eigen::Vector2d&, double) { return q; };
    a.temperatures = {{"left", [T_D](double) { return T_D; }}, {"right", [T_D](double) { return T_D; }}};
    Heat2DProblem b = a;
    b.source_W_m3 = nullptr;
    b.temperatures = {{"right", [T_D](double) { return T_D; }}};
    b.fluxes = {{"left", [q_in](double) { return q_in; }}};
    const Eigen::VectorXd Ta = solve_steady(build_operators(a), 0.0), Tb = solve_steady(build_operators(b), 0.0);
    double ea = 0.0, eb = 0.0;
    for (int i = 0; i < a.mesh.n_nodes(); ++i) {
      const double x = a.mesh.x_m(i, 0);
      ea = std::max(ea, std::abs(Ta(i) - (T_D + q * x * (L - x) / (2 * k))) / (q * L * L / (8 * k)));
      eb = std::max(eb, std::abs(Tb(i) - (T_D + q_in * (L - x) / k)) / (q_in * L / k));
    }
    INFO("T" << t << ": Dirichlet + source " << ea << ", flux " << eb);
    if (t == 6) CHECK(ea < 1e-10);
    CHECK(eb < 1e-10);
  }
}
