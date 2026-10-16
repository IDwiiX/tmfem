#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <stdexcept>

#include "fem/line2.hpp"
#include "fem/quadrature.hpp"
#include "io/timeseries.hpp"
#include "mesh/interval_mesh.hpp"
#include "support/paths.hpp"
#include "thermal/diffusion_analytic.hpp"
#include "thermal/heat1d.hpp"
#include "thermal/theta_integrator.hpp"

// Supporting unit tests for phase 1 (no reference row): they check the building blocks that the
// P1-* validation tests rely on, against closed forms from docs/learning/phase1.md.

using namespace tmfem;

TEST_CASE("element matrices match closed form", "[phase1]") {
  const double k = 3.0, rho_c = 5.0, x1 = 0.2, x2 = 0.7, h = x2 - x1;
  Eigen::Matrix2d K_ref;
  K_ref << 1, -1, -1, 1;
  K_ref *= k / h;
  Eigen::Matrix2d M_ref;
  M_ref << 2, 1, 1, 2;
  M_ref *= rho_c * h / 6.0;
  const Eigen::Matrix2d M_lump_ref = Eigen::Vector2d::Constant(rho_c * h / 2.0).asDiagonal();
  CHECK((element_stiffness(k, x1, x2) - K_ref).norm() < 1e-14);
  CHECK((element_capacity(rho_c, x1, x2, MassType::Consistent) - M_ref).norm() < 1e-14);
  CHECK((element_capacity(rho_c, x1, x2, MassType::Lumped) - M_lump_ref).norm() < 1e-14);
  const Eigen::Vector2d F = element_source([](double) { return 4.0; }, x1, x2);
  CHECK(F(0) == Catch::Approx(4.0 * h / 2.0).epsilon(1e-14));
  CHECK(F(1) == Catch::Approx(4.0 * h / 2.0).epsilon(1e-14));
}

TEST_CASE("Gauss-Legendre rules integrate polynomials of degree 2n-1 exactly", "[phase1]") {
  for (int n = 1; n <= 5; ++n) {
    const QuadRule1D q = gauss_legendre(n);
    const int p = 2 * n - 1;
    for (int d = 0; d <= p; ++d) {
      double s = 0.0;
      for (std::size_t g = 0; g < q.xi.size(); ++g) s += q.w[g] * std::pow(q.xi[g], d);
      const double exact = (d % 2 == 1) ? 0.0 : 2.0 / (d + 1);
      CHECK(std::abs(s - exact) < 1e-14);
    }
  }
}

TEST_CASE("graded mesh spans [0, L] with a constant length ratio", "[phase1]") {
  const IntervalMesh m = graded_interval(2.0, 10, 1.3);
  REQUIRE(m.n_nodes() == 11);
  CHECK(m.x_m.front() == 0.0);
  CHECK(m.x_m.back() == 2.0);
  for (int i = 1; i < 10; ++i) {
    const auto u = static_cast<std::size_t>(i);
    const double ratio = (m.x_m[u + 1] - m.x_m[u]) / (m.x_m[u] - m.x_m[u - 1]);
    CHECK(ratio == Catch::Approx(1.3).epsilon(1e-12));
  }
}

TEST_CASE("time series interpolates linearly and refuses to extrapolate", "[phase1]") {
  const auto path = test::results_dir() / "phase1" / "unit_timeseries.csv";
  std::filesystem::create_directories(path.parent_path());
  std::ofstream(path) << "# comment\nt_s,value\n0,1.0\n10,3.0\n20,-1.0\n";
  const TimeSeries ts = TimeSeries::read_csv(path);
  CHECK(ts(0.0) == 1.0);
  CHECK(ts(5.0) == Catch::Approx(2.0).epsilon(1e-15));
  CHECK(ts(15.0) == Catch::Approx(1.0).epsilon(1e-15));
  CHECK(ts(20.0) == -1.0);
  CHECK_THROWS_AS(ts(20.1), std::out_of_range);
  CHECK_THROWS_AS(ts(-0.1), std::out_of_range);
}

TEST_CASE("erfinv inverts std::erf", "[phase1]") {
  // Forward residual erf(erfinv(y)) - y is at roundoff for every y.
  for (double y : {-0.99999, -0.5, 1e-9, 0.3, 0.8421, 0.999999}) {
    CHECK(std::abs(std::erf(erfinv(y)) - y) <= 4e-16);
  }
  // Inverse error: d z / d y = sqrt(pi)/2 exp(z^2) amplifies the 1e-16 rounding of erf(z),
  // so the attainable accuracy in z is about 1e-16 exp(z^2): 1e-13 needs |z| <= 2.
  for (double z : {-2.0, -0.7, 1e-6, 0.3, 1.0, 2.0}) {
    CHECK(erfinv(std::erf(z)) == Catch::Approx(z).epsilon(1e-13));
  }
}

TEST_CASE("steady flux-in, Dirichlet-out gives the exact linear profile", "[phase1]") {
  // -k T'' = 0, k T'(0) = -q_in (heat enters at x = 0), T(L) = T_L: T = T_L + q_in (L - x) / k.
  const double k = 2.0, q_in = 50.0, L = 0.3, T_L = 300.0;
  Heat1DProblem p{graded_interval(L, 7, 1.2), k, 1.0, FluxBC{[q_in](double) { return q_in; }},
                  DirichletBC{[T_L](double) { return T_L; }}, {}};
  const Eigen::VectorXd T = solve_steady(build_operators(p, MassType::Consistent), 0.0);
  for (int i = 0; i < p.mesh.n_nodes(); ++i) {
    const double x = p.mesh.x_m[static_cast<std::size_t>(i)];
    CHECK(T(i) == Catch::Approx(T_L + q_in * (L - x) / k).epsilon(1e-13));
  }
}

TEST_CASE("steady Robin on both faces matches the resistance network", "[phase1]") {
  // Series resistances 1/h1 + L/k + 1/h2 carry q = (T1 - T2) / R; face temperatures follow.
  const double k = 60.0, L = 0.005, h1 = 22.0, h2 = 9.0, T1 = 340.0, T2 = 299.0;
  Heat1DProblem p{uniform_interval(L, 10), k, 1.0, ConvectiveBC{h1, [T1](double) { return T1; }},
                  ConvectiveBC{h2, [T2](double) { return T2; }}, {}};
  const Eigen::VectorXd T = solve_steady(build_operators(p, MassType::Consistent), 0.0);
  const double q = (T1 - T2) / (1.0 / h1 + L / k + 1.0 / h2);
  // The exact answer is in the FE space, so the only error is roundoff: about eps * cond(A), with
  // cond(A) ~ 1.7e5 here (lambda_max ~ 4k/h_element, lambda_min ~ (h1 + h2)/(n + 1)).
  constexpr double kRoundoff = 1e-11;
  CHECK(T(0) == Catch::Approx(T1 - q / h1).epsilon(kRoundoff));
  CHECK(T(10) == Catch::Approx(T2 + q / h2).epsilon(kRoundoff));
}

TEST_CASE("steady uniform source gives the exact parabola at the nodes", "[phase1]") {
  // -k T'' = q, T(0) = T(L) = 0: T = q x (L - x) / (2k); linear elements are nodally exact in 1D.
  const double k = 1.5, q = 8.0, L = 1.0;
  Heat1DProblem p{graded_interval(L, 9, 0.9), k, 1.0, DirichletBC{[](double) { return 0.0; }},
                  DirichletBC{[](double) { return 0.0; }}, [q](double, double) { return q; }};
  const Eigen::VectorXd T = solve_steady(build_operators(p, MassType::Consistent), 0.0);
  for (int i = 0; i < p.mesh.n_nodes(); ++i) {
    const double x = p.mesh.x_m[static_cast<std::size_t>(i)];
    CHECK(std::abs(T(i) - q * x * (L - x) / (2.0 * k)) < 1e-13);
  }
}
