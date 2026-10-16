#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstdio>
#include <string>

#include "fem/line2.hpp"
#include "fem/quadrature.hpp"
#include "io/csv_writer.hpp"
#include "support/convergence.hpp"
#include "support/paths.hpp"
#include "support/validation_record.hpp"
#include "thermal/heat1d.hpp"
#include "thermal/theta_integrator.hpp"

// P1-m observed orders: BE time 1, CN time 2, space 2, each +/- 0.15 on four refinements
// (docs/REFERENCE_VALUES.md row P1-m; method in docs/learning/phase1.md section 10).
// Case: T_t = alpha T_xx on (0, 1), T = 0 at both ends, T(x, 0) = sin(pi x); exact
// T = exp(-pi^2 alpha t) sin(pi x). On a uniform mesh the nodal sine is an exact discrete eigenvector.

using namespace tmfem;

namespace {

constexpr double kTol = 0.15;  // REFERENCE_VALUES P1-m
constexpr double kL = 1.0;
constexpr double kAlpha = 1.0;
constexpr double kTend = 0.1;
const double kPi = std::acos(-1.0);

Heat1DProblem sine_problem(int n_elements) {
  return {uniform_interval(kL, n_elements), kAlpha, 1.0, DirichletBC{[](double) { return 0.0; }},
          DirichletBC{[](double) { return 0.0; }}, {}};
}

Eigen::VectorXd sine(const IntervalMesh& m) {
  return nodal_values(m, [](double x) { return std::sin(kPi * x / kL); });
}

// Discrete eigenvalue lambda_h = v^T A v / v^T M v (exact Rayleigh quotient for an eigenvector).
double discrete_eigenvalue(const ThermalOperators& ops, const Eigen::VectorXd& v) {
  return v.dot(ops.A * v) / v.dot(ops.M * v);
}

// Time error alone: FE solution vs the semi-discrete exact exp(-lambda_h t) v (zero spatial error).
double time_error(double theta, double dt, MassType mass) {
  const Heat1DProblem p = sine_problem(16);
  const ThetaIntegrator integ(build_operators(p, mass), theta, dt);
  const Eigen::VectorXd v = sine(p.mesh);
  const double lambda_h = discrete_eigenvalue(integ.operators(), v);
  const int n = static_cast<int>(std::lround(kTend / dt));
  const Eigen::VectorXd T = march(integ, v, 0.0, n, {});
  return (T - std::exp(-lambda_h * kTend) * v).cwiseAbs().maxCoeff();
}

// L2 error of the FE function against the exact PDE solution, by 3-point Gauss per element.
double l2_error(const IntervalMesh& m, const Eigen::VectorXd& T, double t) {
  const QuadRule1D q = gauss_legendre(3);
  double sum = 0.0;
  for (int e = 0; e < m.n_elements(); ++e) {
    const auto ue = static_cast<std::size_t>(e);
    const double x1 = m.x_m[ue], x2 = m.x_m[ue + 1];
    for (std::size_t g = 0; g < q.xi.size(); ++g) {
      const Eigen::Vector2d N = line2_shape(q.xi[g]);
      const double x = line2_coordinate(x1, x2, q.xi[g]);
      const double exact = std::exp(-kPi * kPi * kAlpha * t / (kL * kL)) * std::sin(kPi * x / kL);
      const double d = N(0) * T(e) + N(1) * T(e + 1) - exact;
      sum += d * d * q.w[g] * 0.5 * (x2 - x1);
    }
  }
  return std::sqrt(sum);
}

void check_time_order(double theta, double expected, const char* id_quantity, const char* test_name) {
  const std::vector<double> dts{0.01, 0.005, 0.0025, 0.00125};
  std::vector<double> err;
  for (double dt : dts) err.push_back(time_error(theta, dt, MassType::Consistent));
  test::print_orders(id_quantity, dts, err);
  const double p = test::observed_order(dts, err);
  test::record("P1-m", id_quantity, expected, p, "+/- 0.15", test_name);
  CHECK(std::abs(p - expected) <= kTol);
  write_csv(test::results_dir() / (std::string("phase1/convergence_time_theta") + (theta == 1.0 ? "1" : "0.5") + ".csv"),
            {"dt_s", "max_nodal_error"}, {dts, err}, {test_name});
}

void check_space_order(MassType mass, const char* id_quantity, const char* test_name, const char* csv) {
  const std::vector<int> ns{8, 16, 32, 64};
  constexpr double kDt = 1e-4;  // CN; its time error is checked below to be < 1% of the spatial error
  std::vector<double> hs, err;
  for (int n : ns) {
    const Heat1DProblem p = sine_problem(n);
    const ThetaIntegrator cn(build_operators(p, mass), 0.5, kDt);
    const Eigen::VectorXd v = sine(p.mesh);
    const Eigen::VectorXd T = march(cn, v, 0.0, static_cast<int>(std::lround(kTend / kDt)), {});
    const double e = l2_error(p.mesh, T, kTend);
    // Time error of this run, exactly: |G(lambda_h dt)^n - exp(-lambda_h t)| times the L2 norm of the sine
    // interpolant (at most sqrt(L/2)).
    const double lh = discrete_eigenvalue(cn.operators(), v);
    const double G = (1.0 - 0.5 * lh * kDt) / (1.0 + 0.5 * lh * kDt);
    const double e_time = std::abs(std::pow(G, std::lround(kTend / kDt)) - std::exp(-lh * kTend)) * std::sqrt(kL / 2.0);
    INFO("n = " << n << ": L2 error " << e << ", time part " << e_time);
    CHECK(e_time < 0.01 * e);
    hs.push_back(kL / n);
    err.push_back(e);
  }
  test::print_orders(id_quantity, hs, err);
  const double p = test::observed_order(hs, err);
  test::record("P1-m", id_quantity, 2.0, p, "+/- 0.15", test_name);
  CHECK(std::abs(p - 2.0) <= kTol);
  write_csv(test::results_dir() / csv, {"h_m", "l2_error"}, {hs, err}, {test_name});
}

// Forced problem: T = cos t + x^2 on (0, 1) with source, time-dependent Dirichlet at x = 0 and Robin with
// time-dependent T_inf at x = 1. This exercises the forcing path theta F^{n+1} + (1 - theta) F^n and
// g(t_{n+1}) that the unforced sine case never touches. T_t is uniform in x and x^2 is nodally exact, so
// the semi-discrete solution is the nodal interpolant (both mass types): the error is pure time error.
constexpr double kForcedK = 0.7, kForcedRhoC = 1.3, kForcedH = 5.0, kForcedTend = 1.0;

double forced_time_error(double theta, double dt, MassType mass) {
  const Heat1DProblem p{uniform_interval(1.0, 8), kForcedK, kForcedRhoC,
                        DirichletBC{[](double t) { return std::cos(t); }},
                        ConvectiveBC{kForcedH, [](double t) { return std::cos(t) + 1.0 + 2.0 * kForcedK / kForcedH; }},
                        [](double, double t) { return -kForcedRhoC * std::sin(t) - 2.0 * kForcedK; }};
  const ThetaIntegrator integ(build_operators(p, mass), theta, dt);
  const auto exact = [&p](double t) { return nodal_values(p.mesh, [t](double x) { return std::cos(t) + x * x; }); };
  const Eigen::VectorXd T = march(integ, exact(0.0), 0.0, static_cast<int>(std::lround(kForcedTend / dt)), {});
  return (T - exact(kForcedTend)).cwiseAbs().maxCoeff();
}

void check_forced_time_order(double theta, double expected, const char* test_name) {
  const std::vector<double> dts{0.1, 0.05, 0.025, 0.0125};
  for (MassType mass : {MassType::Consistent, MassType::Lumped}) {
    const std::string label = std::string(theta == 1.0 ? "forced time order BE" : "forced time order CN") +
                              (mass == MassType::Consistent ? " consistent mass" : " lumped mass");
    std::vector<double> err;
    for (double dt : dts) err.push_back(forced_time_error(theta, dt, mass));
    test::print_orders(label.c_str(), dts, err);
    const double p = test::observed_order(dts, err);
    test::record("P1-m", label, expected, p, "+/- 0.15", test_name);
    CHECK(std::abs(p - expected) <= kTol);
  }
}

}  // namespace

TEST_CASE("P1-m temporal order with time-dependent forcing, backward Euler", "[phase1][P1-m]") {
  check_forced_time_order(1.0, 1.0, "P1-m temporal order with time-dependent forcing, backward Euler");
}

TEST_CASE("P1-m temporal order with time-dependent forcing, Crank-Nicolson", "[phase1][P1-m]") {
  check_forced_time_order(0.5, 2.0, "P1-m temporal order with time-dependent forcing, Crank-Nicolson");
}

TEST_CASE("P1-m temporal order backward Euler", "[phase1][P1-m]") {
  check_time_order(1.0, 1.0, "time order BE", "P1-m temporal order backward Euler");
}

TEST_CASE("P1-m temporal order Crank-Nicolson", "[phase1][P1-m]") {
  check_time_order(0.5, 2.0, "time order CN", "P1-m temporal order Crank-Nicolson");
}

TEST_CASE("P1-m spatial order consistent mass", "[phase1][P1-m]") {
  check_space_order(MassType::Consistent, "space order consistent mass", "P1-m spatial order consistent mass",
                    "phase1/convergence_space_consistent.csv");
}

TEST_CASE("P1-m spatial order lumped mass", "[phase1][P1-m]") {
  check_space_order(MassType::Lumped, "space order lumped mass", "P1-m spatial order lumped mass",
                    "phase1/convergence_space_lumped.csv");
}
