#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <vector>

#include "fem/quadrature.hpp"
#include "fem/shape_functions.hpp"

// Shape functions and triangle quadrature (docs/learning/phase2.md sections 5 and 6). Not a reference
// row: partition of unity, Kronecker delta, derivative and exactness checks.

using namespace tmfem;

namespace {

const std::vector<Eigen::Vector2d> kTriPoints{{0.2, 0.1}, {0.6, 0.3}, {0.05, 0.85}, {1.0 / 3.0, 1.0 / 3.0}};
const std::vector<Eigen::Vector2d> kLinePoints{{-0.7, 0.0}, {0.1, 0.0}, {0.9, 0.0}};

void check_type(CellType type, const std::vector<Eigen::Vector2d>& points) {
  const int dim = cell_dim(type);
  for (const auto& xi : points) {
    const ShapeEval s = shape(type, xi);
    CHECK(std::abs(s.N.sum() - 1.0) < 1e-14);                     // partition of unity
    CHECK(s.dN_dxi.colwise().sum().cwiseAbs().maxCoeff() < 1e-13);  // gradients of 1 vanish
    for (int d = 0; d < dim; ++d) {                               // central differences
      const double h = 1e-6;
      Eigen::Vector2d xp = xi, xm = xi;
      xp(d) += h;
      xm(d) -= h;
      const Eigen::VectorXd fd = (shape(type, xp).N - shape(type, xm).N) / (2 * h);
      CHECK((fd - s.dN_dxi.col(d)).cwiseAbs().maxCoeff() < 1e-8);
    }
  }
  const Eigen::MatrixXd nodes = reference_nodes(type);
  for (Eigen::Index a = 0; a < nodes.rows(); ++a) {  // Kronecker delta at the nodes
    const Eigen::VectorXd N = shape(type, nodes.row(a).transpose()).N;
    for (Eigen::Index b = 0; b < N.size(); ++b) CHECK(std::abs(N(b) - (a == b ? 1.0 : 0.0)) < 1e-14);
  }
}

double factorial(int n) { return n <= 1 ? 1.0 : n * factorial(n - 1); }

// Max relative error of the rule over all monomials xi^p eta^q with p + q <= degree.
double max_monomial_error(const QuadRuleTri& rule, int degree) {
  double worst = 0.0;
  for (int p = 0; p <= degree; ++p) {
    for (int q = 0; p + q <= degree; ++q) {
      const double exact = factorial(p) * factorial(q) / factorial(p + q + 2);
      double sum = 0.0;
      for (std::size_t i = 0; i < rule.w.size(); ++i) sum += rule.w[i] * std::pow(rule.xi[i](0), p) * std::pow(rule.xi[i](1), q);
      worst = std::max(worst, std::abs(sum - exact) / exact);
    }
  }
  return worst;
}

}  // namespace

TEST_CASE("Shape functions: Line2, Line3, Tri3, Tri6", "[phase2][fem]") {
  check_type(CellType::Line2, kLinePoints);
  check_type(CellType::Line3, kLinePoints);
  check_type(CellType::Tri3, kTriPoints);
  check_type(CellType::Tri6, kTriPoints);
}

TEST_CASE("Triangle quadrature exactness", "[phase2][fem]") {
  CHECK(max_monomial_error(triangle_rule(1), 1) < 1e-14);
  CHECK(max_monomial_error(triangle_rule(3), 2) < 1e-14);
  CHECK(max_monomial_error(triangle_rule(3), 3) > 1e-3);  // degree 2 only: the claim is not loose
  for (int n = 1; n <= 5; ++n) {
    INFO("collapsed Gauss n = " << n);
    CHECK(max_monomial_error(triangle_collapsed_gauss(n), 2 * n - 2) < 1e-13);
  }
}
