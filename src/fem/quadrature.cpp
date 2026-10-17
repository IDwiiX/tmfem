#include "fem/quadrature.hpp"

#include <cmath>
#include <stdexcept>

namespace tmfem {

QuadRule1D gauss_legendre(int n) {
  switch (n) {
    case 1:
      return {{0.0}, {2.0}};
    case 2: {
      const double a = 1.0 / std::sqrt(3.0);
      return {{-a, a}, {1.0, 1.0}};
    }
    case 3: {
      const double a = std::sqrt(3.0 / 5.0);
      return {{-a, 0.0, a}, {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0}};
    }
    case 4: {
      const double a = std::sqrt(3.0 / 7.0 - 2.0 / 7.0 * std::sqrt(6.0 / 5.0));
      const double b = std::sqrt(3.0 / 7.0 + 2.0 / 7.0 * std::sqrt(6.0 / 5.0));
      const double wa = (18.0 + std::sqrt(30.0)) / 36.0;
      const double wb = (18.0 - std::sqrt(30.0)) / 36.0;
      return {{-b, -a, a, b}, {wb, wa, wa, wb}};
    }
    case 5: {
      const double a = std::sqrt(5.0 - 2.0 * std::sqrt(10.0 / 7.0)) / 3.0;
      const double b = std::sqrt(5.0 + 2.0 * std::sqrt(10.0 / 7.0)) / 3.0;
      const double wa = (322.0 + 13.0 * std::sqrt(70.0)) / 900.0;
      const double wb = (322.0 - 13.0 * std::sqrt(70.0)) / 900.0;
      return {{-b, -a, 0.0, a, b}, {wb, wa, 128.0 / 225.0, wa, wb}};
    }
    default:
      throw std::invalid_argument("gauss_legendre: n must be 1..5");
  }
}

QuadRuleTri triangle_rule(int n_points) {
  switch (n_points) {
    case 1:
      return {{Eigen::Vector2d(1.0 / 3.0, 1.0 / 3.0)}, {0.5}};
    case 3:
      return {{Eigen::Vector2d(1.0 / 6.0, 1.0 / 6.0), Eigen::Vector2d(2.0 / 3.0, 1.0 / 6.0),
               Eigen::Vector2d(1.0 / 6.0, 2.0 / 3.0)},
              {1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0}};
    case 6: {
      // Strang-Fix / Dunavant degree-4 rule; weights for unit area, halved for the reference triangle.
      const double a = 0.445948490915965, wa = 0.5 * 0.223381589678011;
      const double b = 0.091576213509771, wb = 0.5 * 0.109951743655322;
      return {{Eigen::Vector2d(a, a), Eigen::Vector2d(1.0 - 2.0 * a, a), Eigen::Vector2d(a, 1.0 - 2.0 * a),
               Eigen::Vector2d(b, b), Eigen::Vector2d(1.0 - 2.0 * b, b), Eigen::Vector2d(b, 1.0 - 2.0 * b)},
              {wa, wa, wa, wb, wb, wb}};
    }
    default:
      throw std::invalid_argument("triangle_rule: n_points must be 1, 3 or 6");
  }
}

QuadRuleTri triangle_collapsed_gauss(int n) {
  const QuadRule1D g = gauss_legendre(n);
  QuadRuleTri rule;
  for (std::size_t i = 0; i < g.xi.size(); ++i) {
    const double u = 0.5 * (1.0 + g.xi[i]);  // [-1, 1] -> [0, 1], weight halved
    for (std::size_t j = 0; j < g.xi.size(); ++j) {
      const double v = 0.5 * (1.0 + g.xi[j]);
      rule.xi.emplace_back(u, v * (1.0 - u));
      rule.w.push_back(0.25 * g.w[i] * g.w[j] * (1.0 - u));
    }
  }
  return rule;
}

}  // namespace tmfem
