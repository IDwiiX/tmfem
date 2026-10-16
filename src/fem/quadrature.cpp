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

}  // namespace tmfem
