#pragma once

#include <vector>

namespace tmfem {

// Quadrature rule on the reference interval [-1, 1].
struct QuadRule1D {
  std::vector<double> xi;
  std::vector<double> w;
};

// n-point Gauss-Legendre rule (n = 1..5), exact for polynomials of degree 2n - 1.
QuadRule1D gauss_legendre(int n);

}  // namespace tmfem
