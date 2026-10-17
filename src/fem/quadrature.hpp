#pragma once

#include <vector>

#include <Eigen/Dense>

namespace tmfem {

// Quadrature rule on the reference interval [-1, 1].
struct QuadRule1D {
  std::vector<double> xi;
  std::vector<double> w;
};

// n-point Gauss-Legendre rule (n = 1..5), exact for polynomials of degree 2n - 1.
QuadRule1D gauss_legendre(int n);

// Quadrature rule on the reference triangle {xi, eta >= 0, xi + eta <= 1}; weights sum to 1/2.
struct QuadRuleTri {
  std::vector<Eigen::Vector2d> xi;
  std::vector<double> w;
};

// Stiffness rules (docs/learning/phase2.md section 6): 1 point (degree 1) or 3 points (degree 2).
QuadRuleTri triangle_rule(int n_points);

// Collapsed (Duffy) Gauss rule with n x n points, xi = u, eta = v (1 - u), built from gauss_legendre(n).
// Exact for polynomials of total degree 2n - 2 on the triangle (the Jacobian 1 - u adds one power of u);
// used for error norms.
QuadRuleTri triangle_collapsed_gauss(int n);

}  // namespace tmfem
