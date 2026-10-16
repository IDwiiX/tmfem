#include "fem/line2.hpp"

#include "fem/quadrature.hpp"

namespace tmfem {

Eigen::Vector2d line2_shape(double xi) {
  return {0.5 * (1.0 - xi), 0.5 * (1.0 + xi)};
}

double line2_coordinate(double x1_m, double x2_m, double xi) {
  const Eigen::Vector2d N = line2_shape(xi);
  return N(0) * x1_m + N(1) * x2_m;
}

Eigen::Matrix2d element_stiffness(double k, double x1_m, double x2_m) {
  const double h = x2_m - x1_m;
  const double detJ = 0.5 * h;
  const Eigen::Vector2d dN_dx(-1.0 / h, 1.0 / h);
  Eigen::Matrix2d K = Eigen::Matrix2d::Zero();
  const QuadRule1D q = gauss_legendre(2);
  for (std::size_t g = 0; g < q.xi.size(); ++g) {
    K += k * dN_dx * dN_dx.transpose() * q.w[g] * detJ;
  }
  return K;
}

Eigen::Matrix2d element_capacity(double rho_c, double x1_m, double x2_m, MassType mass) {
  const double detJ = 0.5 * (x2_m - x1_m);
  Eigen::Matrix2d M = Eigen::Matrix2d::Zero();
  const QuadRule1D q = gauss_legendre(2);
  for (std::size_t g = 0; g < q.xi.size(); ++g) {
    const Eigen::Vector2d N = line2_shape(q.xi[g]);
    M += rho_c * N * N.transpose() * q.w[g] * detJ;
  }
  if (mass == MassType::Lumped) {
    const Eigen::Vector2d row_sums = M.rowwise().sum();
    M = row_sums.asDiagonal();
  }
  return M;
}

Eigen::Vector2d element_source(const std::function<double(double x_m)>& q_fn, double x1_m, double x2_m) {
  const double detJ = 0.5 * (x2_m - x1_m);
  Eigen::Vector2d F = Eigen::Vector2d::Zero();
  const QuadRule1D q = gauss_legendre(3);
  for (std::size_t g = 0; g < q.xi.size(); ++g) {
    const double x = line2_coordinate(x1_m, x2_m, q.xi[g]);
    F += q_fn(x) * line2_shape(q.xi[g]) * q.w[g] * detJ;
  }
  return F;
}

}  // namespace tmfem
