#include "fem/element_geometry.hpp"

#include <stdexcept>

namespace tmfem {

namespace {

constexpr double kInsideSlack = 1e-10;
constexpr int kMaxNewton = 25;

}  // namespace

CellGeometry cell_geometry(const ShapeEval& s, const Eigen::MatrixXd& x_e) {
  const Eigen::Matrix2d J = x_e.transpose() * s.dN_dxi;  // J_ij = dx_i / dxi_j
  const double detJ = J.determinant();
  if (!(detJ > 0.0)) throw std::runtime_error("cell_geometry: det J <= 0 (inverted or degenerate element)");
  return {x_e.transpose() * s.N, detJ, s.dN_dxi * J.inverse()};
}

FacetGeometry facet_geometry(const ShapeEval& s, const Eigen::MatrixXd& x_e) {
  const Eigen::Vector2d dx_dxi = x_e.transpose() * s.dN_dxi.col(0);
  return {x_e.transpose() * s.N, dx_dxi.norm()};
}

Eigen::Vector2d map_to_physical(CellType type, const Eigen::MatrixXd& x_e, const Eigen::Vector2d& xi) {
  return x_e.transpose() * shape(type, xi).N;
}

std::optional<Eigen::Vector2d> inverse_map(CellType type, const Eigen::MatrixXd& x_e, const Eigen::Vector2d& x) {
  // Cheap rejection: outside the node bounding box enlarged by 10% (curved T6 edges can bulge past it).
  const Eigen::Vector2d lo = x_e.colwise().minCoeff();
  const Eigen::Vector2d hi = x_e.colwise().maxCoeff();
  const Eigen::Vector2d pad = 0.1 * (hi - lo);
  if ((x.array() < (lo - pad).array()).any() || (x.array() > (hi + pad).array()).any()) return std::nullopt;

  Eigen::Vector2d xi(1.0 / 3.0, 1.0 / 3.0);
  const double scale = (hi - lo).norm();
  for (int it = 0; it < kMaxNewton; ++it) {
    const ShapeEval s = shape(type, xi);
    const Eigen::Vector2d r = x - x_e.transpose() * s.N;
    const Eigen::Matrix2d J = x_e.transpose() * s.dN_dxi;
    const Eigen::Vector2d dxi = J.lu().solve(r);
    xi += dxi;
    if (dxi.norm() < 1e-14 && r.norm() < 1e-13 * scale) break;
  }
  if ((map_to_physical(type, x_e, xi) - x).norm() > 1e-10 * scale) {
    return std::nullopt;  // Newton did not converge (point far outside a curved cell)
  }
  const double L1 = 1.0 - xi(0) - xi(1);
  if (xi(0) < -kInsideSlack || xi(1) < -kInsideSlack || L1 < -kInsideSlack) return std::nullopt;
  return xi;
}

}  // namespace tmfem
