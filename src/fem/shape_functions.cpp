#include "fem/shape_functions.hpp"

#include <stdexcept>

namespace tmfem {

namespace {

ShapeEval line2(double xi) {
  ShapeEval s{Eigen::VectorXd(2), Eigen::MatrixXd(2, 1)};
  s.N << 0.5 * (1.0 - xi), 0.5 * (1.0 + xi);
  s.dN_dxi << -0.5, 0.5;
  return s;
}

ShapeEval line3(double xi) {
  ShapeEval s{Eigen::VectorXd(3), Eigen::MatrixXd(3, 1)};
  s.N << 0.5 * xi * (xi - 1.0), 0.5 * xi * (xi + 1.0), 1.0 - xi * xi;
  s.dN_dxi << xi - 0.5, xi + 0.5, -2.0 * xi;
  return s;
}

ShapeEval tri3() {
  ShapeEval s{Eigen::VectorXd(3), Eigen::MatrixXd(3, 2)};
  s.dN_dxi << -1.0, -1.0, 1.0, 0.0, 0.0, 1.0;
  return s;
}

// Area coordinates L1 = 1 - xi - eta, L2 = xi, L3 = eta, with dL1 = (-1, -1), dL2 = (1, 0), dL3 = (0, 1).
ShapeEval tri6(double xi, double eta) {
  const double L1 = 1.0 - xi - eta, L2 = xi, L3 = eta;
  ShapeEval s{Eigen::VectorXd(6), Eigen::MatrixXd(6, 2)};
  s.N << L1 * (2 * L1 - 1), L2 * (2 * L2 - 1), L3 * (2 * L3 - 1), 4 * L1 * L2, 4 * L2 * L3, 4 * L3 * L1;
  s.dN_dxi << -(4 * L1 - 1), -(4 * L1 - 1),  //
      4 * L2 - 1, 0.0,                         //
      0.0, 4 * L3 - 1,                         //
      4 * (L1 - L2), -4 * L2,                  //
      4 * L3, 4 * L2,                          //
      -4 * L3, 4 * (L1 - L3);
  return s;
}

}  // namespace

ShapeEval shape(CellType type, const Eigen::Vector2d& xi) {
  switch (type) {
    case CellType::Line2: return line2(xi(0));
    case CellType::Line3: return line3(xi(0));
    case CellType::Tri3: {
      ShapeEval s = tri3();
      s.N << 1.0 - xi(0) - xi(1), xi(0), xi(1);
      return s;
    }
    case CellType::Tri6: return tri6(xi(0), xi(1));
    case CellType::Point1: break;
  }
  throw std::invalid_argument("shape: unsupported cell type");
}

Eigen::MatrixXd reference_nodes(CellType type) {
  Eigen::MatrixXd r;
  switch (type) {
    case CellType::Line2:
      r.resize(2, 2);
      r << -1, 0, 1, 0;
      return r;
    case CellType::Line3:
      r.resize(3, 2);
      r << -1, 0, 1, 0, 0, 0;
      return r;
    case CellType::Tri3:
      r.resize(3, 2);
      r << 0, 0, 1, 0, 0, 1;
      return r;
    case CellType::Tri6:
      r.resize(6, 2);
      r << 0, 0, 1, 0, 0, 1, 0.5, 0, 0.5, 0.5, 0, 0.5;
      return r;
    case CellType::Point1: break;
  }
  throw std::invalid_argument("reference_nodes: unsupported cell type");
}

}  // namespace tmfem
