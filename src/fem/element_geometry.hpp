#pragma once

#include <optional>

#include <Eigen/Dense>

#include "fem/shape_functions.hpp"

namespace tmfem {

// Isoparametric geometry of a triangle at one quadrature point (docs/learning/phase2.md section 5):
// physical point x, det J, and shape-function gradients dN/dx = J^-T dN/dxi (n x 2).
struct CellGeometry {
  Eigen::Vector2d x;
  double detJ;
  Eigen::MatrixXd dN_dx;
};

// x_e: node coordinates (n x 2). Throws if det J <= 0 (inverted or degenerate element).
CellGeometry cell_geometry(const ShapeEval& s, const Eigen::MatrixXd& x_e);

// Geometry of a boundary edge at one quadrature point: x and |dx/dxi| (so ds = ds_dxi dxi).
struct FacetGeometry {
  Eigen::Vector2d x;
  double ds_dxi;
};

FacetGeometry facet_geometry(const ShapeEval& s, const Eigen::MatrixXd& x_e);

// Physical position of reference point xi in a cell.
Eigen::Vector2d map_to_physical(CellType type, const Eigen::MatrixXd& x_e, const Eigen::Vector2d& xi);

// Reference coordinates of x in a triangle, by Newton on x(xi) = x, if x lies inside the cell
// (all area coordinates >= -1e-10); std::nullopt otherwise (section 8).
std::optional<Eigen::Vector2d> inverse_map(CellType type, const Eigen::MatrixXd& x_e, const Eigen::Vector2d& x);

}  // namespace tmfem
