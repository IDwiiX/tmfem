#pragma once

#include <Eigen/Dense>

#include "mesh/mesh.hpp"

namespace tmfem {

// Shape functions at one reference point: N (n) and dN/dxi (n x reference dim; 1 for lines, 2 for
// triangles). Lines live on xi in [-1, 1] (xi(0) is used); triangles on {xi, eta >= 0, xi + eta <= 1}.
// Formulas and node order: docs/learning/phase2.md section 5.
struct ShapeEval {
  Eigen::VectorXd N;
  Eigen::MatrixXd dN_dxi;
};

ShapeEval shape(CellType type, const Eigen::Vector2d& xi);

// Reference coordinates of the nodes (n x 2; second column 0 for lines).
Eigen::MatrixXd reference_nodes(CellType type);

}  // namespace tmfem
