#pragma once

#include <optional>
#include <vector>

#include <Eigen/Dense>

#include "mechanics/elastic2d.hpp"

namespace tmfem {

// sigma = D B(xi) u_e in one cell of a region block (docs/learning/phase2.md section 8). Pa.
Eigen::Vector3d cell_stress(const Elastic2DProblem& p, const Eigen::VectorXd& u_m, const CellBlock& block, int cell,
                            const Eigen::Vector2d& xi);

// Stress at every stiffness quadrature point, cell by cell over the region blocks.
std::vector<Eigen::Vector3d> gauss_point_stresses(const Elastic2DProblem& p, const Eigen::VectorXd& u_m);

// Unweighted mean of D B(xi_node) u_e over the cells sharing each node (n_nodes x 3, Pa).
// For VTU output and plots only; validation uses stress_at_point.
Eigen::MatrixXd nodal_average_stress(const Elastic2DProblem& p, const Eigen::VectorXd& u_m);

// Stress at a physical point: D B(xi) u_e in the cell containing x, averaged over all containing cells
// when x is on a shared edge or node. std::nullopt if no cell contains x.
std::optional<Eigen::Vector3d> stress_at_point(const Elastic2DProblem& p, const Eigen::VectorXd& u_m,
                                               const Eigen::Vector2d& x_m);

}  // namespace tmfem
