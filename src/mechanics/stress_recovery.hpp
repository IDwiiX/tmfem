#pragma once

#include <optional>
#include <vector>

#include <Eigen/Dense>

#include "mechanics/elastic2d.hpp"

namespace tmfem {

// sigma = D B(xi) u_e - beta dT_h(xi) m in one cell of a region block (docs/learning/phase2.md section 8,
// phase4.md section 9). Pa.
Eigen::Vector3d cell_stress(const Elastic2DProblem& p, const Eigen::VectorXd& u_m, const CellBlock& block, int cell,
                            const Eigen::Vector2d& xi);

// Stress at every stiffness quadrature point, cell by cell over the region blocks.
std::vector<Eigen::Vector3d> gauss_point_stresses(const Elastic2DProblem& p, const Eigen::VectorXd& u_m);

// Unweighted mean of the cell stress (D B u_e - beta dT_h m) over the cells sharing each node (n_nodes x 3, Pa).
// For VTU output and plots only; validation uses stress_at_point.
Eigen::MatrixXd nodal_average_stress(const Elastic2DProblem& p, const Eigen::VectorXd& u_m);

// Nodal output fields (plots and VTU): nodal-averaged in-plane stress, sigma_zz with the thermal term
// (docs/learning/phase4.md section 9) and von Mises, all in Pa. dT is continuous at nodes, so sigma_zz from
// the averaged stress and the nodal dT equals the average of the per-cell sigma_zz.
struct NodalStressFields {
  Eigen::MatrixXd sigma_Pa;  // n_nodes x 3, [xx, yy, xy]
  Eigen::VectorXd sigma_zz_Pa;
  Eigen::VectorXd von_mises_Pa;
};
NodalStressFields nodal_stress_fields(const Elastic2DProblem& p, const Eigen::VectorXd& u_m);

// Stress at a physical point: the cell stress (D B u_e - beta dT_h m) in the cell containing x, averaged over all containing cells
// when x is on a shared edge or node. std::nullopt if no cell contains x.
std::optional<Eigen::Vector3d> stress_at_point(const Elastic2DProblem& p, const Eigen::VectorXd& u_m,
                                               const Eigen::Vector2d& x_m);

}  // namespace tmfem
