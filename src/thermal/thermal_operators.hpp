#pragma once

#include <functional>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace tmfem {

using TimeFunction = std::function<double(double t_s)>;

// Semi-discrete system M dT/dt + A T = F(t) with Dirichlet dofs (docs/learning/phase1.md section 3).
// Dimension-agnostic: built from a 1D mesh (heat1d) or a 2D mesh (heat2d, docs/learning/phase4.md section 3).
struct ThermalOperators {
  Eigen::SparseMatrix<double> M;
  Eigen::SparseMatrix<double> A;  // K + H (+ face sink in 2D)
  std::function<Eigen::VectorXd(double t_s)> F;
  std::vector<int> dirichlet_dofs;
  std::function<Eigen::VectorXd(double t_s)> dirichlet_values;
};

}  // namespace tmfem
