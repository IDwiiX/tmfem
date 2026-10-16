#pragma once

#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace tmfem {

using SparseMatrix = Eigen::SparseMatrix<double>;

// Splits dofs into free (f) and constrained (c) sets for Dirichlet elimination:
// S_ff x_f = b_f - S_fc g (docs/learning/phase1.md section 6). Dimension-agnostic.
class DirichletPartition {
 public:
  DirichletPartition(int n_dofs, std::vector<int> constrained);

  int n_dofs() const { return n_dofs_; }
  int n_free() const { return static_cast<int>(free_.size()); }
  const std::vector<int>& free_dofs() const { return free_; }
  const std::vector<int>& constrained_dofs() const { return constrained_; }

  SparseMatrix block_ff(const SparseMatrix& S) const;
  SparseMatrix block_fc(const SparseMatrix& S) const;
  Eigen::VectorXd restrict_free(const Eigen::VectorXd& full) const;
  // Full vector from free values and constrained values.
  Eigen::VectorXd combine(const Eigen::VectorXd& x_f, const Eigen::VectorXd& g) const;

 private:
  SparseMatrix block(const SparseMatrix& S, const std::vector<int>& col_index, int n_cols) const;

  int n_dofs_;
  std::vector<int> free_;
  std::vector<int> constrained_;
  std::vector<int> free_index_;         // dof -> position in free_, or -1
  std::vector<int> constrained_index_;  // dof -> position in constrained_, or -1
};

}  // namespace tmfem
