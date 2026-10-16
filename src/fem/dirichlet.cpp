#include "fem/dirichlet.hpp"

#include <algorithm>
#include <stdexcept>

namespace tmfem {

DirichletPartition::DirichletPartition(int n_dofs, std::vector<int> constrained)
    : n_dofs_(n_dofs),
      constrained_(std::move(constrained)),
      free_index_(static_cast<std::size_t>(n_dofs), -1),
      constrained_index_(static_cast<std::size_t>(n_dofs), -1) {
  for (std::size_t i = 0; i < constrained_.size(); ++i) {
    const int d = constrained_[i];
    if (d < 0 || d >= n_dofs || constrained_index_[static_cast<std::size_t>(d)] != -1) {
      throw std::invalid_argument("DirichletPartition: invalid or repeated constrained dof");
    }
    constrained_index_[static_cast<std::size_t>(d)] = static_cast<int>(i);
  }
  for (int d = 0; d < n_dofs; ++d) {
    if (constrained_index_[static_cast<std::size_t>(d)] == -1) {
      free_index_[static_cast<std::size_t>(d)] = static_cast<int>(free_.size());
      free_.push_back(d);
    }
  }
}

SparseMatrix DirichletPartition::block(const SparseMatrix& S, const std::vector<int>& col_index,
                                       int n_cols) const {
  std::vector<Eigen::Triplet<double>> triplets;
  for (int col = 0; col < S.outerSize(); ++col) {
    const int jc = col_index[static_cast<std::size_t>(col)];
    if (jc < 0) continue;
    for (SparseMatrix::InnerIterator it(S, col); it; ++it) {
      const int ir = free_index_[static_cast<std::size_t>(it.row())];
      if (ir >= 0) triplets.emplace_back(ir, jc, it.value());
    }
  }
  SparseMatrix B(n_free(), n_cols);
  B.setFromTriplets(triplets.begin(), triplets.end());
  return B;
}

SparseMatrix DirichletPartition::block_ff(const SparseMatrix& S) const {
  return block(S, free_index_, n_free());
}

SparseMatrix DirichletPartition::block_fc(const SparseMatrix& S) const {
  return block(S, constrained_index_, static_cast<int>(constrained_.size()));
}

Eigen::VectorXd DirichletPartition::restrict_free(const Eigen::VectorXd& full) const {
  Eigen::VectorXd x_f(n_free());
  for (std::size_t i = 0; i < free_.size(); ++i) x_f(static_cast<Eigen::Index>(i)) = full(free_[i]);
  return x_f;
}

Eigen::VectorXd DirichletPartition::combine(const Eigen::VectorXd& x_f, const Eigen::VectorXd& g) const {
  Eigen::VectorXd full(n_dofs_);
  for (std::size_t i = 0; i < free_.size(); ++i) full(free_[i]) = x_f(static_cast<Eigen::Index>(i));
  for (std::size_t i = 0; i < constrained_.size(); ++i) {
    full(constrained_[i]) = g(static_cast<Eigen::Index>(i));
  }
  return full;
}

}  // namespace tmfem
