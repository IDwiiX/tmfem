#include "fem/assembly.hpp"

#include <stdexcept>

namespace tmfem {

std::vector<int> cell_dofs(const CellBlock& block, int cell, int n_comp) {
  const auto n = static_cast<int>(block.nodes.cols());
  std::vector<int> dofs(static_cast<std::size_t>(n * n_comp));
  for (int a = 0; a < n; ++a) {
    for (int c = 0; c < n_comp; ++c) dofs[static_cast<std::size_t>(n_comp * a + c)] = n_comp * block.nodes(cell, a) + c;
  }
  return dofs;
}

Eigen::SparseMatrix<double> assemble_matrix(const Mesh& mesh, const std::vector<const CellBlock*>& blocks,
                                            int n_comp, const CellMatrixKernel& kernel) {
  std::vector<Eigen::Triplet<double>> triplets;
  for (const CellBlock* b : blocks) {
    const auto n_dof_e = static_cast<std::size_t>(b->nodes.cols() * n_comp);
    triplets.reserve(triplets.size() + static_cast<std::size_t>(b->nodes.rows()) * n_dof_e * n_dof_e);
    for (int e = 0; e < b->nodes.rows(); ++e) {
      const Eigen::MatrixXd Ke = kernel(b->type, mesh.cell_coordinates(*b, e));
      const std::vector<int> dofs = cell_dofs(*b, e, n_comp);
      if (static_cast<std::size_t>(Ke.rows()) != dofs.size() || Ke.rows() != Ke.cols()) {
        throw std::logic_error("assemble_matrix: element matrix size does not match its dofs");
      }
      for (std::size_t i = 0; i < dofs.size(); ++i) {
        for (std::size_t j = 0; j < dofs.size(); ++j) {
          triplets.emplace_back(dofs[i], dofs[j], Ke(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)));
        }
      }
    }
  }
  const int n = n_comp * mesh.n_nodes();
  Eigen::SparseMatrix<double> K(n, n);
  K.setFromTriplets(triplets.begin(), triplets.end());
  return K;
}

void assemble_vector(const Mesh& mesh, const std::vector<const CellBlock*>& blocks, int n_comp,
                     const CellVectorKernel& kernel, Eigen::VectorXd& f) {
  for (const CellBlock* b : blocks) {
    for (int e = 0; e < b->nodes.rows(); ++e) {
      const Eigen::VectorXd fe = kernel(b->type, mesh.cell_coordinates(*b, e));
      const std::vector<int> dofs = cell_dofs(*b, e, n_comp);
      if (static_cast<std::size_t>(fe.size()) != dofs.size()) {
        throw std::logic_error("assemble_vector: element vector size does not match its dofs");
      }
      for (std::size_t i = 0; i < dofs.size(); ++i) f(dofs[i]) += fe(static_cast<Eigen::Index>(i));
    }
  }
}

}  // namespace tmfem
