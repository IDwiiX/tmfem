#pragma once

#include <functional>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "mesh/mesh.hpp"

namespace tmfem {

// Element kernels return the element matrix or vector for one cell, given its type and node
// coordinates (n x 2). Dofs are interleaved: dof = n_comp * node + comp (docs/learning/phase2.md section 3).
using CellMatrixKernel = std::function<Eigen::MatrixXd(CellType type, const Eigen::MatrixXd& x_e)>;
using CellVectorKernel = std::function<Eigen::VectorXd(CellType type, const Eigen::MatrixXd& x_e)>;

// Global dofs of one cell, in element order.
std::vector<int> cell_dofs(const CellBlock& block, int cell, int n_comp);

// Sum of element matrices over the given blocks (triplets pushed in block and cell order, so the
// result is deterministic). Size n_comp * n_nodes.
Eigen::SparseMatrix<double> assemble_matrix(const Mesh& mesh, const std::vector<const CellBlock*>& blocks,
                                            int n_comp, const CellMatrixKernel& kernel);

// Sum of element vectors over the given blocks (cells or boundary edges), added into f.
void assemble_vector(const Mesh& mesh, const std::vector<const CellBlock*>& blocks, int n_comp,
                     const CellVectorKernel& kernel, Eigen::VectorXd& f);

}  // namespace tmfem
