#include "mesh/mesh.hpp"

#include <algorithm>
#include <stdexcept>

namespace tmfem {

int nodes_per_cell(CellType type) {
  switch (type) {
    case CellType::Point1: return 1;
    case CellType::Line2: return 2;
    case CellType::Line3: return 3;
    case CellType::Tri3: return 3;
    case CellType::Tri6: return 6;
  }
  throw std::invalid_argument("nodes_per_cell: unknown cell type");
}

int cell_dim(CellType type) {
  switch (type) {
    case CellType::Point1: return 0;
    case CellType::Line2:
    case CellType::Line3: return 1;
    case CellType::Tri3:
    case CellType::Tri6: return 2;
  }
  throw std::invalid_argument("cell_dim: unknown cell type");
}

std::vector<const CellBlock*> Mesh::blocks_of(const std::string& name) const {
  const auto it = physical.find(name);
  if (it == physical.end()) throw std::invalid_argument("Mesh: unknown physical group '" + name + "'");
  std::vector<const CellBlock*> out;
  for (const auto& b : blocks) {
    if (b.physical_tag == it->second.tag && cell_dim(b.type) == it->second.dim) out.push_back(&b);
  }
  return out;
}

std::vector<int> Mesh::nodes_of(const std::string& name) const {
  std::vector<int> nodes;
  for (const CellBlock* b : blocks_of(name)) {
    nodes.insert(nodes.end(), b->nodes.data(), b->nodes.data() + b->nodes.size());
  }
  std::sort(nodes.begin(), nodes.end());
  nodes.erase(std::unique(nodes.begin(), nodes.end()), nodes.end());
  return nodes;
}

Eigen::MatrixXd Mesh::cell_coordinates(const CellBlock& block, int cell) const {
  const auto n = static_cast<Eigen::Index>(block.nodes.cols());
  Eigen::MatrixXd x_e(n, 2);
  for (Eigen::Index a = 0; a < n; ++a) x_e.row(a) = x_m.row(block.nodes(cell, a));
  return x_e;
}

}  // namespace tmfem
