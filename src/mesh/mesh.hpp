#pragma once

#include <map>
#include <string>
#include <vector>

#include <Eigen/Dense>

namespace tmfem {

// Cell types used in 2D (docs/learning/phase2.md section 5). Node order follows Gmsh:
// Line3 = (end, end, mid); Tri6 = (3 corners, then mid-sides of edges 01, 12, 20).
enum class CellType { Point1, Line2, Line3, Tri3, Tri6 };

int nodes_per_cell(CellType type);
int cell_dim(CellType type);

using NodeMatrix = Eigen::Matrix<double, Eigen::Dynamic, 2, Eigen::RowMajor>;
using ConnMatrix = Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

// Cells of one type belonging to one physical group; `nodes` holds 0-based node indices, one row per cell.
struct CellBlock {
  CellType type;
  int physical_tag;
  ConnMatrix nodes;
};

struct PhysicalGroup {
  int dim;
  int tag;
};

// 2D mesh: node coordinates in m, cell blocks, and physical group names.
struct Mesh {
  NodeMatrix x_m;
  std::vector<CellBlock> blocks;
  std::map<std::string, PhysicalGroup> physical;

  int n_nodes() const { return static_cast<int>(x_m.rows()); }
  // Blocks of the named physical group; throws if the name is unknown.
  std::vector<const CellBlock*> blocks_of(const std::string& name) const;
  // Sorted, unique node indices of the named group.
  std::vector<int> nodes_of(const std::string& name) const;
  // Coordinates of one cell's nodes (n x 2).
  Eigen::MatrixXd cell_coordinates(const CellBlock& block, int cell) const;
};

}  // namespace tmfem
