#pragma once

#include <filesystem>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include "mesh/mesh.hpp"

namespace tmfem {

// Nodal field for VTU output: one row per mesh node, 1 to 3 columns.
struct VtuPointField {
  std::string name;
  Eigen::MatrixXd values;
};

// ASCII VTK UnstructuredGrid for ParaView. Cells: VTK_TRIANGLE (5) for Tri3 and
// VTK_QUADRATIC_TRIANGLE (22) for Tri6, whose node order matches Gmsh. Creates parent directories.
void write_vtu(const std::filesystem::path& path, const Mesh& mesh, const std::vector<const CellBlock*>& cells,
               const std::vector<VtuPointField>& point_data);

}  // namespace tmfem
