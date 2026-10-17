#pragma once

#include <filesystem>

#include "mesh/mesh.hpp"

namespace tmfem {

// Reads a Gmsh ASCII msh 4.1 file: $PhysicalNames, $Entities (entity -> physical tags), $Nodes and
// $Elements. Only elements of entities with a physical tag are kept, one block per (entity, physical
// tag). Node tags are mapped to contiguous 0-based indices in increasing tag order; z is dropped and
// must be 0. Triangles that Gmsh wrote clockwise are reordered to counter-clockwise.
// Supported element types: 15 (point), 1 (line2), 8 (line3), 2 (tri3), 9 (tri6). Anything else throws.
Mesh read_gmsh41(const std::filesystem::path& path);

}  // namespace tmfem
