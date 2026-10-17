#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <fstream>

#include "fem/element_geometry.hpp"
#include "fem/shape_functions.hpp"
#include "mechanics/elasticity.hpp"
#include "mesh/gmsh_reader.hpp"
#include "phase2/kirsch_common.hpp"
#include "support/paths.hpp"

// Gmsh msh 4.1 reader (docs/learning/phase2.md section 5). Not a reference row: these checks make the
// Kirsch results trustworthy (curved hole, orientation, physical groups).

using namespace tmfem;

TEST_CASE("Gmsh reader: Kirsch T6 coarse mesh", "[phase2][mesh]") {
  const Mesh m = read_gmsh41(test::kirsch_mesh(6, 0));
  CHECK(m.n_nodes() == 437);  // header of kirsch_t6_l0.msh: 437 nodes
  Eigen::Index n_cells = 0;
  for (const auto& b : m.blocks) n_cells += b.nodes.rows();
  CHECK(n_cells == 236);  // header: 236 elements (triangles + boundary edges)
  for (const char* name : {"sym_x", "sym_y", "outer_right", "outer_top", "hole", "domain"}) {
    CHECK(m.physical.count(name) == 1);
  }

  for (const CellBlock* b : m.blocks_of("domain")) CHECK(b->type == CellType::Tri6);
  for (const CellBlock* b : m.blocks_of("hole")) CHECK(b->type == CellType::Line3);

  // Every hole node, mid-side nodes included, lies on r = a: the T6 boundary is curved.
  const double a = test::kKirschA_m;
  for (const int node : m.nodes_of("hole")) {
    CHECK(std::abs(m.x_m.row(node).norm() - a) <= 1e-12 * a);
  }
  for (const int node : m.nodes_of("sym_x")) CHECK(m.x_m(node, 0) == 0.0);
  for (const int node : m.nodes_of("sym_y")) CHECK(m.x_m(node, 1) == 0.0);

  // det J > 0 at every stiffness point (cell_geometry throws otherwise).
  int n_checked = 0;
  for (const CellBlock* b : m.blocks_of("domain")) {
    const QuadRuleTri rule = stiffness_rule(b->type);
    for (int e = 0; e < b->nodes.rows(); ++e) {
      for (const auto& xi : rule.xi) {
        CHECK(cell_geometry(shape(b->type, xi), m.cell_coordinates(*b, e)).detJ > 0.0);
        ++n_checked;
      }
    }
  }
  CHECK(n_checked > 0);
}

TEST_CASE("Gmsh reader: hole nodes lie on r = a on every T6 level", "[phase2][mesh]") {
  const double a = test::kKirschA_m;
  for (int level = 0; level < test::kKirschLevels; ++level) {
    const Mesh m = read_gmsh41(test::kirsch_mesh(6, level));
    double worst = 0.0;
    for (const int node : m.nodes_of("hole")) worst = std::max(worst, std::abs(m.x_m.row(node).norm() - a) / a);
    INFO("level " << level);
    CHECK(worst <= 1e-12);
  }
}

TEST_CASE("Gmsh reader: bar mesh with a physical point", "[phase2][mesh]") {
  const Mesh m = read_gmsh41(test::source_dir() / "cases/bar/bar_t3.msh");
  const auto corner = m.nodes_of("corner");
  REQUIRE(corner.size() == 1);
  CHECK(m.x_m(corner[0], 0) == 0.0);
  CHECK(m.x_m(corner[0], 1) == 0.0);
  for (const int node : m.nodes_of("left")) CHECK(m.x_m(node, 0) == 0.0);
  for (const int node : m.nodes_of("right")) CHECK(std::abs(m.x_m(node, 0) - 0.1) < 1e-15);
}

TEST_CASE("Gmsh reader: clockwise triangles are reordered counter-clockwise", "[phase2][mesh]") {
  // One Tri6 written clockwise: corners (0,0), (0,1), (1,0) and mid-sides of edges 01, 12, 20.
  const auto path = test::results_dir() / "phase2" / "clockwise_tri6.msh";
  std::filesystem::create_directories(path.parent_path());
  {
    std::ofstream out(path);
    out << "$MeshFormat\n4.1 0 8\n$EndMeshFormat\n"
        << "$PhysicalNames\n1\n2 1 \"domain\"\n$EndPhysicalNames\n"
        << "$Entities\n0 0 1 0\n1 0 0 0 1 1 0 1 1 0\n$EndEntities\n"
        << "$Nodes\n1 6 1 6\n2 1 0 6\n1\n2\n3\n4\n5\n6\n"
        << "0 0 0\n0 1 0\n1 0 0\n0 0.5 0\n0.5 0.5 0\n0.5 0 0\n$EndNodes\n"
        << "$Elements\n1 1 1 1\n2 1 9 1\n1 1 2 3 4 5 6\n$EndElements\n";
  }
  const Mesh m = read_gmsh41(path);
  const CellBlock& b = *m.blocks_of("domain").at(0);
  const Eigen::MatrixXd x_e = m.cell_coordinates(b, 0);
  CHECK(cell_geometry(shape(CellType::Tri6, {1.0 / 3.0, 1.0 / 3.0}), x_e).detJ > 0.0);
  // After reordering each mid-side node is still the midpoint of its edge (3: 01, 4: 12, 5: 20).
  CHECK((x_e.row(3) - 0.5 * (x_e.row(0) + x_e.row(1))).norm() < 1e-15);
  CHECK((x_e.row(4) - 0.5 * (x_e.row(1) + x_e.row(2))).norm() < 1e-15);
  CHECK((x_e.row(5) - 0.5 * (x_e.row(2) + x_e.row(0))).norm() < 1e-15);
}
