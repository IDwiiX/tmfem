#include <catch2/catch_test_macros.hpp>

#include <fstream>
#include <sstream>
#include <string>

#include "io/vtu_writer.hpp"
#include "mesh/gmsh_reader.hpp"
#include "phase2/kirsch_common.hpp"
#include "support/paths.hpp"

// VTU writer (V2.3 output for ParaView). Not a reference row: checks counts and cell types.

using namespace tmfem;

namespace {

int count_occurrences(const std::string& text, const std::string& token) {
  int n = 0;
  for (auto pos = text.find(token); pos != std::string::npos; pos = text.find(token, pos + 1)) ++n;
  return n;
}

}  // namespace

TEST_CASE("VTU writer: counts and quadratic cell type", "[phase2][io]") {
  const Mesh m = read_gmsh41(test::kirsch_mesh(6, 0));
  const auto cells = m.blocks_of("domain");
  Eigen::Index n_cells = 0;
  for (const CellBlock* b : cells) n_cells += b->nodes.rows();
  const auto path = test::results_dir() / "phase2" / "vtu_writer_test.vtu";
  write_vtu(path, m, cells, {{"node_x", m.x_m.col(0)}});

  std::ifstream in(path);
  std::stringstream ss;
  ss << in.rdbuf();
  const std::string text = ss.str();
  CHECK(text.find("NumberOfPoints=\"" + std::to_string(m.n_nodes()) + "\"") != std::string::npos);
  CHECK(text.find("NumberOfCells=\"" + std::to_string(n_cells) + "\"") != std::string::npos);
  CHECK(count_occurrences(text, "\n22\n") + 1 >= static_cast<int>(n_cells));  // VTK_QUADRATIC_TRIANGLE
  CHECK(count_occurrences(text, "Name=\"node_x\"") == 1);
  CHECK(text.find("</VTKFile>") != std::string::npos);
}
