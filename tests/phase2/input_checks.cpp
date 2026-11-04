#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>

#include "mesh/gmsh_reader.hpp"
#include "support/paths.hpp"

// Malformed input must fail loudly, not produce a quietly wrong result. Not a benchmark: no reference row.

using namespace tmfem;
using Catch::Matchers::ContainsSubstring;

namespace {

std::filesystem::path write_file(const std::string& name, const std::string& text) {
  const auto dir = test::results_dir() / "phase2" / "input_checks";
  std::filesystem::create_directories(dir);
  const auto path = dir / name;
  std::ofstream(path) << text;
  return path;
}

}  // namespace

TEST_CASE("Gmsh reader names a truncated node block", "[phase2][mesh]") {
  const auto path = write_file("truncated.msh", "$MeshFormat\n4.1 0 8\n$EndMeshFormat\n$Nodes\n1 3 1 3\n2 1 0 3\n1\n2\n");
  CHECK_THROWS_WITH(read_gmsh41(path), ContainsSubstring("truncated"));
}

TEST_CASE("Gmsh reader requires a $MeshFormat section", "[phase2][mesh]") {
  const auto path = write_file("no_format.msh", "$Nodes\n1 1 1 1\n0 1 0 1\n1\n0 0 0\n$EndNodes\n");
  CHECK_THROWS_WITH(read_gmsh41(path), ContainsSubstring("$MeshFormat"));
}
