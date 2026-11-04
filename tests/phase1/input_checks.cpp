#include <catch2/catch_test_macros.hpp>

#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "io/csv_reader.hpp"
#include "mesh/interval_mesh.hpp"
#include "support/paths.hpp"

// Malformed input must fail loudly, not produce a quietly wrong result. Not a benchmark: no reference row.

using namespace tmfem;

namespace {

std::filesystem::path write_file(const std::string& name, const std::string& text) {
  const auto dir = test::results_dir() / "phase1" / "input_checks";
  std::filesystem::create_directories(dir);
  const auto path = dir / name;
  std::ofstream(path) << text;
  return path;
}

}  // namespace

TEST_CASE("read_numeric_csv rejects a row with extra columns", "[phase1][io]") {
  const auto path = write_file("extra_column.csv", "a,b\n1,2\n3,4,5\n");
  CHECK_THROWS_AS(read_numeric_csv(path, 2), std::runtime_error);
}

TEST_CASE("read_numeric_csv treats a blank CRLF line as blank", "[phase1][io]") {
  const auto path = write_file("crlf.csv", "a,b\r\n1,2\r\n\r\n3,4\r\n");
  const auto cols = read_numeric_csv(path, 2);
  REQUIRE(cols[0].size() == 2);
  CHECK(cols[1][1] == 4.0);
}

TEST_CASE("interpolate rejects a nodal vector of the wrong size", "[phase1][mesh]") {
  const IntervalMesh mesh = uniform_interval(1.0, 4);
  const std::vector<double> too_short{0.0, 1.0};
  CHECK_THROWS_AS(interpolate(mesh, too_short, 0.9), std::invalid_argument);
}
