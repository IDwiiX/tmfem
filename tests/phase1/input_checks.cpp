#include <catch2/catch_test_macros.hpp>

#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>

#include "io/csv_reader.hpp"
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
