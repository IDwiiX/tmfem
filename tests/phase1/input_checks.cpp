#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "io/case_heat1d.hpp"
#include "io/csv_reader.hpp"
#include "mesh/interval_mesh.hpp"
#include "support/paths.hpp"
#include "thermal/heat1d.hpp"
#include "thermal/theta_integrator.hpp"

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

// Well-posed slab: insulated left, fixed right.
Heat1DProblem slab() {
  Heat1DProblem p;
  p.mesh = uniform_interval(1.0, 4);
  p.k = 1.0;
  p.rho_c = 1.0;
  p.left = FluxBC{[](double) { return 0.0; }};
  p.right = DirichletBC{[](double) { return 300.0; }};
  return p;
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

TEST_CASE("run_heat1d_case rejects a non-positive time grid", "[phase1][io]") {
  Heat1DCase c;
  c.problem = slab();
  c.T0_K = 300.0;
  c.theta = 1.0;
  c.mass = MassType::Consistent;
  c.dt_s = 1.0;
  c.t_end_s = -10.0;  // otherwise a "multiple of dt" that runs zero steps and reports success
  CHECK_THROWS_AS(run_heat1d_case(c), std::invalid_argument);
  c.dt_s = 0.0;
  c.t_end_s = 0.0;  // otherwise t_end / dt = NaN goes into lround
  CHECK_THROWS_AS(run_heat1d_case(c), std::invalid_argument);
}

TEST_CASE("1D build_operators rejects non-positive k and rho_c", "[phase1][thermal]") {
  Heat1DProblem p = slab();
  p.k = 0.0;
  CHECK_THROWS_AS(build_operators(p, MassType::Consistent), std::invalid_argument);
  p = slab();
  p.rho_c = -1.0;
  CHECK_THROWS_AS(build_operators(p, MassType::Consistent), std::invalid_argument);
}

TEST_CASE("ThetaIntegrator rejects NaN theta and a non-finite dt", "[phase1][thermal]") {
  const ThermalOperators ops = build_operators(slab(), MassType::Consistent);
  CHECK_THROWS_AS(ThetaIntegrator(ops, std::nan(""), 1.0), std::invalid_argument);
  CHECK_THROWS_AS(ThetaIntegrator(ops, 1.0, std::nan("")), std::invalid_argument);
  CHECK_THROWS_AS(ThetaIntegrator(ops, 1.0, HUGE_VAL), std::invalid_argument);
}
