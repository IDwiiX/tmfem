#include <catch2/catch_test_macros.hpp>

#include <stdexcept>

#include "mesh/gmsh_reader.hpp"
#include "support/paths.hpp"
#include "thermal/heat2d.hpp"

// Malformed input must fail loudly, not produce a quietly wrong result. Not a benchmark: no reference row.

using namespace tmfem;

TEST_CASE("Face cooling needs a positive thickness", "[phase4][thermal]") {
  Heat2DProblem p;
  p.mesh = read_gmsh41(test::source_dir() / "cases/bar/bar_t3.msh");
  p.k_W_mK = 1.0;
  p.rho_c_J_m3K = 1.0;
  p.face_cooling = FaceCooling{10.0, 0.0, [](double) { return 300.0; }};
  CHECK_THROWS_AS(build_operators(p), std::invalid_argument);
}
