#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>

#include "io/case_thermomech2d.hpp"
#include "mesh/gmsh_reader.hpp"
#include "support/paths.hpp"
#include "thermal/heat2d.hpp"

// Malformed input must fail loudly, not produce a quietly wrong result. Not a benchmark: no reference row.

using namespace tmfem;
using Catch::Matchers::ContainsSubstring;

namespace {

// A valid thermomech2d case file; each test replaces one fragment before loading.
std::string beam_case(const std::string& old_text = "", const std::string& new_text = "") {
  std::string json = R"({
    "case_type": "thermomech2d",
    "mesh_pattern": "block_s{s_mm}.msh",
    "beam_sizes_mm": [1],
    "material": {"name": "m", "rho_kg_m3": 1, "c_J_kgK": 1, "k_W_mK": 1, "alpha_per_K": 1e-6, "E_GPa": 1, "nu": 0.2,
                 "strength": {"type": "brittle", "flexural_MPa": 60, "compressive_MPa": 130}},
    "plane": "plane_strain",
    "beam": {"E_per_length_J_m": 1, "pulse_s": 1e-5, "pulse_steps": 10},
    "cooling": {"physical": "outer", "h_W_m2K": 10, "T_inf_K": 293, "stage_ends_s": [1e-4, 1e-3], "steps_per_stage": 40},
    "T0_K": 293, "T_ref_K": 293, "quarter_model_factor": 4, "section_area_m2": 0.0016,
    "supports": [], "output_dir": "out"
  })";
  if (!old_text.empty()) {
    const auto at = json.find(old_text);
    REQUIRE(at != std::string::npos);
    json.replace(at, old_text.size(), new_text);
  }
  return json;
}

BeamBlockCase load(const std::string& json) {
  const auto dir = test::results_dir() / "phase4" / "input_checks";
  std::filesystem::create_directories(dir);
  const auto path = dir / "case.json";
  std::ofstream(path) << json;
  return load_thermomech2d_case(path, dir);
}

}  // namespace

TEST_CASE("Face cooling needs a positive thickness", "[phase4][thermal]") {
  Heat2DProblem p;
  p.mesh = read_gmsh41(test::source_dir() / "cases/bar/bar_t3.msh");
  p.k_W_mK = 1.0;
  p.rho_c_J_m3K = 1.0;
  p.face_cooling = FaceCooling{10.0, 0.0, [](double) { return 300.0; }};
  CHECK_THROWS_AS(build_operators(p), std::invalid_argument);
}

TEST_CASE("thermomech2d loader accepts the valid base case", "[phase4][io]") { CHECK_NOTHROW(load(beam_case())); }

TEST_CASE("thermomech2d loader rejects an unknown strength type", "[phase4][io]") {
  CHECK_THROWS_WITH(load(beam_case("\"brittle\"", "\"britle\"")), ContainsSubstring("strength type"));
}

TEST_CASE("thermomech2d loader rejects bad cooling stage ends", "[phase4][io]") {
  const std::string stages = "[1e-4, 1e-3]";
  CHECK_THROWS_WITH(load(beam_case(stages, "[]")), ContainsSubstring("stage_ends_s"));  // run_beam_size calls front()
  CHECK_THROWS_WITH(load(beam_case(stages, "[1e-3, 1e-3]")), ContainsSubstring("stage_ends_s"));
  CHECK_THROWS_WITH(load(beam_case(stages, "[1e-6, 1e-3]")), ContainsSubstring("stage_ends_s"));  // before the pulse ends
}
