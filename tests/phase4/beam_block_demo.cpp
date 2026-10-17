#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstdio>

#include "io/case_thermomech2d.hpp"
#include "phase4/phase4_common.hpp"

// Demo D4 smoke test (not a validation, so no REFERENCE_VALUES row): the copper and graphite cases load,
// run for one beam size, and conserve the deposited energy through the adiabatic pulse (stored energy
// 1^T M (T - T0) equals the discrete deposit tau 1^T F_q to round-off; docs/learning/phase4.md section 4).
// The physics numbers of the demo are reported by the demo itself, not asserted here.

using namespace tmfem;

TEST_CASE("D4 beam-block demo runs and conserves the pulse energy", "[phase4][demo]") {
  for (const char* name : {"copper", "graphite"}) {
    BeamBlockCase c = load_thermomech2d_case(test::source_dir() / "cases/beam_block" / (std::string(name) + ".json"),
                                             test::source_dir());
    c.output_dir = test::results_dir() / "phase4" / "demo_smoke" / name;
    const BeamRunResult r = run_beam_size(c, 1.0);
    const double rel = std::abs(r.energy_stored_J_m - r.energy_deposited_J_m) / r.energy_deposited_J_m;
    std::printf("D4 smoke %s s = 1 mm: stored %.9f J/m, deposited %.9f J/m (rel %.2e), %zu snapshots\n", name,
                r.energy_stored_J_m, r.energy_deposited_J_m, rel, r.history.size());
    INFO(name);
    CHECK(rel < 1e-10);
    CHECK(r.history.size() == c.cooling_stage_ends_s.size() + 1);
    for (const BeamSnapshot& s : r.history) CHECK(std::isfinite(s.vm_max_Pa));
  }
}
