#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstdio>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "io/csv_writer.hpp"
#include "mechanics/elastoplastic2d.hpp"
#include "mesh/gmsh_reader.hpp"
#include "phase3/phase3_common.hpp"
#include "support/validation_record.hpp"

// P3-d (V3.2): global Newton converges quadratically, which shows the tangent is consistent.
// Reference: docs/REFERENCE_VALUES.md row P3-d, observed rate >= 1.8 (Ludwik hardening, kickoff decision Q11:
// smooth for alpha > 0; the tabulated law's slope jumps void the theorem, so its rates are printed only).
// Rate p = ln(r+/r) / ln(r/r-) over consecutive Newton iterates in which the set of plastic points is the same at
// all three iterates and r+ >= 1e-11 F_ref (above round-off). At least one triple must qualify, and the last
// qualifying triple must give p >= 1.8.
//
// Set-up change, approved by the user after the first run (docs/learning/phase3.md section 15): the plan
// asserted this on the Kirsch plate (T6 level 2, plane strain, S = 100 -> 140 MPa in one step). There Newton
// converged so fast once the plastic set settled (93 -> 99 -> 100 points) that only two residuals stayed above
// round-off (1.9e-6, 1.4e-10, then the 7e-14 floor), so no triple qualified. The assertion now uses the
// FE plane-strain bar (cases/bar/bar_t6.msh, u_x = 0.01 L, then one step to 0.05 L), where every point is plastic
// from the first correction, so the active set is fixed by construction. The criterion is unchanged. The
// Kirsch sequences are still printed and written as a diagnostic.

using namespace tmfem;

namespace {

constexpr double kMinRate = 1.8;  // P3-d
constexpr double kFloor = 1e-11;  // r+ >= kFloor * F_ref
constexpr double kL_m = 0.100;
constexpr double kEpsEnd = 0.05;

struct RateSummary {
  int qualifying = 0;
  double last_rate = std::numeric_limits<double>::quiet_NaN();
};

NewtonOptions to_stagnation() {
  NewtonOptions opt;
  opt.run_to_stagnation = true;
  opt.max_iterations = 30;
  return opt;
}

// Bar: u_x = 0 on the left, u_y = 0 at the corner, u_x = lambda 0.05 L on the right; steps lambda = 0.2, 1.
StepReport bar_step(std::shared_ptr<const J2Plasticity> model, PlaneCondition plane = PlaneCondition::PlaneStrain) {
  const auto zero = [](const Eigen::Vector2d&) { return 0.0; };
  const auto pull = [](const Eigen::Vector2d&) { return kEpsEnd * kL_m; };
  ElastoplasticSolver solver({read_gmsh41(test::source_dir() / "cases/bar/bar_t6.msh"),
                              std::move(model),
                              plane,
                              "domain",
                              {{"left", 0, zero}, {"corner", 1, zero}, {"right", 0, pull}},
                              {}},
                             to_stagnation());
  solver.step_to(0.2);
  return solver.step_to(1.0);
}

// Kirsch plate (diagnostic): T6 level 2, plane strain, S = 100 MPa (elastic), then 100 -> 140 MPa in one step.
StepReport kirsch_step(std::shared_ptr<const J2Plasticity> model) {
  ElastoplasticSolver solver(test::kirsch_plastic_problem(6, 2, PlaneCondition::PlaneStrain, std::move(model)),
                             to_stagnation());
  solver.step_to(100.0 / 140.0);
  for (const QuadPoint& q : solver.quad_points()) REQUIRE_FALSE(q.plastic);  // S = 100 MPa is below first yield
  return solver.step_to(1.0);
}

RateSummary analyse(const StepReport& rep, const std::string& label, const std::string& description) {
  RateSummary out;
  const auto& it = rep.iterates;
  // One normalisation for the whole sequence (F_ref of the last iterate), so the printed ratios match the rates;
  // the qualification test itself uses each iterate's own F_ref, as fixed in the plan.
  const double F = it.back().reference_force_N;
  std::vector<double> idx, rel, n_plastic, rate;
  std::printf("P3-d %s: %s\n  iter  r / F_ref,final  plastic points  rate (qualifying triple?)\n", label.c_str(),
              description.c_str());
  for (std::size_t i = 0; i < it.size(); ++i) {
    int n = 0;
    for (const char c : it[i].plastic) n += c;
    double p = std::numeric_limits<double>::quiet_NaN();
    bool ok = false;
    if (i >= 2) {
      const NewtonIterate &a = it[i - 2], &b = it[i - 1], &c = it[i];
      p = std::log(c.residual_N / b.residual_N) / std::log(b.residual_N / a.residual_N);
      ok = a.plastic == b.plastic && b.plastic == c.plastic && c.residual_N >= kFloor * c.reference_force_N;
      if (ok) {
        ++out.qualifying;
        out.last_rate = p;
      }
    }
    char rate_text[32] = "";
    if (i >= 2) std::snprintf(rate_text, sizeof rate_text, "%.4f (%s)", p, ok ? "yes" : "no");
    std::printf("  %4zu  %.3e      %8d        %s\n", i, it[i].residual_N / F, n, rate_text);
    idx.push_back(static_cast<double>(i));
    rel.push_back(it[i].residual_N / F);
    n_plastic.push_back(n);
    rate.push_back(p);
  }
  write_csv(test::results_dir() / ("phase3/newton_history_" + label + ".csv"),
            {"iteration", "residual_over_Fref_final", "n_plastic_points", "rate"}, {idx, rel, n_plastic, rate},
            {"P3-d " + description});
  std::printf("  qualifying triples: %d, last qualifying rate: %.4f\n", out.qualifying, out.last_rate);
  return out;
}

}  // namespace

TEST_CASE("P3-d Newton convergence rate on the plane-strain bar", "[phase3][P3-d]") {
  const auto ludwik = [] { return test::j2(test::ludwik_hardening(test::source_ludwik())); };
  const RateSummary bar = analyse(bar_step(ludwik()), "bar_ludwik",
                                  "FE plane-strain bar T6, source Ludwik, step eps_xx = 0.01 -> 0.05");
  REQUIRE(bar.qualifying >= 1);
  CHECK(bar.last_rate >= kMinRate);
  test::record("P3-d", "Newton rate, last qualifying triple, plane-strain bar, Ludwik", 2.0, bar.last_rate, ">= 1.8",
               "P3-d Newton convergence rate on the plane-strain bar");

  // Printed only: the tabulated bar, and the Kirsch plate for both laws (the original, pre-fixed set-up).
  analyse(bar_step(test::j2(test::tabulated_hardening())), "bar_tabulated",
          "FE plane-strain bar T6, tabulated, step eps_xx = 0.01 -> 0.05 (printed only)");
  analyse(kirsch_step(ludwik()), "kirsch_ludwik",
          "Kirsch T6 level 2 plane strain, source Ludwik, S = 100 -> 140 MPa (diagnostic, printed only)");
  analyse(kirsch_step(test::j2(test::tabulated_hardening())), "kirsch_tabulated",
          "Kirsch T6 level 2 plane strain, tabulated, S = 100 -> 140 MPa (diagnostic, printed only)");
}

TEST_CASE("P3-d Newton convergence rate on the plane-stress bar", "[phase3][P3-d]") {
  // Stage B: the same criterion on the same bar in plane stress, where every point runs the local sigma_zz = 0
  // iteration and the statically condensed tangent (docs/learning/phase3.md section 8).
  const RateSummary bar = analyse(bar_step(test::j2(test::ludwik_hardening(test::source_ludwik())), PlaneCondition::PlaneStress),
                                  "bar_plane_stress_ludwik",
                                  "FE plane-stress bar T6, source Ludwik, step eps_xx = 0.01 -> 0.05");
  REQUIRE(bar.qualifying >= 1);
  CHECK(bar.last_rate >= kMinRate);
  test::record("P3-d", "Newton rate, last qualifying triple, plane-stress bar, Ludwik", 2.0, bar.last_rate, ">= 1.8",
               "P3-d Newton convergence rate on the plane-stress bar");
}
