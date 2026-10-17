#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstdio>
#include <set>
#include <string>

#include "mechanics/elastic2d.hpp"
#include "mechanics/kirsch_analytic.hpp"
#include "mechanics/stress_recovery.hpp"
#include "phase2/kirsch_common.hpp"
#include "support/validation_record.hpp"

// P2-c, P2-d, P2-e: Kirsch plate with a hole at S = 140 MPa (V2.3), plane stress, exact Kirsch tractions
// on the outer edges of the quarter model (docs/learning/phase2.md section 9). Point stresses come from
// the containing element(s) (section 8). "Finest" is level 4: levels 0-3 were fixed before any result
// existed, and level 4 was added later for P2-f (see kirsch_convergence.cpp). On level 3 the results
// were Kt = 2.99876 and max errors 0.036 MPa (P2-d) and 0.031 MPa (P2-e); P2-c prints every level.
// Reference: docs/REFERENCE_VALUES.md rows P2-c (Kt = 3 within 1% on the finest T6 mesh), P2-d (rows at
// r/a = 2 and 4 within 2 MPa) and P2-e (rows at r/a = 1.5 within 2 MPa).

using namespace tmfem;

namespace {

constexpr double kTolKt = 0.01;        // P2-c: 1% relative
constexpr double kTolStress_MPa = 2.0;  // P2-d, P2-e: 2 MPa absolute

double kt(const Elastic2DProblem& p, const Elastic2DSolution& sol) {
  // sigma_tt at r = a, theta = 90 deg (from the load axis) is sigma_yy at (a, 0).
  const auto s = stress_at_point(p, sol.u_m, {test::kKirschA_m, 0.0});
  REQUIRE(s.has_value());
  return (*s)(1) / test::kKirschS_Pa;
}

// Compares every CSV row with r/a in `radii` against the finest T6 solution; returns the max abs error (MPa).
double check_rows(const std::set<double>& radii) {
  const Elastic2DProblem p = test::kirsch_problem(6, test::kKirschFinest, PlaneCondition::PlaneStress);
  const Elastic2DSolution sol = solve_elastic2d(p);
  double worst = 0.0;
  std::printf("  r/a   theta | srr FE    CSV      | stt FE    CSV      | srt FE    CSV      (MPa)\n");
  for (const auto& r : test::read_kirsch_csv()) {
    if (radii.count(r.r_over_a) == 0) continue;
    const test::QuarterPoint q = test::map_to_quarter(r.r_over_a * test::kKirschA_m, r.theta_deg);
    const auto s = stress_at_point(p, sol.u_m, q.x_m);
    REQUIRE(s.has_value());
    const Eigen::Vector3d pol = cartesian_to_polar_from_y(*s, q.theta_rad) / test::kPaPerMPa;
    const double srt = q.rt_sign * pol(2);
    std::printf("  %-5.2g %5.0f | %8.3f %8.3f | %8.3f %8.3f | %8.3f %8.3f\n", r.r_over_a, r.theta_deg, pol(0),
                r.srr_MPa, pol(1), r.stt_MPa, srt, r.srt_MPa);
    INFO("r/a = " << r.r_over_a << ", theta = " << r.theta_deg);
    CHECK(std::abs(pol(0) - r.srr_MPa) <= kTolStress_MPa);
    CHECK(std::abs(pol(1) - r.stt_MPa) <= kTolStress_MPa);
    CHECK(std::abs(srt - r.srt_MPa) <= kTolStress_MPa);
    worst = std::max({worst, std::abs(pol(0) - r.srr_MPa), std::abs(pol(1) - r.stt_MPa), std::abs(srt - r.srt_MPa)});
  }
  std::printf("  max abs error %.4f MPa\n", worst);
  return worst;
}

}  // namespace

TEST_CASE("P2-c Kirsch Kt on finest T6 mesh", "[phase2][P2-c]") {
  double kt_finest_t6 = 0.0;
  std::printf("P2-c Kt = sigma_yy(a, 0) / S, plane stress (only the finest T6 value is asserted)\n");
  for (const int t : {3, 6}) {
    for (int level = 0; level < test::kKirschLevels; ++level) {
      const Elastic2DProblem p = test::kirsch_problem(t, level, PlaneCondition::PlaneStress);
      const double k = kt(p, solve_elastic2d(p));
      std::printf("  T%d level %d (%d nodes): Kt = %.5f (error %+.3f%%)\n", t, level, p.mesh.n_nodes(), k,
                  100 * (k - 3) / 3);
      if (t == 6 && level == test::kKirschFinest) kt_finest_t6 = k;
      if (t == 3 && level == test::kKirschFinest) {
        test::record("P2-c", "Kt finest T3 (report only)", 3.0, k, "not asserted", "P2-c Kirsch Kt on finest T6 mesh");
      }
    }
  }
  CHECK(std::abs(kt_finest_t6 - 3.0) / 3.0 <= kTolKt);
  test::record("P2-c", "Kt finest T6", 3.0, kt_finest_t6, "1%", "P2-c Kirsch Kt on finest T6 mesh");
}

TEST_CASE("P2-d Kirsch stresses at r/a = 2 and 4", "[phase2][P2-d]") {
  std::printf("P2-d finest T6, plane stress, CSV rows at r/a = 2 and 4\n");
  const double worst = check_rows({2.0, 4.0});
  test::record("P2-d", "max abs stress error, 20 rows x 3 components [MPa]", 0.0, worst, "2 MPa",
               "P2-d Kirsch stresses at r/a = 2 and 4");
}

TEST_CASE("P2-e Kirsch stresses at r/a = 1.5", "[phase2][P2-e]") {
  // The source table at r/a = 1.5 (e.g. 2100, -3080, 0 at theta = 0) is wrong; the CSV holds the
  // recomputed Kirsch values (25.93, -10.37, 0), which are the target.
  std::printf("P2-e finest T6, plane stress, CSV rows at r/a = 1.5\n");
  const double worst = check_rows({1.5});
  test::record("P2-e", "max abs stress error, 10 rows x 3 components [MPa]", 0.0, worst, "2 MPa",
               "P2-e Kirsch stresses at r/a = 1.5");
}
