#include <catch2/catch_test_macros.hpp>

#include <cmath>

#include "mechanics/elasticity.hpp"
#include "mechanics/kirsch_analytic.hpp"
#include "phase2/kirsch_common.hpp"

// The exact Kirsch field used by P2-c..f (docs/learning/phase2.md section 9). Not a reference row
// itself, but a typo here would fake every Kirsch result, so it is checked against the CSV and
// against D eps(u) with u the exact displacement.

using namespace tmfem;

TEST_CASE("Kirsch polar stress matches kirsch_sigma140MPa.csv", "[phase2][kirsch]") {
  const KirschPlate k = test::kirsch_plate();
  const auto rows = test::read_kirsch_csv();
  REQUIRE(rows.size() == 40);
  for (const auto& r : rows) {
    INFO("r/a = " << r.r_over_a << ", theta = " << r.theta_deg);
    const Eigen::Vector3d s = kirsch_stress_polar(k, r.r_over_a * k.a_m, r.theta_deg * test::kPi / 180.0) / test::kPaPerMPa;
    CHECK(std::abs(s(0) - r.srr_MPa) < 1e-4);  // CSV has 4 decimals
    CHECK(std::abs(s(1) - r.stt_MPa) < 1e-4);
    CHECK(std::abs(s(2) - r.srt_MPa) < 1e-4);
  }
}

TEST_CASE("Kirsch Cartesian stress, reflected into the quarter, matches the polar formulas", "[phase2][kirsch]") {
  const KirschPlate k = test::kirsch_plate();
  for (const auto& r : test::read_kirsch_csv()) {
    INFO("r/a = " << r.r_over_a << ", theta = " << r.theta_deg);
    const test::QuarterPoint q = test::map_to_quarter(r.r_over_a * k.a_m, r.theta_deg);
    const Eigen::Vector3d p = cartesian_to_polar_from_y(kirsch_stress_cartesian(k, q.x_m), q.theta_rad) / test::kPaPerMPa;
    CHECK(std::abs(p(0) - r.srr_MPa) < 1e-4);
    CHECK(std::abs(p(1) - r.stt_MPa) < 1e-4);
    CHECK(std::abs(q.rt_sign * p(2) - r.srt_MPa) < 1e-4);
  }
}

TEST_CASE("Kirsch exact displacement: D eps(u) equals the Kirsch stress", "[phase2][kirsch]") {
  const KirschPlate k = test::kirsch_plate();
  const ElasticMaterial m = test::phase2_material();
  const double a = k.a_m, h = 1e-7 * a;
  for (const auto plane : {PlaneCondition::PlaneStress, PlaneCondition::PlaneStrain}) {
    const Eigen::Matrix3d D = elasticity_matrix(m, plane);
    for (const Eigen::Vector2d& x : {Eigen::Vector2d(1.2 * a, 0.3 * a), Eigen::Vector2d(0.5 * a, 2.0 * a),
                                    Eigen::Vector2d(3.1 * a, 2.7 * a), Eigen::Vector2d(4.9 * a, 0.2 * a)}) {
      const auto u = [&](double dx, double dy) { return kirsch_displacement(k, m, plane, x + Eigen::Vector2d(dx, dy)); };
      const Eigen::Vector2d du_dx = (u(h, 0) - u(-h, 0)) / (2 * h);
      const Eigen::Vector2d du_dy = (u(0, h) - u(0, -h)) / (2 * h);
      const Eigen::Vector3d eps(du_dx(0), du_dy(1), du_dy(0) + du_dx(1));
      const Eigen::Vector3d s_exact = kirsch_stress_cartesian(k, x);
      INFO("x/a = " << x.transpose() / a);
      CHECK((D * eps - s_exact).norm() < 1e-6 * k.S_Pa);
    }
  }
}

TEST_CASE("Kirsch exact field: free hole and symmetry conditions", "[phase2][kirsch]") {
  const KirschPlate k = test::kirsch_plate();
  const ElasticMaterial m = test::phase2_material();
  for (int i = 0; i <= 18; ++i) {
    const double t = i * 5.0 * test::kPi / 180.0;
    const Eigen::Vector2d x(k.a_m * std::sin(t), k.a_m * std::cos(t));
    const Eigen::Vector3d p = cartesian_to_polar_from_y(kirsch_stress_cartesian(k, x), t);
    CHECK(std::abs(p(0)) < 1e-9 * k.S_Pa);  // sigma_rr = 0 on r = a
    CHECK(std::abs(p(2)) < 1e-9 * k.S_Pa);  // sigma_rt = 0 on r = a
  }
  for (const auto plane : {PlaneCondition::PlaneStress, PlaneCondition::PlaneStrain}) {
    for (const double s : {1.0, 2.5, 5.0}) {
      CHECK(std::abs(kirsch_displacement(k, m, plane, {0.0, s * k.a_m})(0)) < 1e-18);  // u_x = 0 on x = 0
      CHECK(std::abs(kirsch_displacement(k, m, plane, {s * k.a_m, 0.0})(1)) < 1e-18);  // u_y = 0 on y = 0
    }
  }
}
