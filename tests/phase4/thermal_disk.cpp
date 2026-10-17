#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "io/csv_writer.hpp"
#include "mechanics/elastic2d.hpp"
#include "mechanics/error_norms.hpp"
#include "mechanics/stress_recovery.hpp"
#include "mesh/gmsh_reader.hpp"
#include "phase4/phase4_common.hpp"
#include "support/convergence.hpp"
#include "support/validation_record.hpp"
#include "thermal/heat2d.hpp"

// P4-c, P4-d thin solid disk in plane stress, T = dT (1 - r^2/b^2), free rim (V4.3; Timoshenko and Goodier):
//   sigma_r / (E alpha dT) = -(1 - r^2/b^2)/4,  sigma_theta / (E alpha dT) = -(1 - 3 r^2/b^2)/4,
//   u_r = alpha dT r/4 [(3 + nu) - (1 + nu) r^2/b^2]   (docs/learning/phase4.md section 10).
// Reference: docs/REFERENCE_VALUES.md rows P4-c, P4-d: at r/b = 0, 0.5, 1 within 1% on the finest T6 mesh,
// absolute 0.0025 (in units of E alpha dT) where the target is 0 (sigma_r at r = b).
// P4-g: L2 displacement and energy-norm orders T3: 2, 1; T6: 3, 2, +/- 0.15 over levels 1-4.
// Quarter disk (cases/disk), u_x = 0 on x = 0, u_y = 0 on y = 0, nodal interpolant of T.
// Points on the ray theta = 30 deg; at r = b the rim node nearest to that ray (it lies on the circle).

using namespace tmfem;

namespace {

constexpr double kTolRel = 0.01;     // P4-c, P4-d
constexpr double kTolAbs = 0.0025;   // P4-c, P4-d where the target is 0
constexpr double kTolOrder = 0.15;   // P4-g

struct Disk {
  ElasticMaterial m = test::thermal_test_material();
  double b = test::kDisk_b_m, dT = test::kdT_K, alpha = test::kAlpha_per_K;
  double T(const Eigen::Vector2d& x) const { return test::kTref_K + dT * (1 - x.squaredNorm() / (b * b)); }
  Eigen::Vector2d u(const Eigen::Vector2d& x) const {
    const double rho2 = x.squaredNorm() / (b * b);
    return alpha * dT / 4 * ((3 + m.nu) - (1 + m.nu) * rho2) * x;  // u_r (x, y)/r
  }
  Eigen::Vector3d sigma(const Eigen::Vector2d& x) const {  // Cartesian, Pa
    const double s = m.E_Pa * alpha * dT, rho2 = x.squaredNorm() / (b * b);
    const double sr = -s * (1 - rho2) / 4, st = -s * (1 - 3 * rho2) / 4;
    const double r = x.norm();
    if (r == 0.0) return {sr, sr, 0.0};
    const double c = x(0) / r, sn = x(1) / r;
    return {sr * c * c + st * sn * sn, sr * sn * sn + st * c * c, (sr - st) * sn * c};
  }
};

Elastic2DProblem disk_problem(int t_nodes, int level, const Disk& d) {
  Mesh mesh = read_gmsh41(test::disk_mesh(t_nodes, level));
  const Eigen::VectorXd T = nodal_values(mesh, [&d](const Eigen::Vector2d& x) { return d.T(x); });
  const auto zero = [](const Eigen::Vector2d&) { return 0.0; };
  return {std::move(mesh), d.m, PlaneCondition::PlaneStress, "domain", {{"sym_x", 0, zero}, {"sym_y", 1, zero}}, {},
          ThermalStrain{d.alpha, test::kTref_K, T}};
}

// sigma_r, sigma_theta at x from a Cartesian stress (at r = 0 both are sigma_xx = sigma_yy).
Eigen::Vector2d polar(const Eigen::Vector3d& s, const Eigen::Vector2d& x) {
  const double r = x.norm();
  if (r == 0.0) return {s(0), s(1)};
  const double c = x(0) / r, sn = x(1) / r;
  return {s(0) * c * c + s(1) * sn * sn + 2 * s(2) * sn * c, s(0) * sn * sn + s(1) * c * c - 2 * s(2) * sn * c};
}

}  // namespace

TEST_CASE("P4-c P4-d thin disk stresses", "[phase4][P4-c][P4-d]") {
  const Disk d;
  const Elastic2DProblem p = disk_problem(6, test::kP4Levels - 1, d);
  const Elastic2DSolution sol = solve_elastic2d(p);
  const double scale = d.m.E_Pa * d.alpha * d.dT;
  const double th = 30.0 * M_PI / 180.0;

  int rim_node = -1;  // rim node nearest to the ray
  double best = 1e30;
  for (const int n : p.mesh.nodes_of("rim")) {
    const double dth = std::abs(std::atan2(p.mesh.x_m(n, 1), p.mesh.x_m(n, 0)) - th);
    if (dth < best) best = dth, rim_node = n;
  }
  const std::vector<std::pair<double, Eigen::Vector2d>> points{
      {0.0, Eigen::Vector2d(0.0, 0.0)},
      {0.5, 0.5 * d.b * Eigen::Vector2d(std::cos(th), std::sin(th))},
      {1.0, p.mesh.x_m.row(rim_node).transpose()}};
  const double sr_ref[] = {-0.25, -0.1875, 0.0}, st_ref[] = {-0.25, -0.0625, 0.5};
  {  // profile along the ray for the figure (scripts/plot_phase4.py)
    std::vector<double> rb, sr, st, sr_x, st_x;
    for (int i = 0; i <= 40; ++i) {
      const Eigen::Vector2d x = (0.999 * d.b * i / 40.0) * Eigen::Vector2d(std::cos(th), std::sin(th));
      const auto s = stress_at_point(p, sol.u_m, x);
      if (!s) continue;
      const Eigen::Vector2d sp = polar(*s, x) / scale, ex = polar(d.sigma(x), x) / scale;
      rb.push_back(x.norm() / d.b), sr.push_back(sp(0)), st.push_back(sp(1)), sr_x.push_back(ex(0)), st_x.push_back(ex(1));
    }
    write_csv(test::results_dir() / "phase4" / "disk_profile.csv",
              {"r_over_b", "sr_fe", "st_fe", "sr_exact", "st_exact"}, {rb, sr, st, sr_x, st_x},
              {"P4-c/d thin disk, T6 finest, ray theta = 30 deg, stresses / (E alpha dT)"});
  }
  for (std::size_t i = 0; i < points.size(); ++i) {
    const auto& [rb, x] = points[i];
    const auto s = stress_at_point(p, sol.u_m, x);
    REQUIRE(s.has_value());
    const Eigen::Vector2d sp = polar(*s, x) / scale;
    const Eigen::Vector2d ex = polar(d.sigma(x), x) / scale;  // exact at the actual point (r = b to 1e-16)
    std::printf("P4-c/d T6 finest, r/b = %.3f: sigma_r/(E a dT) = %.6f (exact %.6f, table %.4f), "
                "sigma_t/(E a dT) = %.6f (exact %.6f, table %.4f)\n",
                x.norm() / d.b, sp(0), ex(0), sr_ref[i], sp(1), ex(1), st_ref[i]);
    INFO("r/b = " << rb);
    CHECK(std::abs(ex(0) - sr_ref[i]) < 1e-12);  // the evaluation point is where the table says
    CHECK(std::abs(ex(1) - st_ref[i]) < 1e-12);
    const auto ok = [](double v, double ref) { return ref == 0.0 ? std::abs(v) <= kTolAbs : std::abs(v - ref) <= kTolRel * std::abs(ref); };
    CHECK(ok(sp(0), sr_ref[i]));
    CHECK(ok(sp(1), st_ref[i]));
    const std::string q = "r/b = " + std::to_string(rb).substr(0, 3);
    test::record("P4-c", "sigma_r over E alpha dT at " + q, sr_ref[i], sp(0),
                 sr_ref[i] == 0.0 ? "0.0025 absolute" : "1%", "P4-c P4-d thin disk stresses");
    test::record("P4-d", "sigma_theta over E alpha dT at " + q, st_ref[i], sp(1), "1%", "P4-c P4-d thin disk stresses");
  }
}

TEST_CASE("P4-g thin disk convergence", "[phase4][P4-g]") {
  const Disk d;
  for (const auto& [t, oL2, oE] : {std::tuple{3, 2.0, 1.0}, std::tuple{6, 3.0, 2.0}}) {
    std::vector<double> h, eL2, eE, lv;
    for (int level = 0; level < test::kP4Levels; ++level) {
      const Elastic2DProblem p = disk_problem(t, level, d);
      const Elastic2DSolution sol = solve_elastic2d(p);
      const ElasticErrors e = elastic_errors(p, sol.u_m, [&d](const Eigen::Vector2d& x) { return d.u(x); },
                                             [&d](const Eigen::Vector2d& x) { return d.sigma(x); });
      h.push_back(e.h_m);
      eL2.push_back(e.L2 / e.L2_exact);
      eE.push_back(e.energy / e.energy_exact);
      lv.push_back(level);
    }
    const std::string label = "T" + std::to_string(t) + " plane stress";
    test::print_orders(("P4-g " + label + ", relative L2 displacement error").c_str(), h, eL2);
    test::print_orders(("P4-g " + label + ", relative energy-norm error").c_str(), h, eE);
    const auto fit = [](const std::vector<double>& v) { return std::vector<double>(v.begin() + test::kP4FitFirstLevel, v.end()); };
    const double pL2 = test::observed_order(fit(h), fit(eL2)), pE = test::observed_order(fit(h), fit(eE));
    std::printf("  asserted orders (levels 1-4): L2 %.4f, energy %.4f\n", pL2, pE);
    INFO(label);
    CHECK(std::abs(pL2 - oL2) <= kTolOrder);
    CHECK(std::abs(pE - oE) <= kTolOrder);
    test::record("P4-g", label + " L2 order (levels 1-4)", oL2, pL2, "+/- 0.15", "P4-g thin disk convergence");
    test::record("P4-g", label + " energy order (levels 1-4)", oE, pE, "+/- 0.15", "P4-g thin disk convergence");
    write_csv(test::results_dir() / "phase4" / ("disk_convergence_t" + std::to_string(t) + ".csv"),
              {"level", "h_m", "L2_rel", "energy_rel"}, {lv, h, eL2, eE}, {"P4-g thin disk, plane stress"});
  }
}
