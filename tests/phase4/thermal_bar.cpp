#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>

#include "coupling/thermal_strain.hpp"
#include "mechanics/elastic2d.hpp"
#include "mechanics/stress_recovery.hpp"
#include "mesh/gmsh_reader.hpp"
#include "phase4/phase4_common.hpp"
#include "support/validation_record.hpp"

// P4-a bar held in x at both ends, free laterally, uniform dT (V4.1).
// Reference: docs/REFERENCE_VALUES.md row P4-a, tolerance 1e-8 relative (scaled by the nonzero target).
// Plane stress: sigma_xx = -E alpha dT, sigma_yy = 0, eps_yy = (1 + nu) alpha dT.
// Plane strain: sigma_xx = sigma_zz = -E alpha dT / (1 - nu), sigma_yy = 0, eps_yy = (1 + nu)/(1 - nu) alpha dT.
// Bar 100 mm x 20 mm (cases/bar), u_x = 0 on left and right, u_y = 0 at the corner (0, 0).

using namespace tmfem;

namespace {

constexpr double kTol = 1e-8;  // P4-a

void run(int t_nodes, PlaneCondition plane) {
  const ElasticMaterial m = test::thermal_test_material();
  const auto zero = [](const Eigen::Vector2d&) { return 0.0; };
  Mesh mesh = read_gmsh41(test::source_dir() / ("cases/bar/bar_t" + std::to_string(t_nodes) + ".msh"));
  const Eigen::VectorXd T = Eigen::VectorXd::Constant(mesh.n_nodes(), test::kTref_K + test::kdT_K);
  const Elastic2DProblem p{std::move(mesh), m, plane, "domain", {{"left", 0, zero}, {"right", 0, zero}, {"corner", 1, zero}},
                           {}, ThermalStrain{test::kAlpha_per_K, test::kTref_K, T}};
  const Elastic2DSolution sol = solve_elastic2d(p);

  const bool ps = plane == PlaneCondition::PlaneStress;
  const double adT = test::kAlpha_per_K * test::kdT_K;
  const double sxx_ref = ps ? -m.E_Pa * adT : -m.E_Pa * adT / (1 - m.nu);
  const double szz_ref = ps ? 0.0 : sxx_ref;
  const double eyy_ref = ps ? (1 + m.nu) * adT : (1 + m.nu) / (1 - m.nu) * adT;
  const double scale = std::abs(sxx_ref);

  double e_sxx = 0, e_syy = 0, e_sxy = 0, e_szz = 0, e_u = 0;
  for (const Eigen::Vector3d& s : gauss_point_stresses(p, sol.u_m)) {
    e_sxx = std::max(e_sxx, std::abs(s(0) - sxx_ref) / scale);
    e_syy = std::max(e_syy, std::abs(s(1)) / scale);
    e_sxy = std::max(e_sxy, std::abs(s(2)) / scale);
    e_szz = std::max(e_szz, std::abs(sigma_zz_thermal(m, plane, s, adT) - szz_ref) / scale);
  }
  const double u_scale = eyy_ref * p.mesh.x_m.col(1).maxCoeff();
  for (int i = 0; i < p.mesh.n_nodes(); ++i) {
    e_u = std::max(e_u, std::hypot(sol.u_m(2 * i), sol.u_m(2 * i + 1) - eyy_ref * p.mesh.x_m(i, 1)) / u_scale);
  }
  const std::string label = "T" + std::to_string(t_nodes) + " " + test::plane_name(plane);
  std::printf("P4-a %s: sxx ref %.6e Pa; max rel error sxx %.2e syy %.2e sxy %.2e szz %.2e u %.2e\n", label.c_str(),
              sxx_ref, e_sxx, e_syy, e_sxy, e_szz, e_u);
  INFO(label);
  CHECK(e_sxx < kTol);
  CHECK(e_syy < kTol);
  CHECK(e_sxy < kTol);
  CHECK(e_szz < kTol);
  CHECK(e_u < kTol);
  test::record("P4-a", label + " max rel error of sigma_xx", 0.0, e_sxx, "1e-8 relative", "P4-a held bar");
}

}  // namespace

TEST_CASE("P4-a held bar", "[phase4][P4-a]") {
  for (const int t : {3, 6}) {
    run(t, PlaneCondition::PlaneStress);
    run(t, PlaneCondition::PlaneStrain);
  }
}
