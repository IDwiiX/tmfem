#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>

#include "coupling/thermal_strain.hpp"
#include "mechanics/elastic2d.hpp"
#include "mechanics/stress_recovery.hpp"
#include "mesh/gmsh_reader.hpp"
#include "phase2/kirsch_common.hpp"
#include "phase4/phase4_common.hpp"
#include "support/validation_record.hpp"

// P4-b unconstrained body, uniform dT (V4.2): zero stress and u = alpha_eff dT x.
// Reference: docs/REFERENCE_VALUES.md row P4-b, tolerance 1e-10 relative. alpha_eff = alpha (plane stress) or
// (1 + nu) alpha (plane strain, where sigma_zz = -E alpha dT). The quarter plate with a hole of
// cases/kirsch (unstructured, curved T6 hole edge) on symmetry rollers: u_x = 0 on x = 0, u_y = 0 on y = 0,
// which does not restrain the expansion about the origin. Section 9 of docs/learning/phase4.md explains
// why this holds to round-off on curved T6.

using namespace tmfem;

namespace {

constexpr double kTol = 1e-10;  // P4-b

void run(int t_nodes, PlaneCondition plane) {
  const ElasticMaterial m = test::thermal_test_material();
  const auto zero = [](const Eigen::Vector2d&) { return 0.0; };
  Mesh mesh = read_gmsh41(test::kirsch_mesh(t_nodes, 1));
  const Eigen::VectorXd T = Eigen::VectorXd::Constant(mesh.n_nodes(), test::kTref_K + test::kdT_K);
  const Elastic2DProblem p{std::move(mesh), m, plane, "domain", {{"sym_x", 0, zero}, {"sym_y", 1, zero}}, {},
                           ThermalStrain{test::kAlpha_per_K, test::kTref_K, T}};
  const Elastic2DSolution sol = solve_elastic2d(p);

  const bool ps = plane == PlaneCondition::PlaneStress;
  const double adT = test::kAlpha_per_K * test::kdT_K;
  const double a_eff_dT = ps ? adT : (1 + m.nu) * adT;
  const double scale = p.beta() * test::kdT_K;  // stress a forgotten thermal term would leave
  const double szz_ref = ps ? 0.0 : -m.E_Pa * adT;

  double e_s = 0, e_szz = 0, e_u = 0, u_max = 0;
  for (const Eigen::Vector3d& s : gauss_point_stresses(p, sol.u_m)) {
    e_s = std::max(e_s, s.cwiseAbs().maxCoeff() / scale);
    e_szz = std::max(e_szz, std::abs(sigma_zz_thermal(m, plane, s, adT) - szz_ref) / (m.E_Pa * adT));
  }
  for (int i = 0; i < p.mesh.n_nodes(); ++i) {
    const Eigen::Vector2d x = p.mesh.x_m.row(i).transpose();
    e_u = std::max(e_u, (sol.u_m.segment<2>(2 * i) - a_eff_dT * x).norm());
    u_max = std::max(u_max, a_eff_dT * x.norm());
  }
  e_u /= u_max;
  const std::string label = "T" + std::to_string(t_nodes) + " " + test::plane_name(plane);
  std::printf("P4-b %s: max |sigma|/(beta dT) %.2e, sigma_zz rel error %.2e, max |u - alpha_eff dT x|/max|u| %.2e\n",
              label.c_str(), e_s, e_szz, e_u);
  INFO(label);
  CHECK(e_s < kTol);
  CHECK(e_szz < kTol);
  CHECK(e_u < kTol);
  test::record("P4-b", label + " max stress over beta dT", 0.0, e_s, "1e-10 relative", "P4-b free body, uniform dT");
  test::record("P4-b", label + " max rel displacement error", 0.0, e_u, "1e-10 relative", "P4-b free body, uniform dT");
}

}  // namespace

TEST_CASE("P4-b free body, uniform dT", "[phase4][P4-b]") {
  for (const int t : {3, 6}) {
    run(t, PlaneCondition::PlaneStress);
    run(t, PlaneCondition::PlaneStrain);
  }
}
