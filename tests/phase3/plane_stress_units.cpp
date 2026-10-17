#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

#include "mechanics/elastoplastic2d.hpp"
#include "mechanics/material_point.hpp"
#include "mechanics/plane_reduction.hpp"
#include "mesh/gmsh_reader.hpp"
#include "phase3/phase3_common.hpp"

// Stage B supporting tests: the plane-stress wrapper (local sigma_zz = 0 Newton and the condensed tangent,
// docs/learning/phase3.md section 8) is checked against the phase 2 elastic D, finite differences and the 3D
// uniaxial material-point driver before the hole plate uses it. No reference row.

using namespace tmfem;

namespace {

constexpr Control S = Control::Stress;
constexpr Control E = Control::Strain;

}  // namespace

TEST_CASE("Plane stress wrapper: elastic limit equals the phase 2 plane-stress D", "[phase3]") {
  const LinearElasticModel m(test::phase2_material());
  const PlaneResponse r = plane_update(m, PlaneCondition::PlaneStress, Eigen::Vector3d(1e-4, -3e-5, 2e-5), {});
  const Eigen::Matrix3d D = elasticity_matrix(m.elastic(), PlaneCondition::PlaneStress);
  CHECK((r.tangent_Pa - D).cwiseAbs().maxCoeff() <= 1e-12 * D(0, 0));
  CHECK((r.stress_Pa - D * Eigen::Vector3d(1e-4, -3e-5, 2e-5)).cwiseAbs().maxCoeff() <= 1e-12 * D(0, 0) * 1e-4);
  CHECK(std::abs(r.sigma_zz_Pa) <= 1e-12 * m.stress_scale_Pa());
}

TEST_CASE("Plane stress wrapper: condensed tangent matches finite differences", "[phase3]") {
  const MaterialState s{Voigt4(1e-3, -6e-4, -4e-4, 8e-4), 0.004, -4e-4};
  const Eigen::Vector3d eps(4e-3, 1e-3, 3e-3);
  for (const auto& [label, h] : {std::make_pair(std::string("Ludwik"), std::shared_ptr<const HardeningLaw>(
                                                                           test::ludwik_hardening(test::source_ludwik()))),
                                 std::make_pair(std::string("tabulated"),
                                                std::shared_ptr<const HardeningLaw>(test::tabulated_hardening()))}) {
    INFO(label);
    const auto m = test::j2(h);
    const PlaneResponse r = plane_update(*m, PlaneCondition::PlaneStress, eps, s);
    REQUIRE(r.plastic);
    CHECK(std::abs(r.sigma_zz_Pa) <= 1e-12 * m->stress_scale_Pa());
    Eigen::Matrix3d fd;
    for (int j = 0; j < 3; ++j) {
      Eigen::Vector3d d = Eigen::Vector3d::Zero();
      d(j) = 1e-8;
      fd.col(j) = (plane_update(*m, PlaneCondition::PlaneStress, eps + d, s).stress_Pa -
                   plane_update(*m, PlaneCondition::PlaneStress, eps - d, s).stress_Pa) /
                  2e-8;
    }
    CHECK((fd - r.tangent_Pa).cwiseAbs().maxCoeff() <= 1e-6 * r.tangent_Pa.cwiseAbs().maxCoeff());
  }
}

TEST_CASE("Plane stress wrapper reproduces the 3D uniaxial driver", "[phase3]") {
  // The 3D driver gives (eps_xx, eps_zz, sigma_yy) for uniaxial stress along y. Feeding its eps_xx and eps_yy to
  // the wrapper must return the same sigma_yy, sigma_xx = 0 and the driver's eps_zz in the state.
  const auto m = test::j2(test::tabulated_hardening());
  const std::vector<double> eps_yy{0.001, 0.002, 0.005, 0.01, 0.02, 0.03158, 0.05};
  const DriverResult drv = MaterialPointDriver(*m).run(uniaxial_strain_steps(eps_yy));
  MaterialState state;
  for (std::size_t k = 0; k < eps_yy.size(); ++k) {
    const Voigt4& e = drv.strain[k];
    const PlaneResponse r = plane_update(*m, PlaneCondition::PlaneStress, Eigen::Vector3d(e(0), e(1), 0.0), state);
    INFO("eps_yy = " << eps_yy[k]);
    CHECK(std::abs(r.stress_Pa(1) - drv.stress_Pa[k](1)) <= 1e-10 * drv.stress_Pa[k](1));
    // The input strains carry the driver's convergence tolerance (lateral stresses to 1e-4 Pa) and the wrapper's
    // (sigma_zz to 1e-12 sigma_y0, 3e-4 Pa), so sigma_xx is zero only to that level: about 1e-12 sigma_yy.
    CHECK(std::abs(r.stress_Pa(0)) <= 1e-12 * drv.stress_Pa[k](1));
    CHECK(std::abs(r.state.eps_zz - e(2)) <= 1e-10 * std::abs(e(2)));
    CHECK(std::abs(r.state.alpha - drv.state[k].alpha) <= 1e-10 * std::max(drv.state[k].alpha, 1e-12));
    state = r.state;
  }
}

TEST_CASE("Plane-stress bar FE equals the uniaxial material-point driver", "[phase3]") {
  // cases/bar T3 and T6, u_x = lambda 0.05 L in 20 steps, plane stress: a uniform uniaxial stress state, so every
  // point must equal the 3D driver with eps_xx prescribed and sigma_yy = sigma_zz = sigma_xy = 0.
  const auto m = test::j2(test::tabulated_hardening());
  std::vector<DriverStep> steps;
  for (int k = 1; k <= 20; ++k) steps.push_back({{E, S, S, S}, Voigt4(0.05 * k / 20.0, 0.0, 0.0, 0.0)});
  const DriverResult drv = MaterialPointDriver(*m).run(steps);
  const auto zero = [](const Eigen::Vector2d&) { return 0.0; };
  const auto pull = [](const Eigen::Vector2d&) { return 0.05 * 0.100; };
  for (const int t : {3, 6}) {
    ElastoplasticSolver solver({read_gmsh41(test::source_dir() / "cases/bar" / ("bar_t" + std::to_string(t) + ".msh")),
                                m,
                                PlaneCondition::PlaneStress,
                                "domain",
                                {{"left", 0, zero}, {"corner", 1, zero}, {"right", 0, pull}},
                                {}});
    double worst = 0.0;
    for (int k = 0; k < 20; ++k) {
      solver.step_to((k + 1) / 20.0);
      const Voigt4& s_ref = drv.stress_Pa[static_cast<std::size_t>(k)];
      const Voigt4& e_ref = drv.strain[static_cast<std::size_t>(k)];
      const double a_ref = drv.state[static_cast<std::size_t>(k)].alpha;
      for (std::size_t q = 0; q < solver.quad_points().size(); ++q) {
        const QuadPoint& pt = solver.quad_points()[q];
        const double errs[] = {std::abs(pt.stress_Pa(0) - s_ref(0)) / s_ref(0), std::abs(pt.stress_Pa(1)) / s_ref(0),
                               std::abs(pt.sigma_zz_Pa) / s_ref(0),
                               std::abs(solver.quad_strain(q)(1) - e_ref(1)) / std::abs(e_ref(1)),
                               std::abs(pt.state.eps_zz - e_ref(2)) / std::abs(e_ref(2)),
                               a_ref > 0.0 ? std::abs(pt.state.alpha - a_ref) / a_ref : pt.state.alpha};
        for (const double err : errs) {
          worst = std::max(worst, err);
          INFO("T" << t << " step " << k + 1 << " point " << q);
          CHECK(err <= 1e-6);
        }
      }
    }
    std::printf("plane-stress bar T%d vs uniaxial driver: max relative difference %.3e\n", t, worst);
  }
}
