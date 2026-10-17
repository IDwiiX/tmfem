#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <memory>
#include <string>

#include "mechanics/j2_plasticity.hpp"
#include "mechanics/material_point.hpp"
#include "phase3/phase3_common.hpp"

// Supporting unit tests for the radial return and its consistent tangent (docs/learning/phase3.md sections
// 5-7). No reference row. The non-proportional state below has shear, so a missing factor 2 on the plastic
// shear strain or a wrong tangent term cannot hide behind uniaxial symmetry.

using namespace tmfem;

namespace {

const Voigt4 kEpsP_n(1e-3, -6e-4, -4e-4, 8e-4);  // deviatoric (trace 0), engineering shear
constexpr double kAlpha_n = 0.004;
const Voigt4 kEps(4e-3, 1e-3, -2e-3, 3e-3);

MaterialState state_n() { return {kEpsP_n, kAlpha_n, 0.0}; }

// Central finite difference of the return map's stress with respect to each strain component.
Eigen::Matrix4d fd_tangent(const ConstitutiveModel& m, const Voigt4& eps, const MaterialState& s) {
  Eigen::Matrix4d T;
  const double h = 1e-8;
  for (int j = 0; j < 4; ++j) {
    Voigt4 d = Voigt4::Zero();
    d(j) = h;
    T.col(j) = (m.update(eps + d, s).stress_Pa - m.update(eps - d, s).stress_Pa) / (2.0 * h);
  }
  return T;
}

void check_plastic_update(const J2Plasticity& m, const Voigt4& eps, const MaterialState& s, const std::string& label) {
  INFO(label);
  const MaterialUpdate u = m.update(eps, s);
  REQUIRE(u.plastic);
  const double sy = m.hardening().yield_stress_Pa(u.state.alpha);
  CHECK(std::abs(equivalent_stress(u.stress_Pa) - sy) <= 1e-12 * sy);  // back on the yield surface
  CHECK(std::abs(u.state.eps_p.head<3>().sum()) <= 1e-18);              // plastic incompressibility
  const Voigt4 eps_p = eps - elastic_tensor4(m.elastic()).inverse() * u.stress_Pa;
  CHECK((eps_p - u.state.eps_p).cwiseAbs().maxCoeff() <= 1e-15);        // engineering shear factor 2
  const Eigen::Matrix4d fd = fd_tangent(m, eps, s);
  const double scale = u.tangent_Pa.cwiseAbs().maxCoeff();
  CHECK((fd - u.tangent_Pa).cwiseAbs().maxCoeff() <= 1e-6 * scale);
  CHECK((u.tangent_Pa - u.tangent_Pa.transpose()).cwiseAbs().maxCoeff() <= 1e-12 * scale);
}

}  // namespace

TEST_CASE("Elastic tensor matches the phase 2 plane-strain D", "[phase3]") {
  const ElasticMaterial m = test::phase2_material();
  const Eigen::Matrix4d C = elastic_tensor4(m);
  const Eigen::Matrix3d D = elasticity_matrix(m, PlaneCondition::PlaneStrain);
  const int idx[3] = {0, 1, 3};
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) CHECK(std::abs(C(idx[i], idx[j]) - D(i, j)) <= 1e-14 * D(0, 0));
  }
}

TEST_CASE("Return mapping: elastic branch returns the elastic tensor", "[phase3]") {
  const J2Plasticity m(test::phase2_material(), std::make_shared<LinearHardening>(300e6, 1e9));
  const Voigt4 eps(1e-4, -5e-5, 2e-5, 1e-4);
  const MaterialUpdate u = m.update(eps, {});
  CHECK_FALSE(u.plastic);
  CHECK(u.tangent_Pa == elastic_tensor4(m.elastic()));
  CHECK((u.stress_Pa - elastic_tensor4(m.elastic()) * eps).norm() == 0.0);
}

TEST_CASE("Return mapping: linear hardening matches the closed-form multiplier", "[phase3]") {
  const double H = 2e9;
  const J2Plasticity m(test::phase2_material(), std::make_shared<LinearHardening>(300e6, H));
  const double G = shear_modulus(m.elastic());
  const Voigt4 trial = elastic_tensor4(m.elastic()) * (kEps - kEpsP_n);
  const double f_trial = equivalent_stress(trial) - (300e6 + H * kAlpha_n);
  const MaterialUpdate u = m.update(kEps, state_n());
  CHECK(std::abs((u.state.alpha - kAlpha_n) - f_trial / (3.0 * G + H)) <= 1e-12 * f_trial / (3.0 * G + H));
  check_plastic_update(m, kEps, state_n(), "linear hardening");
}

TEST_CASE("Return mapping: consistent tangent matches finite differences", "[phase3]") {
  const J2Plasticity ludwik(test::phase2_material(), test::ludwik_hardening(test::source_ludwik()));
  const J2Plasticity tabulated(test::phase2_material(), test::tabulated_hardening());
  check_plastic_update(ludwik, kEps, state_n(), "Ludwik, alpha_n = 0.004");
  check_plastic_update(tabulated, kEps, state_n(), "tabulated, alpha_n = 0.004");
  // First yield from alpha_n = 0, where the Ludwik slope is infinite.
  const Voigt4 eps0(-1e-3, 4e-3, -1e-3, 1e-3);
  check_plastic_update(ludwik, eps0, {}, "Ludwik, first yield from alpha_n = 0");
  check_plastic_update(tabulated, eps0, {}, "tabulated, first yield from alpha_n = 0");
}

TEST_CASE("Return mapping: the multiplier crosses a tabulated knot", "[phase3]") {
  const auto h = test::tabulated_hardening();
  const J2Plasticity m(test::phase2_material(), h);
  const MaterialState s{Voigt4::Zero(), 0.008, 0.0};  // just below the knot at 0.0086197
  const MaterialUpdate u = m.update(Voigt4(-4e-3, 1e-2, -4e-3, 0.0), s);
  REQUIRE(u.plastic);
  CHECK(u.state.alpha > h->alpha()[1]);
  CHECK(std::abs(equivalent_stress(u.stress_Pa) - h->yield_stress_Pa(u.state.alpha)) <= 1e-12 * 4e8);
}

TEST_CASE("Material point: uniaxial stress with linear hardening gives the bilinear curve", "[phase3]") {
  const double sy0 = 300e6, H = 2e9;
  const J2Plasticity m(test::phase2_material(), std::make_shared<LinearHardening>(sy0, H));
  const double E = m.elastic().E_Pa, nu = m.elastic().nu;
  const DriverResult r = MaterialPointDriver(m).run(uniaxial_stress_steps(400e6, 10));
  for (std::size_t k = 0; k < r.strain.size(); ++k) {
    const double s = r.stress_Pa[k](1);
    const double alpha = s > sy0 ? (s - sy0) / H : 0.0;
    CHECK(std::abs(r.strain[k](1) - (s / E + alpha)) <= 1e-12 * (s / E + alpha));
    CHECK(std::abs(r.strain[k](0) - (-nu * s / E - 0.5 * alpha)) <= 1e-12 * (s / E + alpha));
    CHECK(r.stress_Pa[k].cwiseAbs()(0) <= 1e-4);
  }
}
