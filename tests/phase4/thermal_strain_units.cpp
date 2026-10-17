#include <catch2/catch_test_macros.hpp>

#include <cmath>

#include "coupling/thermal_strain.hpp"
#include "fem/element_geometry.hpp"
#include "fem/quadrature.hpp"
#include "fem/shape_functions.hpp"
#include "mechanics/elasticity.hpp"
#include "phase4/phase4_common.hpp"

// Unit checks of the thermal load (docs/learning/phase4.md sections 7 to 9).

using namespace tmfem;

namespace {

// f = int B^T m beta dT_h dA with a high-order collapsed Gauss rule, as an independent reference.
Eigen::VectorXd reference_thermal_load(CellType type, const Eigen::MatrixXd& x_e, const Eigen::VectorXd& dT_e,
                                       double beta) {
  const QuadRuleTri rule = triangle_collapsed_gauss(5);
  Eigen::VectorXd f = Eigen::VectorXd::Zero(2 * x_e.rows());
  for (std::size_t q = 0; q < rule.w.size(); ++q) {
    const ShapeEval s = shape(type, rule.xi[q]);
    const CellGeometry g = cell_geometry(s, x_e);
    const double w = beta * s.N.dot(dT_e) * g.detJ * rule.w[q];
    for (Eigen::Index a = 0; a < x_e.rows(); ++a) f.segment<2>(2 * a) += w * g.dN_dx.row(a).transpose();
  }
  return f;
}

Eigen::MatrixXd straight_t6() {
  Eigen::MatrixXd x(6, 2);
  x << 0.0, 0.0, 2.0, 0.3, 0.4, 1.5, 1.0, 0.15, 1.2, 0.9, 0.2, 0.75;
  return x;
}

Eigen::MatrixXd curved_t6() {  // edge 1-2 bowed outward: its mid-side node is off the chord
  Eigen::MatrixXd x = straight_t6();
  x.row(4) << 1.4, 1.1;
  return x;
}

}  // namespace

TEST_CASE("Thermal stress coefficient beta", "[phase4][units]") {
  const ElasticMaterial m = test::thermal_test_material();
  const double Ea = m.E_Pa * test::kAlpha_per_K;
  CHECK(std::abs(thermal_stress_coefficient(m, test::kAlpha_per_K, PlaneCondition::PlaneStress) - Ea / (1 - m.nu)) <=
        1e-14 * Ea);
  CHECK(std::abs(thermal_stress_coefficient(m, test::kAlpha_per_K, PlaneCondition::PlaneStrain) -
                 Ea / (1 - 2 * m.nu)) <= 1e-14 * Ea);
  // beta m = D m alpha_eff (section 8): D_ps m alpha = beta_ps m and D_pe m (1 + nu) alpha = beta_pe m.
  for (const auto plane : {PlaneCondition::PlaneStress, PlaneCondition::PlaneStrain}) {
    const double a_eff = plane == PlaneCondition::PlaneStress ? test::kAlpha_per_K : (1 + m.nu) * test::kAlpha_per_K;
    const Eigen::Vector3d Dm = elasticity_matrix(m, plane) * Eigen::Vector3d(1, 1, 0) * a_eff;
    const double beta = thermal_stress_coefficient(m, test::kAlpha_per_K, plane);
    CHECK((Dm - beta * Eigen::Vector3d(1, 1, 0)).norm() <= 1e-14 * beta);
  }
}

TEST_CASE("Thermal load quadrature is exact (T3 linear dT, straight T6 quadratic dT)", "[phase4][units]") {
  const double beta = 1.0e6;
  Eigen::MatrixXd x3 = straight_t6().topRows(3);
  const Eigen::Vector3d dT3(10.0, -4.0, 7.0);
  const Eigen::VectorXd f3 = element_thermal_load(CellType::Tri3, x3, dT3, beta);
  CHECK((f3 - reference_thermal_load(CellType::Tri3, x3, dT3, beta)).norm() <= 1e-12 * f3.norm());

  Eigen::VectorXd dT6(6);
  dT6 << 10.0, -4.0, 7.0, 2.0, 12.0, -1.0;  // any nodal values: dT_h is quadratic
  const Eigen::VectorXd f6 = element_thermal_load(CellType::Tri6, straight_t6(), dT6, beta);
  CHECK((f6 - reference_thermal_load(CellType::Tri6, straight_t6(), dT6, beta)).norm() <= 1e-12 * f6.norm());
}

TEST_CASE("Uniform dT: element thermal load equals K_e u* on a curved T6 cell", "[phase4][units]") {
  // Section 9: both integrands are beta dT adj(J)^T dN/dxi, degree 2 in xi even on a curved cell.
  const ElasticMaterial m = test::thermal_test_material();
  const Eigen::MatrixXd x = curved_t6();
  for (const auto plane : {PlaneCondition::PlaneStress, PlaneCondition::PlaneStrain}) {
    const double a_eff = plane == PlaneCondition::PlaneStress ? test::kAlpha_per_K : (1 + m.nu) * test::kAlpha_per_K;
    const double beta = thermal_stress_coefficient(m, test::kAlpha_per_K, plane);
    Eigen::VectorXd u_star(12);
    for (int a = 0; a < 6; ++a) u_star.segment<2>(2 * a) = a_eff * test::kdT_K * x.row(a).transpose();
    const Eigen::VectorXd Ku = element_stiffness(CellType::Tri6, x, elasticity_matrix(m, plane)) * u_star;
    const Eigen::VectorXd f = element_thermal_load(CellType::Tri6, x, Eigen::VectorXd::Constant(6, test::kdT_K), beta);
    INFO(test::plane_name(plane));
    CHECK((f - Ku).norm() <= 1e-12 * f.norm());
  }
}
