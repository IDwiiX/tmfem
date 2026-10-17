#include <catch2/catch_test_macros.hpp>

#include <cmath>

#include <Eigen/Eigenvalues>

#include "mechanics/elasticity.hpp"
#include "phase2/kirsch_common.hpp"

// D matrices, von Mises and element stiffness (docs/learning/phase2.md sections 3, 4 and 7).
// Not a reference row.

using namespace tmfem;

TEST_CASE("Elasticity matrices: symmetry, definiteness, plane strain substitution", "[phase2][mechanics]") {
  const ElasticMaterial m = test::phase2_material();
  for (const auto plane : {PlaneCondition::PlaneStress, PlaneCondition::PlaneStrain}) {
    const Eigen::Matrix3d D = elasticity_matrix(m, plane);
    CHECK((D - D.transpose()).norm() == 0.0);
    CHECK(Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d>(D).eigenvalues().minCoeff() > 0.0);
    CHECK(std::abs(D(2, 2) - shear_modulus(m)) < 1e-12 * shear_modulus(m));
  }
  // Plane strain with (E, nu) = plane stress with E* = E/(1 - nu^2), nu* = nu/(1 - nu).
  const ElasticMaterial star{m.E_Pa / (1 - m.nu * m.nu), m.nu / (1 - m.nu)};
  const Eigen::Matrix3d Dpe = elasticity_matrix(m, PlaneCondition::PlaneStrain);
  CHECK((Dpe - elasticity_matrix(star, PlaneCondition::PlaneStress)).norm() < 1e-12 * Dpe.norm());

  CHECK(std::abs(kappa(m, PlaneCondition::PlaneStress) - (3 - m.nu) / (1 + m.nu)) < 1e-15);
  CHECK(std::abs(kappa(m, PlaneCondition::PlaneStrain) - (3 - 4 * m.nu)) < 1e-15);
}

TEST_CASE("von Mises: uniaxial and pure shear", "[phase2][mechanics]") {
  CHECK(std::abs(von_mises({100.0, 0.0, 0.0}, 0.0) - 100.0) < 1e-12);
  CHECK(std::abs(von_mises({0.0, 0.0, 50.0}, 0.0) - std::sqrt(3.0) * 50.0) < 1e-12);
  CHECK(std::abs(von_mises({70.0, 70.0, 0.0}, 70.0)) < 1e-12);  // hydrostatic
}

TEST_CASE("Element stiffness: symmetric, three rigid-body modes, T3 and curved T6", "[phase2][mechanics]") {
  const Eigen::Matrix3d D = elasticity_matrix(test::phase2_material(), PlaneCondition::PlaneStress);
  Eigen::MatrixXd x3(3, 2);
  x3 << 0.0, 0.0, 0.03, 0.005, 0.01, 0.02;
  Eigen::MatrixXd x6(6, 2);  // same corners, mid-side nodes pushed off the chords (curved edges)
  x6 << 0.0, 0.0, 0.03, 0.005, 0.01, 0.02, 0.015, 0.0005, 0.021, 0.013, 0.004, 0.011;
  for (const auto& [type, x_e] : {std::pair{CellType::Tri3, x3}, std::pair{CellType::Tri6, x6}}) {
    const Eigen::MatrixXd K = element_stiffness(type, x_e, D);
    CHECK((K - K.transpose()).norm() < 1e-12 * K.norm());
    const Eigen::VectorXd ev = Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>(K).eigenvalues();
    const double tol = 1e-10 * ev.maxCoeff();
    int n_zero = 0;
    for (Eigen::Index i = 0; i < ev.size(); ++i) n_zero += std::abs(ev(i)) < tol ? 1 : 0;
    CHECK(n_zero == 3);
    CHECK(ev.minCoeff() > -tol);
    // Translations and the infinitesimal rotation (-y, x) produce no force.
    const auto n = x_e.rows();
    Eigen::VectorXd tx(2 * n), ty(2 * n), rot(2 * n);
    for (Eigen::Index a = 0; a < n; ++a) {
      tx.segment<2>(2 * a) << 1, 0;
      ty.segment<2>(2 * a) << 0, 1;
      rot.segment<2>(2 * a) << -x_e(a, 1), x_e(a, 0);
    }
    CHECK((K * tx).norm() < 1e-10 * K.norm());
    CHECK((K * ty).norm() < 1e-10 * K.norm());
    CHECK((K * rot).norm() < 1e-10 * K.norm() * 0.03);
  }
}
