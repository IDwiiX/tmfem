#include "mechanics/constitutive.hpp"

#include <cmath>

namespace tmfem {

namespace {

const Voigt4 kM(1.0, 1.0, 1.0, 0.0);

}  // namespace

double bulk_modulus(const ElasticMaterial& m) { return m.E_Pa / (3.0 * (1.0 - 2.0 * m.nu)); }

Eigen::Matrix4d deviator_projector() {
  Eigen::Matrix4d P = Eigen::Vector4d(1.0, 1.0, 1.0, 0.5).asDiagonal();
  P -= kM * kM.transpose() / 3.0;
  return P;
}

Eigen::Matrix4d elastic_tensor4(const ElasticMaterial& m) {
  return bulk_modulus(m) * kM * kM.transpose() + 2.0 * shear_modulus(m) * deviator_projector();
}

double equivalent_stress(const Voigt4& s) {
  const double p = (s(0) + s(1) + s(2)) / 3.0;
  const double dx = s(0) - p, dy = s(1) - p, dz = s(2) - p;
  return std::sqrt(1.5 * (dx * dx + dy * dy + dz * dz + 2.0 * s(3) * s(3)));
}

LinearElasticModel::LinearElasticModel(const ElasticMaterial& m) : m_(m), C_(elastic_tensor4(m)) {}

MaterialUpdate LinearElasticModel::update(const Voigt4& eps, const MaterialState& converged) const {
  MaterialUpdate out;
  out.stress_Pa = C_ * eps;
  out.tangent_Pa = C_;
  out.state = converged;
  return out;
}

}  // namespace tmfem
