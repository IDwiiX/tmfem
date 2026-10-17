#pragma once

#include <Eigen/Dense>

#include "mechanics/elasticity.hpp"

namespace tmfem {

// docs/learning/phase3.md section 6. 4-vectors in the order [xx, yy, zz, xy]: strains carry engineering
// shear (gamma_xy = 2 eps_xy), stresses carry the tensor component sigma_xy.
using Voigt4 = Eigen::Vector4d;

// History at one material point, as converged at the end of the last load step.
struct MaterialState {
  Voigt4 eps_p = Voigt4::Zero();  // plastic strain, engineering shear
  double alpha = 0.0;             // equivalent plastic strain
  double eps_zz = 0.0;            // total eps_zz; used only by plane stress (section 8)
};

struct MaterialUpdate {
  Voigt4 stress_Pa = Voigt4::Zero();
  Eigen::Matrix4d tangent_Pa = Eigen::Matrix4d::Zero();  // d stress / d strain (consistent, section 7)
  MaterialState state;
  bool plastic = false;
};

// Constitutive model: total strain at the end of the step plus the converged state -> stress, tangent and
// new state. Equivalent to the strain-increment form, since sigma_trial = C:(eps - eps_p,n) = sigma_n + C:d_eps.
class ConstitutiveModel {
 public:
  virtual ~ConstitutiveModel() = default;
  virtual MaterialUpdate update(const Voigt4& eps, const MaterialState& converged) const = 0;
  virtual const ElasticMaterial& elastic() const = 0;
  // Stress scale for tolerances (Pa): the initial yield stress, or E / 1000 for an elastic model.
  virtual double stress_scale_Pa() const = 0;
};

// Isotropic elasticity in 4-vector form, C = K_b m m^T + 2 G I_dev (section 6).
Eigen::Matrix4d elastic_tensor4(const ElasticMaterial& m);
// I_dev = diag(1, 1, 1, 1/2) - (1/3) m m^T: engineering strain -> tensor deviator.
Eigen::Matrix4d deviator_projector();
double bulk_modulus(const ElasticMaterial& m);

// von Mises equivalent stress q = sqrt(3/2 s:s) of a stress 4-vector (section 2).
double equivalent_stress(const Voigt4& stress);

class LinearElasticModel : public ConstitutiveModel {
 public:
  explicit LinearElasticModel(const ElasticMaterial& m);
  MaterialUpdate update(const Voigt4& eps, const MaterialState& converged) const override;
  const ElasticMaterial& elastic() const override { return m_; }
  double stress_scale_Pa() const override { return 1e-3 * m_.E_Pa; }

 private:
  ElasticMaterial m_;
  Eigen::Matrix4d C_;
};

}  // namespace tmfem
