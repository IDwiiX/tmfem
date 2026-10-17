#pragma once

#include <memory>

#include "mechanics/constitutive.hpp"
#include "mechanics/hardening.hpp"

namespace tmfem {

struct PlasticMultiplier {
  double dgamma;
  int iterations;
};

// Solves the consistency equation r(dg) = q_trial - 3 G dg - sigma_y(alpha_n + dg) = 0 on its bracket
// [0, f_trial / 3G] by Newton, started at the upper end and safeguarded by bisection
// (docs/learning/phase3.md section 5). Requires f_trial = q_trial - sigma_y(alpha_n) > 0.
PlasticMultiplier solve_plastic_multiplier(const HardeningLaw& h, double G_Pa, double q_trial_Pa, double alpha_n);

// Small-strain von Mises plasticity with isotropic hardening: radial return (section 5) and the
// consistent tangent (section 7), on the 4-vectors of section 6.
class J2Plasticity : public ConstitutiveModel {
 public:
  J2Plasticity(const ElasticMaterial& m, std::shared_ptr<const HardeningLaw> hardening);

  MaterialUpdate update(const Voigt4& eps, const MaterialState& converged) const override;
  const ElasticMaterial& elastic() const override { return m_; }
  double stress_scale_Pa() const override { return hardening_->yield_stress_Pa(0.0); }
  const HardeningLaw& hardening() const { return *hardening_; }

 private:
  ElasticMaterial m_;
  std::shared_ptr<const HardeningLaw> hardening_;
  Eigen::Matrix4d C_;
  Eigen::Matrix4d I_dev_;
  double G_, K_b_;
};

}  // namespace tmfem
