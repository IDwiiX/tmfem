#pragma once

#include <Eigen/Dense>

#include "mechanics/constitutive.hpp"
#include "mechanics/elasticity.hpp"

namespace tmfem {

// 2D response of a 3D constitutive model (docs/learning/phase3.md section 8). In-plane Voigt order
// [xx, yy, xy] with engineering shear strain, as in phase 2.
struct PlaneResponse {
  Eigen::Vector3d stress_Pa = Eigen::Vector3d::Zero();
  Eigen::Matrix3d tangent_Pa = Eigen::Matrix3d::Zero();
  double sigma_zz_Pa = 0.0;
  MaterialState state;
  bool plastic = false;
};

// Plane strain: eps_zz = 0, so the 3D model runs exactly and the 2D quantities are its {xx, yy, xy} block.
// Plane stress: a local Newton on eps_zz until sigma_zz = 0 (start: the converged eps_zz, stored back in the
// state), then the statically condensed tangent C_ab - C_az C_zb / C_zz.
PlaneResponse plane_update(const ConstitutiveModel& model, PlaneCondition plane, const Eigen::Vector3d& eps,
                           const MaterialState& converged);

}  // namespace tmfem
