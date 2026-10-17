#include "mechanics/plane_reduction.hpp"

#include <cmath>
#include <stdexcept>

namespace tmfem {

namespace {

constexpr int kInPlane[3] = {0, 1, 3};  // [xx, yy, xy] inside [xx, yy, zz, xy]
// Plane stress: the local Newton stops at |sigma_zz| <= this fraction of the model's stress scale.
constexpr double kZZTolerance = 1e-12;
constexpr int kMaxZZIterations = 50;

PlaneResponse in_plane_block(const MaterialUpdate& u) {
  PlaneResponse r;
  for (int i = 0; i < 3; ++i) {
    r.stress_Pa(i) = u.stress_Pa(kInPlane[i]);
    for (int j = 0; j < 3; ++j) r.tangent_Pa(i, j) = u.tangent_Pa(kInPlane[i], kInPlane[j]);
  }
  r.sigma_zz_Pa = u.stress_Pa(2);
  r.state = u.state;
  r.plastic = u.plastic;
  return r;
}

// sigma_zz = 0: Newton on eps_zz, each evaluation a full 3D update from the converged state, started from
// the converged eps_zz; then the condensed tangent D = C_ab - C_az C_zb / C_zz (docs/learning/phase3.md section 8).
PlaneResponse plane_stress(const ConstitutiveModel& model, const Eigen::Vector3d& eps, const MaterialState& converged) {
  const double tol = kZZTolerance * model.stress_scale_Pa();
  double eps_zz = converged.eps_zz;
  MaterialUpdate u = model.update(Voigt4(eps(0), eps(1), eps_zz, eps(2)), converged);
  for (int it = 0; std::abs(u.stress_Pa(2)) > tol; ++it) {
    if (it == kMaxZZIterations) throw std::runtime_error("plane_update: sigma_zz = 0 iteration did not converge");
    eps_zz -= u.stress_Pa(2) / u.tangent_Pa(2, 2);
    u = model.update(Voigt4(eps(0), eps(1), eps_zz, eps(2)), converged);
  }
  PlaneResponse r = in_plane_block(u);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      r.tangent_Pa(i, j) -= u.tangent_Pa(kInPlane[i], 2) * u.tangent_Pa(2, kInPlane[j]) / u.tangent_Pa(2, 2);
    }
  }
  r.state.eps_zz = eps_zz;
  return r;
}

}  // namespace

PlaneResponse plane_update(const ConstitutiveModel& model, PlaneCondition plane, const Eigen::Vector3d& eps,
                           const MaterialState& converged) {
  if (plane == PlaneCondition::PlaneStrain) {
    return in_plane_block(model.update(Voigt4(eps(0), eps(1), 0.0, eps(2)), converged));
  }
  return plane_stress(model, eps, converged);
}

}  // namespace tmfem
