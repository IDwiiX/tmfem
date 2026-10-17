#include "mechanics/j2_plasticity.hpp"

#include <cmath>
#include <stdexcept>
#include <string>

namespace tmfem {

namespace {

const double kSqrt32 = std::sqrt(1.5);
const Voigt4 kM(1.0, 1.0, 1.0, 0.0);

// A trial state counts as plastic only if f_trial exceeds this fraction of sigma_y: a converged point
// lying exactly on the yield surface must stay elastic when re-evaluated (section 5, step 2).
constexpr double kYieldTolerance = 1e-12;
// Scalar Newton stops at |r| <= this fraction of the initial yield stress (about 3e-4 Pa here).
constexpr double kMultiplierTolerance = 1e-12;
constexpr int kMaxMultiplierIterations = 200;

// Norm of a stress-like deviator: the xy entry appears twice in s:s.
double deviator_norm(const Voigt4& s) { return std::sqrt(s(0) * s(0) + s(1) * s(1) + s(2) * s(2) + 2.0 * s(3) * s(3)); }

}  // namespace

PlasticMultiplier solve_plastic_multiplier(const HardeningLaw& h, double G_Pa, double q_trial_Pa, double alpha_n) {
  const double f_trial = q_trial_Pa - h.yield_stress_Pa(alpha_n);
  if (!(f_trial > 0.0)) throw std::invalid_argument("solve_plastic_multiplier: trial state is not plastic");
  const double tol = kMultiplierTolerance * h.yield_stress_Pa(0.0);
  double lo = 0.0, hi = f_trial / (3.0 * G_Pa);
  double dg = hi;  // H is finite here even for Ludwik at alpha_n = 0
  for (int it = 1; it <= kMaxMultiplierIterations; ++it) {
    const double r = q_trial_Pa - 3.0 * G_Pa * dg - h.yield_stress_Pa(alpha_n + dg);
    if (std::abs(r) <= tol) return {dg, it};
    (r > 0.0 ? lo : hi) = dg;
    const double newton = dg + r / (3.0 * G_Pa + h.slope_Pa(alpha_n + dg));
    dg = (newton > lo && newton < hi) ? newton : 0.5 * (lo + hi);
    if (hi - lo <= 1e-17 * hi) return {dg, it};  // bracket at round-off: r cannot get smaller
  }
  throw std::runtime_error("solve_plastic_multiplier: no convergence in " + std::to_string(kMaxMultiplierIterations) +
                           " iterations");
}

J2Plasticity::J2Plasticity(const ElasticMaterial& m, std::shared_ptr<const HardeningLaw> hardening)
    : m_(m),
      hardening_(std::move(hardening)),
      C_(elastic_tensor4(m)),
      I_dev_(deviator_projector()),
      G_(shear_modulus(m)),
      K_b_(bulk_modulus(m)) {
  if (!hardening_) throw std::invalid_argument("J2Plasticity: no hardening law");
}

MaterialUpdate J2Plasticity::update(const Voigt4& eps, const MaterialState& converged) const {
  MaterialUpdate out;
  out.state = converged;
  out.stress_Pa = C_ * (eps - converged.eps_p);  // elastic trial
  out.tangent_Pa = C_;

  Voigt4 s = out.stress_Pa;
  s.head<3>().array() -= (s(0) + s(1) + s(2)) / 3.0;
  const double s_norm = deviator_norm(s);
  const double q_trial = kSqrt32 * s_norm;
  const double sigma_y_n = hardening_->yield_stress_Pa(converged.alpha);
  if (q_trial - sigma_y_n <= kYieldTolerance * sigma_y_n) return out;

  // Radial return: n_{n+1} = n_trial, q_{n+1} = q_trial - 3 G dg (section 5).
  const double dg = solve_plastic_multiplier(*hardening_, G_, q_trial, converged.alpha).dgamma;
  const Voigt4 n = s / s_norm;
  out.stress_Pa -= 2.0 * G_ * kSqrt32 * dg * n;
  Voigt4 n_eng = n;
  n_eng(3) *= 2.0;  // engineering shear in the plastic strain (section 6)
  out.state.eps_p += kSqrt32 * dg * n_eng;
  out.state.alpha += dg;
  out.plastic = true;

  // Consistent tangent: K_b m m^T + 2 G theta I_dev - 2 G theta_bar n n^T (section 7).
  const double H = hardening_->slope_Pa(out.state.alpha);
  const double ratio = 3.0 * G_ * dg / q_trial;
  const double theta = 1.0 - ratio;
  const double theta_bar = 3.0 * G_ / (3.0 * G_ + H) - ratio;
  out.tangent_Pa = K_b_ * kM * kM.transpose() + 2.0 * G_ * theta * I_dev_ - 2.0 * G_ * theta_bar * n * n.transpose();
  return out;
}

}  // namespace tmfem
