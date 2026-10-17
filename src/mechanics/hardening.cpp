#include "mechanics/hardening.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace tmfem {

LinearHardening::LinearHardening(double sigma_y0_Pa, double H_Pa) : sigma_y0_Pa_(sigma_y0_Pa), H_Pa_(H_Pa) {
  if (sigma_y0_Pa <= 0.0 || H_Pa < 0.0) throw std::invalid_argument("LinearHardening: need sigma_y0 > 0, H >= 0");
}

double LinearHardening::yield_stress_Pa(double alpha) const { return sigma_y0_Pa_ + H_Pa_ * alpha; }

double LinearHardening::slope_Pa(double) const { return H_Pa_; }

LudwikHardening::LudwikHardening(double sigma_y0_Pa, double K_Pa, double m) : sigma_y0_Pa_(sigma_y0_Pa), K_Pa_(K_Pa), m_(m) {
  if (sigma_y0_Pa <= 0.0 || K_Pa < 0.0 || m <= 0.0) {
    throw std::invalid_argument("LudwikHardening: need sigma_y0 > 0, K >= 0, m > 0");
  }
}

double LudwikHardening::yield_stress_Pa(double alpha) const {
  return alpha > 0.0 ? sigma_y0_Pa_ + K_Pa_ * std::pow(alpha, m_) : sigma_y0_Pa_;
}

double LudwikHardening::slope_Pa(double alpha) const {
  if (alpha > 0.0) return m_ * K_Pa_ * std::pow(alpha, m_ - 1.0);
  return m_ < 1.0 ? std::numeric_limits<double>::infinity() : (m_ == 1.0 ? K_Pa_ : 0.0);
}

TabulatedHardening::TabulatedHardening(std::vector<double> alpha, std::vector<double> sigma_Pa)
    : alpha_(std::move(alpha)), sigma_(std::move(sigma_Pa)) {
  if (alpha_.size() < 2 || alpha_.size() != sigma_.size()) {
    throw std::invalid_argument("TabulatedHardening: need at least 2 points of (alpha, sigma)");
  }
  if (alpha_.front() != 0.0) throw std::invalid_argument("TabulatedHardening: alpha must start at 0");
  for (std::size_t i = 1; i < alpha_.size(); ++i) {
    if (!(alpha_[i] > alpha_[i - 1])) throw std::invalid_argument("TabulatedHardening: alpha must increase strictly");
    if (sigma_[i] < sigma_[i - 1]) throw std::invalid_argument("TabulatedHardening: softening is not supported");
  }
}

TabulatedHardening TabulatedHardening::from_true_curve(const std::vector<double>& eps_true,
                                                       const std::vector<double>& sigma_true_Pa, double E_Pa,
                                                       double sigma_start_Pa) {
  std::vector<double> alpha, sigma;
  for (std::size_t i = 0; i < eps_true.size(); ++i) {
    if (sigma_true_Pa[i] < sigma_start_Pa) continue;
    const double a = eps_true[i] - sigma_true_Pa[i] / E_Pa;
    // First plastic point: clipped to alpha = 0. Only data round-off may be clipped (-3.2e-6 for tensile_true.csv);
    // anything more negative means a wrong E or sigma_start, so it is an error, not something to hide.
    if (alpha.empty() && a < -1e-5) {
      throw std::invalid_argument("TabulatedHardening::from_true_curve: alpha < -1e-5 at the first plastic point "
                                  "(wrong E or sigma_start?)");
    }
    alpha.push_back(alpha.empty() ? 0.0 : a);
    sigma.push_back(sigma_true_Pa[i]);
  }
  return TabulatedHardening(std::move(alpha), std::move(sigma));
}

std::size_t TabulatedHardening::segment(double alpha) const {
  const auto it = std::upper_bound(alpha_.begin(), alpha_.end(), alpha);
  const auto i = static_cast<std::size_t>(std::max<std::ptrdiff_t>(it - alpha_.begin() - 1, 0));
  return std::min(i, alpha_.size() - 2);
}

double TabulatedHardening::yield_stress_Pa(double alpha) const {
  const std::size_t i = segment(alpha);
  return sigma_[i] + slope_Pa(alpha) * (alpha - alpha_[i]);
}

double TabulatedHardening::slope_Pa(double alpha) const {
  const std::size_t i = segment(alpha);
  return (sigma_[i + 1] - sigma_[i]) / (alpha_[i + 1] - alpha_[i]);
}

double plastic_strain_at_stress(const HardeningLaw& h, double sigma_Pa) {
  if (sigma_Pa <= h.yield_stress_Pa(0.0)) return 0.0;
  double lo = 0.0, hi = 1e-3;
  while (h.yield_stress_Pa(hi) < sigma_Pa) {
    lo = hi;
    hi *= 2.0;
    if (hi > 1e3) throw std::runtime_error("plastic_strain_at_stress: stress not reached by the hardening law");
  }
  for (int it = 0; it < 200 && hi - lo > 1e-16 * hi; ++it) {
    const double mid = 0.5 * (lo + hi);
    (h.yield_stress_Pa(mid) < sigma_Pa ? lo : hi) = mid;
  }
  return 0.5 * (lo + hi);
}

}  // namespace tmfem
