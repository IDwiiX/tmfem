#include "thermal/diffusion_analytic.hpp"

#include <cmath>
#include <stdexcept>

namespace tmfem {

namespace {
const double kPi = std::acos(-1.0);
}

double erfinv(double y) {
  if (!(std::abs(y) < 1.0)) throw std::domain_error("erfinv: need |y| < 1");
  if (y == 0.0) return 0.0;
  // Bisection to a good bracket, then Newton: d erf / dz = 2/sqrt(pi) exp(-z^2).
  double lo = 0.0;
  double hi = 6.0;
  const double target = std::abs(y);
  for (int i = 0; i < 30; ++i) {
    const double mid = 0.5 * (lo + hi);
    (std::erf(mid) < target ? lo : hi) = mid;
  }
  double z = 0.5 * (lo + hi);
  for (int i = 0; i < 10; ++i) {
    const double dz = (std::erf(z) - target) / (2.0 / std::sqrt(kPi) * std::exp(-z * z));
    z -= dz;
    if (std::abs(dz) < 1e-16 * z) break;
  }
  return y > 0.0 ? z : -z;
}

double erf_profile(double x_m, double t_s, double D_m2_s, double Cs, double C0) {
  return Cs - (Cs - C0) * std::erf(x_m / (2.0 * std::sqrt(D_m2_s * t_s)));
}

double time_to_reach(double C, double x_m, double D_m2_s, double Cs, double C0) {
  const double z = erfinv((Cs - C) / (Cs - C0));
  return x_m * x_m / (4.0 * D_m2_s * z * z);
}

double predeposition_dose(double Cs, double D_m2_s, double t_s) {
  return 2.0 * Cs * std::sqrt(D_m2_s * t_s / kPi);
}

double thin_film_profile(double x_m, double t_s, double D_m2_s, double M) {
  return M / std::sqrt(kPi * D_m2_s * t_s) * std::exp(-x_m * x_m / (4.0 * D_m2_s * t_s));
}

}  // namespace tmfem
