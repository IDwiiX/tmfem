#include "mechanics/kirsch_analytic.hpp"

#include <cmath>

namespace tmfem {

double theta_from_y(const Eigen::Vector2d& x_m) { return std::atan2(x_m(0), x_m(1)); }

Eigen::Vector3d kirsch_stress_polar(const KirschPlate& k, double r_m, double t) {
  const double S = k.S_Pa;
  const double q2 = (k.a_m / r_m) * (k.a_m / r_m), q4 = q2 * q2;
  return {S / 2 * (1 - q2) + S / 2 * (1 - 4 * q2 + 3 * q4) * std::cos(2 * t),
          S / 2 * (1 + q2) - S / 2 * (1 + 3 * q4) * std::cos(2 * t),
          -S / 2 * (1 + 2 * q2 - 3 * q4) * std::sin(2 * t)};
}

Eigen::Vector3d kirsch_stress_cartesian(const KirschPlate& k, const Eigen::Vector2d& x_m) {
  // Textbook form for load along X with theta from X, then relabel X = y, Y = x (a reflection, so the
  // components just swap): sigma_yy = sigma_XX, sigma_xx = sigma_YY, sigma_xy = sigma_XY.
  const double S = k.S_Pa;
  const double r = x_m.norm(), t = theta_from_y(x_m);
  const double q2 = (k.a_m / r) * (k.a_m / r), q4 = q2 * q2;
  const double c2 = std::cos(2 * t), c4 = std::cos(4 * t), s2 = std::sin(2 * t), s4 = std::sin(4 * t);
  const double sXX = S * (1 - q2 * (1.5 * c2 + c4) + 1.5 * q4 * c4);
  const double sYY = S * (-q2 * (0.5 * c2 - c4) - 1.5 * q4 * c4);
  const double sXY = S * (-q2 * (0.5 * s2 + s4) + 1.5 * q4 * s4);
  return {sYY, sXX, sXY};
}

Eigen::Vector2d kirsch_displacement(const KirschPlate& k, const ElasticMaterial& m, PlaneCondition plane,
                                    const Eigen::Vector2d& x_m) {
  const double r = x_m.norm(), t = theta_from_y(x_m);
  const double a = k.a_m, kap = kappa(m, plane), c = k.S_Pa / (4.0 * shear_modulus(m));
  const double a2_r = a * a / r, a4_r3 = a * a * a * a / (r * r * r);
  const double c2 = std::cos(2 * t), s2 = std::sin(2 * t);
  const double u_r = c * (r * ((kap - 1) / 2 + c2) + a2_r * (1 + (1 + kap) * c2) - a4_r3 * c2);
  const double u_t = c * ((1 - kap) * a2_r - r - a4_r3) * s2;
  const double st = std::sin(t), ct = std::cos(t);
  return {u_r * st + u_t * ct, u_r * ct - u_t * st};
}

Eigen::Vector3d cartesian_to_polar_from_y(const Eigen::Vector3d& s, double t) {
  const Eigen::Vector2d er(std::sin(t), std::cos(t)), et(std::cos(t), -std::sin(t));
  Eigen::Matrix2d sig;
  sig << s(0), s(2), s(2), s(1);
  return {er.dot(sig * er), et.dot(sig * et), er.dot(sig * et)};
}

}  // namespace tmfem
