#pragma once

#include <Eigen/Dense>

#include "mechanics/elasticity.hpp"

namespace tmfem {

// Infinite plate with a hole of radius a under remote stress S along y; theta measured from +y toward
// +x, so x = r sin(theta), y = r cos(theta) (docs/learning/phase2.md section 9).
struct KirschPlate {
  double S_Pa;
  double a_m;
};

// Angle from the +y axis toward +x.
double theta_from_y(const Eigen::Vector2d& x_m);

// [sigma_rr, sigma_tt, sigma_rt] at (r, theta). Pa.
Eigen::Vector3d kirsch_stress_polar(const KirschPlate& k, double r_m, double theta_rad);

// [sigma_xx, sigma_yy, sigma_xy] at x. Pa.
Eigen::Vector3d kirsch_stress_cartesian(const KirschPlate& k, const Eigen::Vector2d& x_m);

// Exact displacement (u_x, u_y) in m; depends on the plane condition through kappa.
Eigen::Vector2d kirsch_displacement(const KirschPlate& k, const ElasticMaterial& m, PlaneCondition plane,
                                    const Eigen::Vector2d& x_m);

// Cartesian [xx, yy, xy] -> polar [rr, tt, rt] with e_r = (sin t, cos t), e_t = (cos t, -sin t).
Eigen::Vector3d cartesian_to_polar_from_y(const Eigen::Vector3d& sigma, double theta_rad);

}  // namespace tmfem
