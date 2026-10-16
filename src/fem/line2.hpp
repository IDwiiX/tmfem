#pragma once

#include <functional>

#include <Eigen/Dense>

namespace tmfem {

enum class MassType { Consistent, Lumped };

// Linear 2-node element on [x1, x2]; reference coordinate xi in [-1, 1].
Eigen::Vector2d line2_shape(double xi);
double line2_coordinate(double x1_m, double x2_m, double xi);

// K_e = int k N_i' N_j' dx (docs/learning/phase1.md section 4).
Eigen::Matrix2d element_stiffness(double k, double x1_m, double x2_m);

// M_e = int rho_c N_i N_j dx (consistent) or its row-sum diagonal (lumped).
Eigen::Matrix2d element_capacity(double rho_c, double x1_m, double x2_m, MassType mass);

// F_e = int q(x) N_i dx by 3-point Gauss.
Eigen::Vector2d element_source(const std::function<double(double x_m)>& q, double x1_m, double x2_m);

}  // namespace tmfem
