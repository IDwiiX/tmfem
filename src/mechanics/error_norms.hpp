#pragma once

#include <functional>

#include <Eigen/Dense>

#include "mechanics/elastic2d.hpp"

namespace tmfem {

using DisplacementField = std::function<Eigen::Vector2d(const Eigen::Vector2d& x_m)>;
using StressField = std::function<Eigen::Vector3d(const Eigen::Vector2d& x_m)>;

// Absolute errors and exact-solution norms over the mesh (docs/learning/phase2.md section 10):
//   L2      = ||u - u_h||_L2,                  L2_exact     = ||u||_L2
//   energy  = (int (s - s_h)^T D^-1 (s - s_h))^1/2, energy_exact = (int s^T D^-1 s)^1/2
// integrated with the collapsed Gauss rule of n_gauss^2 points per cell. h_m = sqrt(area / n_cells).
struct ElasticErrors {
  double L2;
  double L2_exact;
  double energy;
  double energy_exact;
  double h_m;
};

ElasticErrors elastic_errors(const Elastic2DProblem& p, const Eigen::VectorXd& u_m, const DisplacementField& u_exact,
                             const StressField& sigma_exact, int n_gauss = 4);

}  // namespace tmfem
