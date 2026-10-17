#pragma once

#include <array>
#include <vector>

#include <Eigen/Dense>

namespace tmfem {

// Ludwik law sigma = sigma_y0 + K alpha^m (docs/learning/phase3.md sections 4 and 11).
struct LudwikParams {
  double sigma_y0_Pa;
  double K_Pa;
  double m;
};

double ludwik_stress_Pa(const LudwikParams& p, double alpha);

// Linear regression of ln(sigma - sigma_y0) on ln(alpha) with sigma_y0 fixed, over the points with
// alpha > 1e-6 (REFERENCE_VALUES.md P3-a).
LudwikParams fit_ludwik_loglog(const std::vector<double>& alpha, const std::vector<double>& sigma_Pa,
                               double sigma_y0_Pa);

// Nonlinear least squares on stress residuals, parameters (sigma_y0, K, m), each free or held at its start
// value. Levenberg-Marquardt with the analytic Jacobian [1, alpha^m, K alpha^m ln alpha]. Uncertainty:
// covariance s^2 (J^T J)^-1 with s^2 = RSS / (n - p), Student-t 95% half-widths, and leave-one-out jackknife
// standard errors. Fixed parameters get zero uncertainty. Deterministic.
struct LudwikFit {
  LudwikParams params;
  std::array<bool, 3> free;
  Eigen::Matrix3d covariance;    // units: Pa, Pa, 1
  Eigen::Vector3d standard_error;
  Eigen::Vector3d ci95_halfwidth;
  Eigen::Vector3d jackknife_se;
  double s_Pa;                   // residual standard deviation sqrt(RSS / (n - p))
  std::vector<double> residuals_Pa;  // data - model
  int n;
  int dof;
  double t95;                    // Student t(0.975, dof)
  double gradient_ratio;         // ||J^T r|| / (||J|| ||r||) at the solution
  int iterations;
};

LudwikFit fit_ludwik_least_squares(const std::vector<double>& alpha, const std::vector<double>& sigma_Pa,
                                   const LudwikParams& start, std::array<bool, 3> free = {true, true, true});

// Half-width of the 95% confidence band of the fitted curve at alpha (delta method, t * sqrt(g^T Cov g)).
double ludwik_band_halfwidth_Pa(const LudwikFit& fit, double alpha);

// Two-sided 95% Student t quantile t(0.975, dof) for dof = 1..10.
double student_t95(int dof);

}  // namespace tmfem
