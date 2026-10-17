#include "mechanics/tensile_fit.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace tmfem {

namespace {

// The fit runs in MPa internally: in Pa the m column of J (K alpha^m ln alpha ~ 1e8) and the sigma_y0
// column (1) differ by 1e8 and J^T J would be needlessly ill-conditioned. Inputs and outputs stay in Pa.
constexpr double kScale = 1e6;
constexpr int kMaxIterations = 500;

using Vec3 = Eigen::Vector3d;

double alpha_pow(double a, double m) { return a > 0.0 ? std::pow(a, m) : 0.0; }

double model(const Vec3& p, double a) { return p(0) + p(1) * alpha_pow(a, p(2)); }

// d model / d (sigma_y0, K, m); alpha^m ln alpha -> 0 as alpha -> 0 (m > 0).
Vec3 gradient(const Vec3& p, double a) {
  const double am = alpha_pow(a, p(2));
  return {1.0, am, a > 0.0 ? p(1) * am * std::log(a) : 0.0};
}

double rss(const Vec3& p, const std::vector<double>& a, const std::vector<double>& s) {
  if (!(p(2) > 0.0)) return std::numeric_limits<double>::infinity();
  double sum = 0.0;
  for (std::size_t i = 0; i < a.size(); ++i) sum += (s[i] - model(p, a[i])) * (s[i] - model(p, a[i]));
  return sum;
}

struct Solution {
  Vec3 p;
  Eigen::MatrixXd J;  // n x n_free
  Eigen::VectorXd r;
  int iterations;
};

std::vector<int> free_indices(const std::array<bool, 3>& free) {
  std::vector<int> idx;
  for (int k = 0; k < 3; ++k) {
    if (free[static_cast<std::size_t>(k)]) idx.push_back(k);
  }
  return idx;
}

void evaluate(const Vec3& p, const std::vector<double>& a, const std::vector<double>& s, const std::vector<int>& idx,
              Eigen::MatrixXd& J, Eigen::VectorXd& r) {
  const auto n = static_cast<Eigen::Index>(a.size());
  J.resize(n, static_cast<Eigen::Index>(idx.size()));
  r.resize(n);
  for (Eigen::Index i = 0; i < n; ++i) {
    const auto ii = static_cast<std::size_t>(i);
    r(i) = s[ii] - model(p, a[ii]);
    const Vec3 g = gradient(p, a[ii]);
    for (std::size_t k = 0; k < idx.size(); ++k) J(i, static_cast<Eigen::Index>(k)) = g(idx[k]);
  }
}

// Levenberg-Marquardt: (J^T J + mu diag(J^T J)) d = J^T r, mu shrinks on success and grows on failure.
Solution levenberg_marquardt(const std::vector<double>& a, const std::vector<double>& s, Vec3 p,
                             const std::vector<int>& idx) {
  Eigen::MatrixXd J;
  Eigen::VectorXd r;
  double mu = 1e-3;
  int it = 0;
  for (; it < kMaxIterations; ++it) {
    evaluate(p, a, s, idx, J, r);
    const Eigen::MatrixXd A = J.transpose() * J;
    const Eigen::VectorXd g = J.transpose() * r;
    const double S = r.squaredNorm();
    bool accepted = false;
    Eigen::VectorXd d;
    while (mu < 1e20) {
      Eigen::MatrixXd M = A;
      M.diagonal() += mu * A.diagonal();
      d = M.ldlt().solve(g);
      Vec3 trial = p;
      for (std::size_t k = 0; k < idx.size(); ++k) trial(idx[k]) += d(static_cast<Eigen::Index>(k));
      if (rss(trial, a, s) < S) {
        p = trial;
        mu = std::max(mu / 3.0, 1e-15);
        accepted = true;
        break;
      }
      mu *= 3.0;
    }
    // Stop when no step lowers the sum of squares or the step is at round-off.
    if (!accepted || d.norm() <= 1e-14 * (1.0 + p.norm())) break;
  }
  evaluate(p, a, s, idx, J, r);
  return {p, J, r, it};
}

Vec3 scaled(const LudwikParams& p) { return {p.sigma_y0_Pa / kScale, p.K_Pa / kScale, p.m}; }
LudwikParams unscaled(const Vec3& p) { return {p(0) * kScale, p(1) * kScale, p(2)}; }
const Vec3 kUnit(kScale, kScale, 1.0);  // scaled -> SI per parameter

}  // namespace

double ludwik_stress_Pa(const LudwikParams& p, double alpha) { return p.sigma_y0_Pa + p.K_Pa * alpha_pow(alpha, p.m); }

double student_t95(int dof) {
  static const double t[] = {12.7062, 4.3027, 3.1824, 2.7764, 2.5706, 2.4469, 2.3646, 2.3060, 2.2622, 2.2281};
  if (dof < 1 || dof > 10) throw std::invalid_argument("student_t95: dof must be 1..10");
  return t[dof - 1];
}

LudwikParams fit_ludwik_loglog(const std::vector<double>& alpha, const std::vector<double>& sigma_Pa,
                               double sigma_y0_Pa) {
  double n = 0, sx = 0, sy = 0, sxx = 0, sxy = 0;
  for (std::size_t i = 0; i < alpha.size(); ++i) {
    if (!(alpha[i] > 1e-6) || !(sigma_Pa[i] > sigma_y0_Pa)) continue;
    const double x = std::log(alpha[i]), y = std::log(sigma_Pa[i] - sigma_y0_Pa);
    n += 1;
    sx += x;
    sy += y;
    sxx += x * x;
    sxy += x * y;
  }
  if (n < 2) throw std::invalid_argument("fit_ludwik_loglog: fewer than 2 plastic points");
  const double m = (n * sxy - sx * sy) / (n * sxx - sx * sx);
  const double lnK = (sy - m * sx) / n;
  return {sigma_y0_Pa, std::exp(lnK), m};
}

LudwikFit fit_ludwik_least_squares(const std::vector<double>& alpha, const std::vector<double>& sigma_Pa,
                                   const LudwikParams& start, std::array<bool, 3> free) {
  const std::vector<int> idx = free_indices(free);
  const int n = static_cast<int>(alpha.size());
  const int p = static_cast<int>(idx.size());
  if (p == 0 || n <= p || sigma_Pa.size() != alpha.size()) {
    throw std::invalid_argument("fit_ludwik_least_squares: need more points than free parameters");
  }
  std::vector<double> s(sigma_Pa.size());
  for (std::size_t i = 0; i < s.size(); ++i) s[i] = sigma_Pa[i] / kScale;

  const Solution sol = levenberg_marquardt(alpha, s, scaled(start), idx);
  LudwikFit fit{};
  fit.params = unscaled(sol.p);
  fit.free = free;
  fit.n = n;
  fit.dof = n - p;
  fit.t95 = student_t95(fit.dof);
  fit.iterations = sol.iterations;
  const double S = sol.r.squaredNorm();
  fit.s_Pa = std::sqrt(S / fit.dof) * kScale;
  for (Eigen::Index i = 0; i < sol.r.size(); ++i) fit.residuals_Pa.push_back(sol.r(i) * kScale);
  fit.gradient_ratio = (sol.J.transpose() * sol.r).norm() / (sol.J.norm() * sol.r.norm());

  // Covariance of the free parameters, embedded in 3 x 3 and converted to SI units.
  const Eigen::MatrixXd cov_free = (S / fit.dof) * (sol.J.transpose() * sol.J).inverse();
  fit.covariance.setZero();
  for (int a = 0; a < p; ++a) {
    for (int b = 0; b < p; ++b) {
      fit.covariance(idx[static_cast<std::size_t>(a)], idx[static_cast<std::size_t>(b)]) =
          cov_free(a, b) * kUnit(idx[static_cast<std::size_t>(a)]) * kUnit(idx[static_cast<std::size_t>(b)]);
    }
  }
  fit.standard_error = fit.covariance.diagonal().cwiseSqrt();
  fit.ci95_halfwidth = fit.t95 * fit.standard_error;

  // Jackknife: refit without each point in turn, starting from the full fit.
  std::vector<Vec3> loo;
  for (int k = 0; k < n; ++k) {
    std::vector<double> a_k, s_k;
    for (int i = 0; i < n; ++i) {
      if (i == k) continue;
      a_k.push_back(alpha[static_cast<std::size_t>(i)]);
      s_k.push_back(s[static_cast<std::size_t>(i)]);
    }
    loo.push_back(levenberg_marquardt(a_k, s_k, sol.p, idx).p);
  }
  Vec3 mean = Vec3::Zero();
  for (const Vec3& q : loo) mean += q / n;
  Vec3 ss = Vec3::Zero();
  for (const Vec3& q : loo) ss += (q - mean).cwiseAbs2();
  fit.jackknife_se = ((n - 1.0) / n * ss).cwiseSqrt().cwiseProduct(kUnit);
  return fit;
}

double ludwik_band_halfwidth_Pa(const LudwikFit& fit, double alpha) {
  const Vec3 p = scaled(fit.params);
  const Vec3 g = gradient(p, alpha).cwiseProduct(Vec3(1.0, 1.0, kScale));  // d sigma[Pa] / d (Pa, Pa, 1)
  return fit.t95 * std::sqrt(g.dot(fit.covariance * g));
}

}  // namespace tmfem
