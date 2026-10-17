#include "mechanics/material_point.hpp"

#include <stdexcept>
#include <string>

namespace tmfem {

namespace {

constexpr double kStressTolerance_Pa = 1e-4;  // about 2e-13 of the stresses here
constexpr int kMaxIterations = 50;

constexpr Control S = Control::Stress;
constexpr Control E = Control::Strain;

}  // namespace

DriverResult MaterialPointDriver::run(const std::vector<DriverStep>& steps, const MaterialState& start) const {
  DriverResult out;
  Voigt4 eps = Voigt4::Zero();
  MaterialState state = start;
  for (std::size_t k = 0; k < steps.size(); ++k) {
    const DriverStep& step = steps[k];
    std::vector<int> s;  // stress-controlled components
    for (int i = 0; i < 4; ++i) {
      if (step.control[static_cast<std::size_t>(i)] == Control::Strain) {
        eps(i) = step.target(i);
      } else {
        s.push_back(i);
      }
    }
    const auto n_s = static_cast<Eigen::Index>(s.size());
    MaterialUpdate u = model_.update(eps, state);
    int it = 0;
    for (;; ++it) {
      Eigen::VectorXd r(n_s);
      Eigen::MatrixXd C_ss(n_s, n_s);
      for (Eigen::Index a = 0; a < n_s; ++a) {
        r(a) = u.stress_Pa(s[static_cast<std::size_t>(a)]) - step.target(s[static_cast<std::size_t>(a)]);
        for (Eigen::Index b = 0; b < n_s; ++b) {
          C_ss(a, b) = u.tangent_Pa(s[static_cast<std::size_t>(a)], s[static_cast<std::size_t>(b)]);
        }
      }
      if (n_s == 0 || r.cwiseAbs().maxCoeff() <= kStressTolerance_Pa) break;
      if (it == kMaxIterations) {
        throw std::runtime_error("MaterialPointDriver: step " + std::to_string(k) + " did not converge");
      }
      const Eigen::VectorXd d = C_ss.ldlt().solve(-r);
      for (Eigen::Index a = 0; a < n_s; ++a) eps(s[static_cast<std::size_t>(a)]) += d(a);
      u = model_.update(eps, state);
    }
    state = u.state;
    out.strain.push_back(eps);
    out.stress_Pa.push_back(u.stress_Pa);
    out.state.push_back(state);
    out.iterations.push_back(it);
  }
  return out;
}

std::vector<DriverStep> uniaxial_stress_steps(double sigma_end_Pa, int n_steps) {
  std::vector<DriverStep> steps;
  for (int k = 1; k <= n_steps; ++k) {
    steps.push_back({{S, S, S, S}, Voigt4(0.0, sigma_end_Pa * k / n_steps, 0.0, 0.0)});
  }
  return steps;
}

std::vector<DriverStep> uniaxial_strain_steps(const std::vector<double>& eps_yy) {
  std::vector<DriverStep> steps;
  for (const double e : eps_yy) steps.push_back({{S, E, S, S}, Voigt4(0.0, e, 0.0, 0.0)});
  return steps;
}

std::vector<DriverStep> plane_strain_bar_steps(const std::vector<double>& eps_xx) {
  std::vector<DriverStep> steps;
  for (const double e : eps_xx) steps.push_back({{E, S, E, S}, Voigt4(e, 0.0, 0.0, 0.0)});
  return steps;
}

}  // namespace tmfem
