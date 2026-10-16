#include "thermal/lumped_plate.hpp"

namespace tmfem {

double biot_number(const ConvectivePlate& p) {
  return p.h_out_W_m2K * p.L_m / p.k_W_mK;
}

double lumped_time_constant_s(const ConvectivePlate& p) {
  return p.rho_c_J_m3K * p.L_m / (p.h_out_W_m2K + p.h_in_W_m2K);
}

double quasi_steady_conduction_fraction(const ConvectivePlate& p) {
  const double R_cond = p.L_m / p.k_W_mK;
  return R_cond / (1.0 / p.h_out_W_m2K + R_cond + 1.0 / p.h_in_W_m2K);
}

double lumped_plate_rate(const ConvectivePlate& p, double T_K, double T_out_K, double T_in_K) {
  return (p.h_out_W_m2K * (T_out_K - T_K) + p.h_in_W_m2K * (T_in_K - T_K)) / (p.rho_c_J_m3K * p.L_m);
}

std::vector<double> rk4_trajectory(const std::function<double(double t_s, double y)>& f, double y0,
                                   double t0_s, double dt_s, int n_steps) {
  std::vector<double> y(static_cast<std::size_t>(n_steps) + 1);
  y[0] = y0;
  for (int n = 0; n < n_steps; ++n) {
    const double t = t0_s + n * dt_s;
    const double yn = y[static_cast<std::size_t>(n)];
    const double k1 = f(t, yn);
    const double k2 = f(t + 0.5 * dt_s, yn + 0.5 * dt_s * k1);
    const double k3 = f(t + 0.5 * dt_s, yn + 0.5 * dt_s * k2);
    const double k4 = f(t + dt_s, yn + dt_s * k3);
    y[static_cast<std::size_t>(n) + 1] = yn + dt_s / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
  }
  return y;
}

}  // namespace tmfem
