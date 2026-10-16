#pragma once

#include <functional>
#include <vector>

namespace tmfem {

// Thin plate with convection on both faces (docs/learning/phase1.md section 9).
struct ConvectivePlate {
  double L_m;
  double k_W_mK;
  double rho_c_J_m3K;
  double h_out_W_m2K;
  double h_in_W_m2K;
};

double biot_number(const ConvectivePlate& p);       // h_out L / k
double lumped_time_constant_s(const ConvectivePlate& p);  // rho c L / (h_out + h_in)
// Quasi-steady T(0) - T(L) per unit (T_out - T_in): (L/k) / (1/h_out + L/k + 1/h_in).
double quasi_steady_conduction_fraction(const ConvectivePlate& p);

// Lumped ODE right-hand side: dT/dt = [h_out (T_out(t) - T) + h_in (T_in(t) - T)] / (rho c L).
double lumped_plate_rate(const ConvectivePlate& p, double T_K, double T_out_K, double T_in_K);

// Classical RK4 for y' = f(t, y); returns y at t0 + n dt for n = 0..n_steps.
std::vector<double> rk4_trajectory(const std::function<double(double t_s, double y)>& f, double y0,
                                   double t0_s, double dt_s, int n_steps);

}  // namespace tmfem
