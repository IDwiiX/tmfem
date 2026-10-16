#include "thermal/arrhenius.hpp"

#include <cmath>

namespace tmfem {

double arrhenius_D(double D0_m2_s, double Q_J_mol, double T_K, double R_J_molK) {
  return D0_m2_s * std::exp(-Q_J_mol / (R_J_molK * T_K));
}

ArrheniusFit fit_arrhenius_two_points(double T1_K, double D1_m2_s, double T2_K, double D2_m2_s,
                                      double R_J_molK) {
  const double Q = R_J_molK * std::log(D2_m2_s / D1_m2_s) / (1.0 / T1_K - 1.0 / T2_K);
  return {Q, D1_m2_s * std::exp(Q / (R_J_molK * T1_K))};
}

double equivalent_time(double t1, double D1_m2_s, double D2_m2_s) {
  return t1 * D1_m2_s / D2_m2_s;
}

double equal_time_depth(double x1, double D1_m2_s, double D2_m2_s) {
  return x1 * std::sqrt(D2_m2_s / D1_m2_s);
}

}  // namespace tmfem
