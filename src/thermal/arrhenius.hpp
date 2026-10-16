#pragma once

namespace tmfem {

// Molar gas constant, J/(mol K) (CODATA 2018).
inline constexpr double kGasConstant_J_molK = 8.314462618;

// D(T) = D0 exp(-Q / (R T)), T in K.
double arrhenius_D(double D0_m2_s, double Q_J_mol, double T_K, double R_J_molK = kGasConstant_J_molK);

struct ArrheniusFit {
  double Q_J_mol;
  double D0_m2_s;
};

// Q = R ln(D2/D1) / (1/T1 - 1/T2), D0 = D1 exp(Q / (R T1)).
ArrheniusFit fit_arrhenius_two_points(double T1_K, double D1_m2_s, double T2_K, double D2_m2_s,
                                      double R_J_molK = kGasConstant_J_molK);

// Equal D t gives equal profiles: time at D2 equivalent to t1 at D1.
double equivalent_time(double t1, double D1_m2_s, double D2_m2_s);

// Depth at D2 with the same concentration as x1 at D1 after equal times: x1 sqrt(D2/D1).
double equal_time_depth(double x1, double D1_m2_s, double D2_m2_s);

}  // namespace tmfem
