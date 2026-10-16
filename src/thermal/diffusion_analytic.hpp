#pragma once

namespace tmfem {

// Closed-form diffusion solutions (docs/learning/phase1.md section 8). SI units.

// Inverse error function for |y| < 1, to machine precision.
double erfinv(double y);

// Semi-infinite body, surface held at Cs, initial C0: C = Cs - (Cs - C0) erf(x / (2 sqrt(D t))).
double erf_profile(double x_m, double t_s, double D_m2_s, double Cs, double C0);

// Time for the erf profile to reach C at depth x: t = x^2 / (4 D z^2), erf(z) = (Cs - C)/(Cs - C0).
double time_to_reach(double C, double x_m, double D_m2_s, double Cs, double C0);

// Dose through the surface for a constant surface concentration: M = 2 Cs sqrt(D t / pi).
double predeposition_dose(double Cs, double D_m2_s, double t_s);

// Thin film of dose M on a zero-flux half-space: C = M / sqrt(pi D t) exp(-x^2 / (4 D t)).
double thin_film_profile(double x_m, double t_s, double D_m2_s, double M);

}  // namespace tmfem
