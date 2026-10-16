#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstdio>

#include "io/json_file.hpp"
#include "support/paths.hpp"
#include "support/validation_record.hpp"
#include "thermal/arrhenius.hpp"
#include "thermal/diffusion_analytic.hpp"
#include "thermal/heat1d.hpp"
#include "thermal/theta_integrator.hpp"

// Diffusion exercises from references/data/diffusion_cases.json (docs/REFERENCE_VALUES.md, Phase 1).
// Inputs in cm are converted to SI here, at the input boundary; results are converted back only
// for comparison with the tabulated targets.

using namespace tmfem;

namespace {

constexpr double kCm2ToM2 = 1e-4;
constexpr double kPerCm3ToPerM3 = 1e6;
constexpr double kPerCm2ToPerM2 = 1e4;
constexpr double kCelsiusToKelvin = 273.15;

nlohmann::json cases() { return read_json(test::source_dir() / "references/data/diffusion_cases.json"); }

double rel_err(double value, double reference) { return std::abs(value - reference) / std::abs(reference); }

void check_rel(const char* id, const char* quantity, double value, double reference, double tol,
               const char* test_name) {
  INFO(id << " " << quantity << ": tmfem " << value << ", reference " << reference << ", rel error "
          << rel_err(value, reference));
  char tol_text[32];
  std::snprintf(tol_text, sizeof tol_text, "%g relative", tol);
  test::record(id, quantity, reference, value, tol_text, test_name);
  CHECK(rel_err(value, reference) <= tol);
}

}  // namespace

// P1-a: steady flux J = -D dC/dx between 5 and 10 mm; target 2.4e-9 kg/(m2 s), 1e-12 relative.
TEST_CASE("P1-a steady flux ex1", "[phase1][P1-a]") {
  constexpr double kTol = 1e-12;  // REFERENCE_VALUES P1-a
  const auto ex = cases().at("ex1_steady_flux");
  const double D = ex.at("D_m2_per_s").get<double>();
  const double C1 = ex.at("C_at_5mm_kg_per_m3").get<double>();
  const double C2 = ex.at("C_at_10mm_kg_per_m3").get<double>();
  const double x1 = 5e-3, x2 = 10e-3;
  const double target = 2.4e-9;

  const double J_formula = -D * (C2 - C1) / (x2 - x1);
  check_rel("P1-a", "J formula", J_formula, target, kTol, "P1-a steady flux ex1");

  // Steady FE solve on [x1, x2] (shifted to [0, x2 - x1]) with Dirichlet ends; diffusion: rho_c = 1, k = D.
  Heat1DProblem p{uniform_interval(x2 - x1, 10), D, 1.0, DirichletBC{[C1](double) { return C1; }},
                  DirichletBC{[C2](double) { return C2; }}, {}};
  const Eigen::VectorXd C = solve_steady(build_operators(p, MassType::Consistent), 0.0);
  for (double J : element_flux(p.mesh, D, C)) {
    CHECK(rel_err(J, target) <= kTol);
  }
  check_rel("P1-a", "J finite element", element_flux(p.mesh, D, C).front(), target, kTol,
            "P1-a steady flux ex1");
}

// P1-c: equivalent time at 500 C for 10 h at 600 C; target 110.4 h, 0.1%.
TEST_CASE("P1-c equivalent time at 500 C", "[phase1][P1-c]") {
  constexpr double kTol = 1e-3;  // REFERENCE_VALUES P1-c
  const auto ex = cases().at("ex3_equivalent_time");
  const double t600_s = ex.at("t_600C_h").get<double>() * 3600.0;
  const double t500_s = equivalent_time(t600_s, ex.at("D_600C_m2_per_s").get<double>(),
                                        ex.at("D_500C_m2_per_s").get<double>());
  check_rel("P1-c", "t at 500 C [h]", t500_s / 3600.0, 110.4, kTol, "P1-c equivalent time at 500 C");
}

// P1-d: D(1000 K), D(800 K) and the depth at 800 K; targets 2.067e-13, 2.140e-15 m2/s, 1.577 mm; 0.5%.
TEST_CASE("P1-d Arrhenius D(1000 K), D(800 K), x at 800 K", "[phase1][P1-d]") {
  constexpr double kTol = 5e-3;  // REFERENCE_VALUES P1-d
  const auto all = cases();
  const double R = all.at("R_J_per_molK").get<double>();
  const auto ex = all.at("ex4_diffusion_couple");
  const double D0 = ex.at("D0_m2_per_s").get<double>();
  const double Q = ex.at("Q_J_per_mol").get<double>();
  const double D1 = arrhenius_D(D0, Q, ex.at("T1_K").get<double>(), R);
  const double D2 = arrhenius_D(D0, Q, ex.at("T2_K").get<double>(), R);
  const double x2_mm = equal_time_depth(ex.at("x1_mm").get<double>(), D1, D2);
  const char* name = "P1-d Arrhenius D(1000 K), D(800 K), x at 800 K";
  check_rel("P1-d", "D(1000 K) [m2/s]", D1, 2.067e-13, kTol, name);
  check_rel("P1-d", "D(800 K) [m2/s]", D2, 2.140e-15, kTol, name);
  check_rel("P1-d", "x at 800 K [mm]", x2_mm, 1.577, kTol, name);
}

namespace {

// Ex. 5 input converted to SI.
struct BoronInSilicon {
  double Cs_per_m3;
  double t_s;
  double D1100_m2_s;
  double D1200_m2_s;
};

BoronInSilicon boron() {
  const auto ex = cases().at("ex5_boron_in_silicon");
  return {ex.at("Cs_at_per_cm3").get<double>() * kPerCm3ToPerM3, ex.at("t_s").get<double>(),
          ex.at("D_1100C_cm2_per_s").get<double>() * kCm2ToM2,
          ex.at("D_1200C_cm2_per_s").get<double>() * kCm2ToM2};
}

// FE predeposition: surface held at Cs, initially clean; returns the dose int C dx at t_end.
double fe_predeposition_dose(const BoronInSilicon& b, int n_elements, double ratio, double dt_s) {
  const double L = 2e-6;  // >= 10 sqrt(D t) = 1.2e-6 m
  Heat1DProblem p{graded_interval(L, n_elements, ratio), b.D1100_m2_s, 1.0,
                  DirichletBC{[Cs = b.Cs_per_m3](double) { return Cs; }}, FluxBC{[](double) { return 0.0; }},
                  {}};
  // Lumped mass: monotone for the step initial condition (docs/learning/phase1.md section 5).
  const ThetaIntegrator be(build_operators(p, MassType::Lumped), 1.0, dt_s);
  Eigen::VectorXd C0 = Eigen::VectorXd::Zero(p.mesh.n_nodes());
  C0(0) = b.Cs_per_m3;
  const int n_steps = static_cast<int>(std::lround(b.t_s / dt_s));
  const Eigen::VectorXd C = march(be, C0, 0.0, n_steps, {});
  return Eigen::VectorXd::Ones(C.size()).dot(be.operators().M * C);  // int C_h dx = 1^T M C
}

}  // namespace

// P1-e: predeposition dose at 1100 C, 6 min; target 4.062e15 at/cm2, 0.5%.
TEST_CASE("P1-e boron predeposition dose", "[phase1][P1-e]") {
  constexpr double kTol = 5e-3;  // REFERENCE_VALUES P1-e
  const BoronInSilicon b = boron();
  const double target = 4.062e15;
  const char* name = "P1-e boron predeposition dose";
  const double M = predeposition_dose(b.Cs_per_m3, b.D1100_m2_s, b.t_s) / kPerCm2ToPerM2;
  check_rel("P1-e", "dose formula [at/cm2]", M, target, kTol, name);

  // FE cross-check; resolution fixed a priori, then halved in h and dt to show it is converged.
  const double M_fe = fe_predeposition_dose(b, 400, 1.01, 0.05) / kPerCm2ToPerM2;
  const double M_fine = fe_predeposition_dose(b, 800, std::sqrt(1.01), 0.025) / kPerCm2ToPerM2;
  INFO("FE dose " << M_fe << ", refined " << M_fine << ", change " << rel_err(M_fe, M_fine));
  check_rel("P1-e", "dose finite element [at/cm2]", M_fe, target, kTol, name);
  CHECK(rel_err(M_fe, M_fine) < kTol / 5.0);
}

// P1-f: surface concentration of the thin film after drive-in at 1200 C; target 6.97e19 at/cm3, 0.5%.
TEST_CASE("P1-f thin-film surface concentration at 1200 C", "[phase1][P1-f]") {
  constexpr double kTol = 5e-3;  // REFERENCE_VALUES P1-f
  const BoronInSilicon b = boron();
  const double M = predeposition_dose(b.Cs_per_m3, b.D1100_m2_s, b.t_s);
  const double Cs = thin_film_profile(0.0, b.t_s, b.D1200_m2_s, M) / kPerCm3ToPerM3;
  check_rel("P1-f", "surface C at 1200 C [at/cm3]", Cs, 6.97e19, kTol,
            "P1-f thin-film surface concentration at 1200 C");
}

// P1-g: Arrhenius fit through D(1100 C), D(1200 C); targets 338.9 kJ/mol, 3.11 cm2/s, 1.74e-11 cm2/s; 0.5%.
TEST_CASE("P1-g Arrhenius fit Q, D0, D(1300 C)", "[phase1][P1-g]") {
  constexpr double kTol = 5e-3;  // REFERENCE_VALUES P1-g
  const double R = cases().at("R_J_per_molK").get<double>();
  const BoronInSilicon b = boron();
  const double T1 = 1100.0 + kCelsiusToKelvin, T2 = 1200.0 + kCelsiusToKelvin, T3 = 1300.0 + kCelsiusToKelvin;
  const ArrheniusFit fit = fit_arrhenius_two_points(T1, b.D1100_m2_s, T2, b.D1200_m2_s, R);
  const double D3 = arrhenius_D(fit.D0_m2_s, fit.Q_J_mol, T3, R);
  const char* name = "P1-g Arrhenius fit Q, D0, D(1300 C)";
  check_rel("P1-g", "Q [kJ/mol]", fit.Q_J_mol / 1e3, 338.9, kTol, name);
  check_rel("P1-g", "D0 [cm2/s]", fit.D0_m2_s / kCm2ToM2, 3.11, kTol, name);
  check_rel("P1-g", "D(1300 C) [cm2/s]", D3 / kCm2ToM2, 1.74e-11, kTol, name);
}
