#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>

#include "io/case_heat1d.hpp"
#include "io/csv_writer.hpp"
#include "io/json_file.hpp"
#include "io/timeseries.hpp"
#include "support/paths.hpp"
#include "support/validation_record.hpp"
#include "thermal/lumped_plate.hpp"

// Steel plate 5 mm, convective both faces (V1.4). Parameters: references/data/steel_plate_case.json.
// Case: cases/steel_plate/synthetic.json (BE, dt = 300 s, synthetic sol-air signal).

using namespace tmfem;

namespace {

constexpr double kCelsiusToKelvin = 273.15;

nlohmann::json reference_case() { return read_json(test::source_dir() / "references/data/steel_plate_case.json"); }

ConvectivePlate plate_from_reference() {
  const auto j = reference_case();
  const auto& m = j.at("material");
  return {j.at("geometry").at("thickness_m").get<double>(), m.at("k_W_per_mK").get<double>(),
          m.at("rho_kg_per_m3").get<double>() * m.at("cp_J_per_kgK").get<double>(),
          j.at("outer_face_x0").at("h_W_per_m2K").get<double>(), j.at("inner_face_xL").at("h_W_per_m2K").get<double>()};
}

Heat1DCase synthetic_case() {
  return load_heat1d_case(test::source_dir() / "cases/steel_plate/synthetic.json", test::source_dir());
}

const TimeFunction& outer_air(const Heat1DCase& c) { return std::get<ConvectiveBC>(c.problem.left).T_inf_K; }

double inner_air_K() {
  return reference_case().at("inner_face_xL").at("T_inf_C").get<double>() + kCelsiusToKelvin;
}

// Lumped ODE by RK4 at 1 s, sampled at the FE output times (multiples of the FE dt).
std::vector<double> lumped_at(const Heat1DCase& c, const std::vector<double>& t_s) {
  const ConvectivePlate p = plate_from_reference();
  const double T_in = inner_air_K();
  const TimeFunction T_out = outer_air(c);
  constexpr double kDt_rk4 = 1.0;
  const int n = static_cast<int>(std::lround(t_s.back() / kDt_rk4));
  const auto y = rk4_trajectory([&](double t, double T) { return lumped_plate_rate(p, T, T_out(t), T_in); }, c.T0_K,
                                0.0, kDt_rk4, n);
  std::vector<double> out;
  for (double t : t_s) out.push_back(y[static_cast<std::size_t>(std::lround(t / kDt_rk4))]);
  return out;
}

double max_abs_diff(const std::vector<double>& a, const std::vector<double>& b) {
  double m = 0.0;
  for (std::size_t i = 0; i < a.size(); ++i) m = std::max(m, std::abs(a[i] - b[i]));
  return m;
}

}  // namespace

// P1-i: Bi = h_e L / k = 0.00183 and tau = rho c L / (h_e + h_i) = 549.5 s ("exact formulas"):
// checked to half a unit in the last quoted digit.
TEST_CASE("P1-i steel plate Biot number and time constant", "[phase1][P1-i]") {
  const ConvectivePlate p = plate_from_reference();
  const double Bi = biot_number(p);
  const double tau = lumped_time_constant_s(p);
  std::printf("P1-i Bi = %.6f, tau = %.3f s\n", Bi, tau);
  test::record("P1-i", "Biot number", 0.00183, Bi, "+/- 5e-6 (last quoted digit)", "P1-i steel plate Biot number and time constant");
  test::record("P1-i", "time constant [s]", 549.5, tau, "+/- 0.05 s (last quoted digit)", "P1-i steel plate Biot number and time constant");
  CHECK(std::abs(Bi - 0.00183) <= 0.5e-5);
  CHECK(std::abs(tau - 549.5) <= 0.05);

  // The case file must carry the same parameters as the reference data.
  const auto ref = reference_case();
  const Heat1DCase c = synthetic_case();
  CHECK(c.problem.mesh.length_m() == p.L_m);
  CHECK(c.problem.k == p.k_W_mK);
  CHECK(c.problem.rho_c == p.rho_c_J_m3K);
  CHECK(std::get<ConvectiveBC>(c.problem.left).h_W_m2K == p.h_out_W_m2K);
  CHECK(std::get<ConvectiveBC>(c.problem.right).h_W_m2K == p.h_in_W_m2K);
  CHECK(std::get<ConvectiveBC>(c.problem.right).T_inf_K(0.0) == inner_air_K());
  CHECK(c.T0_K == ref.at("initial_T_C").get<double>() + kCelsiusToKelvin);
  CHECK(c.dt_s == ref.at("dt_s").get<double>());
  CHECK(c.t_end_s == ref.at("duration_h").get<double>() * 3600.0);
}

// P1-j: 1D FE inner face vs lumped ODE (RK4, 1 s), max 0.05 K.
TEST_CASE("P1-j steel plate 1D FE vs lumped ODE", "[phase1][P1-j]") {
  constexpr double kTol_K = 0.05;  // REFERENCE_VALUES P1-j
  const Heat1DCase c = synthetic_case();
  REQUIRE(c.theta == 1.0);  // SPEC V1.4: backward Euler at dt = 300 s (CN rings there)
  const ProbeHistory h = run_heat1d_case(c);
  const std::vector<double> lumped = lumped_at(c, h.t_s);
  const double diff = max_abs_diff(h.T_K[1], lumped);

  // Convergence check: halving the prescribed dt changes the inner face by < 1/5 of the tolerance.
  Heat1DCase half = synthetic_case();
  half.dt_s /= 2.0;
  const ProbeHistory hh = run_heat1d_case(half);
  std::vector<double> hh_inner_at_coarse;
  for (std::size_t i = 0; i < hh.t_s.size(); i += 2) hh_inner_at_coarse.push_back(hh.T_K[1][i]);
  const double dt_sensitivity = max_abs_diff(h.T_K[1], hh_inner_at_coarse);
  const double diff_half = max_abs_diff(hh.T_K[1], lumped_at(half, hh.t_s));

  const auto peak = std::max_element(h.T_K[1].begin(), h.T_K[1].end());
  std::printf("P1-j max |T_inner FE - T_lumped| = %.4f K (dt = 150 s: %.4f K); inner peak %.3f C; "
              "FE change dt 300 -> 150 s: %.4f K\n",
              diff, diff_half, *peak - kCelsiusToKelvin, dt_sensitivity);
  test::record("P1-j", "max abs inner face FE - lumped [K]", 0.0, diff, "0.05 K", "P1-j steel plate 1D FE vs lumped ODE");
  CHECK(diff <= kTol_K);
  CHECK(dt_sensitivity < kTol_K / 5.0);

  std::vector<double> T_sa;
  for (double t : h.t_s) T_sa.push_back(outer_air(c)(t));
  write_csv(test::results_dir() / "phase1/steel_plate_synthetic.csv",
            {"t_s", "T_sa_K", "T_outer_fe_K", "T_inner_fe_K", "T_lumped_K"}, {h.t_s, T_sa, h.T_K[0], h.T_K[1], lumped},
            {"P1-j steel plate, synthetic T_sa, BE dt = 300 s, 20 elements; lumped RK4 dt = 1 s"});
}

// P1-k: max through-thickness difference vs the quasi-steady bound (L/k)/R_tot * max|T_sa - T_in|
// = 0.0255 K, 5% (row approved and applied at phase 1 start; the old 0.01 K is infeasible).
TEST_CASE("P1-k steel plate through-thickness difference", "[phase1][P1-k]") {
  constexpr double kTol = 0.05;  // REFERENCE_VALUES P1-k
  const Heat1DCase c = synthetic_case();
  const ProbeHistory h = run_heat1d_case(c);
  const double dT_fe = max_abs_diff(h.T_K[0], h.T_K[1]);

  const TimeSeries tsa = TimeSeries::read_csv(test::source_dir() / "cases/steel_plate/synthetic_tsa.csv");
  const double T_in_C = inner_air_K() - kCelsiusToKelvin;
  double swing = 0.0;
  for (double v : tsa.values()) swing = std::max(swing, std::abs(v - T_in_C));
  const double bound = quasi_steady_conduction_fraction(plate_from_reference()) * swing;
  std::printf("P1-k max |T(0) - T(L)|: FE %.5f K, quasi-steady %.5f K (rel %.3e)\n", dT_fe, bound,
              (dT_fe - bound) / bound);
  test::record("P1-k", "max through-thickness difference [K]", bound, dT_fe, "5%",
               "P1-k steel plate through-thickness difference");
  CHECK(std::abs(bound - 0.0255) <= 0.00005);  // the row's quoted value
  CHECK(std::abs(dT_fe - bound) / bound <= kTol);
}

// P1-l: real ttab vs the ANSYS inner-face curve, report only (no tuning). Skipped until the files exist.
TEST_CASE("P1-l steel plate vs ANSYS inner face", "[phase1][P1-l]") {
  const auto ttab = test::source_dir() / "references/data/ttab.csv";
  const auto ansys = test::source_dir() / "references/data/ansys_inner_face.csv";
  if (!std::filesystem::exists(ttab) || !std::filesystem::exists(ansys)) {
    SKIP("references/data/ttab.csv or ansys_inner_face.csv not provided yet");
  }
  Heat1DCase c = synthetic_case();
  const auto series = std::make_shared<TimeSeries>(TimeSeries::read_csv(ttab));
  std::get<ConvectiveBC>(c.problem.left).T_inf_K = [series](double t) { return (*series)(t) + kCelsiusToKelvin; };
  c.t_end_s = std::floor(series->t_last() / c.dt_s) * c.dt_s;
  const ProbeHistory h = run_heat1d_case(c);
  const TimeSeries ref = TimeSeries::read_csv(ansys);
  const std::vector<double> inner(h.T_K[1]);
  const TimeSeries fe(h.t_s, inner);
  double sum2 = 0.0;
  std::size_t n = 0;
  for (std::size_t i = 0; i < ref.times().size(); ++i) {
    const double t = ref.times()[i];
    if (t < fe.t_first() || t > fe.t_last()) continue;
    const double d = fe(t) - (ref.values()[i] + kCelsiusToKelvin);
    sum2 += d * d;
    ++n;
  }
  REQUIRE(n > 0);
  const double rms = std::sqrt(sum2 / static_cast<double>(n));
  std::printf("P1-l RMS(FE - ANSYS) inner face = %.3f K over %zu points\n", rms, n);
  test::record("P1-l", "RMS FE - ANSYS inner face [K]", 0.0, rms, "report only", "P1-l steel plate vs ANSYS inner face");
  const double source_peaks_C[3] = {51.0, 53.0, 60.0};  // "about 51, 53, 60 C" in the source report
  for (int day = 0; day * 86400.0 < h.t_s.back(); ++day) {
    double peak = -1e300;
    for (std::size_t i = 0; i < h.t_s.size(); ++i) {
      if (h.t_s[i] >= day * 86400.0 && h.t_s[i] < (day + 1) * 86400.0) peak = std::max(peak, inner[i]);
    }
    std::printf("P1-l day %d inner-face peak %.2f C (source reports about 51, 53, 60 C)\n", day + 1,
                peak - kCelsiusToKelvin);
    if (day < 3) {
      test::record("P1-l", "inner-face peak day " + std::to_string(day + 1) + " [C]", source_peaks_C[day],
                   peak - kCelsiusToKelvin, "report only", "P1-l steel plate vs ANSYS inner face");
    }
  }
}
