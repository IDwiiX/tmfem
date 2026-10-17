#include "io/case_thermomech2d.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>

#include "io/case_elastic2d.hpp"
#include "io/csv_writer.hpp"
#include "io/json_file.hpp"
#include "io/vtu_writer.hpp"
#include "mechanics/stress_recovery.hpp"
#include "mesh/gmsh_reader.hpp"
#include "thermal/heat2d.hpp"
#include "thermal/theta_integrator.hpp"

namespace tmfem {

namespace {

constexpr double kPaPerMPa = 1e6;
constexpr double kPaPerGPa = 1e9;
constexpr double kPi = 3.14159265358979323846;

std::string format_mm(double s_mm) {
  char buf[32];
  std::snprintf(buf, sizeof buf, "%g", s_mm);
  return buf;
}

std::filesystem::path mesh_path(const BeamBlockCase& c, double s_mm) {
  std::string p = c.mesh_pattern;
  const auto at = p.find("{s_mm}");
  if (at == std::string::npos) throw std::runtime_error("thermomech2d: mesh_pattern must contain {s_mm}");
  return p.replace(at, 6, format_mm(s_mm));
}

BeamMaterial load_material(const nlohmann::json& j) {
  BeamMaterial m{j.at("name").get<std::string>(),
                 j.at("rho_kg_m3").get<double>(),
                 j.at("c_J_kgK").get<double>(),
                 j.at("k_W_mK").get<double>(),
                 j.at("alpha_per_K").get<double>(),
                 j.at("E_GPa").get<double>() * kPaPerGPa,
                 j.at("nu").get<double>(),
                 j.at("strength").at("type").get<std::string>() == "ductile",
                 0.0,
                 0.0,
                 0.0,
                 {}};
  const auto& s = j.at("strength");
  if (m.ductile) {
    m.yield_Pa = s.at("yield_MPa").get<double>() * kPaPerMPa;
  } else {
    m.flexural_Pa = s.at("flexural_MPa").get<double>() * kPaPerMPa;
    m.compressive_Pa = s.at("compressive_MPa").get<double>() * kPaPerMPa;
  }
  m.nu_sensitivity = j.value("nu_sensitivity", std::vector<double>{});
  return m;
}

// Volumetric source of the Gaussian beam during the pulse (section 11), W/m3.
SpaceTimeField beam_source(const BeamBlockCase& c, double s_m) {
  const double q0 = c.E_l_J_m / (2 * kPi * s_m * s_m * c.pulse_s);
  return [q0, s_m](const Eigen::Vector2d& x, double) { return q0 * std::exp(-x.squaredNorm() / (2 * s_m * s_m)); };
}

Heat2DProblem thermal_problem(const BeamBlockCase& c, const Mesh& mesh) {
  Heat2DProblem p;
  p.mesh = mesh;
  p.k_W_mK = c.material.k_W_mK;
  p.rho_c_J_m3K = c.material.rho_kg_m3 * c.material.c_J_kgK;
  return p;
}

Elastic2DProblem mechanical_problem(const BeamBlockCase& c, const Mesh& mesh, const Eigen::VectorXd& T, double nu) {
  return {mesh, {c.material.E_Pa, nu}, c.plane, "domain", c.supports, {},
          ThermalStrain{c.material.alpha_per_K, c.T_ref_K, T}};
}

// In-plane principal stresses combined with sigma_zz (a principal stress: no xz, yz shear).
std::pair<double, double> principal_extremes(const Eigen::Vector3d& s, double szz) {
  const double mid = 0.5 * (s(0) + s(1)), rad = std::hypot(0.5 * (s(0) - s(1)), s(2));
  return {std::max(mid + rad, szz), std::min(mid - rad, szz)};
}

struct Fields {
  BeamSnapshot snap;
  NodalStressFields stress;
  Eigen::MatrixXd s1s3;  // n x 2
  Eigen::VectorXd u_m;
};

Fields solve_snapshot(const BeamBlockCase& c, const Mesh& mesh, const Eigen::VectorXd& T, double t_s, double nu,
                      int centre) {
  const Elastic2DProblem p = mechanical_problem(c, mesh, T, nu);
  const Elastic2DSolution sol = solve_elastic2d(p);
  Fields f{{}, nodal_stress_fields(p, sol.u_m), Eigen::MatrixXd(mesh.n_nodes(), 2), sol.u_m};
  BeamSnapshot& s = f.snap;
  s.t_s = t_s;
  s.T_max_K = T.maxCoeff();
  s.dT_centre_K = T(centre) - c.T0_K;
  s.vm_max_Pa = f.stress.von_mises_Pa.maxCoeff();
  for (int i = 0; i < mesh.n_nodes(); ++i) {
    const auto [s1, s3] = principal_extremes(f.stress.sigma_Pa.row(i).transpose(), f.stress.sigma_zz_Pa(i));
    f.s1s3.row(i) << s1, s3;
  }
  s.s1_max_Pa = f.s1s3.col(0).maxCoeff();
  s.s3_min_Pa = f.s1s3.col(1).minCoeff();
  const auto sc = stress_at_point(p, sol.u_m, mesh.x_m.row(centre).transpose());
  if (!sc) throw std::runtime_error("thermomech2d: centre outside the mesh");
  const double szz = sigma_zz_thermal(p.material, p.plane, *sc, c.material.alpha_per_K * (T(centre) - c.T_ref_K));
  s.sxx_centre_Pa = (*sc)(0);
  s.szz_centre_Pa = szz;
  s.vm_centre_Pa = von_mises(*sc, szz);
  return f;
}

std::vector<double> col(const Eigen::MatrixXd& m, Eigen::Index j, double scale = 1.0) {
  std::vector<double> v(static_cast<std::size_t>(m.rows()));
  for (Eigen::Index i = 0; i < m.rows(); ++i) v[static_cast<std::size_t>(i)] = m(i, j) * scale;
  return v;
}

void write_end_of_pulse(const std::filesystem::path& dir, const Mesh& mesh, const Eigen::VectorXd& T, const Fields& f) {
  const int n = mesh.n_nodes();
  Eigen::MatrixXd u(n, 3);
  for (int i = 0; i < n; ++i) u.row(i) << f.u_m(2 * i), f.u_m(2 * i + 1), 0.0;
  const Eigen::MatrixXd s_MPa = f.stress.sigma_Pa / kPaPerMPa, s13_MPa = f.s1s3 / kPaPerMPa;
  const Eigen::MatrixXd szz_MPa = f.stress.sigma_zz_Pa / kPaPerMPa, vm_MPa = f.stress.von_mises_Pa / kPaPerMPa;
  write_vtu(dir / "end_of_pulse.vtu", mesh, mesh.blocks_of("domain"),
            {{"temperature_K", T}, {"displacement_m", u}, {"stress_xx_yy_xy_MPa", s_MPa}, {"stress_zz_MPa", szz_MPa},
             {"von_mises_MPa", vm_MPa}, {"principal_max_min_MPa", s13_MPa}});
  write_csv(dir / "end_of_pulse_nodal.csv",
            {"x_m", "y_m", "T_K", "sxx_MPa", "syy_MPa", "sxy_MPa", "szz_MPa", "von_mises_MPa", "s1_MPa", "s3_MPa"},
            {col(mesh.x_m, 0), col(mesh.x_m, 1), col(T, 0), col(s_MPa, 0), col(s_MPa, 1), col(s_MPa, 2), col(szz_MPa, 0),
             col(vm_MPa, 0), col(s13_MPa, 0), col(s13_MPa, 1)},
            {"demo D4 end of pulse, quarter model, nodal-averaged stresses"});
}

void write_history(const std::filesystem::path& path, const std::vector<BeamSnapshot>& h) {
  std::vector<std::vector<double>> cols(9);
  for (const BeamSnapshot& s : h) {
    const double v[] = {s.t_s, s.T_max_K, s.dT_centre_K, s.vm_max_Pa / kPaPerMPa, s.s1_max_Pa / kPaPerMPa,
                        s.s3_min_Pa / kPaPerMPa, s.sxx_centre_Pa / kPaPerMPa, s.szz_centre_Pa / kPaPerMPa,
                        s.vm_centre_Pa / kPaPerMPa};
    for (std::size_t j = 0; j < cols.size(); ++j) cols[j].push_back(v[j]);
  }
  write_csv(path, {"t_s", "T_max_K", "dT_centre_K", "vm_max_MPa", "s1_max_MPa", "s3_min_MPa", "sxx_centre_MPa",
                   "szz_centre_MPa", "vm_centre_MPa"},
            cols, {"demo D4 history; first row = end of pulse"});
}

}  // namespace

BeamBlockCase load_thermomech2d_case(const std::filesystem::path& json_path, const std::filesystem::path& base_dir) {
  const nlohmann::json j = read_json(json_path);
  if (j.at("case_type").get<std::string>() != "thermomech2d") {
    throw std::runtime_error("load_thermomech2d_case: case_type is not thermomech2d in " + json_path.string());
  }
  const auto& beam = j.at("beam");
  const auto& cool = j.at("cooling");
  BeamBlockCase c{(base_dir / j.at("mesh_pattern").get<std::string>()).string(),
                  j.at("beam_sizes_mm").get<std::vector<double>>(),
                  load_material(j.at("material")),
                  plane_condition(j.at("plane").get<std::string>()),
                  beam.at("E_per_length_J_m").get<double>(),
                  beam.at("pulse_s").get<double>(),
                  beam.at("pulse_steps").get<int>(),
                  cool.at("physical").get<std::string>(),
                  cool.at("h_W_m2K").get<double>(),
                  cool.at("T_inf_K").get<double>(),
                  j.at("T0_K").get<double>(),
                  j.at("T_ref_K").get<double>(),
                  cool.at("stage_ends_s").get<std::vector<double>>(),
                  cool.at("steps_per_stage").get<int>(),
                  j.at("quarter_model_factor").get<double>(),
                  j.at("section_area_m2").get<double>(),
                  {},
                  base_dir / j.at("output_dir").get<std::string>()};
  for (const auto& s : j.at("supports")) {
    const auto comp = s.at("component").get<std::string>();
    if (comp != "x" && comp != "y") throw std::runtime_error("thermomech2d: support component must be x or y");
    c.supports.push_back({s.at("physical").get<std::string>(), comp == "x" ? 0 : 1, [](const Eigen::Vector2d&) { return 0.0; }});
  }
  return c;
}

BeamRunResult run_beam_size(const BeamBlockCase& c, double s_mm) {
  const double s_m = s_mm * 1e-3;
  const Mesh mesh = read_gmsh41(mesh_path(c, s_mm));
  const std::vector<int> centre_nodes = mesh.nodes_of("centre");
  if (centre_nodes.size() != 1) throw std::runtime_error("thermomech2d: physical point 'centre' must hold one node");
  const int centre = centre_nodes.front();
  const auto dir = c.output_dir / ("s" + format_mm(s_mm));

  // Pulse stage (section 5): adiabatic, source on at every step end, BE.
  Heat2DProblem pulse = thermal_problem(c, mesh);
  pulse.source_W_m3 = beam_source(c, s_m);
  const ThermalOperators pulse_ops = build_operators(pulse);
  const Eigen::VectorXd T0 = Eigen::VectorXd::Constant(mesh.n_nodes(), c.T0_K);
  const double dt_pulse = c.pulse_s / c.pulse_steps;
  Eigen::VectorXd T = march(ThetaIntegrator(pulse_ops, 1.0, dt_pulse), T0, 0.0, c.pulse_steps, nullptr);

  BeamRunResult r{s_mm, 0, 0, c.E_l_J_m, {}, {}};
  const Eigen::VectorXd one = Eigen::VectorXd::Ones(mesh.n_nodes());
  r.energy_stored_J_m = c.quarter_model_factor * one.dot(pulse_ops.M * (T - T0));
  r.energy_deposited_J_m = c.quarter_model_factor * c.pulse_s * one.dot(pulse_ops.F(c.pulse_s));  // BE: sum dt F^{n+1}

  const Fields end = solve_snapshot(c, mesh, T, c.pulse_s, c.material.nu, centre);
  write_end_of_pulse(dir, mesh, T, end);
  r.history.push_back(end.snap);
  for (const double nu : c.material.nu_sensitivity) {
    r.nu_sensitivity.emplace_back(nu, solve_snapshot(c, mesh, T, c.pulse_s, nu, centre).snap);
  }

  // Cooling stages: source off, Robin on the cooled edges, BE with a fixed step per stage.
  Heat2DProblem cooling = thermal_problem(c, mesh);
  const double T_inf = c.T_inf_K;
  cooling.convection.push_back({c.cooled_physical, c.h_W_m2K, [T_inf](double) { return T_inf; }});
  double t = c.pulse_s;
  ThetaIntegrator integ(build_operators(cooling), 1.0, (c.cooling_stage_ends_s.front() - t) / c.cooling_steps);
  for (const double t_end : c.cooling_stage_ends_s) {
    integ.set_dt((t_end - t) / c.cooling_steps);
    T = march(integ, T, t, c.cooling_steps, nullptr);
    t = t_end;
    r.history.push_back(solve_snapshot(c, mesh, T, t, c.material.nu, centre).snap);
  }
  write_history(dir / "history.csv", r.history);
  return r;
}

std::vector<BeamRunResult> run_thermomech2d_case(const BeamBlockCase& c) {
  const BeamMaterial& m = c.material;
  const double rho_c = m.rho_kg_m3 * m.c_J_kgK;
  std::printf("Demo D4 (not a validation): %s, %s, E_l = %g J/m, tau = %g s\n", m.name.c_str(),
              c.plane == PlaneCondition::PlaneStrain ? "plane strain" : "plane stress", c.E_l_J_m, c.pulse_s);
  std::vector<BeamRunResult> results;
  std::vector<std::vector<double>> sweep(16), sens(6);
  for (const double s_mm : c.beam_sizes_mm) {
    const BeamRunResult r = run_beam_size(c, s_mm);
    const BeamSnapshot& e = r.history.front();
    const double s_m = s_mm * 1e-3;
    const double dT0 = c.E_l_J_m / (2 * kPi * s_m * s_m * rho_c);
    const double vm_est = m.E_Pa * m.alpha_per_K * dT0 / (2 * (1 - m.nu));  // infinite plane strain, centre
    const double szz_est = -m.E_Pa * m.alpha_per_K * dT0 / (1 - m.nu);
    const double gps = (1 - m.nu) * 2 * kPi * s_m * s_m / c.section_area_m2;  // generalised plane strain shift
    double peak_vm = 0, peak_s1 = 0, peak_s3 = 0;
    for (const BeamSnapshot& h : r.history) {
      peak_vm = std::max(peak_vm, h.vm_max_Pa);
      peak_s1 = std::max(peak_s1, h.s1_max_Pa);
      peak_s3 = std::min(peak_s3, h.s3_min_Pa);
    }
    const double ratio = m.ductile ? peak_vm / m.yield_Pa : peak_s1 / m.flexural_Pa;
    const double ratio_c = m.ductile ? 0.0 : -peak_s3 / m.compressive_Pa;
    std::printf("  s = %g mm: dT_centre %.2f K (estimate %.2f); centre von Mises %.2f MPa (estimate %.2f), "
                "sigma_zz %.2f MPa (estimate %.2f)\n",
                s_mm, e.dT_centre_K, dT0, e.vm_centre_Pa / kPaPerMPa, vm_est / kPaPerMPa, e.szz_centre_Pa / kPaPerMPa,
                szz_est / kPaPerMPa);
    std::printf("    peak over time: von Mises %.2f MPa, s1 %.2f MPa, s3 %.2f MPa; T_max %.1f K; ", peak_vm / kPaPerMPa,
                peak_s1 / kPaPerMPa, peak_s3 / kPaPerMPa, e.T_max_K);
    if (m.ductile) {
      std::printf("von Mises / yield %.2f%s\n", ratio, ratio > 1 ? " (above yield: elastic result is only an indicator)" : "");
    } else {
      std::printf("s1 / flexural %.3f (optimistic: flexural > tensile), |s3| / compressive %.3f\n", ratio, ratio_c);
    }
    std::printf("    energy: stored %.6f J/m, deposited (discrete) %.6f, nominal %.6f; generalised plane strain "
                "shift %.2f%% of |sigma_zz|\n",
                r.energy_stored_J_m, r.energy_deposited_J_m, r.energy_nominal_J_m, 100 * gps);
    const double row[] = {s_mm, dT0, e.dT_centre_K, vm_est / kPaPerMPa, e.vm_centre_Pa / kPaPerMPa,
                          szz_est / kPaPerMPa, e.szz_centre_Pa / kPaPerMPa, peak_vm / kPaPerMPa, peak_s1 / kPaPerMPa,
                          peak_s3 / kPaPerMPa, e.T_max_K, ratio, ratio_c, 100 * gps, r.energy_stored_J_m,
                          r.energy_deposited_J_m};
    for (std::size_t j = 0; j < sweep.size(); ++j) sweep[j].push_back(row[j]);
    for (const auto& [nu, sn] : r.nu_sensitivity) {
      const double srow[] = {s_mm, nu, sn.vm_max_Pa / kPaPerMPa, sn.s1_max_Pa / kPaPerMPa, sn.s3_min_Pa / kPaPerMPa,
                             sn.szz_centre_Pa / kPaPerMPa};
      for (std::size_t j = 0; j < sens.size(); ++j) sens[j].push_back(srow[j]);
      std::printf("    nu = %.2f: von Mises %.2f MPa (nominal nu %.2f: %.2f), s1 %.2f MPa\n", nu,
                  sn.vm_max_Pa / kPaPerMPa, m.nu, e.vm_max_Pa / kPaPerMPa, sn.s1_max_Pa / kPaPerMPa);
    }
    results.push_back(r);
  }
  write_csv(c.output_dir / "sweep.csv",
            {"s_mm", "dT0_estimate_K", "dT_centre_K", "vm_centre_estimate_MPa", "vm_centre_MPa", "szz_centre_estimate_MPa",
             "szz_centre_MPa", "vm_peak_MPa", "s1_peak_MPa", "s3_peak_MPa", "T_max_end_of_pulse_K", "strength_ratio",
             "compressive_ratio", "gps_shift_percent", "energy_stored_J_m", "energy_deposited_J_m"},
            sweep, {"demo D4 " + m.name + " (not a validation); estimates: infinite plane strain, section 11"});
  if (!m.nu_sensitivity.empty()) {
    write_csv(c.output_dir / "nu_sensitivity.csv",
              {"s_mm", "nu", "vm_max_MPa", "s1_max_MPa", "s3_min_MPa", "szz_centre_MPa"}, sens,
              {"demo D4 " + m.name + ", end of pulse"});
  }
  return results;
}

}  // namespace tmfem
