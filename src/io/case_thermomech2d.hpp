#pragma once

#include <filesystem>
#include <string>
#include <vector>

#include "mechanics/elastic2d.hpp"

namespace tmfem {

// Demo D4, beam-heated block (docs/learning/phase4.md section 11; "case_type": "thermomech2d").
// A Gaussian beam of rms size s deposits E_l (J per m of beam path) during a pulse of length tau in the
// block cross-section; adiabatic pulse, then convective cooling of the outer edges. One-way coupling:
// each output time solves linear elasticity with the thermal strain. The file gives E in GPa and
// strengths in MPa; the loader converts to SI (project convention: convert at input boundaries only).
struct BeamMaterial {
  std::string name;
  double rho_kg_m3, c_J_kgK, k_W_mK, alpha_per_K, E_Pa, nu;
  bool ductile;  // true: von Mises against yield_Pa; false: principal stresses against the strengths below
  double yield_Pa = 0.0, flexural_Pa = 0.0, compressive_Pa = 0.0;
  std::vector<double> nu_sensitivity;  // extra nu values re-solved at the end of the pulse
};

struct BeamBlockCase {
  std::string mesh_pattern;  // "{s_mm}" is replaced by the beam size, formatted with %g
  std::vector<double> beam_sizes_mm;
  BeamMaterial material;
  PlaneCondition plane;
  double E_l_J_m, pulse_s;
  int pulse_steps;
  std::string cooled_physical;
  double h_W_m2K, T_inf_K, T0_K, T_ref_K;
  std::vector<double> cooling_stage_ends_s;  // after the pulse; each stage has cooling_steps BE steps
  int cooling_steps;
  double quarter_model_factor;  // 4: the mesh is a quarter of the section
  double section_area_m2;       // full section, for the generalised plane strain estimate
  std::vector<DisplacementBC> supports;
  std::filesystem::path output_dir;
};

// Relative paths inside the file (meshes, outputs) are resolved against base_dir.
BeamBlockCase load_thermomech2d_case(const std::filesystem::path& json_path, const std::filesystem::path& base_dir);

// One output time. Stresses are nodal-averaged (plots); the centre values come from stress_at_point.
struct BeamSnapshot {
  double t_s, T_max_K, dT_centre_K;
  double vm_max_Pa, s1_max_Pa, s3_min_Pa;  // s1, s3: largest and smallest principal stress, zz included
  double sxx_centre_Pa, szz_centre_Pa, vm_centre_Pa;
};

struct BeamRunResult {
  double s_mm;
  double energy_stored_J_m, energy_deposited_J_m, energy_nominal_J_m;  // whole section (quarter x 4)
  std::vector<BeamSnapshot> history;  // history[0] is the end of the pulse
  std::vector<std::pair<double, BeamSnapshot>> nu_sensitivity;  // (nu, end-of-pulse snapshot)
};

// Runs one beam size: writes <output_dir>/s<s>/end_of_pulse.vtu, end_of_pulse_nodal.csv and history.csv.
BeamRunResult run_beam_size(const BeamBlockCase& c, double s_mm);

// Runs the sweep, writes <output_dir>/sweep.csv and nu_sensitivity.csv, and prints the comparison with the
// closed-form estimates of section 11 (sanity checks, not validation).
std::vector<BeamRunResult> run_thermomech2d_case(const BeamBlockCase& c);

}  // namespace tmfem
