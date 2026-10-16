#pragma once

#include <filesystem>
#include <vector>

#include "thermal/heat1d.hpp"

namespace tmfem {

// A 1D heat case read from cases/*.json ("case_type": "heat1d"). Temperatures are in K inside;
// the file gives them in degrees C and the loader converts (project convention: convert at input boundaries only).
struct Heat1DCase {
  Heat1DProblem problem;
  double T0_K;
  double theta;
  double dt_s;
  double t_end_s;
  MassType mass;
  std::vector<double> probes_m;
  std::filesystem::path output_csv;
};

// Relative paths inside the file (time series, output) are resolved against base_dir.
Heat1DCase load_heat1d_case(const std::filesystem::path& json_path, const std::filesystem::path& base_dir);

// Probe temperatures (K) at t = n dt for n = 0..N.
struct ProbeHistory {
  std::vector<double> t_s;
  std::vector<std::vector<double>> T_K;  // one vector per probe
};

ProbeHistory run_heat1d_case(const Heat1DCase& c);

}  // namespace tmfem
