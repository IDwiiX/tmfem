#pragma once

#include <filesystem>

#include "mechanics/elastoplastic2d.hpp"

namespace tmfem {

// A 2D elastoplastic case read from cases/*.json ("case_type": "elastoplastic2d"; docs/learning/phase3.md).
// Fields as in elastic2d (mesh, material E_GPa / nu, plane, region, boundary_conditions at load factor 1), plus
// "hardening": {"type": "tabulated", "csv": true-curve file (eps_true, sigma_true_MPa), "sigma_start_MPa"} or
// {"type": "ludwik", "sigma_y0_MPa", "K_MPa", "m"}; "load_steps": equal steps up to load factor 1; and the
// outputs "output_vtu" and "output_gauss_csv". Converted to SI on input.
struct Elastoplastic2DCase {
  Elastoplastic2DProblem problem;
  int load_steps;
  std::filesystem::path output_vtu;
  std::filesystem::path output_gauss_csv;
};

// Relative paths inside the file (mesh, hardening csv, outputs) are resolved against base_dir.
Elastoplastic2DCase load_elastoplastic2d_case(const std::filesystem::path& json_path,
                                              const std::filesystem::path& base_dir);

// Runs every load step, then writes the VTU (u in m; cell-mean equivalent plastic strain and von Mises stress
// averaged to the nodes, plots only) and the quadrature-point CSV (stresses in MPa).
void run_elastoplastic2d_case(const Elastoplastic2DCase& c);

}  // namespace tmfem
