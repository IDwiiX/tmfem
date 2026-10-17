#pragma once

#include <filesystem>

#include "mechanics/elastic2d.hpp"

namespace tmfem {

// A 2D linear-elastic case read from cases/*.json ("case_type": "elastic2d"). The file gives E in GPa
// and stresses in MPa; the loader converts to SI (project convention: convert at input boundaries only).
// Boundary condition types: "displacement" (component "x"/"y", value_m), "traction" (traction_MPa [tx, ty]),
// "kirsch_traction" (S_MPa, a_m, outward normal [nx, ny]: t = sigma_Kirsch(x) n, docs/learning/phase2.md
// section 9).
struct Elastic2DCase {
  Elastic2DProblem problem;
  std::filesystem::path output_vtu;
  std::filesystem::path output_nodal_csv;
};

// Relative paths inside the file (mesh, outputs) are resolved against base_dir.
Elastic2DCase load_elastic2d_case(const std::filesystem::path& json_path, const std::filesystem::path& base_dir);

// Solves and writes the VTU (u in m, stresses and von Mises in MPa, nodal averages) and a nodal CSV.
void run_elastic2d_case(const Elastic2DCase& c);

}  // namespace tmfem
