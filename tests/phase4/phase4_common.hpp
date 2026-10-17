#pragma once

#include <filesystem>
#include <string>

#include <Eigen/Dense>

#include "mechanics/elasticity.hpp"
#include "support/paths.hpp"

// Shared data of the phase 4 tests (docs/learning/phase4.md section 10).
namespace tmfem::test {

constexpr int kP4Levels = 5;        // {disk,square}_t{3,6}_l0..l4, -clscale 1, 1/2, ..., 1/16
constexpr int kP4FitFirstLevel = 1;  // P4-g, P4-h: levels 1-4 enter the fit (REFERENCE_VALUES.md)
constexpr double kDisk_b_m = 0.010;  // cases/disk/disk.geo
constexpr double kSquare_L_m = 0.012;  // cases/heat2d/square.geo

// A steel-like material for the mechanical tests: the targets are ratios to E alpha dT.
inline ElasticMaterial thermal_test_material() { return {200e9, 0.3}; }
constexpr double kAlpha_per_K = 1.2e-5;
constexpr double kTref_K = 293.15;
constexpr double kdT_K = 100.0;

inline std::filesystem::path disk_mesh(int t_nodes, int level) {
  return source_dir() / "cases/disk/msh" / ("disk_t" + std::to_string(t_nodes) + "_l" + std::to_string(level) + ".msh");
}
inline std::filesystem::path square_mesh(int t_nodes, int level) {
  return source_dir() / "cases/heat2d/msh" /
         ("square_t" + std::to_string(t_nodes) + "_l" + std::to_string(level) + ".msh");
}

inline const char* plane_name(PlaneCondition p) { return p == PlaneCondition::PlaneStress ? "plane stress" : "plane strain"; }

}  // namespace tmfem::test
