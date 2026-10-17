#pragma once

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "io/json_file.hpp"
#include "mechanics/elastic2d.hpp"
#include "mechanics/kirsch_analytic.hpp"
#include "mesh/gmsh_reader.hpp"
#include "support/paths.hpp"

// Shared set-up for the Kirsch tests (docs/learning/phase2.md section 9): quarter model [0, 5a]^2,
// symmetry on x = 0 and y = 0, exact Kirsch tractions on the outer edges, remote 140 MPa along y.

namespace tmfem::test {

constexpr double kPi = 3.14159265358979323846;
constexpr double kPaPerMPa = 1e6;
constexpr double kKirschS_Pa = 140e6;  // REFERENCE_VALUES.md P2-c..e: S = 140 MPa
constexpr double kKirschA_m = 0.010;   // hole radius in cases/kirsch/kirsch.geo
constexpr int kKirschLevels = 5;       // kirsch_t{3,6}_l0..l4, -clscale 1, 1/2, 1/4, 1/8, 1/16
constexpr int kKirschFinest = kKirschLevels - 1;

// E = 207.423 GPa, nu = 0.3 from material_tensile.json (used for every phase 2 test).
inline ElasticMaterial phase2_material() {
  const auto j = read_json(source_dir() / "references/data/material_tensile.json");
  return {j.at("E_GPa").get<double>() * 1e9, j.at("nu_assumed").get<double>()};
}

inline KirschPlate kirsch_plate() { return {kKirschS_Pa, kKirschA_m}; }

inline std::filesystem::path kirsch_mesh(int t_nodes, int level) {
  return source_dir() / "cases/kirsch/msh" / ("kirsch_t" + std::to_string(t_nodes) + "_l" + std::to_string(level) + ".msh");
}

// t = sigma_Kirsch(x) n on an outer edge with outward normal n.
inline TractionField kirsch_traction(const KirschPlate& k, const Eigen::Vector2d& n) {
  return [k, n](const Eigen::Vector2d& x) {
    const Eigen::Vector3d s = kirsch_stress_cartesian(k, x);
    return Eigen::Vector2d(s(0) * n(0) + s(2) * n(1), s(2) * n(0) + s(1) * n(1));
  };
}

inline Elastic2DProblem kirsch_problem(int t_nodes, int level, PlaneCondition plane) {
  const KirschPlate k = kirsch_plate();
  const auto zero = [](const Eigen::Vector2d&) { return 0.0; };
  return {read_gmsh41(kirsch_mesh(t_nodes, level)),
          phase2_material(),
          plane,
          "domain",
          {{"sym_x", 0, zero}, {"sym_y", 1, zero}},
          {{"outer_right", kirsch_traction(k, {1.0, 0.0})}, {"outer_top", kirsch_traction(k, {0.0, 1.0})}}};
}

// One row of references/data/kirsch_sigma140MPa.csv (stresses in MPa, theta in degrees from y).
struct KirschRow {
  double r_over_a, theta_deg, srr_MPa, stt_MPa, srt_MPa;
};

inline std::vector<KirschRow> read_kirsch_csv() {
  std::ifstream in(source_dir() / "references/data/kirsch_sigma140MPa.csv");
  if (!in) throw std::runtime_error("cannot open kirsch_sigma140MPa.csv");
  std::vector<KirschRow> rows;
  std::string line;
  bool header = true;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    if (header) {
      header = false;
      continue;
    }
    std::replace(line.begin(), line.end(), ',', ' ');
    std::istringstream ss(line);
    KirschRow r{};
    ss >> r.r_over_a >> r.theta_deg >> r.srr_MPa >> r.stt_MPa >> r.srt_MPa;
    rows.push_back(r);
  }
  return rows;
}

// A CSV angle mapped into the quarter by reflection (section 9 table): theta' in [0, 90] degrees,
// and the sign that sigma_rt picks up (-1 for an odd number of reflections).
struct QuarterPoint {
  Eigen::Vector2d x_m;
  double theta_rad;
  double rt_sign;
};

inline QuarterPoint map_to_quarter(double r_m, double theta_deg) {
  double t_deg = theta_deg, sign = 1.0;
  if (theta_deg > 270.0) {
    t_deg = 360.0 - theta_deg;
    sign = -1.0;
  } else if (theta_deg > 180.0) {
    t_deg = theta_deg - 180.0;
  } else if (theta_deg > 90.0) {
    t_deg = 180.0 - theta_deg;
    sign = -1.0;
  }
  const double t = t_deg * kPi / 180.0;
  return {Eigen::Vector2d(r_m * std::sin(t), r_m * std::cos(t)), t, sign};
}

}  // namespace tmfem::test
