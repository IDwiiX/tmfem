#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <string>
#include <vector>

#include "io/case_elastic2d.hpp"
#include "io/csv_reader.hpp"
#include "mesh/gmsh_reader.hpp"
#include "support/paths.hpp"

// Plane-strain thermal stress through the case output (nodal CSV): sigma_zz must carry -E alpha dT
// (docs/learning/phase4.md sections 8 and 9), and von Mises must use that sigma_zz.
// Reference: docs/REFERENCE_VALUES.md rows P4-a (held bar, plane strain: sigma_xx = sigma_zz = -E alpha dT/(1-nu),
// sigma_yy = 0, so von Mises = E alpha dT/(1-nu); 1e-8 relative) and P4-b (free body, plane strain: in-plane
// stress 0, sigma_zz = -E alpha dT, so von Mises = E alpha dT; 1e-10 relative).
// Bar 100 mm x 20 mm (cases/bar/bar.geo), uniform dT. The phase 2 formula sigma_zz = nu (sxx + syy) gives
// -nu E alpha dT/(1-nu) (held bar) and 0 (free body) instead, so this test fails without the thermal term.

using namespace tmfem;

namespace {

constexpr double kTolHeld = 1e-8;   // P4-a: 1e-8 relative
constexpr double kTolFree = 1e-10;  // P4-b: 1e-10 relative
const ElasticMaterial kSteel{200e9, 0.3};
constexpr double kAlpha_per_K = 1.2e-5, kTref_K = 293.15, kdT_K = 100.0;

struct NodalOutput {
  std::vector<double> sxx, syy, sxy, szz, vm;  // MPa
};

NodalOutput run_case(int t_nodes, bool held, const std::string& tag) {
  const auto zero = [](const Eigen::Vector2d&) { return 0.0; };
  std::vector<DisplacementBC> bcs{{"left", 0, zero}, {"corner", 1, zero}};
  if (held) bcs.push_back({"right", 0, zero});
  Mesh mesh = read_gmsh41(test::source_dir() / ("cases/bar/bar_t" + std::to_string(t_nodes) + ".msh"));
  const Eigen::VectorXd T = Eigen::VectorXd::Constant(mesh.n_nodes(), kTref_K + kdT_K);
  const auto dir = test::results_dir() / "phase4";
  std::filesystem::create_directories(dir);
  const Elastic2DCase c{{std::move(mesh), kSteel, PlaneCondition::PlaneStrain, "domain", bcs, {},
                         ThermalStrain{kAlpha_per_K, kTref_K, T}},
                        dir / ("thermal_output_" + tag + ".vtu"),
                        dir / ("thermal_output_" + tag + ".csv")};
  run_elastic2d_case(c);
  const auto cols = read_numeric_csv(c.output_nodal_csv, 9);  // x y ux uy sxx syy sxy szz vm
  return {cols[4], cols[5], cols[6], cols[7], cols[8]};
}

double max_abs_error(const std::vector<double>& v, double target) {
  double worst = 0.0;
  for (const double x : v) worst = std::max(worst, std::abs(x - target));
  return worst;
}

}  // namespace

TEST_CASE("Plane-strain thermal sigma_zz and von Mises in the case output", "[phase4][P4-a][P4-b][output]") {
  const double scale_MPa = kSteel.E_Pa * kAlpha_per_K * kdT_K / 1e6;  // E alpha dT = 240 MPa
  for (const int t : {3, 6}) {
    const std::string name = "T" + std::to_string(t);
    {  // P4-a held bar
      const NodalOutput o = run_case(t, true, "held_t" + std::to_string(t));
      const double s = -scale_MPa / (1.0 - kSteel.nu);
      const double e_sxx = max_abs_error(o.sxx, s), e_syy = max_abs_error(o.syy, 0.0);
      const double e_szz = max_abs_error(o.szz, s), e_vm = max_abs_error(o.vm, std::abs(s));
      std::printf("P4-a %s plane strain output: max |err|/|ref| sxx %.2e syy %.2e szz %.2e von Mises %.2e\n",
                  name.c_str(), e_sxx / std::abs(s), e_syy / std::abs(s), e_szz / std::abs(s), e_vm / std::abs(s));
      INFO(name << " held bar");
      CHECK(e_sxx < kTolHeld * std::abs(s));
      CHECK(e_syy < kTolHeld * std::abs(s));
      CHECK(e_szz < kTolHeld * std::abs(s));
      CHECK(e_vm < kTolHeld * std::abs(s));
    }
    {  // P4-b free body
      const NodalOutput o = run_case(t, false, "free_t" + std::to_string(t));
      const double e_in = std::max({max_abs_error(o.sxx, 0.0), max_abs_error(o.syy, 0.0), max_abs_error(o.sxy, 0.0)});
      const double e_szz = max_abs_error(o.szz, -scale_MPa), e_vm = max_abs_error(o.vm, scale_MPa);
      std::printf("P4-b %s plane strain output: max |err|/(E alpha dT) in-plane %.2e szz %.2e von Mises %.2e\n",
                  name.c_str(), e_in / scale_MPa, e_szz / scale_MPa, e_vm / scale_MPa);
      INFO(name << " free body");
      CHECK(e_in < kTolFree * scale_MPa);
      CHECK(e_szz < kTolFree * scale_MPa);
      CHECK(e_vm < kTolFree * scale_MPa);
    }
  }
}
