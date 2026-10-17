#pragma once

#include <cmath>
#include <memory>
#include <vector>

#include "io/csv_reader.hpp"
#include "io/json_file.hpp"
#include "mechanics/elastoplastic2d.hpp"
#include "mechanics/hardening.hpp"
#include "mechanics/j2_plasticity.hpp"
#include "mechanics/tensile_fit.hpp"
#include "phase2/kirsch_common.hpp"
#include "support/paths.hpp"

// Shared data for the phase 3 tests (docs/learning/phase3.md sections 4 and 11): the true tensile curve,
// the tabulated and source-Ludwik hardening laws, and E = 207.423 GPa, nu = 0.3 (phase2_material()).

namespace tmfem::test {

// The tabulated curve starts at the true-curve point at yield, 310.5611 MPa; its computed eps_p of
// -3.2e-6 is clipped to 0 (kickoff decision, docs/learning/phase3.md section 4).
constexpr double kSigmaStart_Pa = 310.5611e6;

struct TrueCurve {
  std::vector<double> eps;
  std::vector<double> sigma_Pa;
};

inline TrueCurve true_curve() {
  const auto cols = read_numeric_csv(source_dir() / "references/data/tensile_true.csv", 2);
  TrueCurve c{cols[0], cols[1]};
  for (double& s : c.sigma_Pa) s *= kPaPerMPa;
  return c;
}

// The plastic points (sigma >= 310.5611 MPa) as (alpha, sigma), alpha = eps - sigma / E, first clipped to 0.
struct PlasticPoints {
  std::vector<double> eps;
  std::vector<double> alpha;
  std::vector<double> sigma_Pa;
};

inline PlasticPoints plastic_points() {
  const TrueCurve c = true_curve();
  const double E = phase2_material().E_Pa;
  PlasticPoints p;
  for (std::size_t i = 0; i < c.eps.size(); ++i) {
    if (c.sigma_Pa[i] < kSigmaStart_Pa) continue;
    p.eps.push_back(c.eps[i]);
    p.alpha.push_back(p.alpha.empty() ? 0.0 : c.eps[i] - c.sigma_Pa[i] / E);
    p.sigma_Pa.push_back(c.sigma_Pa[i]);
  }
  return p;
}

inline std::shared_ptr<const TabulatedHardening> tabulated_hardening() {
  const TrueCurve c = true_curve();
  return std::make_shared<TabulatedHardening>(
      TabulatedHardening::from_true_curve(c.eps, c.sigma_Pa, phase2_material().E_Pa, kSigmaStart_Pa));
}

// Source Ludwik parameters (material_tensile.json): sigma_y0 = Re = 310.862 MPa, K = 559.27 MPa, m = 0.4077.
inline LudwikParams source_ludwik() {
  const auto j = read_json(source_dir() / "references/data/material_tensile.json");
  return {j.at("Re_MPa").get<double>() * kPaPerMPa, j.at("ludwik_K_MPa").get<double>() * kPaPerMPa,
          j.at("ludwik_m").get<double>()};
}

inline std::shared_ptr<const LudwikHardening> ludwik_hardening(const LudwikParams& p) {
  return std::make_shared<LudwikHardening>(p.sigma_y0_Pa, p.K_Pa, p.m);
}

inline std::shared_ptr<const J2Plasticity> j2(std::shared_ptr<const HardeningLaw> h) {
  return std::make_shared<J2Plasticity>(phase2_material(), std::move(h));
}

// Elastoplastic version of the phase 2 Kirsch quarter model (elastic Kirsch tractions, scaled by lambda).
inline Elastoplastic2DProblem kirsch_plastic_problem(int t_nodes, int level, PlaneCondition plane,
                                                     std::shared_ptr<const ConstitutiveModel> model) {
  Elastic2DProblem e = kirsch_problem(t_nodes, level, plane);
  return {std::move(e.mesh), std::move(model), plane, e.region, std::move(e.displacements), std::move(e.tractions)};
}

// Equivalent stress of a 2D stress plus sigma_zz.
inline double equivalent(const Eigen::Vector3d& s, double szz) { return equivalent_stress(Voigt4(s(0), s(1), szz, s(2))); }

}  // namespace tmfem::test
