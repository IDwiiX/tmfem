#pragma once

#include <array>
#include <vector>

#include "mechanics/constitutive.hpp"

namespace tmfem {

// Mixed control of one material point (docs/learning/phase3.md section 10). Each component of
// [xx, yy, zz, xy] is either strain-controlled (engineering shear) or stress-controlled (Pa).
enum class Control { Strain, Stress };

struct DriverStep {
  std::array<Control, 4> control;
  Voigt4 target;
};

struct DriverResult {
  std::vector<Voigt4> strain;
  std::vector<Voigt4> stress_Pa;
  std::vector<MaterialState> state;
  std::vector<int> iterations;  // Newton corrections per step
};

class MaterialPointDriver {
 public:
  explicit MaterialPointDriver(const ConstitutiveModel& model) : model_(model) {}

  // Runs the steps in order; each starts from the previous converged strain and state. Throws if a step
  // does not reach max |sigma_s - target_s| <= 1e-4 Pa in 50 Newton corrections.
  DriverResult run(const std::vector<DriverStep>& steps, const MaterialState& start = {}) const;

 private:
  const ConstitutiveModel& model_;
};

// Uniaxial stress along y, stress-controlled: sigma_yy ramped linearly to sigma_end in n equal steps,
// all other stresses zero.
std::vector<DriverStep> uniaxial_stress_steps(double sigma_end_Pa, int n_steps);

// Uniaxial stress along y, strain-controlled: eps_yy takes each given value, the other stresses are zero.
std::vector<DriverStep> uniaxial_strain_steps(const std::vector<double>& eps_yy);

// Plane-strain bar along x: eps_xx takes each given value, eps_zz = 0, sigma_yy = sigma_xy = 0.
std::vector<DriverStep> plane_strain_bar_steps(const std::vector<double>& eps_xx);

}  // namespace tmfem
