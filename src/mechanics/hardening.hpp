#pragma once

#include <vector>

namespace tmfem {

// Isotropic hardening: yield stress sigma_y(alpha) and slope H = d sigma_y / d alpha, alpha being the
// equivalent plastic strain (docs/learning/phase3.md section 4).
class HardeningLaw {
 public:
  virtual ~HardeningLaw() = default;
  virtual double yield_stress_Pa(double alpha) const = 0;
  virtual double slope_Pa(double alpha) const = 0;
};

// sigma_y = sigma_y0 + H alpha (unit tests: closed-form answers).
class LinearHardening : public HardeningLaw {
 public:
  LinearHardening(double sigma_y0_Pa, double H_Pa);
  double yield_stress_Pa(double alpha) const override;
  double slope_Pa(double alpha) const override;

 private:
  double sigma_y0_Pa_, H_Pa_;
};

// Ludwik: sigma_y = sigma_y0 + K alpha^m. For m < 1 the slope m K alpha^(m-1) is infinite at alpha = 0.
class LudwikHardening : public HardeningLaw {
 public:
  LudwikHardening(double sigma_y0_Pa, double K_Pa, double m);
  double yield_stress_Pa(double alpha) const override;
  double slope_Pa(double alpha) const override;

 private:
  double sigma_y0_Pa_, K_Pa_, m_;
};

// Piecewise-linear sigma_y(alpha) through (alpha_i, sigma_i), alpha_0 = 0 and alpha strictly increasing.
// The slope at a knot is the right-hand segment's; the last segment's slope is extended beyond the table.
class TabulatedHardening : public HardeningLaw {
 public:
  TabulatedHardening(std::vector<double> alpha, std::vector<double> sigma_Pa);

  // From a true stress-strain curve: keeps the points with sigma >= sigma_start_Pa, sets
  // alpha = eps_true - sigma / E and clips the first alpha to 0 (it is -3.2e-6 for tensile_true.csv
  // because the data's elastic slope is 210 GPa, not E; section 4).
  static TabulatedHardening from_true_curve(const std::vector<double>& eps_true,
                                            const std::vector<double>& sigma_true_Pa, double E_Pa,
                                            double sigma_start_Pa);

  double yield_stress_Pa(double alpha) const override;
  double slope_Pa(double alpha) const override;
  const std::vector<double>& alpha() const { return alpha_; }
  const std::vector<double>& sigma_Pa() const { return sigma_; }

 private:
  std::size_t segment(double alpha) const;  // index i of [alpha_i, alpha_i+1) containing alpha (clamped)

  std::vector<double> alpha_;
  std::vector<double> sigma_;
};

// alpha with sigma_y(alpha) = sigma (0 if sigma is at or below the initial yield stress), by bisection.
double plastic_strain_at_stress(const HardeningLaw& h, double sigma_Pa);

}  // namespace tmfem
