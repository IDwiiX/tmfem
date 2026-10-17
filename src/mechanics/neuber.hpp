#pragma once

#include "mechanics/hardening.hpp"

namespace tmfem {

// Neuber's rule at a notch root with uniaxial stress (docs/learning/phase3.md section 12):
// sigma * eps = (Kt S)^2 / E with eps = sigma / E + alpha(sigma) on the hardening curve.
struct NeuberEstimate {
  double sigma_Pa;
  double eps;
  double alpha;
};

NeuberEstimate neuber_uniaxial(const HardeningLaw& h, double E_Pa, double KtS_Pa);

}  // namespace tmfem
