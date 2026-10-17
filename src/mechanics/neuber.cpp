#include "mechanics/neuber.hpp"

#include <stdexcept>

namespace tmfem {

NeuberEstimate neuber_uniaxial(const HardeningLaw& h, double E_Pa, double KtS_Pa) {
  if (!(KtS_Pa > 0.0) || !(E_Pa > 0.0)) throw std::invalid_argument("neuber_uniaxial: need Kt S > 0 and E > 0");
  const double target = KtS_Pa * KtS_Pa / E_Pa;
  const auto product = [&](double s) { return s * (s / E_Pa + plastic_strain_at_stress(h, s)); };
  // product(s) increases with s; product(0) = 0 < target <= product(Kt S), so the root lies in (0, Kt S].
  double lo = 0.0, hi = KtS_Pa;
  for (int it = 0; it < 200 && hi - lo > 1e-15 * hi; ++it) {
    const double mid = 0.5 * (lo + hi);
    (product(mid) < target ? lo : hi) = mid;
  }
  const double s = 0.5 * (lo + hi);
  const double a = plastic_strain_at_stress(h, s);
  return {s, s / E_Pa + a, a};
}

}  // namespace tmfem
