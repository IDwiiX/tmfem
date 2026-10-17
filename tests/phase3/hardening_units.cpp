#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <stdexcept>

#include "mechanics/hardening.hpp"
#include "phase3/phase3_common.hpp"

// Supporting unit tests for the hardening laws (docs/learning/phase3.md section 4). No reference row.

using namespace tmfem;

TEST_CASE("Tabulated hardening built from the true curve", "[phase3]") {
  const auto h = test::tabulated_hardening();
  const auto& a = h->alpha();
  const auto& s = h->sigma_Pa();
  REQUIRE(a.size() == 8);
  CHECK(a[0] == 0.0);
  CHECK(s[0] == test::kSigmaStart_Pa);
  CHECK(std::abs(a[1] - 0.0086197) < 1e-7);  // 0.010524 - 394.9973 / 207423
  CHECK(std::abs(a[2] - 0.0293766) < 1e-7);  // 0.031581 - 457.2405 / 207423
  for (std::size_t i = 0; i < a.size(); ++i) CHECK(std::abs(h->yield_stress_Pa(a[i]) - s[i]) <= 1e-15 * s[i]);
  for (std::size_t i = 0; i + 1 < a.size(); ++i) {
    const double slope = (s[i + 1] - s[i]) / (a[i + 1] - a[i]);
    CHECK(std::abs(h->slope_Pa(a[i]) - slope) <= 1e-12 * slope);  // right-hand slope at the knot
    const double mid = 0.5 * (a[i] + a[i + 1]);
    CHECK(std::abs(h->yield_stress_Pa(mid) - 0.5 * (s[i] + s[i + 1])) <= 1e-12 * s[i]);
  }
  // Beyond the table the last slope is extended.
  const double last = (s[7] - s[6]) / (a[7] - a[6]);
  CHECK(std::abs(h->yield_stress_Pa(0.2) - (s[7] + last * (0.2 - a[7]))) <= 1e-12 * s[7]);
  CHECK(std::abs(plastic_strain_at_stress(*h, s[2]) - a[2]) <= 1e-12 * a[2]);
  CHECK(plastic_strain_at_stress(*h, 300e6) == 0.0);
  // Clipping is limited to data round-off: with a wrong E of 150 GPa the yield point would have alpha = -5.8e-4.
  const test::TrueCurve c = test::true_curve();
  CHECK_THROWS_AS(TabulatedHardening::from_true_curve(c.eps, c.sigma_Pa, 150e9, test::kSigmaStart_Pa),
                  std::invalid_argument);
}

TEST_CASE("Ludwik hardening slope matches finite differences", "[phase3]") {
  const LudwikHardening h(310.86e6, 559.27e6, 0.4077);
  for (const double a : {1e-4, 1e-3, 1e-2, 0.1}) {
    const double d = 1e-6 * a;
    const double fd = (h.yield_stress_Pa(a + d) - h.yield_stress_Pa(a - d)) / (2.0 * d);
    CHECK(std::abs(fd - h.slope_Pa(a)) <= 1e-8 * h.slope_Pa(a));
  }
  CHECK(std::isinf(h.slope_Pa(0.0)));
  CHECK(h.yield_stress_Pa(0.0) == 310.86e6);
}
