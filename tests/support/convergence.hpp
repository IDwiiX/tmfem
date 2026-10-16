#pragma once

#include <cmath>
#include <cstdio>
#include <vector>

namespace tmfem::test {

// Least-squares slope of log(error) against log(h): the observed order of convergence.
inline double observed_order(const std::vector<double>& h, const std::vector<double>& error) {
  const auto n = static_cast<double>(h.size());
  double sx = 0, sy = 0, sxx = 0, sxy = 0;
  for (std::size_t i = 0; i < h.size(); ++i) {
    const double x = std::log(h[i]);
    const double y = std::log(error[i]);
    sx += x;
    sy += y;
    sxx += x * x;
    sxy += x * y;
  }
  return (n * sxy - sx * sy) / (n * sxx - sx * sx);
}

// Prints h, error and the pairwise order between consecutive levels.
inline void print_orders(const char* label, const std::vector<double>& h, const std::vector<double>& error) {
  std::printf("%s\n  %-12s %-12s %s\n", label, "h", "error", "pairwise order");
  for (std::size_t i = 0; i < h.size(); ++i) {
    if (i == 0) {
      std::printf("  %-12.4e %-12.4e\n", h[i], error[i]);
    } else {
      const double p = std::log(error[i - 1] / error[i]) / std::log(h[i - 1] / h[i]);
      std::printf("  %-12.4e %-12.4e %.4f\n", h[i], error[i], p);
    }
  }
  std::printf("  least-squares order: %.4f\n", observed_order(h, error));
}

}  // namespace tmfem::test
