#include "mesh/interval_mesh.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace tmfem {

IntervalMesh uniform_interval(double L_m, int n_elements) {
  return graded_interval(L_m, n_elements, 1.0);
}

IntervalMesh graded_interval(double L_m, int n_elements, double ratio) {
  if (L_m <= 0.0 || n_elements < 1 || ratio <= 0.0) {
    throw std::invalid_argument("graded_interval: need L > 0, n >= 1, ratio > 0");
  }
  IntervalMesh mesh;
  mesh.x_m.resize(static_cast<std::size_t>(n_elements) + 1);
  if (std::abs(ratio - 1.0) < 1e-14) {
    for (int i = 0; i <= n_elements; ++i) {
      mesh.x_m[static_cast<std::size_t>(i)] = L_m * i / n_elements;
    }
    return mesh;
  }
  // Partial sums of a geometric series: x_i = h0 (ratio^i - 1) / (ratio - 1), x_n = L.
  const double denom = std::pow(ratio, n_elements) - 1.0;
  for (int i = 0; i <= n_elements; ++i) {
    mesh.x_m[static_cast<std::size_t>(i)] = L_m * (std::pow(ratio, i) - 1.0) / denom;
  }
  mesh.x_m.back() = L_m;
  return mesh;
}

int nearest_node(const IntervalMesh& mesh, double x_m) {
  const auto it = std::min_element(mesh.x_m.begin(), mesh.x_m.end(), [x_m](double a, double b) {
    return std::abs(a - x_m) < std::abs(b - x_m);
  });
  return static_cast<int>(it - mesh.x_m.begin());
}

double interpolate(const IntervalMesh& mesh, const std::vector<double>& nodal, double x_m) {
  const auto& x = mesh.x_m;
  if (x_m < x.front() || x_m > x.back()) {
    throw std::out_of_range("interpolate: x outside the mesh");
  }
  auto upper = std::upper_bound(x.begin(), x.end(), x_m);
  if (upper == x.end()) {
    return nodal.back();
  }
  const auto i = static_cast<std::size_t>(upper - x.begin()) - 1;
  const double s = (x_m - x[i]) / (x[i + 1] - x[i]);
  return (1.0 - s) * nodal[i] + s * nodal[i + 1];
}

}  // namespace tmfem
