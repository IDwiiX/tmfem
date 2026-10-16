#pragma once

#include <vector>

namespace tmfem {

// 1D mesh of linear elements: node coordinates in m, strictly increasing from 0 to L.
struct IntervalMesh {
  std::vector<double> x_m;

  int n_nodes() const { return static_cast<int>(x_m.size()); }
  int n_elements() const { return n_nodes() - 1; }
  double length_m() const { return x_m.back() - x_m.front(); }
};

// n_elements equal elements on [0, L].
IntervalMesh uniform_interval(double L_m, int n_elements);

// Geometric grading on [0, L]: element i has length h0 * ratio^i, finest at x = 0.
// ratio = 1 gives the uniform mesh.
IntervalMesh graded_interval(double L_m, int n_elements, double ratio);

// Index of the node closest to x.
int nearest_node(const IntervalMesh& mesh, double x_m);

// Linear interpolation of a nodal field at x (x inside the mesh).
double interpolate(const IntervalMesh& mesh, const std::vector<double>& nodal, double x_m);

}  // namespace tmfem
