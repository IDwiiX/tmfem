#pragma once

#include <functional>
#include <variant>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "fem/line2.hpp"
#include "mesh/interval_mesh.hpp"
#include "thermal/thermal_operators.hpp"

namespace tmfem {

// Boundary conditions at one end of the slab (docs/learning/phase1.md section 1).
struct DirichletBC {
  TimeFunction value;  // T_K (or C for diffusion)
};
struct FluxBC {
  TimeFunction q_in;  // W/m2 entering the body
};
struct ConvectiveBC {
  double h_W_m2K;
  TimeFunction T_inf_K;  // -k dT/dn = h (T - T_inf)
};
using BoundaryCondition = std::variant<DirichletBC, FluxBC, ConvectiveBC>;

// rho_c dT/dt = d/dx(k dT/dx) + q(x, t). Diffusion: rho_c = 1, k = D.
struct Heat1DProblem {
  IntervalMesh mesh;
  double k;      // W/(m K), or D in m2/s
  double rho_c;  // J/(m3 K), or 1
  BoundaryCondition left;
  BoundaryCondition right;
  std::function<double(double x_m, double t_s)> source;  // W/m3; empty means none
};

ThermalOperators build_operators(const Heat1DProblem& problem, MassType mass);

// Nodal interpolant of f on the mesh.
Eigen::VectorXd nodal_values(const IntervalMesh& mesh, const std::function<double(double x_m)>& f);

// Diffusive flux -k dT/dx in each element of a nodal field (constant per linear element).
std::vector<double> element_flux(const IntervalMesh& mesh, double k, const Eigen::VectorXd& T);

}  // namespace tmfem
