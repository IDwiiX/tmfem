#pragma once

#include <functional>
#include <optional>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include "mesh/mesh.hpp"
#include "thermal/thermal_operators.hpp"

namespace tmfem {

// docs/learning/phase4.md sections 1 to 4: rho c dT/dt = div(k grad T) + q - (2 h_f / t_z)(T - T_inf)
// on a 2D mesh of T3 or T6 cells, consistent capacity only (no row-sum lumping: T6 corners would get 0).
using ScalarField = std::function<double(const Eigen::Vector2d& x_m)>;
using SpaceTimeField = std::function<double(const Eigen::Vector2d& x_m, double t_s)>;

struct TemperatureBC2D {  // T = T_K(t) on every node of a physical group
  std::string physical;
  TimeFunction T_K;
};
struct ConvectionBC2D {  // -k dT/dn = h (T - T_inf(t)) on the edges of a physical group
  std::string physical;
  double h_W_m2K;
  TimeFunction T_inf_K;
};
struct FluxBC2D {  // -k dT/dn = -q_in(t): q_in W/m2 entering the body
  std::string physical;
  TimeFunction q_in_W_m2;
};
// Thin plate cooled on both faces (section 1): a sink (2 h / t_z)(T - T_inf) per unit volume.
struct FaceCooling {
  double h_W_m2K;
  double thickness_m;
  TimeFunction T_inf_K;
};

struct Heat2DProblem {
  Mesh mesh;
  std::string region = "domain";
  double k_W_mK = 0.0;
  double rho_c_J_m3K = 0.0;
  SpaceTimeField source_W_m3;  // empty means none
  std::vector<TemperatureBC2D> temperatures;
  std::vector<ConvectionBC2D> convection;
  std::vector<FluxBC2D> fluxes;
  std::optional<FaceCooling> face_cooling;
};

// Element kernels (section 3). Sizes n x n or n, n = nodes per cell (or per edge).
Eigen::MatrixXd element_conductivity(CellType type, const Eigen::MatrixXd& x_e, double k);  // int k gradN.gradN
Eigen::MatrixXd element_capacity(CellType type, const Eigen::MatrixXd& x_e, double c);      // int c N N^T, exact on straight cells
Eigen::VectorXd element_source(CellType type, const Eigen::MatrixXd& x_e, const ScalarField& q);  // 6-point
Eigen::MatrixXd edge_robin_matrix(CellType type, const Eigen::MatrixXd& x_e, double h);     // int_edge h N N^T
Eigen::VectorXd edge_load(CellType type, const Eigen::MatrixXd& x_e, double value);         // int_edge value N

// M, A = K + H + H_f, F(t) and the Dirichlet data. The source is re-integrated at every F(t) call;
// the convection, flux and face terms are precomputed per unit T_inf or q_in.
ThermalOperators build_operators(const Heat2DProblem& p);

// Nodal interpolant of f.
Eigen::VectorXd nodal_values(const Mesh& mesh, const ScalarField& f);

// int_region w(x) T_h(x) dA with the 6-point rule (exact for straight T6 and polynomial w of degree <= 2).
double integrate_field(const Mesh& mesh, const std::string& region, const Eigen::VectorXd& T, const ScalarField& w);

// ||T_h - T_exact||_L2 over the region (collapsed Gauss, n_gauss^2 points per cell) and h = sqrt(area / cells).
struct FieldError {
  double L2;
  double L2_exact;
  double h_m;
};
FieldError field_l2_error(const Mesh& mesh, const std::string& region, const Eigen::VectorXd& T,
                          const ScalarField& exact, int n_gauss = 5);

// T_h at a physical point (mean over the containing cells); std::nullopt if outside.
std::optional<double> field_at_point(const Mesh& mesh, const std::string& region, const Eigen::VectorXd& T,
                                     const Eigen::Vector2d& x_m);

}  // namespace tmfem
