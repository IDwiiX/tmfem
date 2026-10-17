#pragma once

#include <vector>

#include <Eigen/Dense>

#include "fem/quadrature.hpp"
#include "mechanics/elasticity.hpp"
#include "mesh/mesh.hpp"

namespace tmfem {

// One-way thermal strain (docs/learning/phase4.md sections 6 to 9): eps_th = alpha (T - T_ref) I.
// T_K holds nodal temperatures on the mechanical mesh, interpolated with the cell's shape functions.
struct ThermalStrain {
  double alpha_per_K;
  double T_ref_K;
  Eigen::VectorXd T_K;
};

// sigma = D eps - beta dT m with m = [1, 1, 0]: beta = E alpha / (1 - nu) (plane stress),
// E alpha / (1 - 2 nu) (plane strain). Sections 7 and 8.
double thermal_stress_coefficient(const ElasticMaterial& m, double alpha_per_K, PlaneCondition plane);

// Out-of-plane stress with the thermal term: 0 (plane stress), nu (sxx + syy) - E alpha dT (plane strain).
double sigma_zz_thermal(const ElasticMaterial& m, PlaneCondition plane, const Eigen::Vector3d& sigma,
                        double alpha_dT);

// Thermal load quadrature (section 9): 1 point for Tri3 (degree 1), 6 points for Tri6 (degree 3).
QuadRuleTri thermal_load_rule(CellType type);

// f_e = int B^T m beta dT_h dA, dT_e = nodal T - T_ref (size n). Size 2n, interleaved dofs.
Eigen::VectorXd element_thermal_load(CellType type, const Eigen::MatrixXd& x_e, const Eigen::VectorXd& dT_e,
                                     double beta);

// Nodal T - T_ref of one cell.
Eigen::VectorXd cell_temperature_rise(const ThermalStrain& th, const CellBlock& block, int cell);

// Global thermal load over the region cells (size 2 n_nodes).
Eigen::VectorXd assemble_thermal_load(const Mesh& mesh, const std::vector<const CellBlock*>& cells,
                                      const ThermalStrain& th, double beta);

}  // namespace tmfem
