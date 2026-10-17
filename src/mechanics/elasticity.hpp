#pragma once

#include <functional>

#include <Eigen/Dense>

#include "fem/quadrature.hpp"
#include "mesh/mesh.hpp"

namespace tmfem {

// docs/learning/phase2.md sections 3 and 4. Voigt order [xx, yy, xy], engineering shear strain.
enum class PlaneCondition { PlaneStress, PlaneStrain };

struct ElasticMaterial {
  double E_Pa;
  double nu;
};

Eigen::Matrix3d elasticity_matrix(const ElasticMaterial& m, PlaneCondition plane);
double shear_modulus(const ElasticMaterial& m);
// Kolosov constant: (3 - nu)/(1 + nu) in plane stress, 3 - 4 nu in plane strain.
double kappa(const ElasticMaterial& m, PlaneCondition plane);

// Out-of-plane normal stress: 0 in plane stress, nu (sigma_xx + sigma_yy) in plane strain.
double sigma_zz(const ElasticMaterial& m, PlaneCondition plane, const Eigen::Vector3d& sigma);
double von_mises(const Eigen::Vector3d& sigma, double sigma_zz);

// Strain-displacement matrix B (3 x 2n) from shape-function gradients dN/dx (n x 2).
Eigen::MatrixXd strain_displacement(const Eigen::MatrixXd& dN_dx);

// Stiffness quadrature (section 6): 1 point for Tri3, 3 points for Tri6.
QuadRuleTri stiffness_rule(CellType type);

// K_e = sum_q B^T D B det J w (unit thickness). Size 2n x 2n.
Eigen::MatrixXd element_stiffness(CellType type, const Eigen::MatrixXd& x_e, const Eigen::Matrix3d& D);

// f_e = int N^T t ds over a Line2/Line3 edge with 5-point Gauss (section 7); traction in Pa.
using TractionField = std::function<Eigen::Vector2d(const Eigen::Vector2d& x_m)>;
Eigen::VectorXd edge_traction_load(CellType type, const Eigen::MatrixXd& x_e, const TractionField& t);

}  // namespace tmfem
