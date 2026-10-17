#pragma once

#include <functional>
#include <map>
#include <string>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "mechanics/elasticity.hpp"
#include "mesh/mesh.hpp"

namespace tmfem {

using ScalarField = std::function<double(const Eigen::Vector2d& x_m)>;

// u_component = value(x) at every node of a physical group (edge or point); component 0 = x, 1 = y.
struct DisplacementBC {
  std::string physical;
  int component;
  ScalarField value_m;
};

// Traction t(x) in Pa on every edge of a physical group.
struct TractionBC {
  std::string physical;
  TractionField traction_Pa;
};

// Linear elasticity on a 2D mesh, unit thickness (docs/learning/phase2.md sections 2 and 7).
struct Elastic2DProblem {
  Mesh mesh;
  ElasticMaterial material;
  PlaneCondition plane;
  std::string region = "domain";
  std::vector<DisplacementBC> displacements;
  std::vector<TractionBC> tractions;

  Eigen::Matrix3d D() const { return elasticity_matrix(material, plane); }
};

// Nodal displacements (interleaved, m) and reactions (N per m of thickness) at the constrained dofs.
struct Elastic2DSolution {
  Eigen::VectorXd u_m;
  std::vector<int> constrained_dofs;
  Eigen::VectorXd reactions_N;
};

Eigen::SparseMatrix<double> assemble_stiffness(const Elastic2DProblem& p);
Eigen::VectorXd assemble_traction_load(const Elastic2DProblem& p);
Eigen::VectorXd assemble_traction_load(const Mesh& mesh, const std::vector<TractionBC>& tractions);

// Constrained dofs (2 node + component) and their values, in increasing dof order. Throws if a dof
// shared by two conditions gets two different values.
std::map<int, double> collect_dirichlet(const Mesh& mesh, const std::vector<DisplacementBC>& displacements);

// Assembles, eliminates the Dirichlet dofs and solves with SimplicialLDLT. Throws if the
// constraints leave a rigid-body mode (a pivot below 1e-10 of the largest).
Elastic2DSolution solve_elastic2d(const Elastic2DProblem& p);

}  // namespace tmfem
