#include "mechanics/elasticity.hpp"

#include <cmath>
#include <stdexcept>

#include "fem/element_geometry.hpp"
#include "fem/shape_functions.hpp"

namespace tmfem {

Eigen::Matrix3d elasticity_matrix(const ElasticMaterial& m, PlaneCondition plane) {
  const double E = m.E_Pa, nu = m.nu;
  Eigen::Matrix3d D;
  if (plane == PlaneCondition::PlaneStress) {
    const double c = E / (1.0 - nu * nu);
    D << c, c * nu, 0.0,  //
        c * nu, c, 0.0,   //
        0.0, 0.0, c * (1.0 - nu) / 2.0;
  } else {
    const double c = E / ((1.0 + nu) * (1.0 - 2.0 * nu));
    D << c * (1.0 - nu), c * nu, 0.0,  //
        c * nu, c * (1.0 - nu), 0.0,   //
        0.0, 0.0, c * (1.0 - 2.0 * nu) / 2.0;
  }
  return D;
}

double shear_modulus(const ElasticMaterial& m) { return m.E_Pa / (2.0 * (1.0 + m.nu)); }

double kappa(const ElasticMaterial& m, PlaneCondition plane) {
  return plane == PlaneCondition::PlaneStress ? (3.0 - m.nu) / (1.0 + m.nu) : 3.0 - 4.0 * m.nu;
}

double sigma_zz(const ElasticMaterial& m, PlaneCondition plane, const Eigen::Vector3d& sigma) {
  return plane == PlaneCondition::PlaneStress ? 0.0 : m.nu * (sigma(0) + sigma(1));
}

double von_mises(const Eigen::Vector3d& s, double szz) {
  const double dxy = s(0) - s(1), dyz = s(1) - szz, dzx = szz - s(0);
  return std::sqrt(0.5 * (dxy * dxy + dyz * dyz + dzx * dzx) + 3.0 * s(2) * s(2));
}

Eigen::MatrixXd strain_displacement(const Eigen::MatrixXd& dN_dx) {
  const auto n = dN_dx.rows();
  Eigen::MatrixXd B = Eigen::MatrixXd::Zero(3, 2 * n);
  for (Eigen::Index a = 0; a < n; ++a) {
    B(0, 2 * a) = dN_dx(a, 0);
    B(1, 2 * a + 1) = dN_dx(a, 1);
    B(2, 2 * a) = dN_dx(a, 1);  // gamma_xy = du_x/dy + du_y/dx
    B(2, 2 * a + 1) = dN_dx(a, 0);
  }
  return B;
}

QuadRuleTri stiffness_rule(CellType type) {
  switch (type) {
    case CellType::Tri3: return triangle_rule(1);
    case CellType::Tri6: return triangle_rule(3);
    default: throw std::invalid_argument("stiffness_rule: not a triangle");
  }
}

Eigen::MatrixXd element_stiffness(CellType type, const Eigen::MatrixXd& x_e, const Eigen::Matrix3d& D) {
  const QuadRuleTri rule = stiffness_rule(type);
  const auto n_dof = 2 * x_e.rows();
  Eigen::MatrixXd Ke = Eigen::MatrixXd::Zero(n_dof, n_dof);
  for (std::size_t q = 0; q < rule.w.size(); ++q) {
    const CellGeometry g = cell_geometry(shape(type, rule.xi[q]), x_e);
    const Eigen::MatrixXd B = strain_displacement(g.dN_dx);
    Ke += B.transpose() * D * B * (g.detJ * rule.w[q]);
  }
  return Ke;
}

Eigen::VectorXd edge_traction_load(CellType type, const Eigen::MatrixXd& x_e, const TractionField& t) {
  if (type != CellType::Line2 && type != CellType::Line3) {
    throw std::invalid_argument("edge_traction_load: not an edge");
  }
  const QuadRule1D rule = gauss_legendre(5);
  const auto n = x_e.rows();
  Eigen::VectorXd fe = Eigen::VectorXd::Zero(2 * n);
  for (std::size_t q = 0; q < rule.w.size(); ++q) {
    const ShapeEval s = shape(type, Eigen::Vector2d(rule.xi[q], 0.0));
    const FacetGeometry g = facet_geometry(s, x_e);
    const Eigen::Vector2d tq = t(g.x);
    for (Eigen::Index a = 0; a < n; ++a) fe.segment<2>(2 * a) += s.N(a) * tq * (g.ds_dxi * rule.w[q]);
  }
  return fe;
}

}  // namespace tmfem
