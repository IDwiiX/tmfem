#include "coupling/thermal_strain.hpp"

#include <stdexcept>

#include "fem/assembly.hpp"
#include "fem/element_geometry.hpp"
#include "fem/shape_functions.hpp"

namespace tmfem {

double thermal_stress_coefficient(const ElasticMaterial& m, double alpha_per_K, PlaneCondition plane) {
  return plane == PlaneCondition::PlaneStress ? m.E_Pa * alpha_per_K / (1.0 - m.nu)
                                              : m.E_Pa * alpha_per_K / (1.0 - 2.0 * m.nu);
}

double sigma_zz_thermal(const ElasticMaterial& m, PlaneCondition plane, const Eigen::Vector3d& sigma,
                        double alpha_dT) {
  return plane == PlaneCondition::PlaneStress ? 0.0 : m.nu * (sigma(0) + sigma(1)) - m.E_Pa * alpha_dT;
}

QuadRuleTri thermal_load_rule(CellType type) {
  switch (type) {
    case CellType::Tri3: return triangle_rule(1);
    case CellType::Tri6: return triangle_rule(6);
    default: throw std::invalid_argument("thermal_load_rule: not a triangle");
  }
}

Eigen::VectorXd element_thermal_load(CellType type, const Eigen::MatrixXd& x_e, const Eigen::VectorXd& dT_e,
                                     double beta) {
  const QuadRuleTri rule = thermal_load_rule(type);
  const auto n = x_e.rows();
  Eigen::VectorXd fe = Eigen::VectorXd::Zero(2 * n);
  for (std::size_t q = 0; q < rule.w.size(); ++q) {
    const ShapeEval s = shape(type, rule.xi[q]);
    const CellGeometry g = cell_geometry(s, x_e);
    const double w = beta * s.N.dot(dT_e) * g.detJ * rule.w[q];
    // B^T m: row 2a = dN_a/dx, row 2a+1 = dN_a/dy (section 9).
    for (Eigen::Index a = 0; a < n; ++a) fe.segment<2>(2 * a) += w * g.dN_dx.row(a).transpose();
  }
  return fe;
}

Eigen::VectorXd cell_temperature_rise(const ThermalStrain& th, const CellBlock& block, int cell) {
  Eigen::VectorXd dT(block.nodes.cols());
  for (Eigen::Index a = 0; a < dT.size(); ++a) dT(a) = th.T_K(block.nodes(cell, a)) - th.T_ref_K;
  return dT;
}

Eigen::VectorXd assemble_thermal_load(const Mesh& mesh, const std::vector<const CellBlock*>& cells,
                                      const ThermalStrain& th, double beta) {
  if (th.T_K.size() != mesh.n_nodes()) throw std::invalid_argument("ThermalStrain: one temperature per node expected");
  Eigen::VectorXd f = Eigen::VectorXd::Zero(2 * mesh.n_nodes());
  for (const CellBlock* b : cells) {
    for (int e = 0; e < b->nodes.rows(); ++e) {
      const Eigen::VectorXd fe = element_thermal_load(b->type, mesh.cell_coordinates(*b, e),
                                                      cell_temperature_rise(th, *b, e), beta);
      const std::vector<int> dofs = cell_dofs(*b, e, 2);
      for (std::size_t i = 0; i < dofs.size(); ++i) f(dofs[i]) += fe(static_cast<Eigen::Index>(i));
    }
  }
  return f;
}

}  // namespace tmfem
