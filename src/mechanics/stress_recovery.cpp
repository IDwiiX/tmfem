#include "mechanics/stress_recovery.hpp"

#include "fem/assembly.hpp"
#include "fem/element_geometry.hpp"
#include "fem/shape_functions.hpp"

namespace tmfem {

namespace {

Eigen::VectorXd cell_displacements(const CellBlock& block, int cell, const Eigen::VectorXd& u_m) {
  const std::vector<int> dofs = cell_dofs(block, cell, 2);
  Eigen::VectorXd u_e(static_cast<Eigen::Index>(dofs.size()));
  for (std::size_t i = 0; i < dofs.size(); ++i) u_e(static_cast<Eigen::Index>(i)) = u_m(dofs[i]);
  return u_e;
}

// D B u_e - beta dT_h(xi) m (docs/learning/phase4.md section 9; the thermal term is 0 without p.thermal).
Eigen::Vector3d stress_in_cell(const Elastic2DProblem& p, const Eigen::Matrix3d& D, const CellBlock& block, int cell,
                               const Eigen::VectorXd& u_m, const Eigen::Vector2d& xi) {
  const ShapeEval s = shape(block.type, xi);
  const CellGeometry g = cell_geometry(s, p.mesh.cell_coordinates(block, cell));
  Eigen::Vector3d sigma = D * strain_displacement(g.dN_dx) * cell_displacements(block, cell, u_m);
  if (p.thermal) {
    const double dT = s.N.dot(cell_temperature_rise(*p.thermal, block, cell));
    sigma -= p.beta() * dT * Eigen::Vector3d(1.0, 1.0, 0.0);
  }
  return sigma;
}

}  // namespace

Eigen::Vector3d cell_stress(const Elastic2DProblem& p, const Eigen::VectorXd& u_m, const CellBlock& block, int cell,
                            const Eigen::Vector2d& xi) {
  return stress_in_cell(p, p.D(), block, cell, u_m, xi);
}

std::vector<Eigen::Vector3d> gauss_point_stresses(const Elastic2DProblem& p, const Eigen::VectorXd& u_m) {
  const Eigen::Matrix3d D = p.D();
  std::vector<Eigen::Vector3d> out;
  for (const CellBlock* b : p.mesh.blocks_of(p.region)) {
    const QuadRuleTri rule = stiffness_rule(b->type);
    for (int e = 0; e < b->nodes.rows(); ++e) {
      for (const auto& xi : rule.xi) out.push_back(stress_in_cell(p, D, *b, e, u_m, xi));
    }
  }
  return out;
}

Eigen::MatrixXd nodal_average_stress(const Elastic2DProblem& p, const Eigen::VectorXd& u_m) {
  const Eigen::Matrix3d D = p.D();
  Eigen::MatrixXd sum = Eigen::MatrixXd::Zero(p.mesh.n_nodes(), 3);
  Eigen::VectorXd count = Eigen::VectorXd::Zero(p.mesh.n_nodes());
  for (const CellBlock* b : p.mesh.blocks_of(p.region)) {
    const Eigen::MatrixXd xi_nodes = reference_nodes(b->type);
    for (int e = 0; e < b->nodes.rows(); ++e) {
      for (Eigen::Index a = 0; a < xi_nodes.rows(); ++a) {
        const int node = b->nodes(e, a);
        sum.row(node) += stress_in_cell(p, D, *b, e, u_m, xi_nodes.row(a).transpose()).transpose();
        count(node) += 1.0;
      }
    }
  }
  for (Eigen::Index i = 0; i < sum.rows(); ++i) {
    if (count(i) > 0.0) sum.row(i) /= count(i);
  }
  return sum;
}

std::optional<Eigen::Vector3d> stress_at_point(const Elastic2DProblem& p, const Eigen::VectorXd& u_m,
                                               const Eigen::Vector2d& x_m) {
  const Eigen::Matrix3d D = p.D();
  Eigen::Vector3d sum = Eigen::Vector3d::Zero();
  int found = 0;
  for (const CellBlock* b : p.mesh.blocks_of(p.region)) {
    for (int e = 0; e < b->nodes.rows(); ++e) {
      const auto xi = inverse_map(b->type, p.mesh.cell_coordinates(*b, e), x_m);
      if (!xi) continue;
      sum += stress_in_cell(p, D, *b, e, u_m, *xi);
      ++found;
    }
  }
  if (found == 0) return std::nullopt;
  return Eigen::Vector3d(sum / found);
}

NodalStressFields nodal_stress_fields(const Elastic2DProblem& p, const Eigen::VectorXd& u_m) {
  NodalStressFields f{nodal_average_stress(p, u_m), Eigen::VectorXd(p.mesh.n_nodes()), Eigen::VectorXd(p.mesh.n_nodes())};
  for (int i = 0; i < p.mesh.n_nodes(); ++i) {
    const Eigen::Vector3d s = f.sigma_Pa.row(i).transpose();
    const double alpha_dT = p.thermal ? p.thermal->alpha_per_K * (p.thermal->T_K(i) - p.thermal->T_ref_K) : 0.0;
    f.sigma_zz_Pa(i) = sigma_zz_thermal(p.material, p.plane, s, alpha_dT);
    f.von_mises_Pa(i) = von_mises(s, f.sigma_zz_Pa(i));
  }
  return f;
}

}  // namespace tmfem
