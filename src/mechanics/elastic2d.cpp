#include "mechanics/elastic2d.hpp"

#include <cmath>
#include <map>
#include <stdexcept>

#include <Eigen/SparseCholesky>

#include "fem/assembly.hpp"
#include "fem/dirichlet.hpp"

namespace tmfem {

namespace {

// A pivot below this fraction of the largest one means K_ff is singular to round-off. Well-posed
// meshes here stay many orders above it; a missed rigid mode gives about 1e-15.
constexpr double kMinRelativePivot = 1e-10;

}  // namespace

std::map<int, double> collect_dirichlet(const Mesh& mesh, const std::vector<DisplacementBC>& displacements) {
  std::map<int, double> values;
  for (const DisplacementBC& bc : displacements) {
    if (bc.component != 0 && bc.component != 1) throw std::invalid_argument("DisplacementBC: component must be 0 or 1");
    for (const int node : mesh.nodes_of(bc.physical)) {
      const double v = bc.value_m(mesh.x_m.row(node).transpose());
      const int dof = 2 * node + bc.component;
      const auto [it, inserted] = values.emplace(dof, v);
      if (!inserted && std::abs(it->second - v) > 1e-14 * (1.0 + std::abs(v))) {
        throw std::invalid_argument("DisplacementBC: conflicting values at a shared node");
      }
    }
  }
  return values;
}

Eigen::SparseMatrix<double> assemble_stiffness(const Elastic2DProblem& p) {
  const Eigen::Matrix3d D = p.D();
  return assemble_matrix(p.mesh, p.mesh.blocks_of(p.region), 2,
                         [&D](CellType type, const Eigen::MatrixXd& x_e) { return element_stiffness(type, x_e, D); });
}

Eigen::VectorXd assemble_traction_load(const Elastic2DProblem& p) { return assemble_traction_load(p.mesh, p.tractions); }

Eigen::VectorXd assemble_traction_load(const Mesh& mesh, const std::vector<TractionBC>& tractions) {
  Eigen::VectorXd f = Eigen::VectorXd::Zero(2 * mesh.n_nodes());
  for (const TractionBC& bc : tractions) {
    assemble_vector(mesh, mesh.blocks_of(bc.physical), 2,
                    [&bc](CellType type, const Eigen::MatrixXd& x_e) {
                      return edge_traction_load(type, x_e, bc.traction_Pa);
                    },
                    f);
  }
  return f;
}

Elastic2DSolution solve_elastic2d(const Elastic2DProblem& p) {
  const Eigen::SparseMatrix<double> K = assemble_stiffness(p);
  const Eigen::VectorXd f = assemble_traction_load(p);
  const std::map<int, double> bc = collect_dirichlet(p.mesh, p.displacements);

  std::vector<int> constrained;
  Eigen::VectorXd g(static_cast<Eigen::Index>(bc.size()));
  for (const auto& [dof, v] : bc) {
    g(static_cast<Eigen::Index>(constrained.size())) = v;
    constrained.push_back(dof);
  }
  const DirichletPartition part(static_cast<int>(f.size()), constrained);
  const Eigen::SparseMatrix<double> K_ff = part.block_ff(K);
  const Eigen::VectorXd rhs = part.restrict_free(f) - part.block_fc(K) * g;

  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt(K_ff);
  if (ldlt.info() != Eigen::Success) throw std::runtime_error("solve_elastic2d: factorisation failed");
  const Eigen::VectorXd u_f = ldlt.solve(rhs);
  if (ldlt.info() != Eigen::Success) throw std::runtime_error("solve_elastic2d: solve failed");
  // LDLT of a singular matrix (a missed rigid mode) "succeeds": the zero pivot comes out as round-off,
  // about 1e-15 of the largest and possibly positive. So compare the pivots relative to the largest.
  const Eigen::VectorXd pivots = ldlt.vectorD();
  if (pivots.minCoeff() < kMinRelativePivot * pivots.cwiseAbs().maxCoeff()) {
    throw std::runtime_error("solve_elastic2d: K_ff is singular or indefinite (rigid-body mode not constrained?)");
  }

  Elastic2DSolution sol;
  sol.u_m = part.combine(u_f, g);
  sol.constrained_dofs = constrained;
  const Eigen::VectorXd Ku_minus_f = K * sol.u_m - f;
  sol.reactions_N.resize(static_cast<Eigen::Index>(constrained.size()));
  for (std::size_t i = 0; i < constrained.size(); ++i) sol.reactions_N(static_cast<Eigen::Index>(i)) = Ku_minus_f(constrained[i]);
  return sol;
}

}  // namespace tmfem
