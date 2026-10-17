#include "mechanics/elastoplastic2d.hpp"

#include <cmath>
#include <map>
#include <sstream>
#include <stdexcept>

#include "fem/assembly.hpp"
#include "fem/element_geometry.hpp"
#include "fem/shape_functions.hpp"

namespace tmfem {

namespace {

// Same singularity guard as solve_elastic2d (docs/learning/phase2.md section 7).
constexpr double kMinRelativePivot = 1e-10;

std::vector<int> constrained_list(const std::map<int, double>& bc, Eigen::VectorXd& values) {
  std::vector<int> dofs;
  values.resize(static_cast<Eigen::Index>(bc.size()));
  for (const auto& [dof, v] : bc) {
    values(static_cast<Eigen::Index>(dofs.size())) = v;
    dofs.push_back(dof);
  }
  return dofs;
}

Eigen::VectorXd gather(const Eigen::VectorXd& u, const std::vector<int>& dofs) {
  Eigen::VectorXd out(static_cast<Eigen::Index>(dofs.size()));
  for (std::size_t i = 0; i < dofs.size(); ++i) out(static_cast<Eigen::Index>(i)) = u(dofs[i]);
  return out;
}

}  // namespace

QuadRuleTri plastic_rule(CellType type) {
  switch (type) {
    case CellType::Tri3: return triangle_rule(1);
    case CellType::Tri6: return triangle_rule(6);
    default: throw std::invalid_argument("plastic_rule: not a triangle");
  }
}

ElastoplasticSolver::ElastoplasticSolver(Elastoplastic2DProblem problem, NewtonOptions options)
    : p_(std::move(problem)),
      opt_(options),
      f_hat_(assemble_traction_load(p_.mesh, p_.tractions)),
      part_(2 * p_.mesh.n_nodes(), constrained_list(collect_dirichlet(p_.mesh, p_.displacements), g_hat_)),
      u_(Eigen::VectorXd::Zero(2 * p_.mesh.n_nodes())),
      f_int_(Eigen::VectorXd::Zero(2 * p_.mesh.n_nodes())) {
  if (!p_.model) throw std::invalid_argument("ElastoplasticSolver: no constitutive model");
  constrained_ = part_.constrained_dofs();
  for (const CellBlock* b : p_.mesh.blocks_of(p_.region)) {
    const QuadRuleTri rule = plastic_rule(b->type);
    for (int e = 0; e < b->nodes.rows(); ++e) {
      const Eigen::MatrixXd x_e = p_.mesh.cell_coordinates(*b, e);
      cells_.push_back({cell_dofs(*b, e, 2), static_cast<int>(qp_.size()), static_cast<int>(rule.w.size())});
      for (std::size_t q = 0; q < rule.w.size(); ++q) {
        const CellGeometry g = cell_geometry(shape(b->type, rule.xi[q]), x_e);
        QuadPoint pt;
        pt.x_m = g.x;
        pt.weight_m2 = g.detJ * rule.w[q];
        pt.B = strain_displacement(g.dN_dx);
        pt.cell = static_cast<int>(cells_.size()) - 1;
        qp_.push_back(std::move(pt));
      }
    }
  }
}

ElastoplasticSolver::Linearization ElastoplasticSolver::linearize(const Eigen::VectorXd& u_m) const {
  const int n = 2 * p_.mesh.n_nodes();
  Linearization lin;
  lin.f_int_N = Eigen::VectorXd::Zero(n);
  lin.response.resize(qp_.size());
  std::vector<Eigen::Triplet<double>> triplets;
  for (const Cell& c : cells_) {
    const Eigen::VectorXd u_e = gather(u_m, c.dofs);
    const auto nd = static_cast<Eigen::Index>(c.dofs.size());
    Eigen::MatrixXd Ke = Eigen::MatrixXd::Zero(nd, nd);
    Eigen::VectorXd fe = Eigen::VectorXd::Zero(nd);
    for (int q = c.first_qp; q < c.first_qp + c.n_qp; ++q) {
      const QuadPoint& pt = qp_[static_cast<std::size_t>(q)];
      PlaneResponse r = plane_update(*p_.model, p_.plane, pt.B * u_e, pt.state);
      fe += pt.B.transpose() * r.stress_Pa * pt.weight_m2;
      Ke += pt.B.transpose() * r.tangent_Pa * pt.B * pt.weight_m2;
      lin.response[static_cast<std::size_t>(q)] = std::move(r);
    }
    for (Eigen::Index i = 0; i < nd; ++i) {
      lin.f_int_N(c.dofs[static_cast<std::size_t>(i)]) += fe(i);
      for (Eigen::Index j = 0; j < nd; ++j) {
        triplets.emplace_back(c.dofs[static_cast<std::size_t>(i)], c.dofs[static_cast<std::size_t>(j)], Ke(i, j));
      }
    }
  }
  lin.K_T.resize(n, n);
  lin.K_T.setFromTriplets(triplets.begin(), triplets.end());
  return lin;
}

NewtonIterate ElastoplasticSolver::measure(const Linearization& lin, const Eigen::VectorXd& f_ext) const {
  const Eigen::VectorXd R = lin.f_int_N - f_ext;
  NewtonIterate it;
  it.residual_N = part_.restrict_free(R).norm();
  it.reference_force_N = std::max(part_.restrict_free(f_ext).norm(), gather(R, constrained_).norm());
  it.plastic.reserve(lin.response.size());
  for (const PlaneResponse& r : lin.response) it.plastic.push_back(r.plastic ? 1 : 0);
  return it;
}

Eigen::VectorXd ElastoplasticSolver::solve_correction(const Linearization& lin, const Eigen::VectorXd& rhs_f) {
  const Eigen::SparseMatrix<double> K_ff = part_.block_ff(lin.K_T);
  if (!pattern_analysed_) {
    ldlt_.analyzePattern(K_ff);
    pattern_analysed_ = true;
  }
  ldlt_.factorize(K_ff);
  if (ldlt_.info() != Eigen::Success) throw std::runtime_error("ElastoplasticSolver: factorisation failed");
  const Eigen::VectorXd pivots = ldlt_.vectorD();
  if (pivots.minCoeff() < kMinRelativePivot * pivots.cwiseAbs().maxCoeff()) {
    throw std::runtime_error("ElastoplasticSolver: K_T,ff is singular or indefinite (rigid-body mode not constrained?)");
  }
  return ldlt_.solve(rhs_f);
}

StepReport ElastoplasticSolver::step_to(double load_factor) {
  const Eigen::VectorXd f_ext = load_factor * f_hat_;
  const Eigen::VectorXd g = load_factor * g_hat_;
  Eigen::VectorXd u = u_;
  Linearization lin = linearize(u);
  StepReport rep;
  rep.load_factor = load_factor;
  rep.iterates.push_back(measure(lin, f_ext));
  const auto small = [&](const NewtonIterate& it) { return it.residual_N <= opt_.tolerance * it.reference_force_N; };
  rep.converged = small(rep.iterates[0]) && (g - gather(u, constrained_)).norm() == 0.0;

  for (int k = 1; k <= opt_.max_iterations && !(rep.converged && !opt_.run_to_stagnation); ++k) {
    // K_ff du_f = -R_f - K_fc du_c, with du_c nonzero only in the first correction (section 9).
    const Eigen::VectorXd du_c = g - gather(u, constrained_);
    const Eigen::VectorXd rhs = -part_.restrict_free(lin.f_int_N - f_ext) - part_.block_fc(lin.K_T) * du_c;
    const Eigen::VectorXd du_f = solve_correction(lin, rhs);
    u = part_.combine(part_.restrict_free(u) + du_f, g);
    lin = linearize(u);
    rep.iterates.push_back(measure(lin, f_ext));
    const double r = rep.iterates.back().residual_N;
    const double r_prev = rep.iterates[rep.iterates.size() - 2].residual_N;
    if (rep.converged && (r == 0.0 || r > 0.5 * r_prev)) break;  // run_to_stagnation: floor reached
    if (small(rep.iterates.back())) rep.converged = true;
  }
  if (!rep.converged) {
    std::ostringstream msg;
    msg << "ElastoplasticSolver: no convergence at load factor " << load_factor << "; residual / F_ref:";
    for (const NewtonIterate& it : rep.iterates) msg << ' ' << it.residual_N / it.reference_force_N;
    throw std::runtime_error(msg.str());
  }

  u_ = u;
  lambda_ = load_factor;
  f_int_ = lin.f_int_N - f_ext;  // reactions live at the constrained dofs
  for (std::size_t q = 0; q < qp_.size(); ++q) {
    const PlaneResponse& r = lin.response[q];
    qp_[q].state = r.state;
    qp_[q].stress_Pa = r.stress_Pa;
    qp_[q].sigma_zz_Pa = r.sigma_zz_Pa;
    qp_[q].plastic = r.plastic;
  }
  return rep;
}

Eigen::VectorXd ElastoplasticSolver::reactions_N() const { return gather(f_int_, constrained_); }

Eigen::Vector3d ElastoplasticSolver::quad_strain(std::size_t q) const {
  const QuadPoint& pt = qp_.at(q);
  return pt.B * gather(u_, cells_[static_cast<std::size_t>(pt.cell)].dofs);
}

std::optional<Eigen::Vector3d> strain_at_point(const Mesh& mesh, const std::string& region, const Eigen::VectorXd& u_m,
                                               const Eigen::Vector2d& x_m) {
  Eigen::Vector3d sum = Eigen::Vector3d::Zero();
  int found = 0;
  for (const CellBlock* b : mesh.blocks_of(region)) {
    for (int e = 0; e < b->nodes.rows(); ++e) {
      const Eigen::MatrixXd x_e = mesh.cell_coordinates(*b, e);
      const auto xi = inverse_map(b->type, x_e, x_m);
      if (!xi) continue;
      const CellGeometry g = cell_geometry(shape(b->type, *xi), x_e);
      sum += strain_displacement(g.dN_dx) * gather(u_m, cell_dofs(*b, e, 2));
      ++found;
    }
  }
  if (found == 0) return std::nullopt;
  return Eigen::Vector3d(sum / found);
}

}  // namespace tmfem
