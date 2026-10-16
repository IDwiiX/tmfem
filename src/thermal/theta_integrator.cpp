#include "thermal/theta_integrator.hpp"

#include <stdexcept>

namespace tmfem {

ThetaIntegrator::ThetaIntegrator(ThermalOperators ops, double theta, double dt_s)
    : ops_(std::move(ops)),
      theta_(theta),
      dt_(dt_s),
      partition_(static_cast<int>(ops_.M.rows()), ops_.dirichlet_dofs) {
  if (theta < 0.5 || theta > 1.0) {
    throw std::invalid_argument("ThetaIntegrator: theta must be in [0.5, 1] (unconditional stability)");
  }
  factorize();
}

void ThetaIntegrator::set_dt(double dt_s) {
  dt_ = dt_s;
  factorize();
}

void ThetaIntegrator::factorize() {
  if (dt_ <= 0.0) throw std::invalid_argument("ThetaIntegrator: dt must be positive");
  const SparseMatrix S = ops_.M + (theta_ * dt_) * ops_.A;
  R_ = ops_.M - ((1.0 - theta_) * dt_) * ops_.A;
  S_fc_ = partition_.block_fc(S);
  solver_.compute(partition_.block_ff(S));
  if (solver_.info() != Eigen::Success) {
    throw std::runtime_error("ThetaIntegrator: factorisation of M + theta dt A failed");
  }
}

Eigen::VectorXd ThetaIntegrator::step(const Eigen::VectorXd& T_n, double t_n) const {
  const double t_next = t_n + dt_;
  Eigen::VectorXd b = R_ * T_n + dt_ * (theta_ * ops_.F(t_next));
  if (theta_ < 1.0) b += dt_ * (1.0 - theta_) * ops_.F(t_n);
  const Eigen::VectorXd g = ops_.dirichlet_values(t_next);
  Eigen::VectorXd rhs_f = partition_.restrict_free(b);
  if (g.size() > 0) rhs_f -= S_fc_ * g;
  const Eigen::VectorXd T_f = solver_.solve(rhs_f);
  return partition_.combine(T_f, g);
}

Eigen::VectorXd march(const ThetaIntegrator& integrator, Eigen::VectorXd T0, double t0_s, int n_steps,
                      const StepObserver& observer) {
  Eigen::VectorXd T = std::move(T0);
  if (observer) observer(t0_s, T);
  for (int n = 0; n < n_steps; ++n) {
    const double t_n = t0_s + n * integrator.dt();
    T = integrator.step(T, t_n);
    if (observer) observer(t0_s + (n + 1) * integrator.dt(), T);
  }
  return T;
}

Eigen::VectorXd solve_steady(const ThermalOperators& ops, double t_s) {
  const DirichletPartition part(static_cast<int>(ops.A.rows()), ops.dirichlet_dofs);
  const Eigen::VectorXd g = ops.dirichlet_values(t_s);
  Eigen::VectorXd rhs_f = part.restrict_free(ops.F(t_s));
  if (g.size() > 0) rhs_f -= part.block_fc(ops.A) * g;
  Eigen::SimplicialLDLT<SparseMatrix> solver(part.block_ff(ops.A));
  if (solver.info() != Eigen::Success) {
    throw std::runtime_error("solve_steady: singular system (no Dirichlet or Robin condition?)");
  }
  return part.combine(solver.solve(rhs_f), g);
}

}  // namespace tmfem
