#pragma once

#include <functional>

#include <Eigen/Dense>
#include <Eigen/SparseCholesky>

#include "fem/dirichlet.hpp"
#include "thermal/heat1d.hpp"

namespace tmfem {

// Theta method (docs/learning/phase1.md section 5):
// (M + theta dt A) T^{n+1} = (M - (1 - theta) dt A) T^n + dt [theta F^{n+1} + (1 - theta) F^n],
// Dirichlet by elimination. theta = 1: backward Euler; theta = 0.5: Crank-Nicolson.
// Works on assembled matrices only, so it is the same for 1D and 2D.
class ThetaIntegrator {
 public:
  ThetaIntegrator(ThermalOperators ops, double theta, double dt_s);
  ThetaIntegrator(const ThetaIntegrator&) = delete;
  ThetaIntegrator& operator=(const ThetaIntegrator&) = delete;

  // Changes the step size and refactorises.
  void set_dt(double dt_s);
  double dt() const { return dt_; }
  double theta() const { return theta_; }
  const ThermalOperators& operators() const { return ops_; }

  // One step from (t_n, T_n) to t_n + dt.
  Eigen::VectorXd step(const Eigen::VectorXd& T_n, double t_n) const;

 private:
  void factorize();

  ThermalOperators ops_;
  double theta_;
  double dt_;
  DirichletPartition partition_;
  SparseMatrix S_fc_;
  SparseMatrix R_;  // M - (1 - theta) dt A
  Eigen::SimplicialLDLT<SparseMatrix> solver_;
};

// Calls observer(t, T) at t0 and after every step; t_n = t0 + n dt (no accumulation).
using StepObserver = std::function<void(double t_s, const Eigen::VectorXd& T)>;
Eigen::VectorXd march(const ThetaIntegrator& integrator, Eigen::VectorXd T0, double t0_s, int n_steps,
                      const StepObserver& observer);

// Steady state A T = F(t) with Dirichlet elimination.
Eigen::VectorXd solve_steady(const ThermalOperators& ops, double t_s);

}  // namespace tmfem
