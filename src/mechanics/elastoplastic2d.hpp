#pragma once

#include <memory>
#include <optional>
#include <string>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

#include "fem/dirichlet.hpp"
#include "fem/quadrature.hpp"
#include "mechanics/constitutive.hpp"
#include "mechanics/elastic2d.hpp"
#include "mechanics/plane_reduction.hpp"
#include "mesh/mesh.hpp"

namespace tmfem {

// Elastoplastic problem on a 2D mesh, unit thickness (docs/learning/phase3.md section 9). The boundary
// conditions are the phase 2 structs, given at load factor 1 and scaled by the load factor lambda.
struct Elastoplastic2DProblem {
  Mesh mesh;
  std::shared_ptr<const ConstitutiveModel> model;
  PlaneCondition plane;
  std::string region = "domain";
  std::vector<DisplacementBC> displacements;
  std::vector<TractionBC> tractions;
};

struct NewtonOptions {
  double tolerance = 1e-10;  // converged when ||R_free|| <= tolerance * F_ref
  int max_iterations = 25;
  // After convergence, keep iterating until the residual no longer halves (the round-off floor), so the
  // whole residual sequence is recorded. Used by the convergence-rate test P3-d.
  bool run_to_stagnation = false;
};

// One Newton iterate: residual norm on the free dofs, the force scale F_ref = max(||f_ext,free||,
// ||reactions||), and which quadrature points were plastic at this iterate.
struct NewtonIterate {
  double residual_N;
  double reference_force_N;
  std::vector<char> plastic;
};

struct StepReport {
  double load_factor = 0.0;
  bool converged = false;
  std::vector<NewtonIterate> iterates;  // iterates[0] is the start of the step, before any correction
  int corrections() const { return static_cast<int>(iterates.size()) - 1; }
};

// Converged data at one quadrature point.
struct QuadPoint {
  Eigen::Vector2d x_m;
  double weight_m2;     // w det J (unit thickness)
  Eigen::MatrixXd B;    // 3 x 2n strain-displacement matrix
  int cell;             // index into the solver's cell list
  MaterialState state;
  Eigen::Vector3d stress_Pa = Eigen::Vector3d::Zero();
  double sigma_zz_Pa = 0.0;
  bool plastic = false;  // yielded during the last converged step
};

// Quadrature for nonlinear residuals and tangents: 1 point on T3 (constant strain, exact), 6 points
// (degree 4) on T6 (SPEC section 5).
QuadRuleTri plastic_rule(CellType type);

class ElastoplasticSolver {
 public:
  explicit ElastoplasticSolver(Elastoplastic2DProblem problem, NewtonOptions options = {});

  // Newton from the last converged state to load factor lambda. Commits the states on convergence and
  // throws (with the residual history) if max_iterations is reached; there is no hidden step cutting.
  StepReport step_to(double load_factor);

  // Internal force and consistent tangent at displacement u, every point evaluated from its converged state.
  struct Linearization {
    Eigen::SparseMatrix<double> K_T;
    Eigen::VectorXd f_int_N;
    std::vector<PlaneResponse> response;  // one per quadrature point
  };
  Linearization linearize(const Eigen::VectorXd& u_m) const;

  const Elastoplastic2DProblem& problem() const { return p_; }
  double load_factor() const { return lambda_; }
  const Eigen::VectorXd& displacement_m() const { return u_; }
  const std::vector<QuadPoint>& quad_points() const { return qp_; }
  const std::vector<int>& constrained_dofs() const { return constrained_; }
  // Traction load at load factor 1 (N per m of thickness).
  const Eigen::VectorXd& traction_load_N() const { return f_hat_; }
  // Reactions f_int - f_ext at the constrained dofs, converged state.
  Eigen::VectorXd reactions_N() const;
  // Total strain [xx, yy, gamma_xy] at quadrature point q, converged displacement.
  Eigen::Vector3d quad_strain(std::size_t q) const;

 private:
  struct Cell {
    std::vector<int> dofs;
    int first_qp;
    int n_qp;
  };

  NewtonIterate measure(const Linearization& lin, const Eigen::VectorXd& f_ext) const;
  Eigen::VectorXd solve_correction(const Linearization& lin, const Eigen::VectorXd& rhs_f);

  Elastoplastic2DProblem p_;
  NewtonOptions opt_;
  std::vector<Cell> cells_;
  std::vector<QuadPoint> qp_;
  std::vector<int> constrained_;
  Eigen::VectorXd g_hat_;  // prescribed displacements at load factor 1
  Eigen::VectorXd f_hat_;
  DirichletPartition part_;
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt_;
  bool pattern_analysed_ = false;
  double lambda_ = 0.0;
  Eigen::VectorXd u_;
  Eigen::VectorXd f_int_;
};

// Total strain [xx, yy, gamma_xy] at a physical point from the displacement field, averaged over the cells
// of the region that contain it; std::nullopt if none does.
std::optional<Eigen::Vector3d> strain_at_point(const Mesh& mesh, const std::string& region, const Eigen::VectorXd& u_m,
                                               const Eigen::Vector2d& x_m);

}  // namespace tmfem
