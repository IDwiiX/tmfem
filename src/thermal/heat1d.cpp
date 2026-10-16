#include "thermal/heat1d.hpp"

namespace tmfem {

namespace {

template <class... Ts>
struct Overloaded : Ts... {
  using Ts::operator()...;
};
template <class... Ts>
Overloaded(Ts...) -> Overloaded<Ts...>;

// Volume matrices M and K from the element contributions of section 4.
void assemble_volume(const Heat1DProblem& p, MassType mass, Eigen::SparseMatrix<double>& M,
                     Eigen::SparseMatrix<double>& K) {
  const auto& x = p.mesh.x_m;
  const int n = p.mesh.n_nodes();
  std::vector<Eigen::Triplet<double>> tm;
  std::vector<Eigen::Triplet<double>> tk;
  for (int e = 0; e < p.mesh.n_elements(); ++e) {
    const auto ue = static_cast<std::size_t>(e);
    const Eigen::Matrix2d Me = element_capacity(p.rho_c, x[ue], x[ue + 1], mass);
    const Eigen::Matrix2d Ke = element_stiffness(p.k, x[ue], x[ue + 1]);
    for (int a = 0; a < 2; ++a) {
      for (int b = 0; b < 2; ++b) {
        if (Me(a, b) != 0.0) tm.emplace_back(e + a, e + b, Me(a, b));
        tk.emplace_back(e + a, e + b, Ke(a, b));
      }
    }
  }
  M.resize(n, n);
  K.resize(n, n);
  M.setFromTriplets(tm.begin(), tm.end());
  K.setFromTriplets(tk.begin(), tk.end());
}

}  // namespace

ThermalOperators build_operators(const Heat1DProblem& problem, MassType mass) {
  ThermalOperators ops;
  assemble_volume(problem, mass, ops.M, ops.A);
  const int n = problem.mesh.n_nodes();
  const std::pair<int, const BoundaryCondition*> ends[2] = {{0, &problem.left}, {n - 1, &problem.right}};

  // Robin term H: h on the boundary node's diagonal (section 3).
  std::vector<TimeFunction> dirichlet_fns;
  for (const auto& [node, bc] : ends) {
    if (const auto* c = std::get_if<ConvectiveBC>(bc)) {
      ops.A.coeffRef(node, node) += c->h_W_m2K;
    } else if (const auto* d = std::get_if<DirichletBC>(bc)) {
      ops.dirichlet_dofs.push_back(node);
      dirichlet_fns.push_back(d->value);
    }
  }
  ops.A.makeCompressed();

  ops.F = [problem, n](double t_s) {
    Eigen::VectorXd F = Eigen::VectorXd::Zero(n);
    if (problem.source) {
      const auto& x = problem.mesh.x_m;
      const auto q = [&problem, t_s](double x_m) { return problem.source(x_m, t_s); };
      for (int e = 0; e < problem.mesh.n_elements(); ++e) {
        const auto ue = static_cast<std::size_t>(e);
        F.segment<2>(e) += element_source(q, x[ue], x[ue + 1]);
      }
    }
    const std::pair<int, const BoundaryCondition*> bcs[2] = {{0, &problem.left},
                                                             {n - 1, &problem.right}};
    for (const auto& [node, bc] : bcs) {
      std::visit(Overloaded{[](const DirichletBC&) {},
                            [&](const FluxBC& f) { F(node) += f.q_in(t_s); },
                            [&](const ConvectiveBC& c) { F(node) += c.h_W_m2K * c.T_inf_K(t_s); }},
                 *bc);
    }
    return F;
  };

  ops.dirichlet_values = [dirichlet_fns](double t_s) {
    Eigen::VectorXd g(static_cast<Eigen::Index>(dirichlet_fns.size()));
    for (std::size_t i = 0; i < dirichlet_fns.size(); ++i) {
      g(static_cast<Eigen::Index>(i)) = dirichlet_fns[i](t_s);
    }
    return g;
  };
  return ops;
}

Eigen::VectorXd nodal_values(const IntervalMesh& mesh, const std::function<double(double x_m)>& f) {
  Eigen::VectorXd v(mesh.n_nodes());
  for (int i = 0; i < mesh.n_nodes(); ++i) v(i) = f(mesh.x_m[static_cast<std::size_t>(i)]);
  return v;
}

std::vector<double> element_flux(const IntervalMesh& mesh, double k, const Eigen::VectorXd& T) {
  std::vector<double> q(static_cast<std::size_t>(mesh.n_elements()));
  for (int e = 0; e < mesh.n_elements(); ++e) {
    const auto ue = static_cast<std::size_t>(e);
    q[ue] = -k * (T(e + 1) - T(e)) / (mesh.x_m[ue + 1] - mesh.x_m[ue]);
  }
  return q;
}

}  // namespace tmfem
