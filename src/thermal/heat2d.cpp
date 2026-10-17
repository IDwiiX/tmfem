#include "thermal/heat2d.hpp"

#include <cmath>
#include <map>
#include <stdexcept>

#include "fem/assembly.hpp"
#include "fem/element_geometry.hpp"
#include "fem/quadrature.hpp"
#include "fem/shape_functions.hpp"

namespace tmfem {

namespace {

void require_triangle(CellType type, const char* who) {
  if (type != CellType::Tri3 && type != CellType::Tri6) throw std::invalid_argument(std::string(who) + ": not a triangle");
}

void require_edge(CellType type, const char* who) {
  if (type != CellType::Line2 && type != CellType::Line3) throw std::invalid_argument(std::string(who) + ": not an edge");
}

// Exact rules on straight cells (section 3 table): grad.grad has degree 0 (T3) or 2 (T6); N N^T degree 2 or 4.
QuadRuleTri conductivity_rule(CellType type) { return triangle_rule(type == CellType::Tri3 ? 1 : 3); }
QuadRuleTri capacity_rule(CellType type) { return triangle_rule(type == CellType::Tri3 ? 3 : 6); }

// Sum over the edges of the given physical groups of the element vector int N ds (unit value).
Eigen::VectorXd unit_edge_vector(const Mesh& mesh, const std::string& physical, double value) {
  Eigen::VectorXd f = Eigen::VectorXd::Zero(mesh.n_nodes());
  assemble_vector(mesh, mesh.blocks_of(physical), 1,
                  [value](CellType type, const Eigen::MatrixXd& x_e) { return edge_load(type, x_e, value); }, f);
  return f;
}

Eigen::VectorXd source_vector(const Mesh& mesh, const std::vector<const CellBlock*>& cells, const SpaceTimeField& q,
                              double t_s) {
  Eigen::VectorXd f = Eigen::VectorXd::Zero(mesh.n_nodes());
  const ScalarField q_t = [&q, t_s](const Eigen::Vector2d& x) { return q(x, t_s); };
  assemble_vector(mesh, cells, 1,
                  [&q_t](CellType type, const Eigen::MatrixXd& x_e) { return element_source(type, x_e, q_t); }, f);
  return f;
}

}  // namespace

Eigen::MatrixXd element_conductivity(CellType type, const Eigen::MatrixXd& x_e, double k) {
  require_triangle(type, "element_conductivity");
  const QuadRuleTri rule = conductivity_rule(type);
  Eigen::MatrixXd Ke = Eigen::MatrixXd::Zero(x_e.rows(), x_e.rows());
  for (std::size_t q = 0; q < rule.w.size(); ++q) {
    const CellGeometry g = cell_geometry(shape(type, rule.xi[q]), x_e);
    Ke += (k * g.detJ * rule.w[q]) * g.dN_dx * g.dN_dx.transpose();
  }
  return Ke;
}

Eigen::MatrixXd element_capacity(CellType type, const Eigen::MatrixXd& x_e, double c) {
  require_triangle(type, "element_capacity");
  const QuadRuleTri rule = capacity_rule(type);
  Eigen::MatrixXd Me = Eigen::MatrixXd::Zero(x_e.rows(), x_e.rows());
  for (std::size_t q = 0; q < rule.w.size(); ++q) {
    const ShapeEval s = shape(type, rule.xi[q]);
    const CellGeometry g = cell_geometry(s, x_e);
    Me += (c * g.detJ * rule.w[q]) * s.N * s.N.transpose();
  }
  return Me;
}

Eigen::VectorXd element_source(CellType type, const Eigen::MatrixXd& x_e, const ScalarField& q_W_m3) {
  require_triangle(type, "element_source");
  const QuadRuleTri rule = triangle_rule(6);
  Eigen::VectorXd fe = Eigen::VectorXd::Zero(x_e.rows());
  for (std::size_t q = 0; q < rule.w.size(); ++q) {
    const ShapeEval s = shape(type, rule.xi[q]);
    const CellGeometry g = cell_geometry(s, x_e);
    fe += (q_W_m3(g.x) * g.detJ * rule.w[q]) * s.N;
  }
  return fe;
}

Eigen::MatrixXd edge_robin_matrix(CellType type, const Eigen::MatrixXd& x_e, double h) {
  require_edge(type, "edge_robin_matrix");
  const QuadRule1D rule = gauss_legendre(5);
  Eigen::MatrixXd He = Eigen::MatrixXd::Zero(x_e.rows(), x_e.rows());
  for (std::size_t q = 0; q < rule.w.size(); ++q) {
    const ShapeEval s = shape(type, Eigen::Vector2d(rule.xi[q], 0.0));
    const FacetGeometry g = facet_geometry(s, x_e);
    He += (h * g.ds_dxi * rule.w[q]) * s.N * s.N.transpose();
  }
  return He;
}

Eigen::VectorXd edge_load(CellType type, const Eigen::MatrixXd& x_e, double value) {
  require_edge(type, "edge_load");
  const QuadRule1D rule = gauss_legendre(5);
  Eigen::VectorXd fe = Eigen::VectorXd::Zero(x_e.rows());
  for (std::size_t q = 0; q < rule.w.size(); ++q) {
    const ShapeEval s = shape(type, Eigen::Vector2d(rule.xi[q], 0.0));
    fe += (value * facet_geometry(s, x_e).ds_dxi * rule.w[q]) * s.N;
  }
  return fe;
}

ThermalOperators build_operators(const Heat2DProblem& p) {
  if (p.k_W_mK <= 0.0 || p.rho_c_J_m3K <= 0.0) throw std::invalid_argument("Heat2DProblem: k and rho c must be positive");
  const Mesh& mesh = p.mesh;
  const std::vector<const CellBlock*> cells = mesh.blocks_of(p.region);
  ThermalOperators ops;
  const double rho_c = p.rho_c_J_m3K, k = p.k_W_mK;
  ops.M = assemble_matrix(mesh, cells, 1,
                          [rho_c](CellType t, const Eigen::MatrixXd& x) { return element_capacity(t, x, rho_c); });
  ops.A = assemble_matrix(mesh, cells, 1,
                          [k](CellType t, const Eigen::MatrixXd& x) { return element_conductivity(t, x, k); });

  // Time-dependent loads as (unit vector, time function) pairs: F_conv = sum T_inf(t) int h N ds, etc.
  std::vector<std::pair<Eigen::VectorXd, TimeFunction>> loads;
  for (const ConvectionBC2D& c : p.convection) {
    const double h = c.h_W_m2K;
    ops.A += assemble_matrix(mesh, mesh.blocks_of(c.physical), 1,
                             [h](CellType t, const Eigen::MatrixXd& x) { return edge_robin_matrix(t, x, h); });
    loads.emplace_back(unit_edge_vector(mesh, c.physical, h), c.T_inf_K);
  }
  for (const FluxBC2D& f : p.fluxes) loads.emplace_back(unit_edge_vector(mesh, f.physical, 1.0), f.q_in_W_m2);
  if (p.face_cooling) {
    const double c = 2.0 * p.face_cooling->h_W_m2K / p.face_cooling->thickness_m;
    ops.A += assemble_matrix(mesh, cells, 1,
                             [c](CellType t, const Eigen::MatrixXd& x) { return element_capacity(t, x, c); });
    Eigen::VectorXd unit = Eigen::VectorXd::Zero(mesh.n_nodes());
    assemble_vector(mesh, cells, 1,
                    [c](CellType t, const Eigen::MatrixXd& x) {
                      return element_source(t, x, [c](const Eigen::Vector2d&) { return c; });
                    },
                    unit);
    loads.emplace_back(unit, p.face_cooling->T_inf_K);
  }
  ops.A.makeCompressed();

  ops.F = [mesh = p.mesh, region = p.region, source = p.source_W_m3, loads](double t_s) {
    Eigen::VectorXd F = Eigen::VectorXd::Zero(mesh.n_nodes());
    if (source) F += source_vector(mesh, mesh.blocks_of(region), source, t_s);
    for (const auto& [unit, fn] : loads) F += fn(t_s) * unit;
    return F;
  };

  std::map<int, TimeFunction> dirichlet;  // node -> T(t), increasing node order
  for (const TemperatureBC2D& bc : p.temperatures) {
    for (const int node : mesh.nodes_of(bc.physical)) dirichlet.emplace(node, bc.T_K);
  }
  std::vector<TimeFunction> fns;
  for (const auto& [node, fn] : dirichlet) {
    ops.dirichlet_dofs.push_back(node);
    fns.push_back(fn);
  }
  ops.dirichlet_values = [fns](double t_s) {
    Eigen::VectorXd g(static_cast<Eigen::Index>(fns.size()));
    for (std::size_t i = 0; i < fns.size(); ++i) g(static_cast<Eigen::Index>(i)) = fns[i](t_s);
    return g;
  };
  return ops;
}

Eigen::VectorXd nodal_values(const Mesh& mesh, const ScalarField& f) {
  Eigen::VectorXd v(mesh.n_nodes());
  for (int i = 0; i < mesh.n_nodes(); ++i) v(i) = f(mesh.x_m.row(i).transpose());
  return v;
}

double integrate_field(const Mesh& mesh, const std::string& region, const Eigen::VectorXd& T, const ScalarField& w) {
  const QuadRuleTri rule = triangle_rule(6);
  double sum = 0.0;
  for (const CellBlock* b : mesh.blocks_of(region)) {
    for (int e = 0; e < b->nodes.rows(); ++e) {
      const Eigen::MatrixXd x_e = mesh.cell_coordinates(*b, e);
      for (std::size_t q = 0; q < rule.w.size(); ++q) {
        const ShapeEval s = shape(b->type, rule.xi[q]);
        const CellGeometry g = cell_geometry(s, x_e);
        double T_h = 0.0;
        for (Eigen::Index a = 0; a < s.N.size(); ++a) T_h += s.N(a) * T(b->nodes(e, a));
        sum += w(g.x) * T_h * g.detJ * rule.w[q];
      }
    }
  }
  return sum;
}

FieldError field_l2_error(const Mesh& mesh, const std::string& region, const Eigen::VectorXd& T,
                          const ScalarField& exact, int n_gauss) {
  const QuadRuleTri rule = triangle_collapsed_gauss(n_gauss);
  double e2 = 0.0, n2 = 0.0, area = 0.0;
  int n_cells = 0;
  for (const CellBlock* b : mesh.blocks_of(region)) {
    for (int e = 0; e < b->nodes.rows(); ++e) {
      const Eigen::MatrixXd x_e = mesh.cell_coordinates(*b, e);
      for (std::size_t q = 0; q < rule.w.size(); ++q) {
        const ShapeEval s = shape(b->type, rule.xi[q]);
        const CellGeometry g = cell_geometry(s, x_e);
        double T_h = 0.0;
        for (Eigen::Index a = 0; a < s.N.size(); ++a) T_h += s.N(a) * T(b->nodes(e, a));
        const double dA = g.detJ * rule.w[q], Tx = exact(g.x);
        e2 += (Tx - T_h) * (Tx - T_h) * dA;
        n2 += Tx * Tx * dA;
        area += dA;
      }
      ++n_cells;
    }
  }
  return {std::sqrt(e2), std::sqrt(n2), std::sqrt(area / n_cells)};
}

std::optional<double> field_at_point(const Mesh& mesh, const std::string& region, const Eigen::VectorXd& T,
                                     const Eigen::Vector2d& x_m) {
  double sum = 0.0;
  int found = 0;
  for (const CellBlock* b : mesh.blocks_of(region)) {
    for (int e = 0; e < b->nodes.rows(); ++e) {
      const auto xi = inverse_map(b->type, mesh.cell_coordinates(*b, e), x_m);
      if (!xi) continue;
      const ShapeEval s = shape(b->type, *xi);
      for (Eigen::Index a = 0; a < s.N.size(); ++a) sum += s.N(a) * T(b->nodes(e, a));
      ++found;
    }
  }
  if (found == 0) return std::nullopt;
  return sum / found;
}

}  // namespace tmfem
