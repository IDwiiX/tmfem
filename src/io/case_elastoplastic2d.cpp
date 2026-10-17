#include "io/case_elastoplastic2d.hpp"

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "io/case_elastic2d.hpp"
#include "io/csv_reader.hpp"
#include "io/csv_writer.hpp"
#include "io/json_file.hpp"
#include "io/vtu_writer.hpp"
#include "mechanics/hardening.hpp"
#include "mechanics/j2_plasticity.hpp"
#include "mesh/gmsh_reader.hpp"

namespace tmfem {

namespace {

constexpr double kPaPerMPa = 1e6;
constexpr double kPaPerGPa = 1e9;

std::shared_ptr<const HardeningLaw> hardening_law(const nlohmann::json& j, double E_Pa,
                                                  const std::filesystem::path& base_dir) {
  const auto type = j.at("type").get<std::string>();
  if (type == "tabulated") {
    auto cols = read_numeric_csv(base_dir / j.at("csv").get<std::string>(), 2);
    for (double& s : cols[1]) s *= kPaPerMPa;
    return std::make_shared<TabulatedHardening>(
        TabulatedHardening::from_true_curve(cols[0], cols[1], E_Pa, j.at("sigma_start_MPa").get<double>() * kPaPerMPa));
  }
  if (type == "ludwik") {
    return std::make_shared<LudwikHardening>(j.at("sigma_y0_MPa").get<double>() * kPaPerMPa,
                                             j.at("K_MPa").get<double>() * kPaPerMPa, j.at("m").get<double>());
  }
  throw std::runtime_error("elastoplastic2d case: hardening type must be tabulated or ludwik, not '" + type + "'");
}

}  // namespace

Elastoplastic2DCase load_elastoplastic2d_case(const std::filesystem::path& json_path,
                                              const std::filesystem::path& base_dir) {
  const nlohmann::json j = read_json(json_path);
  if (j.at("case_type").get<std::string>() != "elastoplastic2d") {
    throw std::runtime_error("load_elastoplastic2d_case: case_type is not elastoplastic2d in " + json_path.string());
  }
  const ElasticMaterial m{j.at("material").at("E_GPa").get<double>() * kPaPerGPa, j.at("material").at("nu").get<double>()};
  Elastoplastic2DCase c{{read_gmsh41(base_dir / j.at("mesh").get<std::string>()),
                         std::make_shared<J2Plasticity>(m, hardening_law(j.at("hardening"), m.E_Pa, base_dir)),
                         plane_condition(j.at("plane").get<std::string>()),
                         j.value("region", std::string("domain")),
                         {},
                         {}},
                        j.at("load_steps").get<int>(),
                        base_dir / j.at("output_vtu").get<std::string>(),
                        base_dir / j.at("output_gauss_csv").get<std::string>()};
  if (c.load_steps < 1) throw std::runtime_error("elastoplastic2d case: load_steps must be >= 1");
  for (const auto& bc : j.at("boundary_conditions")) {
    add_boundary_condition(bc, c.problem.displacements, c.problem.tractions);
  }
  return c;
}

void run_elastoplastic2d_case(const Elastoplastic2DCase& c) {
  ElastoplasticSolver solver(c.problem);
  for (int k = 1; k <= c.load_steps; ++k) solver.step_to(static_cast<double>(k) / c.load_steps);

  // Quadrature-point CSV, and cell means of alpha and von Mises for the VTU.
  const auto& qp = solver.quad_points();
  std::vector<double> x, y, alpha, vm, sxx, syy, sxy, szz;
  std::vector<double> cell_alpha, cell_vm, cell_n;
  for (const QuadPoint& q : qp) {
    const double q_eq = equivalent_stress(Voigt4(q.stress_Pa(0), q.stress_Pa(1), q.sigma_zz_Pa, q.stress_Pa(2)));
    x.push_back(q.x_m(0));
    y.push_back(q.x_m(1));
    alpha.push_back(q.state.alpha);
    vm.push_back(q_eq / kPaPerMPa);
    sxx.push_back(q.stress_Pa(0) / kPaPerMPa);
    syy.push_back(q.stress_Pa(1) / kPaPerMPa);
    sxy.push_back(q.stress_Pa(2) / kPaPerMPa);
    szz.push_back(q.sigma_zz_Pa / kPaPerMPa);
    const auto e = static_cast<std::size_t>(q.cell);
    if (e >= cell_n.size()) {
      cell_alpha.resize(e + 1, 0.0);
      cell_vm.resize(e + 1, 0.0);
      cell_n.resize(e + 1, 0.0);
    }
    cell_alpha[e] += q.state.alpha;
    cell_vm[e] += q_eq / kPaPerMPa;
    cell_n[e] += 1.0;
  }
  write_csv(c.output_gauss_csv, {"x_m", "y_m", "alpha", "von_mises_MPa", "sxx_MPa", "syy_MPa", "sxy_MPa", "szz_MPa"},
            {x, y, alpha, vm, sxx, syy, sxy, szz}, {"quadrature-point values at load factor 1; unit thickness"});

  // Node values: mean over the cells sharing the node of each cell's mean (cells in solver order).
  const Mesh& mesh = c.problem.mesh;
  const int n = mesh.n_nodes();
  Eigen::MatrixXd u(n, 3), node_alpha = Eigen::MatrixXd::Zero(n, 1), node_vm = Eigen::MatrixXd::Zero(n, 1);
  Eigen::VectorXd count = Eigen::VectorXd::Zero(n);
  std::size_t e = 0;
  for (const CellBlock* b : mesh.blocks_of(c.problem.region)) {
    for (int k = 0; k < b->nodes.rows(); ++k, ++e) {
      for (Eigen::Index a = 0; a < b->nodes.cols(); ++a) {
        const int node = b->nodes(k, a);
        node_alpha(node, 0) += cell_alpha[e] / cell_n[e];
        node_vm(node, 0) += cell_vm[e] / cell_n[e];
        count(node) += 1.0;
      }
    }
  }
  for (int i = 0; i < n; ++i) {
    u.row(i) << solver.displacement_m()(2 * i), solver.displacement_m()(2 * i + 1), 0.0;
    if (count(i) > 0.0) {
      node_alpha(i, 0) /= count(i);
      node_vm(i, 0) /= count(i);
    }
  }
  write_vtu(c.output_vtu, mesh, mesh.blocks_of(c.problem.region),
            {{"displacement_m", u}, {"equivalent_plastic_strain", node_alpha}, {"von_mises_MPa", node_vm}});
}

}  // namespace tmfem
