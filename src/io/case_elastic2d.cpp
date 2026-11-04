#include "io/case_elastic2d.hpp"

#include <stdexcept>
#include <string>
#include <vector>

#include "io/csv_writer.hpp"
#include "io/json_file.hpp"
#include "io/units.hpp"
#include "io/vtu_writer.hpp"
#include "mechanics/kirsch_analytic.hpp"
#include "mechanics/stress_recovery.hpp"
#include "mesh/gmsh_reader.hpp"

namespace tmfem {

namespace {

Eigen::Vector2d vector2(const nlohmann::json& j) { return {j.at(0).get<double>(), j.at(1).get<double>()}; }

}  // namespace

PlaneCondition plane_condition(const std::string& s) {
  if (s == "plane_stress") return PlaneCondition::PlaneStress;
  if (s == "plane_strain") return PlaneCondition::PlaneStrain;
  throw std::runtime_error("case: plane must be plane_stress or plane_strain, not '" + s + "'");
}

int displacement_component(const std::string& s) {
  if (s == "x") return 0;
  if (s == "y") return 1;
  throw std::runtime_error("case: component must be x or y, not '" + s + "'");
}

void add_boundary_condition(const nlohmann::json& j, std::vector<DisplacementBC>& displacements,
                            std::vector<TractionBC>& tractions) {
  const auto type = j.at("type").get<std::string>();
  const auto physical = j.at("physical").get<std::string>();
  if (type == "displacement") {
    const int component = displacement_component(j.at("component").get<std::string>());
    const double v = j.at("value_m").get<double>();
    displacements.push_back({physical, component, [v](const Eigen::Vector2d&) { return v; }});
  } else if (type == "traction") {
    const Eigen::Vector2d t = vector2(j.at("traction_MPa")) * kPaPerMPa;
    tractions.push_back({physical, [t](const Eigen::Vector2d&) { return t; }});
  } else if (type == "kirsch_traction") {
    const KirschPlate k{j.at("S_MPa").get<double>() * kPaPerMPa, j.at("a_m").get<double>()};
    const Eigen::Vector2d n = vector2(j.at("normal"));
    tractions.push_back({physical, [k, n](const Eigen::Vector2d& x) {
                           const Eigen::Vector3d s = kirsch_stress_cartesian(k, x);
                           return Eigen::Vector2d(s(0) * n(0) + s(2) * n(1), s(2) * n(0) + s(1) * n(1));
                         }});
  } else {
    throw std::runtime_error("case: unknown boundary condition type '" + type + "'");
  }
}

Elastic2DCase load_elastic2d_case(const std::filesystem::path& json_path, const std::filesystem::path& base_dir) {
  const nlohmann::json j = read_case_json(json_path, "elastic2d");
  Elastic2DCase c{{read_gmsh41(base_dir / j.at("mesh").get<std::string>()),
                   {j.at("material").at("E_GPa").get<double>() * kPaPerGPa, j.at("material").at("nu").get<double>()},
                   plane_condition(j.at("plane").get<std::string>()),
                   j.value("region", std::string("domain")),
                   {},
                   {}},
                  base_dir / j.at("output_vtu").get<std::string>(),
                  base_dir / j.at("output_nodal_csv").get<std::string>()};
  for (const auto& bc : j.at("boundary_conditions")) add_boundary_condition(bc, c.problem.displacements, c.problem.tractions);
  return c;
}

void run_elastic2d_case(const Elastic2DCase& c) {
  const Elastic2DProblem& p = c.problem;
  const Elastic2DSolution sol = solve_elastic2d(p);
  const NodalStressFields f = nodal_stress_fields(p, sol.u_m);  // sigma_zz with the thermal term
  const int n = p.mesh.n_nodes();
  Eigen::MatrixXd u(n, 3);
  for (int i = 0; i < n; ++i) u.row(i) << sol.u_m(2 * i), sol.u_m(2 * i + 1), 0.0;
  const Eigen::MatrixXd s_MPa = f.sigma_Pa / kPaPerMPa;
  const Eigen::MatrixXd szz_MPa = f.sigma_zz_Pa / kPaPerMPa, vm_MPa = f.von_mises_Pa / kPaPerMPa;
  write_vtu(c.output_vtu, p.mesh, p.mesh.blocks_of(p.region),
            {{"displacement_m", u}, {"stress_xx_yy_xy_MPa", s_MPa}, {"stress_zz_MPa", szz_MPa}, {"von_mises_MPa", vm_MPa}});
  write_csv(c.output_nodal_csv,
            {"x_m", "y_m", "ux_m", "uy_m", "sxx_MPa", "syy_MPa", "sxy_MPa", "szz_MPa", "von_mises_MPa"},
            {matrix_column(p.mesh.x_m, 0), matrix_column(p.mesh.x_m, 1), matrix_column(u, 0), matrix_column(u, 1),
             matrix_column(s_MPa, 0), matrix_column(s_MPa, 1), matrix_column(s_MPa, 2), matrix_column(szz_MPa, 0),
             matrix_column(vm_MPa, 0)},
            {"nodal-averaged stresses (plots only); unit thickness"});
}

}  // namespace tmfem
