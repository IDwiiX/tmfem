#include "io/case_heat1d.hpp"

#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>

#include "io/json_file.hpp"
#include "io/timeseries.hpp"
#include "thermal/theta_integrator.hpp"

namespace tmfem {

namespace {

constexpr double kCelsiusToKelvin = 273.15;

TimeFunction constant(double v) {
  return [v](double) { return v; };
}

// T(t) in K from either "<key>_C" (constant) or "<key>_file" (CSV t_s,value in C).
TimeFunction temperature_input(const nlohmann::json& j, const std::string& key,
                               const std::filesystem::path& base_dir) {
  if (j.contains(key + "_C")) return constant(j.at(key + "_C").get<double>() + kCelsiusToKelvin);
  if (j.contains(key + "_file")) {
    const auto series = std::make_shared<TimeSeries>(
        TimeSeries::read_csv(base_dir / j.at(key + "_file").get<std::string>()));
    return [series](double t_s) { return (*series)(t_s) + kCelsiusToKelvin; };
  }
  throw std::runtime_error("heat1d case: boundary needs '" + key + "_C' or '" + key + "_file'");
}

BoundaryCondition boundary(const nlohmann::json& j, const std::filesystem::path& base_dir) {
  const auto type = j.at("type").get<std::string>();
  if (type == "dirichlet") return DirichletBC{temperature_input(j, "T", base_dir)};
  if (type == "flux") return FluxBC{constant(j.at("q_in_W_per_m2").get<double>())};
  if (type == "convective") {
    return ConvectiveBC{j.at("h_W_per_m2K").get<double>(), temperature_input(j, "T_inf", base_dir)};
  }
  throw std::runtime_error("heat1d case: unknown boundary type '" + type + "'");
}

}  // namespace

Heat1DCase load_heat1d_case(const std::filesystem::path& json_path, const std::filesystem::path& base_dir) {
  const nlohmann::json j = read_json(json_path);
  if (j.at("case_type").get<std::string>() != "heat1d") {
    throw std::runtime_error("load_heat1d_case: case_type is not heat1d in " + json_path.string());
  }
  const auto& mesh = j.at("mesh");
  const auto& mat = j.at("material");
  const auto& time = j.at("time");
  Heat1DCase c;
  c.problem.mesh = graded_interval(mesh.at("length_m").get<double>(), mesh.at("n_elements").get<int>(),
                                   mesh.value("grading_ratio", 1.0));
  c.problem.k = mat.at("k_W_per_mK").get<double>();
  c.problem.rho_c = mat.at("rho_kg_per_m3").get<double>() * mat.at("cp_J_per_kgK").get<double>();
  c.problem.left = boundary(j.at("left"), base_dir);
  c.problem.right = boundary(j.at("right"), base_dir);
  c.T0_K = j.at("initial_T_C").get<double>() + kCelsiusToKelvin;
  c.theta = time.at("theta").get<double>();
  c.dt_s = time.at("dt_s").get<double>();
  c.t_end_s = time.at("t_end_s").get<double>();
  const auto mass = time.at("mass").get<std::string>();
  if (mass != "consistent" && mass != "lumped") throw std::runtime_error("heat1d case: mass must be consistent or lumped");
  c.mass = mass == "lumped" ? MassType::Lumped : MassType::Consistent;
  c.probes_m = j.at("output").at("probes_m").get<std::vector<double>>();
  c.output_csv = base_dir / j.at("output").at("csv").get<std::string>();
  return c;
}

ProbeHistory run_heat1d_case(const Heat1DCase& c) {
  const int n_steps = static_cast<int>(std::lround(c.t_end_s / c.dt_s));
  if (std::abs(n_steps * c.dt_s - c.t_end_s) > 1e-9 * c.t_end_s) {
    throw std::runtime_error("run_heat1d_case: t_end must be a multiple of dt");
  }
  const ThetaIntegrator integrator(build_operators(c.problem, c.mass), c.theta, c.dt_s);
  ProbeHistory h;
  h.T_K.resize(c.probes_m.size());
  const Eigen::VectorXd T0 = Eigen::VectorXd::Constant(c.problem.mesh.n_nodes(), c.T0_K);
  march(integrator, T0, 0.0, n_steps, [&](double t_s, const Eigen::VectorXd& T) {
    h.t_s.push_back(t_s);
    const std::vector<double> nodal(T.data(), T.data() + T.size());
    for (std::size_t p = 0; p < c.probes_m.size(); ++p) {
      h.T_K[p].push_back(interpolate(c.problem.mesh, nodal, c.probes_m[p]));
    }
  });
  return h;
}

}  // namespace tmfem
