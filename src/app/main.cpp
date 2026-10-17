#include <exception>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#include "core/version.hpp"
#include "io/case_elastic2d.hpp"
#include "io/case_elastoplastic2d.hpp"
#include "io/case_heat1d.hpp"
#include "io/csv_writer.hpp"
#include "io/json_file.hpp"

namespace {

// Runs a heat1d case and writes the probe temperatures (K) to the case's output CSV.
int run_heat1d(const std::filesystem::path& case_path) {
  const tmfem::Heat1DCase c = tmfem::load_heat1d_case(case_path, std::filesystem::current_path());
  const tmfem::ProbeHistory h = tmfem::run_heat1d_case(c);
  std::vector<std::string> header{"t_s"};
  std::vector<std::vector<double>> columns{h.t_s};
  for (std::size_t p = 0; p < c.probes_m.size(); ++p) {
    header.push_back("T_K_at_x" + std::to_string(p));
    columns.push_back(h.T_K[p]);
  }
  std::vector<std::string> comments{"case " + case_path.string()};
  for (std::size_t p = 0; p < c.probes_m.size(); ++p) {
    comments.push_back("x" + std::to_string(p) + " = " + std::to_string(c.probes_m[p]) + " m");
  }
  tmfem::write_csv(c.output_csv, header, columns, comments);
  std::cout << "wrote " << c.output_csv.string() << " (" << h.t_s.size() << " rows)\n";
  return 0;
}

// Runs an elastic2d case: VTU and nodal CSV as named in the case file.
int run_elastic2d(const std::filesystem::path& case_path) {
  const tmfem::Elastic2DCase c = tmfem::load_elastic2d_case(case_path, std::filesystem::current_path());
  tmfem::run_elastic2d_case(c);
  std::cout << "wrote " << c.output_vtu.string() << " and " << c.output_nodal_csv.string() << '\n';
  return 0;
}

// Runs an elastoplastic2d case: VTU and quadrature-point CSV as named in the case file.
int run_elastoplastic2d(const std::filesystem::path& case_path) {
  const tmfem::Elastoplastic2DCase c = tmfem::load_elastoplastic2d_case(case_path, std::filesystem::current_path());
  tmfem::run_elastoplastic2d_case(c);
  std::cout << "wrote " << c.output_vtu.string() << " and " << c.output_gauss_csv.string() << '\n';
  return 0;
}

}  // namespace

// Entry point: `tmfem <case.json>` dispatches on the case's "case_type".
int main(int argc, char** argv) {
  if (argc == 2 && std::string(argv[1]) == "--version") {
    std::cout << "tmfem " << tmfem::version() << '\n';
    return 0;
  }
  if (argc != 2) {
    std::cerr << "usage: tmfem <case.json> | tmfem --version\n";
    return 1;
  }
  try {
    const std::filesystem::path case_path(argv[1]);
    const auto type = tmfem::read_json(case_path).at("case_type").get<std::string>();
    if (type == "heat1d") return run_heat1d(case_path);
    if (type == "elastic2d") return run_elastic2d(case_path);
    if (type == "elastoplastic2d") return run_elastoplastic2d(case_path);
    std::cerr << "tmfem: unsupported case_type '" << type << "'\n";
    return 1;
  } catch (const std::exception& e) {
    std::cerr << "tmfem: " << e.what() << '\n';
    return 1;
  }
}
