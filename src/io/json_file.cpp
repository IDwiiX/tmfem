#include "io/json_file.hpp"

#include <fstream>
#include <stdexcept>

namespace tmfem {

nlohmann::json read_json(const std::filesystem::path& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("read_json: cannot open " + path.string());
  try {
    return nlohmann::json::parse(in);
  } catch (const nlohmann::json::exception& e) {
    throw std::runtime_error("read_json: " + path.string() + ": " + e.what());
  }
}

nlohmann::json read_case_json(const std::filesystem::path& path, const std::string& expected_type) {
  nlohmann::json j = read_json(path);
  if (j.at("case_type").get<std::string>() != expected_type) {
    throw std::runtime_error("read_case_json: case_type is not " + expected_type + " in " + path.string());
  }
  return j;
}

}  // namespace tmfem
