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

}  // namespace tmfem
