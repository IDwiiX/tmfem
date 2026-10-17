#include "io/csv_reader.hpp"

#include <algorithm>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace tmfem {

std::vector<std::vector<double>> read_numeric_csv(const std::filesystem::path& path, int n_cols) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("read_numeric_csv: cannot open " + path.string());
  std::vector<std::vector<double>> cols(static_cast<std::size_t>(n_cols));
  std::string line;
  bool header_seen = false;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    if (!header_seen) {
      header_seen = true;
      continue;
    }
    std::replace(line.begin(), line.end(), ',', ' ');
    std::istringstream row(line);
    for (auto& c : cols) {
      double v = 0.0;
      if (!(row >> v)) throw std::runtime_error("read_numeric_csv: bad row '" + line + "' in " + path.string());
      c.push_back(v);
    }
  }
  return cols;
}

}  // namespace tmfem
