#include "io/csv_writer.hpp"

#include <cstdio>
#include <fstream>
#include <stdexcept>

namespace tmfem {

void write_csv(const std::filesystem::path& path, const std::vector<std::string>& header,
               const std::vector<std::vector<double>>& columns, const std::vector<std::string>& comments) {
  if (header.size() != columns.size()) throw std::invalid_argument("write_csv: header/column count mismatch");
  for (const auto& c : columns) {
    if (c.size() != columns.front().size()) throw std::invalid_argument("write_csv: columns differ in length");
  }
  if (path.has_parent_path()) std::filesystem::create_directories(path.parent_path());
  std::ofstream out(path);
  if (!out) throw std::runtime_error("write_csv: cannot open " + path.string());
  for (const auto& c : comments) out << "# " << c << '\n';
  for (std::size_t j = 0; j < header.size(); ++j) out << (j ? "," : "") << header[j];
  out << '\n';
  const std::size_t n_rows = columns.empty() ? 0 : columns.front().size();
  char buf[32];
  for (std::size_t i = 0; i < n_rows; ++i) {
    for (std::size_t j = 0; j < columns.size(); ++j) {
      std::snprintf(buf, sizeof buf, "%.17g", columns[j][i]);
      out << (j ? "," : "") << buf;
    }
    out << '\n';
  }
}

}  // namespace tmfem
