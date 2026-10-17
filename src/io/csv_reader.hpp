#pragma once

#include <filesystem>
#include <vector>

namespace tmfem {

// Numeric CSV: lines starting with '#' are comments, the first other line is a header, and every
// following line holds n_cols numbers. Returns the columns. Throws with the path on a bad row.
std::vector<std::vector<double>> read_numeric_csv(const std::filesystem::path& path, int n_cols);

}  // namespace tmfem
