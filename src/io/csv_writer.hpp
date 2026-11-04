#pragma once

#include <filesystem>
#include <string>
#include <vector>

#include <Eigen/Dense>

namespace tmfem {

// Writes equal-length numeric columns with a header row, full precision (%.17g).
// Creates parent directories. Lines in `comments` are written first, prefixed with '#'.
void write_csv(const std::filesystem::path& path, const std::vector<std::string>& header,
               const std::vector<std::vector<double>>& columns,
               const std::vector<std::string>& comments = {});

// One column of a matrix, as write_csv takes it.
std::vector<double> matrix_column(const Eigen::MatrixXd& m, Eigen::Index c);

}  // namespace tmfem
