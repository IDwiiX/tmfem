#pragma once

#include <filesystem>
#include <string>
#include <vector>

namespace tmfem {

// Writes equal-length numeric columns with a header row, full precision (%.17g).
// Creates parent directories. Lines in `comments` are written first, prefixed with '#'.
void write_csv(const std::filesystem::path& path, const std::vector<std::string>& header,
               const std::vector<std::vector<double>>& columns,
               const std::vector<std::string>& comments = {});

}  // namespace tmfem
