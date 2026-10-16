#pragma once

#include <filesystem>

namespace tmfem::test {

// Repository root, from the TMFEM_SOURCE_DIR compile definition.
inline std::filesystem::path source_dir() { return TMFEM_SOURCE_DIR; }

// Generated data (gitignored).
inline std::filesystem::path results_dir() { return source_dir() / "results"; }

}  // namespace tmfem::test
