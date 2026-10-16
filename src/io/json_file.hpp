#pragma once

#include <filesystem>

#include <nlohmann/json.hpp>

namespace tmfem {

// Parses a JSON file; throws with the path in the message on failure.
nlohmann::json read_json(const std::filesystem::path& path);

}  // namespace tmfem
