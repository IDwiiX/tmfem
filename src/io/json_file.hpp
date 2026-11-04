#pragma once

#include <filesystem>
#include <string>

#include <nlohmann/json.hpp>

namespace tmfem {

// Parses a JSON file; throws with the path in the message on failure.
nlohmann::json read_json(const std::filesystem::path& path);

// read_json for a case file; also throws unless its "case_type" is the expected one.
nlohmann::json read_case_json(const std::filesystem::path& path, const std::string& expected_type);

}  // namespace tmfem
