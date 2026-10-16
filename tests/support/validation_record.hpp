#pragma once

#include <cctype>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <string>

#include "support/paths.hpp"

namespace tmfem::test {

// Writes one validated number to results/validation/<id>__<quantity>.csv so the README table
// can be built from what the tests actually measured. Columns: id, quantity, reference, tmfem,
// relative error, tolerance, test.
inline void record(const std::string& id, const std::string& quantity, double reference, double value,
                   const std::string& tolerance, const std::string& test_name) {
  std::string slug = quantity;
  for (char& ch : slug) {
    if (!std::isalnum(static_cast<unsigned char>(ch))) ch = '_';
  }
  const auto dir = results_dir() / "validation";
  std::filesystem::create_directories(dir);
  std::ofstream out(dir / (id + "__" + slug + ".csv"));
  const double rel = reference != 0.0 ? (value - reference) / std::abs(reference) : value - reference;
  char buf[256];
  std::snprintf(buf, sizeof buf, "%.10g,%.10g,%.3e", reference, value, rel);
  out << "id,quantity,reference,tmfem,rel_error,tolerance,test\n"
      << id << ",\"" << quantity << "\"," << buf << ",\"" << tolerance << "\",\"" << test_name << "\"\n";
}

}  // namespace tmfem::test
