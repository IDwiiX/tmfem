#include <iostream>
#include <string>

#include "core/version.hpp"

// Entry point. Case dispatch (`tmfem cases/<case>.json`) arrives with phase 1.
int main(int argc, char** argv) {
  if (argc == 2 && std::string(argv[1]) == "--version") {
    std::cout << "tmfem " << tmfem::version() << '\n';
    return 0;
  }
  std::cerr << "usage: tmfem --version\n"
               "(case files are not supported yet)\n";
  return 1;
}
