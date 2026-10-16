cmake_minimum_required(VERSION 3.16)

project(tmfem
  VERSION 0.1.0
  DESCRIPTION "Validated 2D thermo-mechanical finite-element solver"
  LANGUAGES CXX)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(TMFEM_BUILD_TESTS "Build the validation test suite" ON)
option(TMFEM_WARNINGS_AS_ERRORS "Treat compiler warnings as errors (CI)" OFF)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(Warnings)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
include(Json)

add_subdirectory(src)

if(TMFEM_BUILD_TESTS)
  enable_testing()
  include(Catch2)
  add_subdirectory(tests)
endif()
