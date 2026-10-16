add_library(tmfem_core STATIC
  core/version.cpp)
target_include_directories(tmfem_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(tmfem_core PUBLIC Eigen3::Eigen)
target_compile_definitions(tmfem_core PRIVATE TMFEM_VERSION_STRING="${PROJECT_VERSION}")
tmfem_set_warnings(tmfem_core)

add_executable(tmfem app/main.cpp)
target_link_libraries(tmfem PRIVATE tmfem_core)
tmfem_set_warnings(tmfem)
