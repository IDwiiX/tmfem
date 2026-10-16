add_library(tmfem_core STATIC
  core/version.cpp
  mesh/interval_mesh.cpp
  fem/quadrature.cpp
  fem/line2.cpp
  fem/dirichlet.cpp
  thermal/heat1d.cpp
  thermal/theta_integrator.cpp
  thermal/diffusion_analytic.cpp
  thermal/arrhenius.cpp
  thermal/lumped_plate.cpp
  io/timeseries.cpp
  io/csv_writer.cpp
  io/json_file.cpp
  io/case_heat1d.cpp)
target_include_directories(tmfem_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(tmfem_core PUBLIC Eigen3::Eigen nlohmann_json::nlohmann_json)
target_compile_definitions(tmfem_core PRIVATE TMFEM_VERSION_STRING="${PROJECT_VERSION}")
tmfem_set_warnings(tmfem_core)

add_executable(tmfem app/main.cpp)
target_link_libraries(tmfem PRIVATE tmfem_core)
tmfem_set_warnings(tmfem)
