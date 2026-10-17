# One test executable per phase. Every Catch2 TEST_CASE becomes a CTest test labelled with
# its phase (`ctest -L phase2`) and with each of its tags (`ctest -L P2-c`).
# SKIP() in Catch2 exits with code 4, which CTest reports as skipped rather than passed.
function(tmfem_add_phase_tests phase)
  cmake_parse_arguments(ARG "" "" "SOURCES" ${ARGN})
  set(exe tests_${phase})
  add_executable(${exe} ${ARG_SOURCES})
  target_link_libraries(${exe} PRIVATE tmfem_core Catch2::Catch2WithMain)
  add_dependencies(${exe} meshes)  # cmake/Meshes.cmake
  target_include_directories(${exe} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  target_compile_definitions(${exe} PRIVATE TMFEM_SOURCE_DIR="${PROJECT_SOURCE_DIR}")
  tmfem_set_warnings(${exe})
  catch_discover_tests(${exe}
    ADD_TAGS_AS_LABELS
    SKIP_RETURN_CODE 4
    PROPERTIES LABELS ${phase})
endfunction()

tmfem_add_phase_tests(phase0 SOURCES phase0/smoke.cpp)
tmfem_add_phase_tests(phase1 SOURCES
  phase1/heat1d_units.cpp
  phase1/diffusion_exercises.cpp
  phase1/carburizing.cpp
  phase1/thin_film.cpp
  phase1/steel_plate.cpp
  phase1/heat1d_convergence.cpp)
tmfem_add_phase_tests(phase2 SOURCES
  phase2/gmsh_reader.cpp
  phase2/shape_functions.cpp
  phase2/elasticity_units.cpp
  phase2/kirsch_analytic.cpp
  phase2/vtu_writer.cpp
  phase2/patch.cpp
  phase2/bar_elastic.cpp
  phase2/kirsch.cpp
  phase2/kirsch_convergence.cpp)
tmfem_add_phase_tests(phase3 SOURCES
  phase3/hardening_units.cpp
  phase3/return_mapping.cpp
  phase3/elastoplastic_units.cpp
  phase3/tensile_fit.cpp
  phase3/material_point.cpp
  phase3/plastic_bar.cpp
  phase3/plane_stress_units.cpp
  phase3/newton_rate.cpp
  phase3/kirsch_plastic.cpp)
tmfem_add_phase_tests(phase4 SOURCES
  phase4/thermal_strain_units.cpp
  phase4/heat2d_units.cpp
  phase4/thermal_output.cpp
  phase4/thermal_bar.cpp
  phase4/thermal_free_body.cpp
  phase4/thermal_disk.cpp
  phase4/heat2d_exactness.cpp
  phase4/heat2d_gaussian.cpp
  phase4/beam_block_demo.cpp)
