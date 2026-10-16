# One test executable per phase. Every Catch2 TEST_CASE becomes a CTest test labelled with
# its phase (`ctest -L phase2`) and with each of its tags (`ctest -L P2-c`).
# SKIP() in Catch2 exits with code 4, which CTest reports as skipped rather than passed.
function(tmfem_add_phase_tests phase)
  cmake_parse_arguments(ARG "" "" "SOURCES" ${ARGN})
  set(exe tests_${phase})
  add_executable(${exe} ${ARG_SOURCES})
  target_link_libraries(${exe} PRIVATE tmfem_core Catch2::Catch2WithMain)
  target_include_directories(${exe} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  target_compile_definitions(${exe} PRIVATE TMFEM_SOURCE_DIR="${PROJECT_SOURCE_DIR}")
  tmfem_set_warnings(${exe})
  catch_discover_tests(${exe}
    ADD_TAGS_AS_LABELS
    SKIP_RETURN_CODE 4
    PROPERTIES LABELS ${phase})
endfunction()

tmfem_add_phase_tests(phase0 SOURCES phase0/smoke.cpp)
