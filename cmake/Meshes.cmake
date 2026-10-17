# Meshes are not committed: every cases/*/*.geo is meshed by scripts/make_meshes.sh (Gmsh 4.12.1, which the
# script checks) into cases/*, where the tests and case files read them. The `meshes` target is part of ALL
# and every test executable depends on it. It reruns when the script or any .geo file changes; to force a
# rerun, delete <build>/meshes.stamp or run the script by hand.
find_program(TMFEM_GMSH gmsh)
if(NOT TMFEM_GMSH)
  message(FATAL_ERROR "Gmsh 4.12.1 is required to generate the meshes (Ubuntu 24.04: sudo apt install gmsh)")
endif()

file(GLOB TMFEM_GEO_FILES CONFIGURE_DEPENDS "${PROJECT_SOURCE_DIR}/cases/*/*.geo")
set(TMFEM_MESH_STAMP "${PROJECT_BINARY_DIR}/meshes.stamp")
add_custom_command(
  OUTPUT "${TMFEM_MESH_STAMP}"
  COMMAND ${CMAKE_COMMAND} -E env "GMSH=${TMFEM_GMSH}" bash "${PROJECT_SOURCE_DIR}/scripts/make_meshes.sh"
  COMMAND ${CMAKE_COMMAND} -E touch "${TMFEM_MESH_STAMP}"
  DEPENDS "${PROJECT_SOURCE_DIR}/scripts/make_meshes.sh" ${TMFEM_GEO_FILES}
  COMMENT "Generating meshes from cases/*/*.geo (scripts/make_meshes.sh)"
  VERBATIM)
add_custom_target(meshes ALL DEPENDS "${TMFEM_MESH_STAMP}")
