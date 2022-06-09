### Basic settings
cmake_minimum_required (VERSION 2.8.12)
enable_testing()

### Define the project
project (driftorbit)
enable_language(Fortran)

### Define paths to external libraries (load external file)
include(${CMAKE_SOURCE_DIR}/ProjectConfig.cmake.in)

### Specify paths
set(CMAKE_BINARY_DIR ${PROJECT_BINARY_DIR}/OBJS)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR})

### Define the name of the executable
set(PROJECT_EXE_NAME "neo_rt.x")                  # Define the name

### Set compiler flags
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/OBJS)

### Find libraries for NEO-2 field routines
#find_package(SuiteSparse REQUIRED)
#find_package(SuperLU REQUIRED)

### Source files (load external file)
include (CMakeSources.txt)

### External modules
#add_subdirectory("../BOOZER_MAGFIE" "BOOZER_MAGFIE")
add_subdirectory("../spline" "spline")
add_subdirectory("../contrib/quadpack" "quadpack")
#add_subdirectory("../contrib/vode" "vode")
add_subdirectory(POTATO/LATEST)  # contains also vode

add_compile_options(-cpp -g -fbacktrace -Wfatal-errors)
add_compile_options(-O3 -march=native -mtune=native)

### Define executable
add_executable(${PROJECT_EXE_NAME}
  main.f90
)

add_executable(orbit_test
  orbit_test.f90
)

add_executable(test_frequencies
  test_frequencies.f90
)


add_library("neo_rt" SHARED
  orbit.f90
  attenuation_factor.f90
  collis_nbi.f90
  polylag_3.f
  common.f90
  driftorbit.f90
  do_magfie_standalone.f90
# For NEO-2 field routines
#  do_magfie_neo.f90
#  ${SUITESPARSE_SRC_FILES}
#  ${SUPERLU_SRC_FILES}
)

#add_dependencies(neo_rt boozer_magfie)
add_dependencies(neo_rt spline)
add_dependencies(neo_rt quadpack)
add_dependencies(neo_rt vode)

find_package(BLAS)
find_package(LAPACK)

target_link_libraries(
  neo_rt
  vode
  spline
  quadpack
  ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES}
)

### Linking
target_link_libraries(
  ${PROJECT_EXE_NAME}
  neo_rt
#  boozer_magfie
#  ${SuperLU_lib}
#  ${Umfpack_lib} ${Amd_lib} ${Cholmod_lib} ${Colamd_lib} ${Camd_lib} ${Metis_lib} ${Ccolamd_lib}
#  ${SuiteSparseConfig_lib}
#  pthread
)


target_link_libraries(
  test_frequencies
  neo_rt
)

### Tests
add_test( test_magfie neo_rt )
