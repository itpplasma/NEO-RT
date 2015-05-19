### Basic settings
cmake_minimum_required (VERSION 2.6)
enable_testing()

### Define the project
project (driftorbit)
enable_language(Fortran)

### Specify paths
set(CMAKE_BINARY_DIR ./OBJS)
set(EXECUTABLE_OUTPUT_PATH .)
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

### Define the name of the executable
set(PROJECT_EXE_NAME "driftorbit_test")                  # Define the name

### Set compiler flags
set(CMAKE_Fortran_MODULE_DIRECTORY ./OBJS)

### External modules
add_subdirectory("../BOOZER_MAGFIE" "BOOZER_MAGFIE")
add_subdirectory("../contrib/rkf45" "rkf45")
add_subdirectory("../spline" "spline")

set(CMAKE_Fortran_FLAGS "-g -Wall -Wno-unused-dummy-argument -O1 -I BOOZER_MAGFIE/OBJS -I spline/OBJS -I rkf45/OBJS")

### Define executable
add_executable(${PROJECT_EXE_NAME}
  driftorbit_test.f90
)

add_library("driftorbit"
  common.f90
  driftorbit.f90
  do_magfie_neo.f90
)

add_dependencies(driftorbit boozer_magfie)
add_dependencies(driftorbit rkf45)
add_dependencies(driftorbit spline)

### Linking
target_link_libraries(
  ${PROJECT_EXE_NAME}
  driftorbit
  lapack
  blas
  boozer_magfie
  rkf45
  spline
)

### For make install
# install (TARGETS ${PROJECT_EXE_NAME} DESTINATION /temp/gernot_k/)
# set_property(TARGET ${PROJECT_EXE_NAME} PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)

### Tests
add_test( test_magfie driftorbit_test )
