### Basic settings
cmake_minimum_required (VERSION 2.8.6)
#project(fortranc Fortran C)
#include(FortranCInterface)
#FortranCInterface_HEADER(FC.h MACRO_NAMESPACE "FC_")
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
#add_subdirectory("../contrib/rkf45" "rkf45")
add_subdirectory("../spline" "spline")
add_subdirectory("../contrib/quadpack" "quadpack")
#add_subdirectory("../contrib/odepack" "odepack")
add_subdirectory("../contrib/vode" "vode")

set(CMAKE_Fortran_FLAGS "-g -Wall -Wno-unused-dummy-argument -O2 -I BOOZER_MAGFIE/OBJS -I spline/OBJS -I quadpack/OBJS -I vode/OBJS")

### Define executable
add_executable(${PROJECT_EXE_NAME}
  driftorbit_test.f90
)

add_library("driftorbit"
  common.f90
  driftorbit.f90
  #do_magfie_standalone.f90
  do_magfie_neo.f90
)

add_dependencies(driftorbit boozer_magfie)
#add_dependencies(driftorbit rkf45)
add_dependencies(driftorbit spline)
add_dependencies(driftorbit quadpack)
#add_dependencies(driftorbit odepack)
add_dependencies(driftorbit vode)

### Linking
target_link_libraries(
  ${PROJECT_EXE_NAME}
  driftorbit
  boozer_magfie
  #rkf45
  #odepack
  vode
  spline
  quadpack
  openblas
  pthread
)

### Copy files
configure_file(neo.in neo.in COPYONLY)
configure_file(tok-synch2-n0.bc tok-synch2-n0.bc COPYONLY)


### For make install
# install (TARGETS ${PROJECT_EXE_NAME} DESTINATION /temp/gernot_k/)
# set_property(TARGET ${PROJECT_EXE_NAME} PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)

### Tests
add_test( test_magfie driftorbit_test )
