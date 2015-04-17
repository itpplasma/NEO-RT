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
set(CMAKE_Fortran_FLAGS "-I BOOZER_MAGFIE/OBJS")

add_subdirectory("../contrib/odepack" "odepack")

### Define executable
add_executable(${PROJECT_EXE_NAME}
  driftorbit_test.f90
)

add_library("driftorbit"
  common.f90
  driftorbit.f90
)

add_dependencies(driftorbit boozer_magfie)
add_dependencies(driftorbit odepack)

### Linking
target_link_libraries(
  ${PROJECT_EXE_NAME}
  driftorbit
  lapack
  blas
  boozer_magfie
  odepack
)

### For make install
# install (TARGETS ${PROJECT_EXE_NAME} DESTINATION /temp/gernot_k/)
# set_property(TARGET ${PROJECT_EXE_NAME} PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)

### Tests
add_test( test_magfie driftorbit_test )
