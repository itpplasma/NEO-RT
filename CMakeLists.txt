### Basic settings
cmake_minimum_required (VERSION 2.8.6)
#project(fortranc Fortran C)
#include(FortranCInterface)
#FortranCInterface_HEADER(FC.h MACRO_NAMESPACE "FC_")
enable_testing()

### Define the project
project (driftorbit)
enable_language(Fortran)

### Define paths to external libraries (load external file)
include(${CMAKE_SOURCE_DIR}/ProjectConfig.cmake.in)

### Specify paths
set(CMAKE_BINARY_DIR ${PROJECT_BINARY_DIR}/OBJS)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

### Define the name of the executable
set(PROJECT_EXE_NAME "driftorbit_test")                  # Define the name

### Set compiler flags
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/OBJS)

### Find libraries
set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
#find_library(SuperLU_lib superlu_4.1 ${SUPERLU_DIR}lib NO_DEFAULT_PATH)
#find_library(Umfpack_lib umfpack ${SUITESPARSE_DIR}lib NO_DEFAULT_PATH)
#find_library(Amd_lib amd ${SUITESPARSE_DIR}lib NO_DEFAULT_PATH)
#find_library(Cholmod_lib cholmod ${SUITESPARSE_DIR}lib NO_DEFAULT_PATH)
#find_library(Colamd_lib colamd ${SUITESPARSE_DIR}lib NO_DEFAULT_PATH)
#find_library(Camd_lib camd ${SUITESPARSE_DIR}lib NO_DEFAULT_PATH)
#find_library(Metis_lib metis ${SUITESPARSE_DIR}lib ${METIS_PATH} NO_DEFAULT_PATH)
#find_library(Ccolamd_lib ccolamd ${SUITESPARSE_DIR}lib NO_DEFAULT_PATH)
#find_library(SuiteSparseConfig_lib suitesparseconfig ${SUITESPARSE_DIR}lib NO_DEFAULT_PATH)

### Source files (load external file)
include (CMakeSources.txt)

### External modules
#add_subdirectory("../BOOZER_MAGFIE" "BOOZER_MAGFIE")
#add_subdirectory("../contrib/rkf45" "rkf45")
add_subdirectory("../spline" "spline")
add_subdirectory("../contrib/quadpack" "quadpack")
#add_subdirectory("../contrib/odepack" "odepack")
add_subdirectory("../contrib/vode" "vode")

#set(CMAKE_Fortran_FLAGS "-O3 -g -traceback -implicitnone -I spline/OBJS -I vode/OBJS") # GNU
#set(CMAKE_Fortran_FLAGS " -O2 -g -I BOOZER_MAGFIE/OBJS -I spline/OBJS -I vode/OBJS") # INTEL
#set(CMAKE_Fortran_FLAGS "-mkl -O3 -g -traceback -implicitnone -I BOOZER_MAGFIE/OBJS -I spline/OBJS -I vode/OBJS")
#set(CMAKE_Fortran_FLAGS "-mkl -O3 -g -traceback -implicitnone -I spline/OBJS -I vode/OBJS")

IF(WIN32)
  set(CMAKE_Fortran_FLAGS "/Qmkl /O3")
ELSE()
#  set(CMAKE_Fortran_FLAGS "-mkl -O3 -init=snan,arrays -g -traceback -implicitnone -I spline/OBJS -I vode/OBJS")
   #set(CMAKE_Fortran_FLAGS "-mkl -O3 -implicitnone -I spline/OBJS -I vode/OBJS")
   set(CMAKE_Fortran_FLAGS "-O3") #GNU
ENDIF()

### Define executable
add_executable(${PROJECT_EXE_NAME}
  driftorbit_test.f90
)

add_executable(orbit_test
  orbit_test.f90
)


add_library("driftorbit"
  orbit.f90
  attenuation_factor.f90
  collis_nbi.f90
  polylag_3.f
  common.f90
  driftorbit.f90
  do_magfie_standalone.f90
#  do_magfie_neo.f90
#  ${SUITESPARSE_SRC_FILES}
#  ${SUPERLU_SRC_FILES}
)

#add_dependencies(driftorbit boozer_magfie)
#add_dependencies(driftorbit rkf45)
add_dependencies(driftorbit spline)
#add_dependencies(driftorbit quadpack)
#add_dependencies(driftorbit odepack)
add_dependencies(driftorbit vode)

### Linking
target_link_libraries(
  ${PROJECT_EXE_NAME}
  driftorbit
#  boozer_magfie
#  ${SuperLU_lib}
#  ${Umfpack_lib} ${Amd_lib} ${Cholmod_lib} ${Colamd_lib} ${Camd_lib} ${Metis_lib} ${Ccolamd_lib}
#  ${SuiteSparseConfig_lib}
  #rkf45
  #odepack
  vode
  spline
  quadpack
  openblas
#  blas
#  lapack
#  pthread
)

### Copy files
#configure_file(neo.in neo.in COPYONLY)
#configure_file(tok-synch2-n0.bc tok-synch2-n0.bc COPYONLY)


### For make install
# install (TARGETS ${PROJECT_EXE_NAME} DESTINATION /temp/gernot_k/)
# set_property(TARGET ${PROJECT_EXE_NAME} PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)

### Tests
add_test( test_magfie driftorbit_test )
