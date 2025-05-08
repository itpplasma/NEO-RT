### Basic settings
cmake_minimum_required (VERSION 3.10)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")
enable_testing()

### Define the project
project (driftorbit)
enable_language(Fortran C)
set(PROJECT_EXE_NAME "neo_rt.x")

### Specify paths
set(CMAKE_BINARY_DIR ${PROJECT_BINARY_DIR}/OBJS)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR})

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

add_compile_options(-fPIC -g)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fbacktrace -fopenmp -fmax-errors=1")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wimplicit-interface -Wimplicit-procedure")
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  add_compile_options(-O0 -Wextra)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fbounds-check -finit-real=nan -fcheck=all")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffpe-trap=zero,overflow,invalid,underflow")
  add_compile_options()
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
  add_compile_options(-O2)
  if(APPLE)
    add_compile_options(-march=native)
  else()
    add_compile_options(-march=x86-64-v2 -mtune=generic)
  endif()
elseif(CMAKE_BUILD_TYPE STREQUAL "Fast")
  add_compile_options(-O3 -ffast-math)
  if(APPLE)
    add_compile_options(-march=native)
  else()
    add_compile_options(-march=native -mtune=native)
  endif()
endif()

include(Util)

find_or_fetch(NEO-2)
find_or_fetch(spline)
find_or_fetch(vode)

### NetCDF
find_program(NF_CONFIG "nf-config")

if (NF_CONFIG)
execute_process(COMMAND nf-config --includedir
                OUTPUT_VARIABLE NETCDFINCLUDE_DIR)
execute_process(COMMAND nc-config --libdir
				OUTPUT_VARIABLE NETCDFLIB_DIR)
execute_process(COMMAND nf-config --flibs
                OUTPUT_VARIABLE NETCDF_FLIBS)
else()
message(SEND_ERROR "nf-config not found. Please install libnetcdff-dev")
endif()

string(STRIP ${NETCDFINCLUDE_DIR} NETCDFINCLUDE_DIR)
string(STRIP ${NETCDFLIB_DIR} NETCDFLIB_DIR)
string(STRIP ${NETCDF_FLIBS} NETCDF_FLIBS)

message(STATUS "NetCDF include path: " ${NETCDFINCLUDE_DIR})
message(STATUS "NetCDF lib path: " ${NETCDFLIB_DIR})
message(STATUS "NetCDF Fortran libs: " ${NETCDF_FLIBS})

# Replace space by semicolon in the Fortran libs
string(REPLACE " " ";" NETCDF_FLIBS ${NETCDF_FLIBS})

include_directories(${NETCDFINCLUDE_DIR})
link_directories(${NETCDFLIB_DIR})
add_link_options(${NETCDF_FLIBS})

add_executable(${PROJECT_EXE_NAME}
  src/main.f90
)

set(NEO_RT_SOURCES
  src/attenuation_factor.f90
  src/collis_nbi.f90
  src/polylag_3.f
  src/util.f90
  src/magfie.f90
  src/profiles.f90
  src/driftorbit.f90
  src/orbit.f90
  src/resonance.f90
  src/nonlin.f90
  src/transport.f90
  src/neort.f90
)

option(USE_STANDALONE "Use do_magfie_standalone" ON)

if(USE_STANDALONE)
  add_library("neo_rt" STATIC
    ${NEO_RT_SOURCES}
    src/do_magfie_standalone.f90
  )
else()

  add_library("neo_rt" STATIC
    ${NEO_RT_SOURCES}
    src/do_magfie_neo.f90
  )
  add_dependencies(neo_rt neo2_ql)
  target_include_directories(neo_rt PUBLIC "${PROJECT_BINARY_DIR}/NEO-2/NEO-2-QL/OBJS")
  target_link_libraries(neo_rt PUBLIC
    neo2_ql
  )
endif()

add_dependencies(neo_rt spline)
target_include_directories(neo_rt PUBLIC "${PROJECT_BINARY_DIR}/spline")
add_dependencies(neo_rt vode)
target_include_directories(neo_rt PUBLIC "${PROJECT_BINARY_DIR}/vode")

find_package(BLAS)
find_package(LAPACK)



target_link_libraries(neo_rt PUBLIC
  vode
  spline
  BLAS::BLAS
  LAPACK::LAPACK
)

target_link_libraries(${PROJECT_EXE_NAME} neo_rt neo2_ql)


add_subdirectory(test)
