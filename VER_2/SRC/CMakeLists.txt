cmake_minimum_required(VERSION 3.13)

project(NEO-RT Fortran)

set(CMAKE_Fortran_FLAGS "-cpp -g -fbacktrace -std=gnu")
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/OBJS)
set(ORBIT_MODE "full" CACHE STRING "Orbit approximation (full/thin)")

message(STATUS "CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")
message(STATUS "ORBIT_MODE=${ORBIT_MODE}")
message(STATUS
"-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
List config using 'cmake -L .'
and change using 'cmake -DVAR_NAME=VALUE .'
or edit CMakeCache.txt
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --"
)

add_subdirectory("../../../contrib/vode" "vode")

add_executable (neo-rt.x main.f90 common.f90)

if (ORBIT_MODE STREQUAL full)
  target_sources(neo-rt.x PRIVATE orbit_full.f90 transport_full.f90)
  add_subdirectory(../../POTATO/VER_3/SRC potato)
endif()

target_link_libraries(neo-rt.x
  vode
  lapack
)
