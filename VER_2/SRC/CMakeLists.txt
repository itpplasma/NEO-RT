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
link_libraries(
  vode
  lapack
)

add_library(neo-rt-full common.f90 orbit_full.f90 transport_full.f90)
add_subdirectory(../../POTATO/VER_3/SRC potato)

add_library(neo-rt-thin common.f90 orbit_thin.f90)

add_executable (neo-rt.x main.f90)
if (ORBIT_MODE STREQUAL full)
  target_link_libraries(neo-rt.x neo-rt-full)
elseif (ORBIT_MODE STREQUAL thin)
  target_link_libraries(neo-rt.x neo-rt-thin)
endif()

add_executable (test_orbit_full.x test_orbit.f90)
target_link_libraries(test_orbit_full.x neo-rt-full)
add_executable (test_orbit_thin.x test_orbit.f90)
target_link_libraries(test_orbit_thin.x neo-rt-thin)
