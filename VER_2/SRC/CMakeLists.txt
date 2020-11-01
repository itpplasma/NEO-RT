cmake_minimum_required(VERSION 3.13)

project(NEO-RT Fortran)

set(CMAKE_Fortran_FLAGS "-cpp -g -fbacktrace")
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/OBJS)
set(ORBIT_MODE "full" CACHE STRING "Orbit approximation (full/thin)")

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Orbit type: ${ORBIT_MODE}")

add_subdirectory("../../../contrib/vode" "vode")

add_executable (neo-rt.x main.f90 common.f90)

if (ORBIT_MODE STREQUAL full)
  add_subdirectory(../../POTATO/VER_1/SRC potato)
  target_sources(neo-rt.x PRIVATE orbit_full.f90 transport_full.f90)
endif()

target_link_libraries(neo-rt.x
  vode
  lapack
)
