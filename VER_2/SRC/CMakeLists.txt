cmake_minimum_required (VERSION 2.8.12)

project(NEO-RT)
enable_language(Fortran)

set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR})
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

add_compile_options(-Og -cpp -g -fbacktrace -Wfatal-errors -std=f2008)

find_package(BLAS)
find_package(LAPACK)

add_subdirectory(../../../contrib/vode vode)
add_subdirectory(../../POTATO/LATEST/SRC potato)
add_library(neo-rt-util util.f90)

add_dependencies(neo-rt-util vode)


link_libraries(
  neo-rt-util
  potato
  vode
  ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES}
)

add_library(neo-rt-full orbit_full.f90 transport_full.f90)
# add_library(neo-rt-thin orbit_thin.f90 transport_thin.f90)

add_executable (neo-rt.x main.f90)
if (ORBIT_MODE STREQUAL full)
  target_link_libraries(neo-rt.x neo-rt-full)
elseif (ORBIT_MODE STREQUAL thin)
  target_link_libraries(neo-rt.x neo-rt-thin)
endif()

add_subdirectory("../../../spline" "spline")
add_executable (test_orbit_full.x test_orbit.f90
  ../../do_magfie_standalone.f90 ../../common.f90)
target_link_libraries(test_orbit_full.x neo-rt-full spline)
#add_executable (test_orbit_thin.x test_orbit.f90)
#target_link_libraries(test_orbit_thin.x neo-rt-thin)
