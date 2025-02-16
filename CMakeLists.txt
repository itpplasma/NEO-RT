### Basic settings
cmake_minimum_required (VERSION 2.8.12)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")
enable_testing()

### Define the project
project (driftorbit)
enable_language(Fortran)

### Specify paths
set(CMAKE_BINARY_DIR ${PROJECT_BINARY_DIR}/OBJS)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR})

### Define the name of the executable
set(PROJECT_EXE_NAME "neo_rt.x")                  # Define the name

### Set compiler flags
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/OBJS)

# Set default build type to Release if not specified
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

add_compile_options(-fPIC -cpp -g)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fbacktrace")
# Set compiler options based on build type
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  add_compile_options(-O0 -ffpe-trap=zero,overflow,invalid,underflow)
  add_compile_options(-fbounds-check -finit-real=nan -fcheck=all -Wextra)
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
  add_compile_options(-O3 -ffast-math)
  if(APPLE)
    add_compile_options(-march=native)
  else()
    add_compile_options(-march=haswell -mtune=generic)
  endif()
elseif(CMAKE_BUILD_TYPE STREQUAL "Profile")
  add_compile_options(-O0 -pg)
  add_link_options(-pg)
endif()

include(Util)

find_or_fetch(BOOZER_MAGFIE)
find_or_fetch(spline)
find_or_fetch(quadpack)
find_or_fetch(vode)

### External modules
# add_subdirectory("../BOOZER_MAGFIE" "boozer_magfie")
# add_subdirectory("../spline" "spline")
# add_subdirectory("../contrib/quadpack" "quadpack")
# add_subdirectory("../contrib/vode" "vode")
#add_subdirectory(POTATO/LATEST)  # contains also vode

### Define executable
add_executable(${PROJECT_EXE_NAME}
  SRC/main.f90
)

add_executable(orbit_test
  TEST/test_orbit.f90
)

set(NEO_RT_SOURCES
  SRC/orbit.f90
  SRC/attenuation_factor.f90
  SRC/collis_nbi.f90
  SRC/polylag_3.f
  SRC/util.f90
  SRC/driftorbit.f90
)

# Set a configuration flag
option(USE_STANDALONE "Use do_magfie_standalone" ON)

# Add source files based on the configuration flag
if(USE_STANDALONE)
  add_library("neo_rt" SHARED
    ${NEO_RT_SOURCES}
    SRC/do_magfie_standalone.f90
  )
else()
  # For NEO-2 field routines
  include_directories(${PROJECT_BINARY_DIR}/boozer_magfie/OBJS/)

  add_library("neo_rt" SHARED
    ${NEO_RT_SOURCES}
    SRC/do_magfie_neo.f90
  )
  add_dependencies(neo_rt boozer_magfie)
  target_link_libraries(neo_rt
    boozer_magfie
  )
endif()

add_dependencies(neo_rt spline)
add_dependencies(neo_rt quadpack)
add_dependencies(neo_rt vode)

find_package(BLAS)
find_package(LAPACK)

target_link_libraries(neo_rt
  vode
  spline
  quadpack
  ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES}
)

### Linking
target_link_libraries(${PROJECT_EXE_NAME}
  neo_rt
)

add_executable(test_frequencies.x
  TEST/test_frequencies.f90
)

target_link_libraries(test_frequencies.x
neo_rt
)

add_executable(test_reslines.x TEST/test_reslines.f90)
target_link_libraries(test_reslines.x neo_rt)

add_executable(test_bounce.x TEST/test_bounce.f90)
target_link_libraries(test_bounce.x neo_rt )

### Tests
add_test(test_magfie neo_rt)
