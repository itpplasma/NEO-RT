cmake_minimum_required (VERSION 3.22)

project(POTATO)
enable_language(Fortran)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/../cmake")

set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR})
add_compile_options(-g -fbacktrace)
add_compile_options(-O3 -march=native -mtune=native)

include(Util)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(HDF5 REQUIRED COMPONENTS Fortran)
# OpenMP parallelizes the find_bounce-heavy resonance grid build
# (SRC/sample_matrix.f90) and the per-mode resonance root search
# (SRC/resonant_int.f90).  Linked to potato_base so -fopenmp reaches POTATO's own
# sources and the runtime is linked.
find_package(OpenMP REQUIRED)

include(FetchContent)
FetchContent_Declare(
    fortplot
    GIT_REPOSITORY https://github.com/lazy-fortran/fortplot.git
    GIT_TAG main
)
FetchContent_MakeAvailable(fortplot)

# libneo supplies the field, magfie and spline code POTATO used to carry as
# relative symlinks into a sibling libneo checkout. Link the upstream library so
# fixes (e.g. the stretch_coords convex-wall load) propagate instead of drifting.
find_or_fetch(libneo)

# Hold POTATO's own sources to f2018, stop at the first diagnostic, and run the C
# preprocessor. Applied after the dependency fetches above so the upstream
# libraries (libneo, fortplot) build under their own standards instead of
# inheriting these; -cpp in particular makes CMake mis-parse libneo's .f90.
add_compile_options(-cpp -Wfatal-errors -std=f2018)

add_subdirectory(SRC)

add_executable (potato.x
    SRC/tt.f90
)
target_link_libraries(potato.x
potato
)


add_executable (test_profile.x
   SRC/test_profile.f90
)
target_link_libraries(test_profile.x
potato
)

add_executable (test_profile_fixed.x
   SRC/test_profile.f90
)
target_link_libraries(test_profile_fixed.x
potato_fixed
)

add_executable(import_neo2_profiles.x
   SRC/import_neo2_profiles.f90
)
target_link_libraries(import_neo2_profiles.x
   potato
)
