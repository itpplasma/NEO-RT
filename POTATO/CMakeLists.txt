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
# OpenMP parallelizes the find_bounce-heavy node-fill loops of the resonance
# grid build (SRC/sample_matrix.f90).  Linked to potato_base below so the
# -fopenmp compile flag reaches POTATO's own sources and the runtime is linked.
find_package(OpenMP REQUIRED)

include(FetchContent)
FetchContent_Declare(
    fortplot
    GIT_REPOSITORY https://github.com/lazy-fortran/fortplot.git
    GIT_TAG main
)
FetchContent_MakeAvailable(fortplot)

# fortnum provides the variable-order Adams ODE integrator with event/root
# detection used by find_bounce to locate the bounce return (replacing the
# fixed-step Poincare-cut crossing search). Same pin as the NEO-RT migration.
if(NOT TARGET fortnum)
    FetchContent_Declare(
        fortnum
        GIT_REPOSITORY https://github.com/lazy-fortran/fortnum.git
        GIT_TAG 92de6e949a772cfffc73bb5295fe5e2b056b9c18
    )
    FetchContent_MakeAvailable(fortnum)
endif()
# fortnum relies on a GNU implicit-typing extension; pin it to the GNU standard.
if(TARGET fortnum)
    target_compile_options(fortnum PRIVATE -std=gnu)
endif()

# libneo supplies the field, magfie and spline code POTATO used to carry as
# relative symlinks into a sibling libneo checkout. Link the upstream library so
# fixes (e.g. the stretch_coords convex-wall load) propagate instead of drifting.
find_or_fetch(libneo)

# Hold POTATO's own sources to f2018, stop at the first diagnostic, and run the C
# preprocessor. Applied after the dependency fetches above so the upstream
# libraries (libneo, fortplot, fortnum) build under their own standards instead
# of inheriting these; -cpp in particular makes CMake mis-parse libneo's .f90.
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

# Public-safe golden-record gate for the OpenMP resonance solver. Runs potato.x
# on a synthetic concentric-circle equilibrium (no AUG #30835 data) and checks
# the per-energy-slice resonance-line counts against a committed golden record,
# plus OMP=1 vs OMP=4 determinism. Catches the get_rescond weight-corruption
# regression. Registered only when Python 3 and numpy are available so the gate
# is opt-in for environments without the Python test deps.
enable_testing()
find_package(Python3 COMPONENTS Interpreter)
if(Python3_Interpreter_FOUND)
    execute_process(
        COMMAND ${Python3_EXECUTABLE} -c "import numpy"
        RESULT_VARIABLE _have_numpy
        OUTPUT_QUIET ERROR_QUIET)
    if(_have_numpy EQUAL 0)
        add_test(
            NAME resonance_golden
            COMMAND ${Python3_EXECUTABLE} -m pytest -v
                ${CMAKE_CURRENT_SOURCE_DIR}/test/golden_record_resonance/test_resonance_golden.py
        )
        set_tests_properties(resonance_golden PROPERTIES
            ENVIRONMENT "POTATO_EXE=$<TARGET_FILE:potato.x>"
            TIMEOUT 300)
    endif()
endif()
