cmake_minimum_required (VERSION 3.22)

project(POTATO)
enable_language(Fortran)
enable_testing()
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/../cmake")

set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR})
add_compile_options(-g -fbacktrace)
if (CMAKE_SYSTEM_PROCESSOR MATCHES "arm64|aarch64")
    add_compile_options(-O3 -mcpu=native)
else ()
    add_compile_options(-O3 -march=native -mtune=native)
endif ()

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

add_executable(test_find_all_roots_bracketed.x
   SRC/test_find_all_roots_bracketed.f90
)
target_link_libraries(test_find_all_roots_bracketed.x
   potato_base
)
add_test(NAME test_find_all_roots_bracketed
   COMMAND test_find_all_roots_bracketed.x
)

add_executable(test_wall_loss.x
   SRC/test_wall_loss.f90
)
target_link_libraries(test_wall_loss.x
   potato_base
)
add_test(NAME test_wall_loss
   COMMAND test_wall_loss.x
)

add_executable(test_wall_profile.x
   SRC/test_wall_profile.f90
)
target_link_libraries(test_wall_profile.x
   potato_base
)

add_executable(import_neo2_profiles.x
   SRC/import_neo2_profiles.f90
)
target_link_libraries(import_neo2_profiles.x
   potato
)

add_executable(potato_resonance_contour.x
   SRC/potato_resonance_contour.f90
)
target_link_libraries(potato_resonance_contour.x
   potato
)

add_executable(potato_resonance_probe.x
   SRC/potato_resonance_probe.f90
)
target_link_libraries(potato_resonance_probe.x
   potato
)

find_package(Python3 COMPONENTS Interpreter)
if(Python3_Interpreter_FOUND)
    add_test(NAME potato_golden_displacement
        COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/test/golden_record_displacement/test_potato_golden.py
            $<TARGET_FILE:potato.x>
            ${CMAKE_CURRENT_SOURCE_DIR}/test/golden_record_displacement/golden_zero_mid.json
    )
    set_tests_properties(potato_golden_displacement PROPERTIES
        SKIP_RETURN_CODE 77
        TIMEOUT 240
    )
    add_test(NAME potato_resonance_probe_clip
        COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/test/resonance_probe_clip/test_probe_clip.py
            $<TARGET_FILE:potato_resonance_probe.x>
    )
    set_tests_properties(potato_resonance_probe_clip PROPERTIES
        SKIP_RETURN_CODE 77
        TIMEOUT 240
    )
endif()
