cmake_minimum_required(VERSION 3.22)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

option(USE_STANDALONE "Use do_magfie_standalone" ON)
if(NOT USE_STANDALONE)
    message(FATAL_ERROR
        "USE_STANDALONE=OFF is currently unsupported: the native NEO-2 adapter "
        "requires compatible NEO-2 field/perturbation APIs and CMake dependency "
        "integration that are not available on main. Use -DUSE_STANDALONE=ON "
        "(the default). See README.md, Building, for details.")
endif()

project(driftorbit LANGUAGES Fortran C)
set(PROJECT_EXE_NAME "neo_rt.x")

# Loaded first, so a consumer's copy (if any) stays the active one.
include(${PROJECT_SOURCE_DIR}/cmake/ItpDeps.cmake)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR})

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

add_compile_options(-fPIC -g)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fbacktrace -fopenmp -fmax-errors=1")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wimplicit-interface -Wimplicit-procedure")

# gfortran on Apple silicon rejects -march=native next to the -mcpu=native that
# dependencies add ("unknown value 'apple-m3' for '-march'"); use -mcpu there.
if(CMAKE_SYSTEM_PROCESSOR MATCHES "arm64|aarch64")
    set(APPLE_NATIVE_CPU_FLAG -mcpu=native)
else()
    set(APPLE_NATIVE_CPU_FLAG -march=native)
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_options(-Og -Wextra)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fbounds-check -finit-real=nan -fcheck=all")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffpe-trap=zero,overflow,invalid,underflow")
    add_compile_options()
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
    add_compile_options(-O2)

    if(APPLE)
        add_compile_options(${APPLE_NATIVE_CPU_FLAG})
    else()
        add_compile_options(-march=x86-64-v2 -mtune=generic)
    endif()
elseif(CMAKE_BUILD_TYPE STREQUAL "Fast")
    add_compile_options(-O3 -ffast-math)

    if(APPLE)
        add_compile_options(${APPLE_NATIVE_CPU_FLAG})
    else()
        add_compile_options(-march=native -mtune=native)
    endif()
endif()

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" AND UNIX AND NOT APPLE AND NOT DEFINED BLA_VENDOR)
    find_package(PkgConfig QUIET)
    if(PKG_CONFIG_FOUND)
        pkg_check_modules(OPENBLAS QUIET openblas)
        if(OPENBLAS_FOUND)
            set(BLA_VENDOR OpenBLAS)
            message(STATUS "Preferring OpenBLAS for BLAS/LAPACK on GNU/Linux")
        endif()
    endif()
endif()

# Consumers get NEO-RT's tests only on request: -DNEORT_BUILD_TESTING=ON plus
# include(CTest) in the consumer.
option(NEORT_BUILD_TESTING "Build tests for NEO-RT" ${PROJECT_IS_TOP_LEVEL})
if(PROJECT_IS_TOP_LEVEL)
    include(CTest)
endif()
if(NEORT_BUILD_TESTING AND NOT BUILD_TESTING)
    message(STATUS "NEORT_BUILD_TESTING is ON, but BUILD_TESTING is OFF; "
        "NEO-RT tests will not be built.")
    set(NEORT_BUILD_TESTING OFF)
endif()

# Symbolic generators are consumer-owned and optional.  The generated kernels
# are committed and are all that a normal NEO-RT build needs.  Set
# FORTSYM_SOURCE_DIR to the local fortsym checkout when regenerating them.
option(NEORT_BUILD_SYMBOLICS "Build NEO-RT fortsym generators" OFF)
set(FORTSYM_SOURCE_DIR "" CACHE PATH "fortsym source tree for NEO-RT generators")

itp_add_dependency(spline REF 1c40553a2cc6c930804cc8ea9cc9393106345a6e
    TARGET spline)
itp_add_dependency(fortnum REF bf0b0c5a859d7ed78e4ae17ee3a9011c8de5b088
    GIT_REPOSITORY https://github.com/lazy-fortran/fortnum.git
    TARGET fortnum)
include(${PROJECT_SOURCE_DIR}/cmake/Fortplot.cmake)

add_executable(${PROJECT_EXE_NAME} src/program.f90)

set(NEO_RT_SOURCES
    src/attenuation_factor.f90
    src/collis_nbi.f90
    src/collisional_layer.f90
    src/config.f90
    src/krook.f90
    src/collision_rates.f90
    src/mixture_collisions.f90
    src/datatypes.f90
    src/driftorbit.f90
    src/kinetic_energy_selection.f90
    src/freq.f90
    src/lib.f90
    src/logging.f90
    src/magfie.f90
    src/main.f90
    src/neort.f90
    src/nonlin.f90
    src/orbit.f90
    src/orbit_classes.f90
    src/output.f90
    src/polylag_3.f
    src/phi1.f90
    src/profiles.f90
    src/resonance.f90
    src/resonance_cutout.f90
    src/generated/neort_separatrix_frequency_symbolic.f90
    src/separatrix.f90
    src/shaing.f90
    src/species.f90
    src/transport.f90
    src/util.f90
)

# Fast builds must retain explicit rejection of nonfinite scientific inputs.
set_source_files_properties(src/phi1.f90 src/separatrix.f90 src/orbit_classes.f90 src/resonance.f90 src/krook.f90 src/profiles.f90 src/species.f90 src/collision_rates.f90 src/mixture_collisions.f90 src/do_magfie_standalone.f90 src/collisional_layer.f90 src/resonance_cutout.f90 PROPERTIES
    COMPILE_OPTIONS "-fno-finite-math-only")

if(USE_STANDALONE)
    itp_add_dependency(libneo REF b8f80b23ad8657dd6c60e1bff469335e347b14a8
        TARGET LIBNEO::boozer)

    add_library("neo_rt" STATIC
        ${NEO_RT_SOURCES}
        src/do_magfie_standalone.f90
    )

    if(NEORT_BUILD_TESTING)
        add_library("util_for_test_${PROJECT_NAME}" STATIC test/util_for_test.f90)
    endif()
else()
    itp_add_dependency(neo_2 REF 881c2336bac918d3e3c79ec2b282fc29f67cfaa3
        GIT_REPOSITORY https://github.com/itpplasma/NEO-2.git
        TARGET neo2_ql)
    FetchContent_GetProperties(neo_2)
    add_library("neo_rt" STATIC
        ${NEO_RT_SOURCES}
        src/do_magfie_neo.f90
    )
    add_dependencies(neo_rt neo2_ql)
    target_include_directories(neo_rt PUBLIC "${neo_2_BINARY_DIR}/NEO-2-QL/OBJS")
    target_link_libraries(neo_rt PUBLIC
        neo2_ql
    )
endif()

target_include_directories(neo_rt PUBLIC "${PROJECT_BINARY_DIR}")
# spline does not export its Fortran module folder.
get_target_property(_neort_spline_modules spline Fortran_MODULE_DIRECTORY)
if(NOT _neort_spline_modules)
    get_target_property(_neort_spline_modules spline BINARY_DIR)
endif()
target_include_directories(neo_rt PUBLIC "${_neort_spline_modules}")

find_package(BLAS)
find_package(LAPACK)

target_link_libraries(neo_rt PUBLIC
    fortnum
    spline
    BLAS::BLAS
    LAPACK::LAPACK
)

target_compile_options(neo_rt PRIVATE
    $<$<COMPILE_LANGUAGE:Fortran>:-cpp>
)

if(USE_STANDALONE)
    target_link_libraries(neo_rt PUBLIC
        LIBNEO::boozer
    )
endif()

# HDF5 output. libneo's hdf5_tools forwards to fortio, which implements the
# required HDF5 subset without linking the system HDF5 library.
if(NOT TARGET LIBNEO::hdf5_tools)
    message(FATAL_ERROR "LIBNEO::hdf5_tools not found; it is required for HDF5 output")
endif()

target_link_libraries(neo_rt PUBLIC LIBNEO::hdf5_tools)

target_link_libraries(${PROJECT_EXE_NAME} neo_rt)

if(NEORT_BUILD_TESTING)
    add_subdirectory(test)
endif()

# Diagnostics
add_subdirectory(src/diag)

if(NEORT_BUILD_SYMBOLICS)
    add_subdirectory(tools/symbolics)
endif()

add_executable(neo_rt_diag.x
    src/diag/neort_diag.f90
)

target_link_libraries(neo_rt_diag.x PRIVATE
    neort_diaglib
    neo_rt
    fortplot::fortplot
)

# Field-only readback; does not initialize profiles, orbits or transport.
add_executable(neo_rt_field_probe.x src/diag/neort_field_probe.f90)
target_link_libraries(neo_rt_field_probe.x PRIVATE neo_rt)

add_executable(neo_rt_profile_probe.x src/diag/neort_profile_probe.f90)
target_link_libraries(neo_rt_profile_probe.x PRIVATE neo_rt)
