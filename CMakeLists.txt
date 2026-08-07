cmake_minimum_required(VERSION 3.22)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

project(driftorbit LANGUAGES Fortran C)
set(PROJECT_EXE_NAME "neo_rt.x")

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
    add_compile_options(-Og -Wextra)
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

if(NOT DEFINED LIBNEO_REF)
    set(LIBNEO_REF 72a05682e51a50ac716cfe70b4917338ce0e08e1)
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

# Testing configuration
set(NEORT_IS_STANDALONE OFF)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    set(NEORT_IS_STANDALONE ON)
endif()

option(NEORT_BUILD_TESTING "Build tests for NEO-RT" ${NEORT_IS_STANDALONE})

if(DEFINED BUILD_TESTING)
    set(NEORT_BUILD_TESTING ${BUILD_TESTING})
endif()

# Disable tests for dependencies
set(BUILD_TESTING OFF)

find_or_fetch(spline)

include(FetchContent)

if(NOT TARGET fortnum)
    FetchContent_Declare(
        fortnum
        GIT_REPOSITORY https://github.com/lazy-fortran/fortnum.git
        GIT_TAG b4e2942240ad8a8a4f191c3d473eea98364cbf6f
    )
    FetchContent_MakeAvailable(fortnum)
endif()

FetchContent_Declare(
    fortplot
    GIT_REPOSITORY https://github.com/lazy-fortran/fortplot.git
    GIT_TAG main
)
FetchContent_MakeAvailable(fortplot)

add_executable(${PROJECT_EXE_NAME} src/program.f90)

set(NEO_RT_SOURCES
    src/attenuation_factor.f90
    src/collis_nbi.f90
    src/boozer_angle_map.f90
    src/wall_polygon_io.f90
    src/config.f90
    src/datatypes.f90
    src/driftorbit.f90
    src/gc_wall_context.f90
    src/gc_perpendicular_invariant.f90
    src/gc_outward_interval.f90
    src/generated/neort_generated_certificate_registry.f90
    src/generated/neort_cylindrical_geometry_symbolic.f90
    src/generated/neort_cylindrical_bstar_symbolic.f90
    src/generated/neort_cylindrical_littlejohn_symbolic.f90
    src/generated/neort_axisymmetric_pphi_symbolic.f90
    src/generated/neort_axisymmetric_noether_symbolic.f90
    src/generated/neort_cylindrical_crossing_symbolic.f90
    src/generated/neort_buchholz_cut_symbolic.f90
    src/generated/neort_buchholz_cdot_symbolic.f90
    src/generated/neort_buchholz_boundary_symbolic.f90
    src/generated/neort_buchholz_action_symbolic.f90
    src/generated/neort_buchholz_energy_symbolic.f90
    src/generated/neort_buchholz_specific_energy_symbolic.f90
    src/generated/neort_physical_mu_symbolic.f90
    src/generated/neort_potato_mu_symbolic.f90
    src/generated/neort_potato_velocity_symbolic.f90
    src/generated/neort_full_fow_action_symbolic.f90
    src/generated/neort_full_fow_perturbation_symbolic.f90
    src/generated/neort_full_fow_harmonic_symbolic.f90
    src/generated/neort_full_fow_resonance_symbolic.f90
    src/generated/neort_full_fow_frequency_contribution_symbolic.f90
    src/generated/neort_full_fow_frequency_identity_symbolic.f90
    src/generated/neort_full_fow_simple_root_symbolic.f90
    src/generated/neort_full_fow_sign_symmetry_symbolic.f90
    src/generated/neort_full_fow_normalization_symbolic.f90
    src/generated/neort_full_fow_quadrature_map_symbolic.f90
    src/generated/neort_polynomial_cell_enclosure_symbolic.f90
    src/generated/neort_full_fow_section_reversal_symbolic.f90
    src/generated/neort_full_fow_eq17_symbolic.f90
    src/generated/neort_full_fow_eq17_outer_symbolic.f90
    src/generated/neort_full_fow_refinement_symbolic.f90
    src/generated/neort_cylindrical_bilinear_complex_symbolic.f90
    src/generated/neort_profile_endpoint_symbolic.f90
    src/generated/neort_profile_potential_segment_symbolic.f90
    src/generated/neort_eqdsk_quintic_cell_jet_symbolic.f90
    src/generated/neort_eqdsk_quintic_cell_jet_interval_symbolic.f90
    src/generated/neort_eqdsk_quintic_profile_jet_symbolic.f90
    src/generated/neort_eqdsk_quintic_profile_jet_interval_symbolic.f90
    src/generated/neort_eqdsk_cut_jet_symbolic.f90
    src/generated/neort_eqdsk_cut_numerator_symbolic.f90
    src/generated/neort_eqdsk_cut_numerator_interval_symbolic.f90
    src/generated/neort_eqdsk_cut_r_chart_symbolic.f90
    src/generated/neort_eqdsk_cut_z_chart_symbolic.f90
    src/generated/neort_eqdsk_cut_r_flux_chart_symbolic.f90
    src/generated/neort_eqdsk_cut_r_flux_chart_interval_symbolic.f90
    src/generated/neort_eqdsk_cut_mean_value_interval_symbolic.f90
    src/generated/neort_cylindrical_hamiltonian_symbolic.f90
    src/generated/neort_cylindrical_canonical_symbolic.f90
    src/generated/neort_cylindrical_vparallel_symbolic.f90
    src/generated/neort_cylindrical_launch_symbolic.f90
    src/gc_dynamics.f90
    src/gc_cylindrical_model.f90
    src/gc_cylindrical_dynamics.f90
    src/gc_cylindrical_physical_return.f90
    src/gc_cylindrical_topology.f90
    src/gc_certified_interval_roots.f90
    src/gc_poincare_cut_atlas.f90
    src/gc_callback_context.f90
    src/gc_cylindrical_class_adapter.f90
    src/gc_cylindrical_orbit.f90
    src/gc_eqdsk_adapter.f90
    src/gc_eqdsk_cut_jet.f90
    src/gc_eqdsk_cut_interval.f90
    src/gc_eqdsk_cut_graph_atlas.f90
    src/gc_frequency_provider.f90
    src/gc_full_resonance.f90
    src/gc_nonlocal_resonance_types.f90
    src/gc_nonlocal_resonance_integral.f90
    src/gc_cylindrical_nonlocal_provider.f90
    src/gc_nonlocal_transport_types.f90
    src/gc_nonlocal_transport_integral.f90
    src/gc_cylindrical_transport_provider.f90
    src/gc_full_fow_normalization_runtime.f90
    src/gc_eqdsk_nonlocal_transport.f90
    src/gc_full_fow_runtime_metadata.f90
    src/gc_full_fow_runtime_delivery.f90
    src/gc_model2_transport_dispatch.f90
    src/gc_frequency_splines.f90
    src/gc_models.f90
    src/gc_orbit_integrator.f90
    src/freq.f90
    src/lib.f90
    src/logging.f90
    src/magfie.f90
    src/main.f90
    src/neort.f90
    src/nonlin.f90
    src/orbit.f90
    src/polylag_3.f
    src/profiles.f90
    src/resonance.f90
    src/shaing.f90
    src/transport.f90
    src/thin_orbit_limit.f90
    src/util.f90
)

option(USE_STANDALONE "Use do_magfie_standalone" ON)

if(USE_STANDALONE)
    list(APPEND NEO_RT_SOURCES
        src/gc_eqdsk_cylindrical_adapter.f90
        src/gc_cylindrical_backend.f90
    )
endif()

add_custom_target(check_neort_full_fow_generated
    COMMAND ${CMAKE_COMMAND} -DNEORT_SOURCE_DIR=${CMAKE_SOURCE_DIR}
        -P ${CMAKE_SOURCE_DIR}/tools/gc_symbolics/check_full_fow_generated.cmake
    USES_TERMINAL
    COMMENT "Regenerating and comparing Fortsym full-FOW kernels")

if(USE_STANDALONE)
    find_or_fetch(libneo)

    add_library("neo_rt" STATIC
        ${NEO_RT_SOURCES}
        src/do_magfie_standalone.f90
    )

    if(NEORT_BUILD_TESTING)
        add_library("util_for_test_${PROJECT_NAME}" STATIC test/util_for_test.f90)
    endif()
else()
    find_or_fetch(NEO-2)
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

target_include_directories(neo_rt PUBLIC "${PROJECT_BINARY_DIR}")
add_dependencies(neo_rt spline)
target_include_directories(neo_rt PUBLIC "${PROJECT_BINARY_DIR}/spline")

find_package(BLAS)
find_package(LAPACK)

target_link_libraries(neo_rt PUBLIC
    fortnum
    spline
    BLAS::BLAS
    LAPACK::LAPACK
)

if(USE_STANDALONE)
    target_compile_definitions(neo_rt PRIVATE NEO_RT_USE_STANDALONE)
    target_link_libraries(neo_rt PUBLIC
        LIBNEO::boozer
        LIBNEO::magfie
    )
endif()

target_link_libraries(${PROJECT_EXE_NAME} neo_rt)

if(NEORT_BUILD_TESTING)
    set(BUILD_TESTING ON)
    include(CTest)
    add_subdirectory(test)
endif()

# Diagnostics
add_subdirectory(src/diag)

add_executable(neo_rt_diag.x
    src/diag/neort_diag.f90
)

target_link_libraries(neo_rt_diag.x PRIVATE
    neort_diaglib
    neo_rt
    fortplot::fortplot
)
