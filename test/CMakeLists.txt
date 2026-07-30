file(GLOB TEST_SOURCES "test_*.f90")

foreach(test_src ${TEST_SOURCES})
    get_filename_component(test_name ${test_src} NAME_WE)
    add_executable(${test_name}.x ${test_src})
    target_link_libraries(${test_name}.x neo_rt "util_for_test_${PROJECT_NAME}")
endforeach()

add_executable(probe_neort_perturbation_fourier_kernel.x
    probe_neort_perturbation_fourier_kernel.f90)
target_link_libraries(probe_neort_perturbation_fourier_kernel.x neo_rt)

add_test(NAME test_timestep
    COMMAND test_timestep.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_test(NAME test_parallel
    COMMAND ${CMAKE_BINARY_DIR}/test_parallel.x
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/examples/base
)

add_test(NAME test_neort_lib
    COMMAND ${CMAKE_BINARY_DIR}/test_neort_lib.x
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/examples/base
)

# test_chartmap_input: Boozer chartmap field I/O test (inp_swi=10).
#
# The reference chartmap fixture (test/fixtures/circ_chartmap.nc) was generated
# once from examples/circ.bc via the libneo converters (bc_to_booz_xform +
# booz_xform_to_boozer_chartmap, nrho=50 ntheta=64 nzeta=4) and committed into
# the repo.  No converter toolchain is required at test time.
#
# Phase 1: sign/sanity check (bmod>0, sqrtg>0, iota>0, Bthcov>0, Bphcov>0).
# Phase 2: cross-check vs the .bc path (relative tolerance 1e-4).
# Both phases run unconditionally in CI.
add_test(NAME test_chartmap_input
    COMMAND ${CMAKE_BINARY_DIR}/test_chartmap_input.x
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/examples/base
)

add_test(NAME test_boozer_handedness
    COMMAND ${CMAKE_BINARY_DIR}/test_boozer_handedness.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_test(NAME test_profile_collision_species
    COMMAND ${CMAKE_BINARY_DIR}/test_profile_collision_species.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_test(NAME test_nonlinear_contract
    COMMAND ${CMAKE_BINARY_DIR}/test_nonlinear_contract.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_test(NAME test_attenuation_asymptotes
    COMMAND ${CMAKE_BINARY_DIR}/test_attenuation_asymptotes.x
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/examples/nonlin
)

add_test(NAME test_harmonic_bounds
    COMMAND test_harmonic_bounds.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_test(NAME test_response_jxb
    COMMAND test_response_jxb.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_test(NAME test_boozer_ampere
    COMMAND test_boozer_ampere.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_test(NAME test_boozer_jxb
    COMMAND test_boozer_jxb.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_test(NAME test_boozer_chartmap_geometry
    COMMAND test_boozer_chartmap_geometry.x
        ${CMAKE_CURRENT_SOURCE_DIR}/fixtures/boozer_torus_geometry.nc
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_test(NAME test_boozer_vector_io
    COMMAND test_boozer_vector_io.x
        ${CMAKE_CURRENT_SOURCE_DIR}/fixtures/boozer_vector_manufactured.nc
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_test(NAME test_boozer_current_program
    COMMAND ${CMAKE_BINARY_DIR}/neo_rt_boozer_current.x
        ${CMAKE_CURRENT_SOURCE_DIR}/fixtures/boozer_vector_manufactured.nc
        ${CMAKE_CURRENT_BINARY_DIR}/boozer_weighted_current.out
        2.5
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
set_tests_properties(test_boozer_current_program PROPERTIES
    PASS_REGULAR_EXPRESSION "computed weighted current on 4 surfaces for 2"
)

add_test(NAME test_boozer_jxb_program
    COMMAND ${CMAKE_BINARY_DIR}/neo_rt_boozer_jxb.x
        ${CMAKE_CURRENT_SOURCE_DIR}/fixtures/boozer_torus_geometry.nc
        ${CMAKE_CURRENT_SOURCE_DIR}/fixtures/boozer_vector_manufactured.nc
        ${CMAKE_CURRENT_BINARY_DIR}/boozer_jxb_profile.out
        2.5
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
set_tests_properties(test_boozer_jxb_program PROPERTIES
    PASS_REGULAR_EXPRESSION "computed Boozer jxb torque density on 4 surfaces"
)

add_test(NAME test_boozer_jxb_profile
    COMMAND test_boozer_jxb_profile.x
        ${CMAKE_CURRENT_SOURCE_DIR}/fixtures/boozer_torus_geometry.nc
        ${CMAKE_CURRENT_SOURCE_DIR}/fixtures/boozer_vector_manufactured.nc
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_test(NAME test_mars_response
    COMMAND test_mars_response.x ${CMAKE_CURRENT_SOURCE_DIR}/fixtures
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_test(NAME test_mars_response_program
    COMMAND ${CMAKE_BINARY_DIR}/neo_rt_jxb_from_mars.x
        ${CMAKE_CURRENT_SOURCE_DIR}/fixtures/jxb_linear_edges.out
        ${CMAKE_CURRENT_SOURCE_DIR}/fixtures/mars_bplasma_small.out
        ${CMAKE_CURRENT_SOURCE_DIR}/fixtures/mars_jplasma_small.out
        ${CMAKE_CURRENT_BINARY_DIR}/mars_jxb_small.out
        2 0 1.0
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
set_tests_properties(test_mars_response_program PROPERTIES
    PASS_REGULAR_EXPRESSION "processed native total:"
)

add_test(NAME test_jxb_profile_program
    COMMAND ${CMAKE_BINARY_DIR}/neo_rt_jxb_profile.x
        ${CMAKE_CURRENT_SOURCE_DIR}/fixtures/jxb_linear_edges.out
        ${CMAKE_CURRENT_SOURCE_DIR}/fixtures/jxb_linear_profile.out
        2.5
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
set_tests_properties(test_jxb_profile_program PROPERTIES
    PASS_REGULAR_EXPRESSION "1.579136704174[0-9]+E\\+02"
)

add_test(NAME test_response_jxb_rejects_shape
    COMMAND test_response_jxb_invalid.x shape
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_test(NAME test_response_jxb_rejects_profile
    COMMAND test_response_jxb_invalid.x profile
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
set_tests_properties(
    test_response_jxb_rejects_shape
    test_response_jxb_rejects_profile
    PROPERTIES WILL_FAIL TRUE
)

add_test(NAME test_vmax_over_vth
    COMMAND test_vmax_over_vth.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_test(NAME test_resonance_bracket
    COMMAND ${CMAKE_BINARY_DIR}/test_resonance_bracket.x
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/examples/base
)
add_test(NAME test_perturbation_fourier_kernel
    COMMAND ${CMAKE_BINARY_DIR}/probe_neort_perturbation_fourier_kernel.x
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
set_tests_properties(test_chartmap_input PROPERTIES
    ENVIRONMENT
        "CHARTMAP_FILE=${CMAKE_CURRENT_SOURCE_DIR}/fixtures/circ_chartmap.nc;BC_FILE=${CMAKE_SOURCE_DIR}/examples/circ.bc"
)
