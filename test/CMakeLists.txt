file(GLOB TEST_SOURCES "test_*.f90")

foreach(test_src ${TEST_SOURCES})
    get_filename_component(test_name ${test_src} NAME_WE)
    add_executable(${test_name}.x ${test_src})
    target_link_libraries(${test_name}.x neo_rt "util_for_test_${PROJECT_NAME}")
endforeach()

add_executable(probe_neort_perturbation_fourier_kernel.x
    probe_neort_perturbation_fourier_kernel.f90)
target_link_libraries(probe_neort_perturbation_fourier_kernel.x neo_rt)

# test_timestep: runs from test/ripple_plateau so the test can read the
# same in_file and driftorbit_test.in configuration as the ripple-plateau test.
add_test(NAME test_timestep
    COMMAND test_timestep.x
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/test/ripple_plateau
)

add_test(NAME test_parallel
    COMMAND ${CMAKE_BINARY_DIR}/test_parallel.x
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/examples/base
)

add_test(NAME test_neort_lib
    COMMAND ${CMAKE_BINARY_DIR}/test_neort_lib.x
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/examples/base
)

add_test(NAME test_parallel_transport
    COMMAND ${CMAKE_BINARY_DIR}/test_parallel_transport.x
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/examples/base
)

set_tests_properties(test_parallel_transport PROPERTIES
    ENVIRONMENT
        "OMP_DYNAMIC=FALSE;OMP_NUM_THREADS=4;OMP_THREAD_LIMIT=4"
    PROCESSORS 4
    TIMEOUT 600
)

# test_chartmap_input: Boozer chartmap field I/O test (inp_swi=10).
#
# The reference chartmap fixture (test/fixtures/circ_chartmap.nc) was generated
# once from examples/circ.bc via the libneo converters (bc_to_booz_xform +
# booz_xform_to_boozer_chartmap, nrho=50 ntheta=64 nzeta=4) and committed into
# the repo.  No converter toolchain is required at test time.  After the
# circ.bc header flux sign fix (#180) only torflux and A_phi change sign; the
# fixture was updated accordingly (identical to a regeneration to 1e-15).
#
# Phase 1: sign/sanity check (bmod>0, sqrtg<0, iota>0, Bthcov>0, Bphcov>0).
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

# Load-time Boozer orientation guard (#180): the consistent circular file and the
# ASDEX-type golden-record file load; flux-flipped copies of both are rejected.
function(add_boozer_orientation_tests name file swi)
    add_test(NAME test_boozer_orientation_accept_${name}
        COMMAND ${CMAKE_BINARY_DIR}/test_boozer_orientation_guard.x accept ${file} ${swi}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
    add_test(NAME test_boozer_orientation_reject_${name}
        COMMAND ${CMAKE_BINARY_DIR}/test_boozer_orientation_guard.x flip ${file} ${swi}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
    set_tests_properties(test_boozer_orientation_reject_${name} PROPERTIES
        PASS_REGULAR_EXPRESSION "Inconsistent Boozer orientation"
    )
    add_test(NAME test_boozer_orientation_magnitude_${name}
        COMMAND ${CMAKE_BINARY_DIR}/test_boozer_orientation_guard.x scale ${file} ${swi}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/orientation_scale_${name}
    )
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/orientation_scale_${name})
    set_tests_properties(test_boozer_orientation_magnitude_${name} PROPERTIES
        PASS_REGULAR_EXPRESSION "Inconsistent Boozer orientation"
    )
endfunction()
add_boozer_orientation_tests(circ ${CMAKE_SOURCE_DIR}/examples/circ.bc 8)
add_boozer_orientation_tests(golden ${CMAKE_SOURCE_DIR}/test/golden_record/input/in_file 9)
add_test(NAME test_boozer_orientation_drift_circ
    COMMAND ${CMAKE_BINARY_DIR}/test_boozer_orientation_guard.x drift
        ${CMAKE_SOURCE_DIR}/examples/circ.bc 8
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

add_test(NAME test_passing_drift_derivative
    COMMAND test_passing_drift_derivative.x
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/examples/base
)

add_test(NAME test_vmax_over_vth
    COMMAND test_vmax_over_vth.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_test(NAME test_resonance_bracket
    COMMAND ${CMAKE_BINARY_DIR}/test_resonance_bracket.x
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/examples/base
)
add_test(NAME test_resonance_oracle
    COMMAND ${CMAKE_BINARY_DIR}/test_resonance_oracle.x
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
