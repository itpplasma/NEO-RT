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

# test_multi_init: repeated do_magfie_init calls must not corrupt or change the
# loaded field data (deallocate-before-reallocate).  The oracle is nmode/nflux
# compared against the values from the first init, and the program stops with a
# nonzero status on any mismatch.  Runs from examples/base, which supplies the
# driftorbit.in and in_file it reads.
add_test(NAME test_multi_init
    COMMAND ${CMAKE_BINARY_DIR}/test_multi_init.x
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/examples/base
)

# test_bounce: the trapped bounce time, against a Gauss-Legendre quadrature of
# the same orbit that does not integrate it, and then the two ODE paths against
# each other.  The absolute oracle runs first because it is the only one of the
# three that can see a defect the two integrations share -- which is how they
# came to agree with each other to 0.26% while both were 0.4% wrong.  Runs from
# examples/base, which supplies the driftorbit.in and in_file it reads.
add_test(NAME test_bounce
    COMMAND ${CMAKE_BINARY_DIR}/test_bounce.x
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

# test_hcurl_axisymmetric: curl of the unit field direction in the standalone
# Boozer path. The oracle is a central difference of hcovar over separate
# do_magfie calls, so it is independent of the analytic formula it checks.
add_test(NAME test_hcurl_axisymmetric
    COMMAND ${CMAKE_BINARY_DIR}/test_hcurl_axisymmetric.x
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

# ---------------------------------------------------------------------------
# Built but deliberately NOT registered as ctest tests
# ---------------------------------------------------------------------------
#
# The glob at the top builds an executable for every test_*.f90 in this
# directory.  Being built is not the same as being run.  The programs listed
# below are built and are intentionally left out of ctest, for two different
# reasons.
#
# (a) Self-checking, but currently failing or missing an input file.
#     Registering them as they stand would turn the suite red; weakening their
#     assertions to make them green would destroy the only oracle they have.
#     They stay unregistered until the underlying problem is fixed.
#
#   test_bounce.f90 is now registered above.  It was left here because
#       bounce_time() and bounce() disagreed by 2.6e-3 against a 1e-8
#       tolerance.  That was real: the pinned fortnum moved its clock onto a
#       located event root without rebuilding the history it interpolates from,
#       so every trapped bounce time -- which continues past the half-bounce
#       crossing by construction -- came out short.  Fixed upstream and
#       repinned; both paths now match a quadrature of the same orbit to
#       1.7e-9.
#
#   test_torque.f90
#       Needs the run configuration driftorbit64.new.in next to in_file,
#       plasma.in and profile.in in the working directory.  driftorbit64.new.in
#       is not in this repository, so the program aborts inside
#       read_and_set_config before it reaches either of its two assertions.
#       Commit that input (or point the program at an existing configuration)
#       before registering it here.
#
# (b) Diagnostics, not tests.  These programs print numbers or write .dat files
#     for plotting and contain no assertion at all, so they exit 0 no matter
#     what the physics does.  Registering them would add permanently green
#     tests that can never detect a regression, which is worse than not running
#     them.  Keep them out of ctest unless an oracle is added to them.
#
#   test_bounce_integrand.f90
#   test_freq_scan.f90
#   test_frequencies.f90  -- writes <runname>_ux*_freq_vs_eta_*.dat and exits 0
#                            for any runname; its only nonzero exit is the
#                            "usage:" guard when no argument is passed.
#   test_misc.f90
#   test_omega_prime.f90
#   test_reslines.f90
