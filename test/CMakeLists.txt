file(GLOB TEST_SOURCES "test_*.f90")

if(NOT USE_STANDALONE)
    list(FILTER TEST_SOURCES EXCLUDE REGEX "test_gc_cylindrical_eqdsk\\.f90$")
endif()

foreach(test_src ${TEST_SOURCES})
    get_filename_component(test_name ${test_src} NAME_WE)
    add_executable(${test_name}.x ${test_src})
    target_link_libraries(${test_name}.x neo_rt "util_for_test_${PROJECT_NAME}")
endforeach()

target_link_libraries(test_frequency_inventory_format.x neort_diaglib)

add_test(NAME test_frequency_inventory_format
    COMMAND test_frequency_inventory_format.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

target_sources(test_gc_orbit_limit.x PRIVATE
    generated_manufactured_torus.f90
    manufactured_gc_field.f90
)

add_executable(probe_neort_perturbation_fourier_kernel.x
    probe_neort_perturbation_fourier_kernel.f90)
target_link_libraries(probe_neort_perturbation_fourier_kernel.x neo_rt)

add_executable(probe_frequency_comparison.x probe_frequency_comparison.f90)
target_link_libraries(probe_frequency_comparison.x neo_rt)

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

add_test(NAME test_eqdsk_axis
    COMMAND ${CMAKE_BINARY_DIR}/test_eqdsk_axis.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_test(NAME test_profile_collision_species
    COMMAND ${CMAKE_BINARY_DIR}/test_profile_collision_species.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_test(NAME test_input_switches
    COMMAND ${CMAKE_BINARY_DIR}/test_input_switches.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_test(NAME test_input_switches_reject_boozer_rz
    COMMAND ${CMAKE_BINARY_DIR}/test_input_switches.x reject_boozer_rz
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
set_tests_properties(test_input_switches_reject_boozer_rz PROPERTIES WILL_FAIL TRUE)
add_test(NAME test_input_switches_reject_full_boozer
    COMMAND ${CMAKE_BINARY_DIR}/test_input_switches.x reject_full_boozer
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
set_tests_properties(test_input_switches_reject_full_boozer PROPERTIES WILL_FAIL TRUE)
add_test(NAME test_input_switches_reject_full_supban
    COMMAND ${CMAKE_BINARY_DIR}/test_input_switches.x reject_full_supban
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
set_tests_properties(test_input_switches_reject_full_supban PROPERTIES WILL_FAIL TRUE)
add_test(NAME test_input_switches_reject_full_nonlin
    COMMAND ${CMAKE_BINARY_DIR}/test_input_switches.x reject_full_nonlin
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
set_tests_properties(test_input_switches_reject_full_nonlin PROPERTIES WILL_FAIL TRUE)
add_test(NAME test_input_switches_reject_model1
    COMMAND ${CMAKE_BINARY_DIR}/test_input_switches.x reject_model1
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
set_tests_properties(test_input_switches_reject_model1 PROPERTIES WILL_FAIL TRUE)

add_test(NAME test_wall_polygon_io
    COMMAND ${CMAKE_BINARY_DIR}/test_wall_polygon_io.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_test(NAME test_wall_config
    COMMAND ${CMAKE_BINARY_DIR}/test_wall_config.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_test(NAME test_wall_config_rejects_missing
    COMMAND ${CMAKE_BINARY_DIR}/test_wall_config.x missing
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
set_tests_properties(test_wall_config_rejects_missing PROPERTIES WILL_FAIL TRUE)
add_test(NAME test_wall_config_rejects_malformed
    COMMAND ${CMAKE_BINARY_DIR}/test_wall_config.x malformed
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
set_tests_properties(test_wall_config_rejects_malformed PROPERTIES WILL_FAIL TRUE)
add_test(NAME test_wall_config_rejects_namelist_missing
    COMMAND ${CMAKE_BINARY_DIR}/test_wall_config.x namelist_missing
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
set_tests_properties(test_wall_config_rejects_namelist_missing PROPERTIES WILL_FAIL TRUE)
add_test(NAME test_wall_config_rejects_model1
    COMMAND ${CMAKE_BINARY_DIR}/test_wall_config.x model1
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
set_tests_properties(test_wall_config_rejects_model1 PROPERTIES WILL_FAIL TRUE)

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
add_test(NAME test_perturbation_fourier_kernel
    COMMAND ${CMAKE_BINARY_DIR}/probe_neort_perturbation_fourier_kernel.x
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
set_tests_properties(test_chartmap_input PROPERTIES
    ENVIRONMENT
        "CHARTMAP_FILE=${CMAKE_CURRENT_SOURCE_DIR}/fixtures/circ_chartmap.nc;BC_FILE=${CMAKE_SOURCE_DIR}/examples/circ.bc"
)
set_tests_properties(test_eqdsk_axis PROPERTIES
    ENVIRONMENT
        "EQDSK_FILE=${CMAKE_SOURCE_DIR}/POTATO/test/golden_record_resonance/circ.eqdsk;RZ_PERT_FILE=${CMAKE_SOURCE_DIR}/POTATO/test/golden_record_resonance/bmod_n.dat"
)

add_test(NAME test_gc_eqdsk_frequency
    COMMAND ${CMAKE_BINARY_DIR}/test_gc_eqdsk_frequency.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_test(NAME test_gc_eqdsk_cut_chart_symbolic
    COMMAND ${CMAKE_BINARY_DIR}/test_gc_eqdsk_cut_chart_symbolic.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_test(NAME test_gc_eqdsk_cut_interval_symbolic
    COMMAND ${CMAKE_BINARY_DIR}/test_gc_eqdsk_cut_interval_symbolic.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_test(NAME test_gc_eqdsk_cut_axis_limit
    COMMAND ${CMAKE_BINARY_DIR}/test_gc_eqdsk_cut_axis_limit.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_test(NAME test_gc_circular_flux_continuation
    COMMAND ${CMAKE_BINARY_DIR}/test_gc_circular_flux_continuation.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_test(NAME test_gc_eqdsk_cut_endpoint_krawczyk
    COMMAND ${CMAKE_BINARY_DIR}/test_gc_eqdsk_cut_endpoint_krawczyk.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_test(NAME test_gc_eqdsk_cut_endpoint_certificate
    COMMAND ${CMAKE_BINARY_DIR}/test_gc_eqdsk_cut_endpoint_certificate.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
set_tests_properties(test_gc_eqdsk_cut_endpoint_certificate PROPERTIES
    ENVIRONMENT
        "EQDSK_FILE=${CMAKE_SOURCE_DIR}/POTATO/test/golden_record_resonance/circ.eqdsk"
    TIMEOUT 180
)

add_test(NAME test_gc_eqdsk_composite_cut_atlas
    COMMAND ${CMAKE_BINARY_DIR}/test_gc_eqdsk_composite_cut_atlas.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
set_tests_properties(test_gc_eqdsk_composite_cut_atlas PROPERTIES
    ENVIRONMENT
        "EQDSK_FILE=${CMAKE_SOURCE_DIR}/POTATO/test/golden_record_resonance/circ.eqdsk"
    TIMEOUT 180
)

add_test(NAME test_gc_eqdsk_cut_flux_curvature
    COMMAND ${CMAKE_BINARY_DIR}/test_gc_eqdsk_cut_flux_curvature.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_test(NAME test_gc_generated_certificates
    COMMAND ${CMAKE_BINARY_DIR}/test_gc_generated_certificates.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_test(NAME test_gc_eqdsk_cut_interval_runtime
    COMMAND ${CMAKE_BINARY_DIR}/test_gc_eqdsk_cut_interval_runtime.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
set_tests_properties(test_gc_eqdsk_cut_interval_runtime PROPERTIES
    ENVIRONMENT
        "EQDSK_FILE=${CMAKE_SOURCE_DIR}/POTATO/test/golden_record_resonance/circ.eqdsk"
)

add_test(NAME test_gc_eqdsk_cut_graph_atlas
    COMMAND ${CMAKE_BINARY_DIR}/test_gc_eqdsk_cut_graph_atlas.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
set_tests_properties(test_gc_eqdsk_cut_graph_atlas PROPERTIES
    ENVIRONMENT
        "EQDSK_FILE=${CMAKE_SOURCE_DIR}/POTATO/test/golden_record_resonance/circ.eqdsk"
    TIMEOUT 180
)

add_test(NAME test_gc_full_orbit_frequency
    COMMAND ${CMAKE_BINARY_DIR}/test_gc_full_orbit_frequency.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
set_tests_properties(test_gc_full_orbit_frequency PROPERTIES
    ENVIRONMENT
        "EQDSK_FILE=${CMAKE_SOURCE_DIR}/POTATO/test/golden_record_resonance/circ.eqdsk"
    TIMEOUT 180
)

set_tests_properties(test_gc_eqdsk_frequency PROPERTIES
    ENVIRONMENT
        "EQDSK_FILE=${CMAKE_SOURCE_DIR}/POTATO/test/golden_record_resonance/circ.eqdsk"
)

add_test(NAME test_boozer_angle_map
    COMMAND ${CMAKE_BINARY_DIR}/test_boozer_angle_map.x
            $ENV{NEORT_TEST_BOOZER_EQUILIBRIUM}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_test(NAME test_trapping_parameter
    COMMAND ${CMAKE_BINARY_DIR}/test_trapping_parameter.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_test(NAME test_gc_dynamics
    COMMAND ${CMAKE_BINARY_DIR}/test_gc_dynamics.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_test(NAME test_gc_cylindrical_dynamics
    COMMAND ${CMAKE_BINARY_DIR}/test_gc_cylindrical_dynamics.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_test(NAME test_gc_perpendicular_invariant
    COMMAND ${CMAKE_BINARY_DIR}/test_gc_perpendicular_invariant.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_test(NAME test_gc_cylindrical_topology
    COMMAND ${CMAKE_BINARY_DIR}/test_gc_cylindrical_topology.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_test(NAME test_gc_cylindrical_orbit
    COMMAND ${CMAKE_BINARY_DIR}/test_gc_cylindrical_orbit.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_test(NAME test_gc_cylindrical_physical_return
    COMMAND ${CMAKE_BINARY_DIR}/test_gc_cylindrical_physical_return.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

if(USE_STANDALONE)
    add_test(NAME test_gc_cylindrical_backend_entry
        COMMAND ${CMAKE_BINARY_DIR}/test_gc_cylindrical_backend_entry.x
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
    set_tests_properties(test_gc_cylindrical_backend_entry PROPERTIES
        ENVIRONMENT
            "EQDSK_FILE=${CMAKE_SOURCE_DIR}/POTATO/test/golden_record_resonance/circ.eqdsk"
        TIMEOUT 180
    )

    add_test(NAME test_gc_cylindrical_eqdsk
        COMMAND ${CMAKE_BINARY_DIR}/test_gc_cylindrical_eqdsk.x
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
    set_tests_properties(test_gc_cylindrical_eqdsk PROPERTIES
        ENVIRONMENT
            "EQDSK_FILE=${CMAKE_SOURCE_DIR}/POTATO/test/golden_record_resonance/circ.eqdsk"
    )
endif()

add_test(NAME test_thin_orbit_limit
    COMMAND ${CMAKE_BINARY_DIR}/test_thin_orbit_limit.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_test(NAME test_gc_orbit_limit
    COMMAND ${CMAKE_BINARY_DIR}/test_gc_orbit_limit.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_test(NAME test_gc_transport_orbit
    COMMAND ${CMAKE_BINARY_DIR}/test_gc_transport_orbit.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_test(NAME test_gc_full_resonance
    COMMAND ${CMAKE_BINARY_DIR}/test_gc_full_resonance.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_test(NAME test_gc_nonlocal_resonance_integral
    COMMAND ${CMAKE_BINARY_DIR}/test_gc_nonlocal_resonance_integral.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_test(NAME test_gc_cylindrical_nonlocal_provider
    COMMAND ${CMAKE_BINARY_DIR}/test_gc_cylindrical_nonlocal_provider.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_test(NAME test_gc_nonlocal_transport_integral
    COMMAND ${CMAKE_BINARY_DIR}/test_gc_nonlocal_transport_integral.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_test(NAME test_gc_cylindrical_class_adapter
    COMMAND ${CMAKE_BINARY_DIR}/test_gc_cylindrical_class_adapter.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_test(NAME test_gc_cylindrical_transport_provider
    COMMAND ${CMAKE_BINARY_DIR}/test_gc_cylindrical_transport_provider.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_test(NAME test_gc_eqdsk_nonlocal_transport
    COMMAND ${CMAKE_BINARY_DIR}/test_gc_eqdsk_nonlocal_transport.x
        ${CMAKE_SOURCE_DIR}/POTATO/test/golden_record_resonance/circ.eqdsk
        ${CMAKE_SOURCE_DIR}/POTATO/test/golden_record_resonance/convexwall.dat
        ${CMAKE_SOURCE_DIR}/POTATO/test/golden_record_resonance/bmod_n.dat
        ${CMAKE_SOURCE_DIR}/examples/orbit/plasma.in
        ${CMAKE_SOURCE_DIR}/examples/orbit/profile.in
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
set_tests_properties(test_gc_eqdsk_nonlocal_transport PROPERTIES TIMEOUT 180)
add_test(NAME test_gc_full_fow_runtime_metadata
    COMMAND ${CMAKE_BINARY_DIR}/test_gc_full_fow_runtime_metadata.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_test(NAME test_gc_full_fow_runtime_delivery
    COMMAND ${CMAKE_BINARY_DIR}/test_gc_full_fow_runtime_delivery.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_test(NAME test_gc_full_fow_symbolic_kernel
    COMMAND ${CMAKE_BINARY_DIR}/test_gc_full_fow_symbolic_kernel.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_test(NAME test_gc_unit_quadrature
    COMMAND ${CMAKE_BINARY_DIR}/test_gc_unit_quadrature.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_test(NAME test_gc_full_fow_simple_root_force
    COMMAND ${CMAKE_BINARY_DIR}/test_gc_full_fow_simple_root_force.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_test(NAME test_gc_full_fow_scales_boundary_cut
    COMMAND ${CMAKE_BINARY_DIR}/test_gc_full_fow_scales_boundary_cut.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_test(NAME test_gc_cylindrical_littlejohn_generated
    COMMAND ${CMAKE_BINARY_DIR}/test_gc_cylindrical_littlejohn_generated.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_test(NAME test_gc_cylindrical_independent_behavior
    COMMAND ${CMAKE_BINARY_DIR}/test_gc_cylindrical_independent_behavior.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_test(NAME test_gc_generated_geometry_interp_endpoints
    COMMAND ${CMAKE_BINARY_DIR}/test_gc_generated_geometry_interp_endpoints.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_test(NAME test_gc_generated_profile_potential_segment
    COMMAND ${CMAKE_BINARY_DIR}/test_gc_generated_profile_potential_segment.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_test(NAME test_gc_profile_potential_map
    COMMAND ${CMAKE_BINARY_DIR}/test_gc_profile_potential_map.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_test(NAME test_gc_eqdsk_allowed_region_interval
    COMMAND ${CMAKE_BINARY_DIR}/test_gc_eqdsk_allowed_region_interval.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_test(NAME test_gc_eqdsk_allowed_region_cut_box
    COMMAND ${CMAKE_BINARY_DIR}/test_gc_eqdsk_allowed_region_cut_box.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
set_tests_properties(test_gc_eqdsk_allowed_region_cut_box PROPERTIES
    ENVIRONMENT
        "EQDSK_FILE=${CMAKE_SOURCE_DIR}/POTATO/test/golden_record_resonance/circ.eqdsk"
    TIMEOUT 180
)
add_test(NAME test_gc_eqdsk_allowed_region_signs
    COMMAND ${CMAKE_BINARY_DIR}/test_gc_eqdsk_allowed_region_signs.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_test(NAME test_gc_eqdsk_allowed_rho_chain
    COMMAND ${CMAKE_BINARY_DIR}/test_gc_eqdsk_allowed_rho_chain.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_test(NAME test_gc_eqdsk_flux_profile_map
    COMMAND ${CMAKE_BINARY_DIR}/test_gc_eqdsk_flux_profile_map.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_test(NAME test_gc_generated_resonance_eq17_seams
    COMMAND ${CMAKE_BINARY_DIR}/test_gc_generated_resonance_eq17_seams.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_test(NAME test_gc_full_fow_normalization_runtime
    COMMAND ${CMAKE_BINARY_DIR}/test_gc_full_fow_normalization_runtime.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_test(NAME test_gc_generated_invariant_conversions
    COMMAND ${CMAKE_BINARY_DIR}/test_gc_generated_invariant_conversions.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_test(NAME test_gc_sign_symmetry_behavior
    COMMAND ${CMAKE_BINARY_DIR}/test_gc_sign_symmetry_behavior.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_test(NAME test_gc_full_fow_canonical_behavior
    COMMAND ${CMAKE_BINARY_DIR}/test_gc_full_fow_canonical_behavior.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_test(NAME test_gc_homoclinic_limit
    COMMAND ${CMAKE_BINARY_DIR}/test_gc_homoclinic_limit.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_test(NAME test_gc_certified_interval_roots
    COMMAND ${CMAKE_BINARY_DIR}/test_gc_certified_interval_roots.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_test(NAME test_gc_poincare_cut_atlas
    COMMAND ${CMAKE_BINARY_DIR}/test_gc_poincare_cut_atlas.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_test(NAME test_gc_transport_status
    COMMAND ${CMAKE_BINARY_DIR}/test_gc_transport_status.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_test(NAME test_gc_model2_transport_dispatch
    COMMAND ${CMAKE_BINARY_DIR}/test_gc_model2_transport_dispatch.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
