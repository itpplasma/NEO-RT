if(NOT DEFINED NEORT_SOURCE_DIR)
    message(FATAL_ERROR "NEORT_SOURCE_DIR is required")
endif()

if(DEFINED ENV{NEORT_TEST_TMPDIR})
    set(regeneration_dir "$ENV{NEORT_TEST_TMPDIR}/neort-full-fow-regeneration")
else()
    set(regeneration_dir "${NEORT_SOURCE_DIR}/.neort-full-fow-regeneration")
endif()

file(REMOVE_RECURSE "${regeneration_dir}")
file(MAKE_DIRECTORY "${regeneration_dir}")
execute_process(
    COMMAND fo exec gen_full_fow_physics "${regeneration_dir}"
    WORKING_DIRECTORY "${NEORT_SOURCE_DIR}/tools/gc_symbolics"
    RESULT_VARIABLE generation_status
)
if(NOT generation_status EQUAL 0)
    message(FATAL_ERROR "Fortsym full-FOW regeneration failed: ${generation_status}")
endif()

set(generated_files
    neort_generated_certificate_registry.f90
    neort_cylindrical_geometry_symbolic.f90
    neort_cylindrical_bstar_symbolic.f90
    neort_cylindrical_littlejohn_symbolic.f90
    neort_axisymmetric_pphi_symbolic.f90
    neort_axisymmetric_noether_symbolic.f90
    neort_cylindrical_crossing_symbolic.f90
    neort_buchholz_cut_symbolic.f90
    neort_buchholz_cdot_symbolic.f90
    neort_buchholz_boundary_symbolic.f90
    neort_buchholz_action_symbolic.f90
    neort_buchholz_energy_symbolic.f90
    neort_buchholz_specific_energy_symbolic.f90
    neort_physical_mu_symbolic.f90
    neort_potato_mu_symbolic.f90
    neort_potato_velocity_symbolic.f90
    neort_full_fow_action_symbolic.f90
    neort_full_fow_perturbation_symbolic.f90
    neort_full_fow_harmonic_symbolic.f90
    neort_full_fow_resonance_symbolic.f90
    neort_full_fow_frequency_contribution_symbolic.f90
    neort_full_fow_frequency_identity_symbolic.f90
    neort_full_fow_simple_root_symbolic.f90
    neort_full_fow_sign_symmetry_symbolic.f90
    neort_full_fow_normalization_symbolic.f90
    neort_full_fow_quadrature_map_symbolic.f90
    neort_polynomial_cell_enclosure_symbolic.f90
    neort_full_fow_section_reversal_symbolic.f90
    neort_full_fow_eq17_symbolic.f90
    neort_full_fow_eq17_outer_symbolic.f90
    neort_full_fow_refinement_symbolic.f90
    neort_cylindrical_bilinear_complex_symbolic.f90
    neort_profile_endpoint_symbolic.f90
    neort_profile_potential_segment_symbolic.f90
    neort_eqdsk_flux_profile_segment_symbolic.f90
    neort_eqdsk_quintic_cell_jet_symbolic.f90
    neort_eqdsk_quintic_cell_jet_interval_symbolic.f90
    neort_eqdsk_quintic_cell_fourth_jet_symbolic.f90
    neort_eqdsk_quintic_cell_fourth_jet_interval_symbolic.f90
    neort_eqdsk_quintic_profile_jet_symbolic.f90
    neort_eqdsk_quintic_profile_jet_interval_symbolic.f90
    neort_eqdsk_cut_jet_symbolic.f90
    neort_eqdsk_cut_numerator_symbolic.f90
    neort_eqdsk_cut_numerator_interval_symbolic.f90
    neort_eqdsk_cut_numerator_hessian_symbolic.f90
    neort_eqdsk_cut_numerator_hessian_interval_symbolic.f90
    neort_eqdsk_cut_r_chart_symbolic.f90
    neort_eqdsk_cut_z_chart_symbolic.f90
    neort_eqdsk_cut_r_flux_chart_symbolic.f90
    neort_eqdsk_cut_r_flux_chart_interval_symbolic.f90
    neort_eqdsk_cut_r_flux_curvature_symbolic.f90
    neort_eqdsk_cut_r_flux_curvature_interval_symbolic.f90
    neort_eqdsk_cut_mean_value_interval_symbolic.f90
    neort_eqdsk_cut_axis_curvature_interval_symbolic.f90
    neort_eqdsk_cut_axis_limit_symbolic.f90
    neort_eqdsk_rho_tor_map_symbolic.f90
    neort_eqdsk_cut_flux_coordinate_symbolic.f90
    neort_eqdsk_cut_axis_rho_limit_symbolic.f90
    neort_cylindrical_hamiltonian_symbolic.f90
    neort_cylindrical_canonical_symbolic.f90
    neort_cylindrical_vparallel_symbolic.f90
    neort_cylindrical_launch_symbolic.f90
)
foreach(generated_name IN LISTS generated_files)
    execute_process(
        COMMAND ${CMAKE_COMMAND} -E compare_files
            "${regeneration_dir}/${generated_name}"
            "${NEORT_SOURCE_DIR}/src/generated/${generated_name}"
        RESULT_VARIABLE compare_status
    )
    if(NOT compare_status EQUAL 0)
        message(FATAL_ERROR "Generated Fortsym file differs: ${generated_name}")
    endif()
endforeach()
message(STATUS "Fortsym full-FOW generated outputs match")
