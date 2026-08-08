program gen_full_fow_physics
    !! Derive and emit the direct real-space full-FOW scalar kernel.
    !!
    !! This is the sole source of the algebra used by the future cylindrical
    !! dynamics seam.  The generated routine is deliberately scalar and
    !! side-effect free: vector components are flattened at the interface so
    !! the consumer can use it from ordinary Fortran, OpenMP, or OpenACC
    !! without a CAS dependency at run time.
    !!
    !! The proofs below use only exact symbolic identities.  In particular,
    !! no random probe is accepted for a conservation or normalization claim.
    !! Canonical p_phi is generated from the axisymmetric phase-space one-form
    !! under the explicit flux convention A_phi_cov=psi and
    !! b_phi_cov=R*b_phi.  The scalar Lagrangian proof closes the Noether
    !! identity; the field/flux convention is checked independently below.
    use, intrinsic :: iso_fortran_env, only: int64, output_unit
    use fortsym_arena, only: arena_t
    use fortsym_check, only: suite_t, suite_begin, suite_end, check_identity
    use fortsym_diff, only: diff
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_native, only: make_native_engine, native_engine_t
    use fortsym_engine_symengine, only: make_symengine_engine, symengine_engine_t
    use fortsym_expr, only: abs, cos, exp, expr_t, func, log, num, rat, operator(+), &
        operator(-), operator(*), operator(/), operator(**), sin, sqrt, sym
    use fortsym_expr, only: pi_expr
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE
    use fortsym_string, only: chars, str, str_t
    use fortsym_subs, only: subs
    implicit none

    character(*), parameter :: DEFAULT_OUTPUT = "../../src/generated"
    character(*), parameter :: FORTSYM_REVISION = &
        "fortsym@545788453a204d58705f735b519c3863c2f734c8"
    character(len=64), parameter :: EQDSK_CELL_COEFFICIENT_NAMES(36) = &
        [character(len=64) :: &
        "coefficient_0_0", "coefficient_0_1", "coefficient_0_2", &
        "coefficient_0_3", "coefficient_0_4", "coefficient_0_5", &
        "coefficient_1_0", "coefficient_1_1", "coefficient_1_2", &
        "coefficient_1_3", "coefficient_1_4", "coefficient_1_5", &
        "coefficient_2_0", "coefficient_2_1", "coefficient_2_2", &
        "coefficient_2_3", "coefficient_2_4", "coefficient_2_5", &
        "coefficient_3_0", "coefficient_3_1", "coefficient_3_2", &
        "coefficient_3_3", "coefficient_3_4", "coefficient_3_5", &
        "coefficient_4_0", "coefficient_4_1", "coefficient_4_2", &
        "coefficient_4_3", "coefficient_4_4", "coefficient_4_5", &
        "coefficient_5_0", "coefficient_5_1", "coefficient_5_2", &
        "coefficient_5_3", "coefficient_5_4", "coefficient_5_5"]
    character(len=64), parameter :: EQDSK_PROFILE_COEFFICIENT_NAMES(6) = &
        [character(len=64) :: "coefficient_0", "coefficient_1", &
        "coefficient_2", "coefficient_3", "coefficient_4", &
        "coefficient_5"]

    type(arena_t), target :: arena
    type(symengine_engine_t) :: proof_engine
    type(native_engine_t) :: simplify_engine
    type(suite_t) :: proofs
    type(expr_t) :: roots(48), action_roots(11), perturbation_roots(6)
    type(expr_t) :: harmonic_integrand_roots(8)
    type(expr_t) :: resonance_roots(3), eq17_roots(2), cylindrical_roots(22)
    type(expr_t) :: noether_roots(4)
    type(expr_t) :: refinement_roots(14), boundary_roots(13), cut_roots(5)
    type(expr_t) :: cutdot_roots(3), section_roots(9), crossing_roots(5)
    type(expr_t) :: sign_symmetry_roots(14)
    type(expr_t) :: normalization_roots(13), quadrature_map_roots(8)
    type(expr_t) :: gauss_interval_map_roots(2)
    type(expr_t) :: polynomial_enclosure_roots(3)
    type(expr_t) :: cylindrical_bstar_roots(5)
    type(expr_t) :: physical_mu_roots(3), buchholz_action_roots(2)
    type(expr_t) :: buchholz_energy_roots(3), potato_velocity_roots(2)
    type(expr_t) :: potato_mu_roots(2), cylindrical_hamiltonian_roots(1)
    type(expr_t) :: cylindrical_canonical_roots(2)
    type(expr_t) :: cylindrical_launch_roots(3)
    type(expr_t) :: cylindrical_vparallel_roots(1)
    type(expr_t) :: buchholz_specific_energy_roots(1)
    type(expr_t) :: mass, charge, c_light, mu, bmod, h, electrostatic_potential
    type(expr_t) :: p_phi
    type(expr_t) :: b1, b2, b3, bhat1, bhat2, bhat3
    type(expr_t) :: curl1, curl2, curl3, gradb1, gradb2, gradb3
    type(expr_t) :: gradphi1, gradphi2, gradphi3, radius, p_parallel
    type(expr_t) :: n_mode, tau, g_prime, temperature, a1, a2, phi_eff
    type(expr_t) :: a_real, a_imag, m_mode, theta, phi
    type(expr_t) :: phase_half_real, phase_half_imag, phase_half_norm
    type(expr_t) :: q_abs, zero, one, jk, omega_c, omega_c_var, jk_omega_c, jk_candidate
    type(expr_t) :: q_phi_energy
    type(expr_t) :: jk_max, d_jk_d_h, d_jk_d_phi, d_jk_d_omega, eq4_coeff
    type(expr_t) :: psi_star, dpsi_dp_phi, f_prime, freq_weight, phase_weight
    type(expr_t) :: resonance_x, resonance_xr, resonance_delta
    type(expr_t) :: g_prime_x, tau_prime_x, g_local, tau_local, f_local
    type(expr_t) :: g_at_root, f_prime_from_series
    type(expr_t) :: force_bracket, bstar1, bstar2, bstar3
    type(expr_t) :: bparallel_star, force1, force2, force3, cross1, cross2
    type(expr_t) :: cross3, v1, v2, v3, dot_r, dot_z, dot_phi, dot_p
    type(expr_t) :: dot_mu, cylindrical_measure, hamiltonian_dot
    type(expr_t) :: energy_on_shell, eta_on_shell, c_on_shell
    type(expr_t) :: perturbation_ratio, boozer_ratio, chi
    type(expr_t) :: harmonic_delta_b_real, harmonic_delta_b_imag
    type(expr_t) :: harmonic_phi, harmonic_phi_launch, harmonic_time
    type(expr_t) :: harmonic_omega_b, harmonic_omega_phi
    type(expr_t) :: harmonic_phase_argument, harmonic_phase_real
    type(expr_t) :: harmonic_phase_imag, harmonic_factor
    type(expr_t) :: harmonic_integrand_real, harmonic_integrand_imag
    type(expr_t) :: harmonic_complete_phase, harmonic_complete_phase_residual
    type(expr_t) :: harmonic_modulus_residual, harmonic_shell_residual
    type(expr_t) :: harmonic_shell_identity
    type(expr_t) :: field_original, field_phase_shifted, field_sign_sum
    type(expr_t) :: phase_cos, phase_sin, chi_cos, chi_sin
    type(expr_t) :: shifted_chi_cos, shifted_chi_sin
    type(expr_t) :: phase_unit_residual, phase_relabel_residual
    type(expr_t) :: field_reversal, field_fixed_conjugate
    type(expr_t) :: fixed_conjugation_difference
    type(expr_t) :: amp_square, amp_square_sign, amp_square_conj
    type(expr_t) :: amplitude_phase, amp_rotated_real, amp_rotated_imag
    type(expr_t) :: amp_square_rotated
    type(expr_t) :: sign_amplitude_residual, phase_amplitude_residual
    type(expr_t) :: pair_field_residual, conjugate_modulus_residual
    type(expr_t) :: n_abs, tau_pos, g_prime_abs, f_prime_abs
    type(expr_t) :: root_weight_abs, phase_weight_abs, h_axisymmetric
    type(expr_t) :: psi, b_phi_cov, cylindrical_bhat_phi
    type(expr_t) :: phi_coordinate, phi_dot_symbol
    type(expr_t) :: psi_r, psi_z, b_r_from_psi, b_z_from_psi
    type(expr_t) :: canonical_p_phi, canonical_p_phi_cylindrical
    type(expr_t) :: phase_space_lagrangian
    type(expr_t) :: d_lagrangian_d_phi, d_lagrangian_d_phi_dot
    type(expr_t) :: canonical_residual, flux_b_r_residual, flux_b_z_residual
    type(expr_t) :: state_r, state_z, state_phi, state_p_parallel, state_mu
    type(expr_t) :: obs_tau, obs_omega_b, obs_omega_phi, hm_real, hm_imag
    type(expr_t) :: root_derivative, torque_density, torque, hm_squared
    type(expr_t) :: r_ref, z_ref, phi_ref, p_parallel_ref, mu_ref
    type(expr_t) :: tau_ref, omega_b_ref, omega_phi_ref, hm_real_ref
    type(expr_t) :: hm_imag_ref, hm_squared_ref, root_ref
    type(expr_t) :: torque_density_ref, torque_ref
    type(expr_t) :: scale_r, scale_z, scale_phi, scale_p_parallel, scale_mu
    type(expr_t) :: scale_tau, scale_omega_b, scale_omega_phi, scale_hm_real
    type(expr_t) :: scale_hm_imag, scale_hm_squared, scale_root
    type(expr_t) :: scale_torque_density, scale_torque
    type(expr_t) :: boundary_distance, boundary_reference_scale, boundary_delta
    type(expr_t) :: regular_coefficient
    type(expr_t) :: reflecting_coefficient, separatrix_coefficient
    type(expr_t) :: lambda_positive, c_tau, regular_form, reflecting_form
    type(expr_t) :: reflecting_derivative, separatrix_form, separatrix_tau
    type(expr_t) :: xpoint_form, xpoint_tau, regular_classification_residual
    type(expr_t) :: reflecting_classification_residual
    type(expr_t) :: separatrix_classification_residual
    type(expr_t) :: xpoint_pair_coefficient_residual
    type(expr_t) :: xpoint_pair_tau_residual
    type(expr_t) :: grad_psi1, grad_psi2, grad_psi3
    type(expr_t) :: grad_phi_coordinate1, grad_phi_coordinate2
    type(expr_t) :: grad_phi_coordinate3, cut_cross1, cut_cross2, cut_cross3
    type(expr_t) :: cut_c
    type(expr_t) :: d_cut_c_d_r, d_cut_c_d_arc_phi, d_cut_c_d_z
    type(expr_t) :: dot_cut_r, dot_cut_arc_phi, dot_cut_z, cut_cdot
    type(expr_t) :: cutdot_transversality, cutdot_orientation
    type(expr_t) :: dpsi_dx_section, f_prime_section
    type(expr_t) :: dpsi_dx_reversed, f_prime_reversed, cdot_section
    type(expr_t) :: cdot_reversed, signed_measure_section
    type(expr_t) :: root_weight_section, root_weight_reversed
    type(expr_t) :: positive_crossing_density, positive_crossing_reversed
    type(expr_t) :: period_section, frequency_section
    type(expr_t) :: period_reversed, frequency_reversed
    type(expr_t) :: crossing_signed_measure, crossing_cdot
    type(expr_t) :: crossing_signed_density, crossing_positive_density
    type(expr_t) :: crossing_abs_measure, crossing_abs_cdot
    type(expr_t) :: geometry_roots(19)
    type(expr_t) :: geom_b(3), geom_d_b(3,3), geom_bhat(3)
    type(expr_t) :: geom_grad_b(3), geom_dbhat(3,3), geom_curl(3)
    type(expr_t) :: geom_bmod, geom_radius, geom_arc_phi, geom_z
    type(expr_t) :: oracle_b(3), oracle_d_b(3,3), oracle_bhat(3)
    type(expr_t) :: oracle_grad_b(3), oracle_dbhat(3,3), oracle_curl(3)
    type(expr_t) :: oracle_bmod, oracle_coordinate(3)
    type(expr_t) :: geom_b_r, geom_b_phi, geom_b_z
    type(expr_t) :: geom_dbr_dr, geom_dbphi_dr, geom_dbz_dr
    type(expr_t) :: geom_dbr_darc_phi, geom_dbphi_darc_phi
    type(expr_t) :: geom_dbz_darc_phi, geom_dbr_dz, geom_dbphi_dz
    type(expr_t) :: geom_dbz_dz
    type(expr_t) :: interpolation_roots(9), endpoint_roots(8)
    type(expr_t) :: profile_potential_roots(3)
    type(expr_t) :: profile_potential_map_roots(3)
    type(expr_t) :: eqdsk_flux_profile_segment_roots(3)
    type(expr_t) :: eqdsk_scaled_flux_normalization_roots(1)
    type(expr_t) :: eqdsk_physical_flux_normalization_roots(1)
    type(expr_t) :: eqdsk_physical_flux_map_roots(2)
    type(expr_t) :: eqdsk_allowed_energy_roots(17)
    type(expr_t) :: eqdsk_canonical_cut_roots(4)
    type(expr_t) :: eqdsk_turning_chart_roots(5)
    type(expr_t) :: eqdsk_cell_jet_roots(10), eqdsk_profile_jet_roots(4)
    type(expr_t) :: eqdsk_cell_fourth_jet_roots(5)
    type(expr_t) :: eqdsk_cut_jet_roots(7)
    type(expr_t) :: eqdsk_cut_numerator_roots(3)
    type(expr_t) :: eqdsk_cut_interval_roots(4)
    type(expr_t) :: canonical_numerator_roots(3)
    type(expr_t) :: canonical_turning_roots(2)
    type(expr_t) :: canonical_symmetry_roots(21)
    type(expr_t) :: eqdsk_cut_r_chart_roots(2), eqdsk_cut_z_chart_roots(2)
    type(expr_t) :: eqdsk_cut_r_flux_chart_roots(2)
    type(expr_t) :: eqdsk_cut_mean_value_roots(1)
    type(expr_t) :: eqdsk_cut_numerator_hessian_roots(3)
    type(expr_t) :: eqdsk_cut_r_flux_curvature_roots(4)
    type(expr_t) :: eqdsk_cut_axis_curvature_roots(2)
    type(expr_t) :: eqdsk_cut_axis_limit_roots(3)
    type(expr_t) :: eqdsk_rho_tor_map_roots(2)
    type(expr_t) :: eqdsk_s_tor_to_rho_roots(1)
    type(expr_t) :: eqdsk_flux_profile_rho_chain_roots(1)
    type(expr_t) :: eqdsk_cut_flux_coordinate_roots(3)
    type(expr_t) :: eqdsk_cut_axis_rho_limit_roots(3)
    type(expr_t) :: eqdsk_allowed_rho_chain_roots(7)
    type(expr_t) :: eqdsk_allowed_axis_rho_chain_roots(2)
    type(expr_t) :: eqdsk_cut_endpoint_system_roots(6)
    type(expr_t) :: eqdsk_cut_endpoint_newton_roots(7)
    type(expr_t) :: eqdsk_cut_endpoint_krawczyk_roots(2)
    type(expr_t) :: eqdsk_axis_stationarity_system_roots(6)
    type(expr_t) :: eqdsk_axis_stationarity_newton_roots(7)
    type(expr_t) :: eqdsk_axis_stationarity_krawczyk_roots(2)
    type(expr_t) :: eq17_outer_roots(1)
    type(expr_t) :: axisymmetric_pphi_roots(3)
    type(expr_t) :: frequency_contribution_roots(2)
    type(expr_t) :: interp_u, interp_v, interp_w00, interp_w10
    type(expr_t) :: interp_w01, interp_w11
    type(expr_t) :: interp_r00, interp_r10, interp_r01, interp_r11
    type(expr_t) :: interp_i00, interp_i10, interp_i01, interp_i11
    type(expr_t) :: interp_value_real, interp_value_imag
    type(expr_t) :: interp_partition_residual
    type(expr_t) :: interp_corner00, interp_corner10, interp_corner01
    type(expr_t) :: interp_corner11
    type(expr_t) :: interp_r, interp_z, interp_r0, interp_r1
    type(expr_t) :: interp_z0, interp_z1
    type(expr_t) :: interp_amplitude_scale
    type(expr_t) :: interp_u0_residual, interp_u1_residual
    type(expr_t) :: interp_v0_residual, interp_v1_residual
    type(expr_t) :: endpoint_s0, endpoint_s1, endpoint_f0, endpoint_f1
    type(expr_t) :: endpoint_slope, endpoint_intercept
    type(expr_t) :: endpoint_value_zero, endpoint_value_one
    type(expr_t) :: endpoint_derivative_zero, endpoint_derivative_one
    type(expr_t) :: endpoint_rho, endpoint_dfd_rho
    type(expr_t) :: endpoint_axis_residual, affine_a, affine_b
    type(expr_t) :: affine_node0, affine_node1, affine_endpoint_zero
    type(expr_t) :: affine_endpoint_one, affine_endpoint_slope
    type(expr_t) :: potential_psi0, potential_psi1, potential_omega0
    type(expr_t) :: potential_omega1, potential_c, delta_phi_segment
    type(expr_t) :: delta_phi_reversed, omega_constant, delta_phi_constant
    type(expr_t) :: potential_psi, potential_phi0, potential_delta_psi
    type(expr_t) :: potential_omega_slope, potential_phi
    type(expr_t) :: potential_dphi_dpsi, potential_d2phi_dpsi2
    type(expr_t) :: potential_phi1, potential_phi_reversed
    type(expr_t) :: potential_reversal_residual
    type(expr_t) :: flux_segment_s, flux_segment_s0, flux_segment_s1
    type(expr_t) :: flux_segment_psi0, flux_segment_psi1
    type(expr_t) :: flux_segment_field_scale, flux_segment_psi_sep
    type(expr_t) :: flux_segment_weight, flux_segment_psihat
    type(expr_t) :: flux_segment_dpsihat_dstor
    type(expr_t) :: flux_segment_target_psihat, flux_segment_inverse_s
    type(expr_t) :: flux_segment_axis_residual, flux_segment_edge_residual
    type(expr_t) :: flux_segment_derivative_residual
    type(expr_t) :: flux_segment_forward_inverse_residual
    type(expr_t) :: flux_segment_inverse_forward_residual
    type(expr_t) :: scaled_flux_value, scaled_flux_psihat
    type(expr_t) :: scaled_flux_round_trip_residual
    type(expr_t) :: physical_flux_normalized_value, physical_flux_psihat
    type(expr_t) :: physical_flux_normalization_residual
    type(expr_t) :: physical_flux_raw_psi, physical_flux_raw_psi_r
    type(expr_t) :: physical_flux_raw_psi_z, physical_flux_dz_dr
    type(expr_t) :: physical_flux_raw_dpsi_dr
    type(expr_t) :: physical_flux_field_scale, physical_flux_psi
    type(expr_t) :: physical_flux_dpsi_dr
    type(expr_t) :: allowed_field_scale, allowed_dz_dr, allowed_d2z_dr2
    type(expr_t) :: allowed_equilibrium_psi_physical
    type(expr_t) :: allowed_equilibrium_dpsi_physical_dr
    type(expr_t) :: allowed_g, allowed_sqrt_g, allowed_dpsi_dr
    type(expr_t) :: allowed_d2psi_dr2, allowed_dpsi_r_dr
    type(expr_t) :: allowed_dpsi_z_dr, allowed_d2psi_r_dr2
    type(expr_t) :: allowed_d2psi_z_dr2, allowed_df_dr, allowed_d2f_dr2
    type(expr_t) :: allowed_dg_dr, allowed_d2g_dr2
    type(expr_t) :: allowed_dsqrt_g_dr, allowed_d2sqrt_g_dr2
    type(expr_t) :: allowed_bmod, allowed_dbmod_dr, allowed_d2bmod_dr2
    type(expr_t) :: allowed_potential, allowed_dphi_dpsi
    type(expr_t) :: allowed_d2phi_dpsi2
    type(expr_t) :: allowed_jk
    type(expr_t) :: allowed_omega_c, allowed_domega_c_dr
    type(expr_t) :: allowed_d2omega_c_dr2
    type(expr_t) :: allowed_energy_margin, allowed_denergy_dr
    type(expr_t) :: allowed_d2energy_dr2
    type(expr_t) :: allowed_vparallel_squared, allowed_dvparallel_squared_dr
    type(expr_t) :: allowed_d2vparallel_squared_dr2
    type(expr_t) :: allowed_bphi_covariant, allowed_dbphi_covariant_dr
    type(expr_t) :: allowed_vparallel_squared_input
    type(expr_t) :: allowed_dvparallel_squared_dr_input
    type(expr_t) :: allowed_bphi_covariant_input
    type(expr_t) :: allowed_dbphi_covariant_dr_input
    type(expr_t) :: allowed_x, allowed_radius_local
    type(expr_t) :: allowed_psi_r_local, allowed_psi_z_local
    type(expr_t) :: allowed_f_local, allowed_g_local
    type(expr_t) :: allowed_bmod_local, allowed_potential_local
    type(expr_t) :: allowed_energy_local
    type(expr_t) :: allowed_g_second_series
    type(expr_t) :: allowed_potential_second_series
    type(expr_t) :: allowed_parallel_sign, allowed_vparallel
    type(expr_t) :: allowed_dvparallel_dr, allowed_psi_physical
    type(expr_t) :: allowed_dpsi_physical_dr, allowed_psi_star
    type(expr_t) :: allowed_dpsi_star_dr
    type(expr_t) :: turning_y, turning_direction, turning_dv2_dx
    type(expr_t) :: turning_psi_root, turning_dpsi_dx
    type(expr_t) :: turning_ratio_root, turning_dratio_dx
    type(expr_t) :: turning_delta_x, turning_vparallel
    type(expr_t) :: turning_psi_local, turning_ratio_local
    type(expr_t) :: turning_psi_star, turning_dpsi_star_dy
    type(expr_t) :: turning_root_derivative
    type(expr_t) :: canonical_x, canonical_sigma, canonical_kappa
    type(expr_t) :: canonical_f0, canonical_f1, canonical_f2, canonical_f3
    type(expr_t) :: canonical_g0, canonical_g1, canonical_g2, canonical_g3
    type(expr_t) :: canonical_p0, canonical_p1, canonical_p2
    type(expr_t) :: canonical_f_local, canonical_g_local, canonical_p_local
    type(expr_t) :: canonical_psi_local, canonical_psi_star_local
    type(expr_t) :: canonical_n_local, canonical_n, canonical_n_r
    type(expr_t) :: canonical_n_rr
    type(expr_t) :: canonical_turning_f_local, canonical_turning_g_local
    type(expr_t) :: canonical_turning_n_local
    type(expr_t) :: canonical_turning_n_root, canonical_turning_n_y_root
    type(expr_t) :: canonical_turning_n_root_residual
    type(expr_t) :: canonical_turning_n_y_root_residual
    type(expr_t) :: canonical_identity_residual
    type(expr_t) :: canonical_n_coordinate, canonical_n_r_coordinate
    type(expr_t) :: canonical_n_rr_coordinate
    type(expr_t) :: canonical_time_phase, canonical_time_phase_reversed
    type(expr_t) :: canonical_time_resonance
    type(expr_t) :: canonical_time_resonance_reversed
    type(expr_t) :: canonical_fourier_pair_phase
    type(expr_t) :: canonical_fourier_pair_resonance
    type(expr_t) :: canonical_combined_resonance
    type(expr_t) :: canonical_toroidal_phase
    type(expr_t) :: canonical_toroidal_phase_relabelled
    type(expr_t) :: canonical_p_phi_covariant
    type(expr_t) :: canonical_p_phi_covariant_relabelled
    type(expr_t) :: canonical_sigma_difference
    type(expr_t) :: canonical_sigma_identity_residual
    type(expr_t) :: canonical_coordinate_residual
    type(expr_t) :: canonical_coordinate_r_residual
    type(expr_t) :: canonical_coordinate_rr_residual
    type(expr_t) :: canonical_time_phase_residual
    type(expr_t) :: canonical_time_resonance_residual
    type(expr_t) :: canonical_fourier_pair_phase_residual
    type(expr_t) :: canonical_fourier_pair_resonance_residual
    type(expr_t) :: canonical_combined_resonance_residual
    type(expr_t) :: canonical_toroidal_phase_residual
    type(expr_t) :: canonical_toroidal_p_phi_residual
    type(expr_t) :: canonical_endpoint_delta_01, canonical_endpoint_delta_12
    type(expr_t) :: canonical_endpoint_variation
    type(expr_t) :: canonical_endpoint_variation_reversed
    type(expr_t) :: canonical_endpoint_variation_residual
    type(expr_t) :: cell_coefficient(6,6), cell_delta_r, cell_delta_z
    type(expr_t) :: cell_psi, cell_psi_r, cell_psi_z, cell_psi_rr
    type(expr_t) :: cell_psi_rz, cell_psi_zz, cell_psi_rrr
    type(expr_t) :: cell_psi_rrz, cell_psi_rzz, cell_psi_zzz
    type(expr_t) :: cell_psi_rrrr, cell_psi_rrrz, cell_psi_rrzz
    type(expr_t) :: cell_psi_rzzz, cell_psi_zzzz
    type(expr_t) :: profile_coefficient(0:5), profile_delta
    type(expr_t) :: profile_value, profile_first, profile_second
    type(expr_t) :: profile_btf, profile_rtf, profile_vacuum_f
    type(expr_t) :: eqcut_radius, eqcut_field_scale, eqcut_orientation
    type(expr_t) :: eqcut_psi0, eqcut_psi_r, eqcut_psi_z
    type(expr_t) :: eqcut_psi_rr, eqcut_psi_rz, eqcut_psi_zz
    type(expr_t) :: eqcut_psi_rrr, eqcut_psi_rrz
    type(expr_t) :: eqcut_psi_rzz, eqcut_psi_zzz
    type(expr_t) :: eqcut_psi_rrrr, eqcut_psi_rrrz, eqcut_psi_rrzz
    type(expr_t) :: eqcut_psi_rzzz, eqcut_psi_zzzz
    type(expr_t) :: eqcut_f0, eqcut_f_hat_first, eqcut_f_hat_second
    type(expr_t) :: eqcut_psi_sep, eqcut_dr, eqcut_darc_phi, eqcut_dz
    type(expr_t) :: eqcut_local_radius, eqcut_local_psi, eqcut_delta_psi
    type(expr_t) :: eqcut_local_f, eqcut_local_psi_r, eqcut_local_psi_z
    type(expr_t) :: eqcut_local_psi_rr, eqcut_local_psi_rz
    type(expr_t) :: eqcut_local_psi_zz, eqcut_local_g, eqcut_local_s
    type(expr_t) :: eqcut_local_k, eqcut_local_n
    type(expr_t) :: eqcut_n, eqcut_n_r, eqcut_n_z
    type(expr_t) :: eqcut_n_rr, eqcut_n_rz, eqcut_n_zz
    type(expr_t) :: eqcut_f_reversal_n_residual
    type(expr_t) :: eqcut_f_reversal_n_rr_residual
    type(expr_t) :: eqcut_f_reversal_n_rz_residual
    type(expr_t) :: eqcut_f_reversal_n_zz_residual
    type(expr_t) :: eqcut_direct_grad_b_r, eqcut_direct_grad_b_z
    type(expr_t) :: eqcut_direct_c
    type(expr_t) :: eqcut_local_c, eqcut_c, eqcut_dc_dr
    type(expr_t) :: eqcut_dc_darc_phi, eqcut_dc_dz
    type(expr_t) :: eqcut_dot_r, eqcut_dot_arc_phi, eqcut_dot_z
    type(expr_t) :: eqcut_cdot, eqcut_abs_cdot, eqcut_orientation_scalar
    type(expr_t) :: eqcut_g, eqcut_s, eqcut_k, eqcut_compact_c
    type(expr_t) :: eqcut_compact_residual, eqcut_emitted_residual
    type(expr_t) :: eqcut_f1_c_residual
    type(expr_t) :: eqcut_n_compact, eqcut_n_residual
    type(expr_t) :: eqcut_c_numerator_residual
    type(expr_t) :: eqcut_f1_n_residual
    type(expr_t) :: eqcut_f2_n_r_residual, eqcut_f2_n_z_residual
    type(expr_t) :: eqcut_f2_c_residual, eqcut_f2_cdot_residual
    type(expr_t) :: eqcut_reversed, eqcut_cdot_reversed
    type(expr_t) :: eqcut_scale_multiplier, eqcut_scaled
    type(expr_t) :: eqcut_midplane
    type(expr_t) :: eqcut_r_chart_slope, eqcut_r_chart_ds
    type(expr_t) :: eqcut_z_chart_slope, eqcut_z_chart_ds
    type(expr_t) :: eqcut_chart_n_r, eqcut_chart_n_z
    type(expr_t) :: eqcut_chart_psi_r, eqcut_chart_psi_z
    type(expr_t) :: eqcut_chart_psi_sep, eqcut_r_chart_dpsihat
    type(expr_t) :: eqcut_r_chart_tangent_residual
    type(expr_t) :: eqcut_z_chart_tangent_residual
    type(expr_t) :: eqcut_r_chart_flux_residual
    type(expr_t) :: eqcut_chart_n_rr, eqcut_chart_n_rz, eqcut_chart_n_zz
    type(expr_t) :: eqcut_chart_psi_rr, eqcut_chart_psi_rz
    type(expr_t) :: eqcut_chart_psi_zz, eqcut_r_chart_second
    type(expr_t) :: eqcut_r_chart_d2psihat
    type(expr_t) :: eqcut_r_chart_second_tangent_residual
    type(expr_t) :: eqcut_r_chart_flux_second_residual
    type(expr_t) :: eqcut_f2_hessian_contraction_residual
    type(expr_t) :: eqcut_r_chart_slope_reversed
    type(expr_t) :: eqcut_r_chart_second_reversed
    type(expr_t) :: eqcut_mean_n_mid, eqcut_mean_n_r, eqcut_mean_n_z
    type(expr_t) :: eqcut_mean_delta_r, eqcut_mean_delta_z
    type(expr_t) :: eqcut_mean_enclosure, eqcut_mean_center_residual
    type(expr_t) :: eqcut_axis_slope, eqcut_axis_psi_rr
    type(expr_t) :: eqcut_axis_psi_rz, eqcut_axis_psi_zz
    type(expr_t) :: eqcut_axis_psi_sep, eqcut_axis_curvature
    type(expr_t) :: eqcut_axis_hessian_determinant
    type(expr_t) :: eqcut_axis_delta_psihat
    type(expr_t) :: eqcut_axis_branch_sign, eqcut_axis_abs_dpsihat
    type(expr_t) :: eqcut_axis_inverse_delta_r
    type(expr_t) :: eqcut_axis_derivative_limit_residual
    type(expr_t) :: eqcut_axis_inverse_limit_residual
    type(expr_t) :: eqcut_rho, eqcut_s_tor, eqcut_dpsihat_dstor
    type(expr_t) :: eqcut_rho_dstor_drho, eqcut_coordinate_dstor_drho
    type(expr_t) :: eqcut_dpsihat_drho, eqcut_dR_drho, eqcut_dZ_drho
    type(expr_t) :: eqcut_coordinate_dpsihat_dR
    type(expr_t) :: eqcut_coordinate_dZ_dR
    type(expr_t) :: eqcut_coordinate_d2psihat_dstor2
    type(expr_t) :: eqcut_coordinate_d2psihat_dR2
    type(expr_t) :: eqcut_d2stor_drho2, eqcut_d2psihat_drho2
    type(expr_t) :: eqcut_d2R_drho2
    type(expr_t) :: eqcut_allowed_g_R, eqcut_allowed_g_RR
    type(expr_t) :: eqcut_allowed_g_rho, eqcut_allowed_g_rhorho
    type(expr_t) :: eqcut_allowed_psistar_R
    type(expr_t) :: eqcut_allowed_psistar_rho
    type(expr_t) :: eqcut_allowed_axis_dR_drho
    type(expr_t) :: eqcut_allowed_axis_g_rho
    type(expr_t) :: eqcut_allowed_axis_psistar_rho
    type(expr_t) :: eqcut_flux_coordinate_chain_residual
    type(expr_t) :: eqcut_flux_inverse_chain_residual
    type(expr_t) :: eqcut_flux_inverse_second_chain_residual
    type(expr_t) :: eqcut_allowed_first_chain_residual
    type(expr_t) :: eqcut_allowed_second_chain_residual
    type(expr_t) :: eqcut_allowed_psistar_chain_residual
    type(expr_t) :: eqcut_allowed_axis_g_chain_residual
    type(expr_t) :: eqcut_allowed_axis_psistar_chain_residual
    type(expr_t) :: eqcut_rho_map_derivative_residual
    type(expr_t) :: eqcut_s_tor_input, eqcut_rho_from_s
    type(expr_t) :: eqcut_s_tor_inverse_residual
    type(expr_t) :: eqcut_axis_rho_dR, eqcut_axis_rho_dZ
    type(expr_t) :: eqcut_axis_rho_limit_residual
    type(expr_t) :: endpoint_n, endpoint_n_r, endpoint_n_z
    type(expr_t) :: endpoint_psi, endpoint_psi_r, endpoint_psi_z
    type(expr_t) :: endpoint_psi_sep, endpoint_target_psihat
    type(expr_t) :: endpoint_flux_residual, endpoint_flux_r, endpoint_flux_z
    type(expr_t) :: endpoint_r0, endpoint_z0, endpoint_n0, endpoint_g0
    type(expr_t) :: endpoint_nr0, endpoint_nz0, endpoint_gr0, endpoint_gz0
    type(expr_t) :: endpoint_determinant, endpoint_y11, endpoint_y12
    type(expr_t) :: endpoint_y21, endpoint_y22
    type(expr_t) :: endpoint_newton_r, endpoint_newton_z
    type(expr_t) :: endpoint_nr_box, endpoint_nz_box
    type(expr_t) :: endpoint_gr_box, endpoint_gz_box
    type(expr_t) :: endpoint_delta_r_box, endpoint_delta_z_box
    type(expr_t) :: endpoint_krawczyk_r, endpoint_krawczyk_z
    type(expr_t) :: endpoint_oracle_delta_r, endpoint_oracle_delta_z
    type(expr_t) :: endpoint_oracle_n0, endpoint_oracle_g0
    type(expr_t) :: endpoint_newton_r_oracle, endpoint_newton_z_oracle
    type(expr_t) :: endpoint_krawczyk_r_oracle, endpoint_krawczyk_z_oracle
    type(expr_t) :: axis_stat_psi_r, axis_stat_psi_z
    type(expr_t) :: axis_stat_psi_rr, axis_stat_psi_rz, axis_stat_psi_zz
    type(expr_t) :: axis_stat_r0, axis_stat_z0
    type(expr_t) :: axis_stat_psi_r0, axis_stat_psi_z0
    type(expr_t) :: axis_stat_psi_rr0, axis_stat_psi_rz0, axis_stat_psi_zz0
    type(expr_t) :: axis_stat_determinant
    type(expr_t) :: axis_stat_y11, axis_stat_y12, axis_stat_y21, axis_stat_y22
    type(expr_t) :: axis_stat_newton_r, axis_stat_newton_z
    type(expr_t) :: axis_stat_psi_rr_box, axis_stat_psi_rz_box
    type(expr_t) :: axis_stat_psi_zz_box
    type(expr_t) :: axis_stat_delta_r_box, axis_stat_delta_z_box
    type(expr_t) :: axis_stat_krawczyk_r, axis_stat_krawczyk_z
    type(expr_t) :: axis_stat_oracle_delta_r, axis_stat_oracle_delta_z
    type(expr_t) :: axis_stat_oracle_psi_r0, axis_stat_oracle_psi_z0
    type(expr_t) :: axis_stat_newton_r_oracle, axis_stat_newton_z_oracle
    type(expr_t) :: axis_stat_krawczyk_r_oracle
    type(expr_t) :: axis_stat_krawczyk_z_oracle
    character(len=64) :: eqdsk_cell_arg_names(38)
    character(len=64) :: eqdsk_profile_arg_names(9)
    type(expr_t) :: axis_b_r, axis_b_phi, axis_b_z, axis_bhat_r
    type(expr_t) :: axis_bmod
    type(expr_t) :: axis_psi
    type(expr_t) :: axis_bhat_phi, axis_bhat_z, axis_dbphi_d_r
    type(expr_t) :: axis_dbphi_d_z, axis_force_r, axis_force_z
    type(expr_t) :: axis_curl_phi
    type(expr_t) :: axis_grad_b_r, axis_grad_b_z, axis_phi_psi
    type(expr_t) :: axis_bparallel_star, axis_dot_r, axis_dot_z
    type(expr_t) :: axis_dot_p, axis_pphi, axis_pphi_dot, axis_pphi_residual
    type(expr_t) :: axis_bstar_r, axis_bstar_phi, axis_bstar_z
    type(expr_t) :: dpsi_star_dx_contribution, abs_hm_squared
    type(expr_t) :: tau_b, frequency_contribution, phase_contribution
    type(expr_t) :: f_prime_contribution
    type(expr_t) :: dpsi_abs, hm_squared_positive
    type(expr_t) :: frequency_contribution_positive
    type(expr_t) :: phase_contribution_positive
    type(expr_t) :: frequency_identity_roots(2)
    type(expr_t) :: simple_root_force_roots(3)
    type(expr_t) :: n_squared_frequency_contribution
    type(expr_t) :: phase_contribution_identity
    type(expr_t) :: frequency_residual_derivative, root_force_value
    type(expr_t) :: root_hm_squared, simple_root_weight
    type(expr_t) :: simple_root_force_contribution
    type(expr_t) :: simple_root_phase_identity
    type(expr_t) :: e_ref, n0, residence, eq17_outer_factor
    type(expr_t) :: eq17_outer_unit_phi
    type(expr_t) :: gauge_constant, eq17_outer_shift, eq17_gauge_residual
    type(expr_t) :: eq17_shifted_energy
    type(expr_t) :: gauge_outer_ledger_residual
    type(expr_t) :: gauge_force_shift, gauge_force_residual
    type(expr_t) :: section_root_residual, crossing_orientation_residual
    type(expr_t) :: sign_omega_b, sign_omega_phi, sign_torque_phi
    type(expr_t) :: sign_rotation_phi, toroidal_resonance
    type(expr_t) :: toroidal_resonance_relabelled
    type(expr_t) :: toroidal_resonance_residual
    type(expr_t) :: toroidal_power, toroidal_power_relabelled
    type(expr_t) :: toroidal_power_residual, torque_component_difference
    type(expr_t) :: charge_reversal_energy_difference
    type(expr_t) :: potential_reversal_energy_difference
    type(expr_t) :: norm_velocity, norm_jk, norm_psi_star
    type(expr_t) :: norm_dpsi_star_dx, norm_tau_b, norm_omega_b
    type(expr_t) :: norm_omega_phi, norm_domega_b_dx, norm_domega_phi_dx
    type(expr_t) :: norm_hm_real, norm_hm_imag
    type(expr_t) :: norm_phi_eff, norm_jk_scale, norm_h_hat, norm_jk_hat
    type(expr_t) :: norm_psi_hat, norm_dpsi_hat_dx, norm_tau_hat
    type(expr_t) :: norm_omega_b_hat, norm_omega_phi_hat
    type(expr_t) :: norm_domega_b_hat_dx, norm_domega_phi_hat_dx
    type(expr_t) :: norm_hm_real_hat, norm_hm_imag_hat
    type(expr_t) :: quad_t, quad_h_weight, quad_j_unit, quad_j_weight
    type(expr_t) :: quad_qphi_min, quad_jk_max, quad_h_hat, quad_h_phys
    type(expr_t) :: quad_dh_hat_dt, quad_weight_h, quad_jk_hat_max
    type(expr_t) :: quad_jk_hat, quad_jk_phys, quad_weight_j, quad_weight
    type(expr_t) :: gauss_standard_node, gauss_standard_weight
    type(expr_t) :: gauss_lower_bound, gauss_upper_bound
    type(expr_t) :: gauss_half_span, gauss_midpoint
    type(expr_t) :: gauss_mapped_node, gauss_mapped_weight
    type(expr_t) :: poly_c0, poly_c1, poly_c2, poly_c3, poly_c4, poly_c5
    type(expr_t) :: poly_width, poly_tail_bound, poly_lower, poly_upper
    type(expr_t) :: conversion_vperp, conversion_mu, conversion_jk
    type(expr_t) :: conversion_jpotato, conversion_v0
    type(expr_t) :: mu_from_vperp, vperp_squared_from_mu
    type(expr_t) :: specific_energy_from_mu, jk_from_mu, mu_from_jk
    type(expr_t) :: conversion_omega_c, vperp_squared_from_jk
    type(expr_t) :: specific_energy_from_jk, jpotato_from_vperp
    type(expr_t) :: vperp_squared_from_jpotato, jpotato_from_mu
    type(expr_t) :: mu_from_jpotato
    type(expr_t) :: invariant_h, invariant_pphi, invariant_psistar
    type(expr_t) :: invariant_vparallel_squared, invariant_ppar_from_pphi
    type(expr_t) :: invariant_ppar_squared, invariant_launch_residual
    type(engine_result_t) :: simplified
    type(engine_result_t) :: resonance_series
    type(engine_result_t) :: allowed_g_second_series_result
    type(engine_result_t) :: allowed_potential_second_series_result
    character(2048) :: output_path
    integer :: argument_status, output_length, k, i, j

    call get_command_argument(1, output_path, length=output_length, &
        status=argument_status)
    if (argument_status /= 0 .or. output_length == 0) then
        output_path = DEFAULT_OUTPUT
    else
        output_path = output_path(:output_length)
    end if

    call arena%init()
    proof_engine = make_symengine_engine(arena)
    simplify_engine = make_native_engine(arena)

    mass = sym(arena, "mass")
    charge = sym(arena, "charge")
    c_light = sym(arena, "c_light")
    mu = sym(arena, "mu")
    bmod = sym(arena, "bmod")
    h = sym(arena, "h")
    electrostatic_potential = sym(arena, "electrostatic_potential")
    p_phi = sym(arena, "p_phi")
    psi = sym(arena, "psi")
    b1 = sym(arena, "b1")
    b2 = sym(arena, "b2")
    b3 = sym(arena, "b3")
    bhat1 = sym(arena, "bhat1")
    bhat2 = sym(arena, "bhat2")
    bhat3 = sym(arena, "bhat3")
    cylindrical_bhat_phi = sym(arena, "bhat_phi")
    curl1 = sym(arena, "curl_bhat1")
    curl2 = sym(arena, "curl_bhat2")
    curl3 = sym(arena, "curl_bhat3")
    gradb1 = sym(arena, "grad_b1")
    gradb2 = sym(arena, "grad_b2")
    gradb3 = sym(arena, "grad_b3")
    gradphi1 = sym(arena, "grad_phi1")
    gradphi2 = sym(arena, "grad_phi2")
    gradphi3 = sym(arena, "grad_phi3")
    radius = sym(arena, "radius")
    p_parallel = sym(arena, "p_parallel")
    n_mode = sym(arena, "n_mode")
    tau = sym(arena, "tau")
    g_prime = sym(arena, "g_prime")
    temperature = sym(arena, "temperature")
    a1 = sym(arena, "a1")
    a2 = sym(arena, "a2")
    phi_eff = sym(arena, "phi_eff")
    a_real = sym(arena, "a_real")
    a_imag = sym(arena, "a_imag")
    m_mode = sym(arena, "m_mode")
    theta = sym(arena, "theta")
    phi = sym(arena, "phi")
    phase_half_real = sym(arena, "phase_half_real")
    phase_half_imag = sym(arena, "phase_half_imag")
    zero = num(arena, 0)
    one = num(arena, 1)

    ! ------------------------------------------------------------------
    ! Exact libneo tensor-quintic cell jet.  The runtime chooses the owning
    ! cell and passes its 6x6 s2dcut coefficients; every interpolation and
    ! derivative operation is generated here.  Powers correspond exactly to
    ! spl(i,j)*delta_R**(i-1)*delta_Z**(j-1).
    cell_delta_r = sym(arena, "delta_R")
    cell_delta_z = sym(arena, "delta_Z")
    eqdsk_cell_arg_names(1:2) = [character(len=64) :: "delta_R", "delta_Z"]
    cell_psi = zero
    do i = 1, 6
        do j = 1, 6
            k = 2+(i-1)*6+j
            eqdsk_cell_arg_names(k) = &
                EQDSK_CELL_COEFFICIENT_NAMES((i-1)*6+j)
            cell_coefficient(i,j) = &
                sym(arena, trim(eqdsk_cell_arg_names(k)))
            cell_psi = cell_psi + cell_coefficient(i,j)* &
                cell_delta_r**(i-1)*cell_delta_z**(j-1)
        end do
    end do
    cell_psi_r = diff(cell_psi, cell_delta_r)
    cell_psi_z = diff(cell_psi, cell_delta_z)
    cell_psi_rr = diff(cell_psi_r, cell_delta_r)
    cell_psi_rz = diff(cell_psi_r, cell_delta_z)
    cell_psi_zz = diff(cell_psi_z, cell_delta_z)
    cell_psi_rrr = diff(cell_psi_rr, cell_delta_r)
    cell_psi_rrz = diff(cell_psi_rr, cell_delta_z)
    cell_psi_rzz = diff(cell_psi_rz, cell_delta_z)
    cell_psi_zzz = diff(cell_psi_zz, cell_delta_z)
    cell_psi_rrrr = diff(cell_psi_rrr, cell_delta_r)
    cell_psi_rrrz = diff(cell_psi_rrr, cell_delta_z)
    cell_psi_rrzz = diff(cell_psi_rrz, cell_delta_z)
    cell_psi_rzzz = diff(cell_psi_rzz, cell_delta_z)
    cell_psi_zzzz = diff(cell_psi_zzz, cell_delta_z)

    ! Exact quintic F(psihat) profile jet in the same coefficient convention.
    profile_delta = sym(arena, "profile_delta")
    eqdsk_profile_arg_names(1) = "profile_delta"
    profile_value = zero
    do i = 0, 5
        eqdsk_profile_arg_names(i+2) = EQDSK_PROFILE_COEFFICIENT_NAMES(i+1)
        profile_coefficient(i) = &
            sym(arena, trim(eqdsk_profile_arg_names(i+2)))
        profile_value = profile_value + &
            profile_coefficient(i)*profile_delta**i
    end do
    profile_first = diff(profile_value, profile_delta)
    profile_second = diff(profile_first, profile_delta)
    profile_btf = sym(arena, "vacuum_Bphi_reference")
    profile_rtf = sym(arena, "vacuum_radius_reference")
    profile_vacuum_f = profile_btf*profile_rtf
    eqdsk_profile_arg_names(8:9) = &
        [character(len=64) :: "vacuum_Bphi_reference", &
        "vacuum_radius_reference"]

    ! Exact axisymmetric Eq. 13 cut jet from a fourth-order local psi jet.  For
    ! R>0, write |B|=field_scale*sqrt(G)/R with
    ! G=psi_R**2+psi_Z**2+F**2.  Reducing
    ! (grad(B) cross grad(psi)).grad(phi) before differentiating gives
    !
    ! C = orientation*field_scale**2 *
    !     (K/(R**2*sqrt(G)) + psi_Z*sqrt(G)/R**3),
    ! K = (psi_R**2-psi_Z**2)*psi_RZ
    !     + psi_R*psi_Z*(psi_ZZ-psi_RR).
    !
    ! The independent direct contraction below proves this reduction and the
    ! cancellation of dF/dpsi from C.  Fortsym then differentiates the reduced
    ! local expression to generate grad(C); d2F/dpsi2 is retained in the local
    ! profile solely to prove that it cancels from the emitted first jet.
    eqcut_radius = sym(arena, "radius")
    eqcut_field_scale = sym(arena, "field_scale")
    eqcut_orientation = sym(arena, "cut_orientation")
    eqcut_psi0 = sym(arena, "psi")
    eqcut_psi_r = sym(arena, "psi_R")
    eqcut_psi_z = sym(arena, "psi_Z")
    eqcut_psi_rr = sym(arena, "psi_RR")
    eqcut_psi_rz = sym(arena, "psi_RZ")
    eqcut_psi_zz = sym(arena, "psi_ZZ")
    eqcut_psi_rrr = sym(arena, "psi_RRR")
    eqcut_psi_rrz = sym(arena, "psi_RRZ")
    eqcut_psi_rzz = sym(arena, "psi_RZZ")
    eqcut_psi_zzz = sym(arena, "psi_ZZZ")
    eqcut_psi_rrrr = sym(arena, "psi_RRRR")
    eqcut_psi_rrrz = sym(arena, "psi_RRRZ")
    eqcut_psi_rrzz = sym(arena, "psi_RRZZ")
    eqcut_psi_rzzz = sym(arena, "psi_RZZZ")
    eqcut_psi_zzzz = sym(arena, "psi_ZZZZ")
    eqcut_f0 = sym(arena, "F")
    eqcut_f_hat_first = sym(arena, "dF_dpsihat")
    eqcut_f_hat_second = sym(arena, "d2F_dpsihat2")
    eqcut_psi_sep = sym(arena, "psi_sep")
    eqcut_dot_r = sym(arena, "dot_R")
    ! The orbit state carries angular dot_phi.  Convert it to the physical
    ! orthonormal arc velocity R*dot_phi inside the generated contract.
    eqcut_dot_arc_phi = eqcut_radius*sym(arena, "dot_phi")
    eqcut_dot_z = sym(arena, "dot_Z")
    eqcut_dr = sym(arena, "local_delta_R")
    eqcut_darc_phi = sym(arena, "local_delta_arc_phi")
    eqcut_dz = sym(arena, "local_delta_Z")
    eqcut_local_radius = eqcut_radius+eqcut_dr
    eqcut_local_psi = eqcut_psi0 + eqcut_psi_r*eqcut_dr + &
        eqcut_psi_z*eqcut_dz + &
        (eqcut_psi_rr*eqcut_dr**2 + &
        2*eqcut_psi_rz*eqcut_dr*eqcut_dz + &
        eqcut_psi_zz*eqcut_dz**2)/2 + &
        (eqcut_psi_rrr*eqcut_dr**3 + &
        3*eqcut_psi_rrz*eqcut_dr**2*eqcut_dz + &
        3*eqcut_psi_rzz*eqcut_dr*eqcut_dz**2 + &
        eqcut_psi_zzz*eqcut_dz**3)/6 + &
        (eqcut_psi_rrrr*eqcut_dr**4 + &
        4*eqcut_psi_rrrz*eqcut_dr**3*eqcut_dz + &
        6*eqcut_psi_rrzz*eqcut_dr**2*eqcut_dz**2 + &
        4*eqcut_psi_rzzz*eqcut_dr*eqcut_dz**3 + &
        eqcut_psi_zzzz*eqcut_dz**4)/24
    eqcut_delta_psi = eqcut_local_psi-eqcut_psi0
    eqcut_local_f = eqcut_f0 + &
        eqcut_f_hat_first/eqcut_psi_sep*eqcut_delta_psi + &
        eqcut_f_hat_second/(2*eqcut_psi_sep**2)*eqcut_delta_psi**2
    eqcut_local_psi_r = exact_derivative(eqcut_local_psi, eqcut_dr, &
        "local EQDSK psi_R")
    eqcut_local_psi_z = exact_derivative(eqcut_local_psi, eqcut_dz, &
        "local EQDSK psi_Z")
    eqcut_local_psi_rr = exact_derivative(eqcut_local_psi_r, eqcut_dr, &
        "local EQDSK psi_RR")
    eqcut_local_psi_rz = exact_derivative(eqcut_local_psi_r, eqcut_dz, &
        "local EQDSK psi_RZ")
    eqcut_local_psi_zz = exact_derivative(eqcut_local_psi_z, eqcut_dz, &
        "local EQDSK psi_ZZ")
    eqcut_local_g = eqcut_local_psi_r**2 + eqcut_local_psi_z**2 + &
        eqcut_local_f**2
    eqcut_local_s = sqrt(eqcut_local_g)
    eqcut_local_k = &
        (eqcut_local_psi_r**2-eqcut_local_psi_z**2)*eqcut_local_psi_rz + &
        eqcut_local_psi_r*eqcut_local_psi_z* &
        (eqcut_local_psi_zz-eqcut_local_psi_rr)
    eqcut_local_c = eqcut_orientation*eqcut_field_scale**2* &
        (eqcut_local_k/(eqcut_local_radius**2*eqcut_local_s) + &
        eqcut_local_psi_z*eqcut_local_s/eqcut_local_radius**3)
    eqcut_c = exact_simplify(subs(subs(eqcut_local_c, eqcut_dr, zero), &
        eqcut_dz, zero), "EQDSK cut value")
    ! This is the denominator-free numerator of the local Eq. 13 field
    ! expression.  Fortsym differentiates it before the point substitution,
    ! so its R/Z jets use the supplied third-order psi jet.
    eqcut_local_n = eqcut_local_radius*eqcut_local_k + &
        eqcut_local_psi_z*eqcut_local_g
    eqcut_n = exact_simplify(subs(subs(eqcut_local_n, eqcut_dr, zero), &
        eqcut_dz, zero), "EQDSK cut numerator")
    eqcut_n_r = exact_simplify(subs(subs(exact_derivative(eqcut_local_n, &
        eqcut_dr, "EQDSK cut numerator radial derivative"), &
        eqcut_dr, zero), eqcut_dz, zero), &
        "EQDSK cut numerator radial jet")
    eqcut_n_z = exact_simplify(subs(subs(exact_derivative(eqcut_local_n, &
        eqcut_dz, "EQDSK cut numerator vertical derivative"), &
        eqcut_dr, zero), eqcut_dz, zero), &
        "EQDSK cut numerator vertical jet")
    eqcut_n_rr = exact_simplify(subs(subs(exact_derivative( &
        exact_derivative(eqcut_local_n, eqcut_dr, &
        "EQDSK cut numerator radial derivative for Hessian"), eqcut_dr, &
        "EQDSK cut numerator RR derivative"), eqcut_dr, zero), &
        eqcut_dz, zero), "EQDSK cut numerator RR jet")
    eqcut_n_rz = exact_simplify(subs(subs(exact_derivative( &
        exact_derivative(eqcut_local_n, eqcut_dr, &
        "EQDSK cut numerator radial derivative for mixed Hessian"), &
        eqcut_dz, "EQDSK cut numerator RZ derivative"), eqcut_dr, zero), &
        eqcut_dz, zero), "EQDSK cut numerator RZ jet")
    eqcut_n_zz = exact_simplify(subs(subs(exact_derivative( &
        exact_derivative(eqcut_local_n, eqcut_dz, &
        "EQDSK cut numerator vertical derivative for Hessian"), eqcut_dz, &
        "EQDSK cut numerator ZZ derivative"), eqcut_dr, zero), &
        eqcut_dz, zero), "EQDSK cut numerator ZZ jet")
    eqcut_f_reversal_n_residual = subs(subs(subs(eqcut_n, &
        eqcut_f0, -eqcut_f0), eqcut_f_hat_first, -eqcut_f_hat_first), &
        eqcut_f_hat_second, -eqcut_f_hat_second)-eqcut_n
    eqcut_f_reversal_n_rr_residual = subs(subs(subs(eqcut_n_rr, &
        eqcut_f0, -eqcut_f0), eqcut_f_hat_first, -eqcut_f_hat_first), &
        eqcut_f_hat_second, -eqcut_f_hat_second)-eqcut_n_rr
    eqcut_f_reversal_n_rz_residual = subs(subs(subs(eqcut_n_rz, &
        eqcut_f0, -eqcut_f0), eqcut_f_hat_first, -eqcut_f_hat_first), &
        eqcut_f_hat_second, -eqcut_f_hat_second)-eqcut_n_rz
    eqcut_f_reversal_n_zz_residual = subs(subs(subs(eqcut_n_zz, &
        eqcut_f0, -eqcut_f0), eqcut_f_hat_first, -eqcut_f_hat_first), &
        eqcut_f_hat_second, -eqcut_f_hat_second)-eqcut_n_zz
    eqcut_f1_n_residual = diff(eqcut_n, eqcut_f_hat_first)
    eqcut_f2_n_r_residual = diff(eqcut_n_r, eqcut_f_hat_second)
    eqcut_f2_n_z_residual = diff(eqcut_n_z, eqcut_f_hat_second)
    eqcut_n = subs(eqcut_n, eqcut_f_hat_second, zero)
    eqcut_n_r = subs(eqcut_n_r, eqcut_f_hat_second, zero)
    eqcut_n_z = subs(eqcut_n_z, eqcut_f_hat_second, zero)
    eqcut_g = eqcut_psi_r**2+eqcut_psi_z**2+eqcut_f0**2
    eqcut_s = sqrt(eqcut_g)
    eqcut_k = (eqcut_psi_r**2-eqcut_psi_z**2)*eqcut_psi_rz + &
        eqcut_psi_r*eqcut_psi_z*(eqcut_psi_zz-eqcut_psi_rr)
    eqcut_n_compact = eqcut_radius*eqcut_k + eqcut_psi_z*eqcut_g
    eqcut_n_residual = eqcut_n-eqcut_n_compact
    eqcut_compact_c = eqcut_orientation*eqcut_field_scale**2* &
        (eqcut_k/(eqcut_radius**2*eqcut_s) + &
        eqcut_psi_z*eqcut_s/eqcut_radius**3)
    ! This is the unreduced definition with the axisymmetric field components
    ! B_R=-psi_Z/R, B_phi=F/R, B_Z=psi_R/R already substituted.  Its F'
    ! terms cancel only in the antisymmetric contraction; that cancellation is
    ! therefore proved independently instead of assumed by the emitted form.
    eqcut_direct_grad_b_r = eqcut_field_scale* &
        ((eqcut_psi_r*eqcut_psi_rr + &
        eqcut_psi_z*eqcut_psi_rz + &
        eqcut_f0*eqcut_f_hat_first/eqcut_psi_sep*eqcut_psi_r)/ &
        (eqcut_radius*eqcut_s) - eqcut_s/eqcut_radius**2)
    eqcut_direct_grad_b_z = eqcut_field_scale* &
        (eqcut_psi_r*eqcut_psi_rz + &
        eqcut_psi_z*eqcut_psi_zz + &
        eqcut_f0*eqcut_f_hat_first/eqcut_psi_sep*eqcut_psi_z)/ &
        (eqcut_radius*eqcut_s)
    eqcut_direct_c = eqcut_orientation* &
        (eqcut_direct_grad_b_z*eqcut_field_scale*eqcut_psi_r - &
        eqcut_direct_grad_b_r*eqcut_field_scale*eqcut_psi_z)/eqcut_radius
    eqcut_compact_residual = eqcut_direct_c-eqcut_compact_c
    eqcut_emitted_residual = eqcut_c-eqcut_compact_c
    eqcut_c_numerator_residual = eqcut_c*eqcut_radius**3*eqcut_s - &
        eqcut_orientation*eqcut_field_scale**2*eqcut_n
    eqcut_f1_c_residual = diff(eqcut_direct_c, eqcut_f_hat_first)
    eqcut_dc_dr = exact_simplify(subs(subs(exact_derivative( &
        eqcut_local_c, eqcut_dr, "EQDSK cut radial derivative"), &
        eqcut_dr, zero), eqcut_dz, zero), "EQDSK cut radial jet")
    eqcut_dc_darc_phi = exact_simplify(subs(subs(subs( &
        exact_derivative(eqcut_local_c, eqcut_darc_phi, &
        "EQDSK cut toroidal derivative"), eqcut_dr, zero), &
        eqcut_darc_phi, zero), eqcut_dz, zero), &
        "EQDSK cut toroidal jet")
    eqcut_dc_dz = exact_simplify(subs(subs(exact_derivative( &
        eqcut_local_c, eqcut_dz, "EQDSK cut vertical derivative"), &
        eqcut_dr, zero), eqcut_dz, zero), "EQDSK cut vertical jet")
    eqcut_cdot = eqcut_dc_dr*eqcut_dot_r + &
        eqcut_dc_darc_phi*eqcut_dot_arc_phi + eqcut_dc_dz*eqcut_dot_z
    eqcut_f2_c_residual = diff(eqcut_c, eqcut_f_hat_second)
    eqcut_f2_cdot_residual = diff(eqcut_cdot, eqcut_f_hat_second)
    ! The exact proofs below certify these residuals as zero.  Remove the
    ! mathematically irrelevant symbol explicitly so the runtime interface is
    ! minimal even if a backend simplifier retains a cancelled subexpression.
    eqcut_c = subs(eqcut_c, eqcut_f_hat_second, zero)
    eqcut_dc_dr = subs(eqcut_dc_dr, eqcut_f_hat_second, zero)
    eqcut_dc_darc_phi = subs(eqcut_dc_darc_phi, eqcut_f_hat_second, zero)
    eqcut_dc_dz = subs(eqcut_dc_dz, eqcut_f_hat_second, zero)
    eqcut_cdot = eqcut_dc_dr*eqcut_dot_r + &
        eqcut_dc_darc_phi*eqcut_dot_arc_phi + eqcut_dc_dz*eqcut_dot_z
    eqcut_abs_cdot = abs(eqcut_cdot)
    eqcut_orientation_scalar = eqcut_cdot
    eqcut_reversed = subs(eqcut_c, eqcut_orientation, -eqcut_orientation)
    eqcut_cdot_reversed = subs(eqcut_cdot, eqcut_orientation, &
        -eqcut_orientation)
    eqcut_scale_multiplier = sym(arena, "field_scale_multiplier")
    eqcut_scaled = subs(eqcut_c, eqcut_field_scale, &
        eqcut_scale_multiplier*eqcut_field_scale)
    eqcut_midplane = subs(subs(eqcut_c, eqcut_psi_z, zero), &
        eqcut_psi_rz, zero)

    ! These pure chart kernels are guarded by their callers: the R chart is
    ! used only when N_Z is nonzero, and the Z chart only when N_R is
    ! nonzero.  Keeping the divisions in separate routines prevents a
    ! tangent point in one chart from evaluating the other denominator.
    eqcut_chart_n_r = sym(arena, "N_R")
    eqcut_chart_n_z = sym(arena, "N_Z")
    eqcut_chart_psi_r = sym(arena, "psi_R")
    eqcut_chart_psi_z = sym(arena, "psi_Z")
    eqcut_chart_psi_sep = sym(arena, "psi_sep")
    eqcut_r_chart_slope = -eqcut_chart_n_r/eqcut_chart_n_z
    eqcut_r_chart_ds = sqrt(one+eqcut_r_chart_slope**2)
    eqcut_r_chart_dpsihat = (eqcut_chart_psi_r + &
        eqcut_chart_psi_z*eqcut_r_chart_slope)/eqcut_chart_psi_sep
    eqcut_z_chart_slope = -eqcut_chart_n_z/eqcut_chart_n_r
    eqcut_z_chart_ds = sqrt(one+eqcut_z_chart_slope**2)
    eqcut_r_chart_tangent_residual = eqcut_chart_n_r + &
        eqcut_chart_n_z*eqcut_r_chart_slope
    eqcut_z_chart_tangent_residual = &
        eqcut_chart_n_r*eqcut_z_chart_slope+eqcut_chart_n_z
    eqcut_r_chart_flux_residual = eqcut_chart_psi_sep* &
        eqcut_r_chart_dpsihat - (eqcut_chart_psi_r + &
        eqcut_chart_psi_z*eqcut_r_chart_slope)
    eqcut_chart_n_rr = sym(arena, "N_RR")
    eqcut_chart_n_rz = sym(arena, "N_RZ")
    eqcut_chart_n_zz = sym(arena, "N_ZZ")
    eqcut_chart_psi_rr = sym(arena, "psi_RR")
    eqcut_chart_psi_rz = sym(arena, "psi_RZ")
    eqcut_chart_psi_zz = sym(arena, "psi_ZZ")
    eqcut_r_chart_second = -(eqcut_chart_n_rr + &
        2*eqcut_chart_n_rz*eqcut_r_chart_slope + &
        eqcut_chart_n_zz*eqcut_r_chart_slope**2)/eqcut_chart_n_z
    eqcut_r_chart_d2psihat = (eqcut_chart_psi_rr + &
        2*eqcut_chart_psi_rz*eqcut_r_chart_slope + &
        eqcut_chart_psi_zz*eqcut_r_chart_slope**2 + &
        eqcut_chart_psi_z*eqcut_r_chart_second)/eqcut_chart_psi_sep
    eqcut_r_chart_second_tangent_residual = eqcut_chart_n_rr + &
        2*eqcut_chart_n_rz*eqcut_r_chart_slope + &
        eqcut_chart_n_zz*eqcut_r_chart_slope**2 + &
        eqcut_chart_n_z*eqcut_r_chart_second
    eqcut_r_chart_flux_second_residual = eqcut_chart_psi_sep &
        *eqcut_r_chart_d2psihat-(eqcut_chart_psi_rr + &
        2*eqcut_chart_psi_rz*eqcut_r_chart_slope + &
        eqcut_chart_psi_zz*eqcut_r_chart_slope**2 + &
        eqcut_chart_psi_z*eqcut_r_chart_second)
    eqcut_f2_hessian_contraction_residual = diff(eqcut_n_rr + &
        2*eqcut_n_rz*eqcut_r_chart_slope + &
        eqcut_n_zz*eqcut_r_chart_slope**2, eqcut_f_hat_second) &
        -2*eqcut_psi_z*eqcut_f0/eqcut_psi_sep**2 &
        *(eqcut_psi_r+eqcut_psi_z*eqcut_r_chart_slope)**2
    eqcut_r_chart_slope_reversed = subs(subs(eqcut_r_chart_slope, &
        eqcut_chart_n_r, -eqcut_chart_n_r), &
        eqcut_chart_n_z, -eqcut_chart_n_z)
    eqcut_r_chart_second_reversed = subs(subs(subs(subs(subs( &
        eqcut_r_chart_second, eqcut_chart_n_r, -eqcut_chart_n_r), &
        eqcut_chart_n_z, -eqcut_chart_n_z), &
        eqcut_chart_n_rr, -eqcut_chart_n_rr), &
        eqcut_chart_n_rz, -eqcut_chart_n_rz), &
        eqcut_chart_n_zz, -eqcut_chart_n_zz)

    ! Mean-value inclusion extension used by the certified atlas.  The
    ! runtime supplies outward enclosures of N_R and N_Z over the complete
    ! box and centered coordinate intervals, so the multivariate mean-value
    ! theorem makes this affine expression an enclosure of N over that box.
    eqcut_mean_n_mid = sym(arena, "N_mid")
    eqcut_mean_n_r = sym(arena, "N_R_box")
    eqcut_mean_n_z = sym(arena, "N_Z_box")
    eqcut_mean_delta_r = sym(arena, "delta_R")
    eqcut_mean_delta_z = sym(arena, "delta_Z")
    eqcut_mean_enclosure = eqcut_mean_n_mid + &
        eqcut_mean_n_r*eqcut_mean_delta_r + &
        eqcut_mean_n_z*eqcut_mean_delta_z
    eqcut_mean_center_residual = subs(subs(eqcut_mean_enclosure, &
        eqcut_mean_delta_r, zero), eqcut_mean_delta_z, zero) - &
        eqcut_mean_n_mid

    ! Nondegenerate magnetic-axis limit along the certified R-chart.  At the
    ! stationary flux point grad(psi)=0, the graph curvature term containing
    ! d2Z/dR2 vanishes, leaving the Hessian contraction below.  The inverse
    ! square-root map is the explicit limiting coordinate required when the
    ! two monotone cut branches meet at the axis.
    eqcut_axis_slope = sym(arena, "dZ_dR")
    eqcut_axis_psi_rr = sym(arena, "psi_RR")
    eqcut_axis_psi_rz = sym(arena, "psi_RZ")
    eqcut_axis_psi_zz = sym(arena, "psi_ZZ")
    eqcut_axis_psi_sep = sym(arena, "psi_sep")
    eqcut_axis_delta_psihat = sym(arena, "delta_psihat")
    eqcut_axis_branch_sign = sym(arena, "branch_sign")
    eqcut_axis_curvature = (eqcut_axis_psi_rr + &
        2*eqcut_axis_psi_rz*eqcut_axis_slope + &
        eqcut_axis_psi_zz*eqcut_axis_slope**2)/eqcut_axis_psi_sep
    eqcut_axis_hessian_determinant = eqcut_axis_psi_rr &
        *eqcut_axis_psi_zz-eqcut_axis_psi_rz**2
    eqcut_axis_abs_dpsihat = sqrt(2*eqcut_axis_curvature &
        *eqcut_axis_delta_psihat)
    eqcut_axis_inverse_delta_r = eqcut_axis_branch_sign &
        *sqrt(2*eqcut_axis_delta_psihat/eqcut_axis_curvature)
    eqcut_axis_derivative_limit_residual = &
        eqcut_axis_abs_dpsihat**2 &
        - 2*eqcut_axis_curvature*eqcut_axis_delta_psihat
    eqcut_axis_inverse_limit_residual = &
        eqcut_axis_curvature*eqcut_axis_inverse_delta_r**2/2 &
        - eqcut_axis_branch_sign**2*eqcut_axis_delta_psihat

    ! Axis-safe class coordinate rho_tor=sqrt(s_tor).  Away from the axis,
    ! the generated chain converts the certified cut R-chart derivative to
    ! d(position)/d(rho_tor).  At rho_tor=0 the two factors proportional to
    ! rho and 1/sqrt(psihat) cancel, leaving the finite branch limit below.
    eqcut_rho = sym(arena, "rho_tor")
    eqcut_dpsihat_dstor = sym(arena, "dpsihat_dstor")
    eqcut_coordinate_dpsihat_dR = sym(arena, "dpsihat_dR")
    eqcut_coordinate_dZ_dR = sym(arena, "dZ_dR")
    eqcut_coordinate_dstor_drho = sym(arena, "dstor_drho_tor")
    eqcut_coordinate_d2psihat_dstor2 = sym(arena, "d2psihat_dstor2")
    eqcut_coordinate_d2psihat_dR2 = sym(arena, "d2psihat_dR2")
    eqcut_allowed_g_R = sym(arena, "dvparallel_squared_dR")
    eqcut_allowed_g_RR = sym(arena, "d2vparallel_squared_dR2")
    eqcut_allowed_psistar_R = sym(arena, "dpsi_star_dR")
    eqcut_allowed_axis_dR_drho = sym(arena, "dR_drho_tor_limit")
    eqcut_s_tor = eqcut_rho**2
    eqcut_rho_dstor_drho = 2*eqcut_rho
    eqcut_d2stor_drho2 = 2*one
    eqcut_dpsihat_drho = eqcut_coordinate_dstor_drho*eqcut_dpsihat_dstor
    eqcut_d2psihat_drho2 = eqcut_coordinate_d2psihat_dstor2* &
        eqcut_coordinate_dstor_drho**2 + &
        eqcut_dpsihat_dstor*eqcut_d2stor_drho2
    eqcut_dR_drho = eqcut_dpsihat_drho/eqcut_coordinate_dpsihat_dR
    eqcut_d2R_drho2 = (eqcut_d2psihat_drho2 - &
        eqcut_coordinate_d2psihat_dR2*eqcut_dR_drho**2)/ &
        eqcut_coordinate_dpsihat_dR
    eqcut_dZ_drho = eqcut_coordinate_dZ_dR*eqcut_dR_drho
    eqcut_allowed_g_rho = eqcut_allowed_g_R*eqcut_dR_drho
    eqcut_allowed_g_rhorho = eqcut_allowed_g_RR*eqcut_dR_drho**2 + &
        eqcut_allowed_g_R*eqcut_d2R_drho2
    eqcut_allowed_psistar_rho = eqcut_allowed_psistar_R*eqcut_dR_drho
    eqcut_flux_coordinate_chain_residual = eqcut_dpsihat_drho &
        -eqcut_coordinate_dstor_drho*eqcut_dpsihat_dstor
    eqcut_flux_inverse_chain_residual = eqcut_coordinate_dpsihat_dR &
        *eqcut_dR_drho-eqcut_dpsihat_drho
    eqcut_flux_inverse_second_chain_residual = &
        eqcut_coordinate_dpsihat_dR*eqcut_d2R_drho2 + &
        eqcut_coordinate_d2psihat_dR2*eqcut_dR_drho**2 - &
        eqcut_d2psihat_drho2
    eqcut_allowed_first_chain_residual = eqcut_allowed_g_rho - &
        eqcut_allowed_g_R*eqcut_dR_drho
    eqcut_allowed_second_chain_residual = eqcut_allowed_g_rhorho - &
        (eqcut_allowed_g_RR*eqcut_dR_drho**2 + &
        eqcut_allowed_g_R*eqcut_d2R_drho2)
    eqcut_allowed_psistar_chain_residual = eqcut_allowed_psistar_rho - &
        eqcut_allowed_psistar_R*eqcut_dR_drho
    eqcut_axis_rho_dR = eqcut_axis_branch_sign &
        *sqrt(2*eqcut_dpsihat_dstor/eqcut_axis_curvature)
    eqcut_axis_rho_dZ = eqcut_axis_slope*eqcut_axis_rho_dR
    eqcut_allowed_axis_g_rho = &
        eqcut_allowed_g_R*eqcut_allowed_axis_dR_drho
    eqcut_allowed_axis_psistar_rho = &
        eqcut_allowed_psistar_R*eqcut_allowed_axis_dR_drho
    eqcut_allowed_axis_g_chain_residual = eqcut_allowed_axis_g_rho - &
        eqcut_allowed_g_R*eqcut_allowed_axis_dR_drho
    eqcut_allowed_axis_psistar_chain_residual = &
        eqcut_allowed_axis_psistar_rho - &
        eqcut_allowed_psistar_R*eqcut_allowed_axis_dR_drho
    eqcut_axis_rho_limit_residual = eqcut_axis_curvature &
        *eqcut_axis_rho_dR**2 &
        -2*eqcut_axis_branch_sign**2*eqcut_dpsihat_dstor
    eqcut_rho_map_derivative_residual = diff(eqcut_s_tor, eqcut_rho) &
        -eqcut_rho_dstor_drho
    eqcut_s_tor_input = sym(arena, "s_tor")
    eqcut_rho_from_s = sqrt(eqcut_s_tor_input)
    eqcut_s_tor_inverse_residual = eqcut_rho_from_s**2 &
        -eqcut_s_tor_input

    ! Certified endpoint chart for the coupled system N=0 and
    ! psi/psi_sep=target.  A generated Newton map supplies the point inverse;
    ! its interval Krawczyk inclusion is the existence-and-uniqueness gate.
    endpoint_n = sym(arena, "cut_numerator")
    endpoint_n_r = sym(arena, "d_cut_numerator_d_R")
    endpoint_n_z = sym(arena, "d_cut_numerator_d_Z")
    endpoint_psi = sym(arena, "psi")
    endpoint_psi_r = sym(arena, "psi_R")
    endpoint_psi_z = sym(arena, "psi_Z")
    endpoint_psi_sep = sym(arena, "psi_sep")
    endpoint_target_psihat = sym(arena, "target_psihat")
    endpoint_flux_residual = endpoint_psi/endpoint_psi_sep &
        -endpoint_target_psihat
    endpoint_flux_r = endpoint_psi_r/endpoint_psi_sep
    endpoint_flux_z = endpoint_psi_z/endpoint_psi_sep
    endpoint_r0 = sym(arena, "R0")
    endpoint_z0 = sym(arena, "Z0")
    endpoint_n0 = sym(arena, "cut_numerator0")
    endpoint_g0 = sym(arena, "flux_residual0")
    endpoint_nr0 = sym(arena, "d_cut_numerator_d_R0")
    endpoint_nz0 = sym(arena, "d_cut_numerator_d_Z0")
    endpoint_gr0 = sym(arena, "d_flux_residual_d_R0")
    endpoint_gz0 = sym(arena, "d_flux_residual_d_Z0")
    endpoint_determinant = endpoint_nr0*endpoint_gz0 &
        -endpoint_nz0*endpoint_gr0
    endpoint_y11 = endpoint_gz0/endpoint_determinant
    endpoint_y12 = -endpoint_nz0/endpoint_determinant
    endpoint_y21 = -endpoint_gr0/endpoint_determinant
    endpoint_y22 = endpoint_nr0/endpoint_determinant
    endpoint_newton_r = endpoint_r0-endpoint_y11*endpoint_n0 &
        -endpoint_y12*endpoint_g0
    endpoint_newton_z = endpoint_z0-endpoint_y21*endpoint_n0 &
        -endpoint_y22*endpoint_g0
    endpoint_nr_box = sym(arena, "d_cut_numerator_d_R_box")
    endpoint_nz_box = sym(arena, "d_cut_numerator_d_Z_box")
    endpoint_gr_box = sym(arena, "d_flux_residual_d_R_box")
    endpoint_gz_box = sym(arena, "d_flux_residual_d_Z_box")
    endpoint_delta_r_box = sym(arena, "delta_R_box")
    endpoint_delta_z_box = sym(arena, "delta_Z_box")
    endpoint_krawczyk_r = endpoint_newton_r &
        +(one-endpoint_y11*endpoint_nr_box &
        -endpoint_y12*endpoint_gr_box)*endpoint_delta_r_box &
        +(-endpoint_y11*endpoint_nz_box &
        -endpoint_y12*endpoint_gz_box)*endpoint_delta_z_box
    endpoint_krawczyk_z = endpoint_newton_z &
        +(-endpoint_y21*endpoint_nr_box &
        -endpoint_y22*endpoint_gr_box)*endpoint_delta_r_box &
        +(one-endpoint_y21*endpoint_nz_box &
        -endpoint_y22*endpoint_gz_box)*endpoint_delta_z_box
    endpoint_oracle_delta_r = sym(arena, "oracle_delta_R")
    endpoint_oracle_delta_z = sym(arena, "oracle_delta_Z")
    endpoint_oracle_n0 = endpoint_nr0*endpoint_oracle_delta_r &
        +endpoint_nz0*endpoint_oracle_delta_z
    endpoint_oracle_g0 = endpoint_gr0*endpoint_oracle_delta_r &
        +endpoint_gz0*endpoint_oracle_delta_z
    endpoint_newton_r_oracle = subs(subs(endpoint_newton_r, endpoint_n0, &
        endpoint_oracle_n0), endpoint_g0, endpoint_oracle_g0)
    endpoint_newton_z_oracle = subs(subs(endpoint_newton_z, endpoint_n0, &
        endpoint_oracle_n0), endpoint_g0, endpoint_oracle_g0)
    endpoint_krawczyk_r_oracle = endpoint_krawczyk_r
    endpoint_krawczyk_z_oracle = endpoint_krawczyk_z
    endpoint_krawczyk_r_oracle = subs(endpoint_krawczyk_r_oracle, &
        endpoint_n0, endpoint_oracle_n0)
    endpoint_krawczyk_r_oracle = subs(endpoint_krawczyk_r_oracle, &
        endpoint_g0, endpoint_oracle_g0)
    endpoint_krawczyk_r_oracle = subs(endpoint_krawczyk_r_oracle, &
        endpoint_nr_box, endpoint_nr0)
    endpoint_krawczyk_r_oracle = subs(endpoint_krawczyk_r_oracle, &
        endpoint_nz_box, endpoint_nz0)
    endpoint_krawczyk_r_oracle = subs(endpoint_krawczyk_r_oracle, &
        endpoint_gr_box, endpoint_gr0)
    endpoint_krawczyk_r_oracle = subs(endpoint_krawczyk_r_oracle, &
        endpoint_gz_box, endpoint_gz0)
    endpoint_krawczyk_z_oracle = subs(endpoint_krawczyk_z_oracle, &
        endpoint_n0, endpoint_oracle_n0)
    endpoint_krawczyk_z_oracle = subs(endpoint_krawczyk_z_oracle, &
        endpoint_g0, endpoint_oracle_g0)
    endpoint_krawczyk_z_oracle = subs(endpoint_krawczyk_z_oracle, &
        endpoint_nr_box, endpoint_nr0)
    endpoint_krawczyk_z_oracle = subs(endpoint_krawczyk_z_oracle, &
        endpoint_nz_box, endpoint_nz0)
    endpoint_krawczyk_z_oracle = subs(endpoint_krawczyk_z_oracle, &
        endpoint_gr_box, endpoint_gr0)
    endpoint_krawczyk_z_oracle = subs(endpoint_krawczyk_z_oracle, &
        endpoint_gz_box, endpoint_gz0)

    ! Magnetic-axis point chart: grad(psi)=0.  The generated Hessian inverse
    ! and Krawczyk map turn the existing convex axis box into a numerical root
    ! enclosure without promoting a sampled minimum to a proof.
    axis_stat_psi_r = sym(arena, "psi_R")
    axis_stat_psi_z = sym(arena, "psi_Z")
    axis_stat_psi_rr = sym(arena, "psi_RR")
    axis_stat_psi_rz = sym(arena, "psi_RZ")
    axis_stat_psi_zz = sym(arena, "psi_ZZ")
    axis_stat_r0 = sym(arena, "R0")
    axis_stat_z0 = sym(arena, "Z0")
    axis_stat_psi_r0 = sym(arena, "psi_R0")
    axis_stat_psi_z0 = sym(arena, "psi_Z0")
    axis_stat_psi_rr0 = sym(arena, "psi_RR0")
    axis_stat_psi_rz0 = sym(arena, "psi_RZ0")
    axis_stat_psi_zz0 = sym(arena, "psi_ZZ0")
    axis_stat_determinant = axis_stat_psi_rr0*axis_stat_psi_zz0 &
        -axis_stat_psi_rz0**2
    axis_stat_y11 = axis_stat_psi_zz0/axis_stat_determinant
    axis_stat_y12 = -axis_stat_psi_rz0/axis_stat_determinant
    axis_stat_y21 = axis_stat_y12
    axis_stat_y22 = axis_stat_psi_rr0/axis_stat_determinant
    axis_stat_newton_r = axis_stat_r0-axis_stat_y11*axis_stat_psi_r0 &
        -axis_stat_y12*axis_stat_psi_z0
    axis_stat_newton_z = axis_stat_z0-axis_stat_y21*axis_stat_psi_r0 &
        -axis_stat_y22*axis_stat_psi_z0
    axis_stat_psi_rr_box = sym(arena, "psi_RR_box")
    axis_stat_psi_rz_box = sym(arena, "psi_RZ_box")
    axis_stat_psi_zz_box = sym(arena, "psi_ZZ_box")
    axis_stat_delta_r_box = sym(arena, "delta_R_box")
    axis_stat_delta_z_box = sym(arena, "delta_Z_box")
    axis_stat_krawczyk_r = axis_stat_newton_r &
        +(one-axis_stat_y11*axis_stat_psi_rr_box &
        -axis_stat_y12*axis_stat_psi_rz_box)*axis_stat_delta_r_box &
        +(-axis_stat_y11*axis_stat_psi_rz_box &
        -axis_stat_y12*axis_stat_psi_zz_box)*axis_stat_delta_z_box
    axis_stat_krawczyk_z = axis_stat_newton_z &
        +(-axis_stat_y21*axis_stat_psi_rr_box &
        -axis_stat_y22*axis_stat_psi_rz_box)*axis_stat_delta_r_box &
        +(one-axis_stat_y21*axis_stat_psi_rz_box &
        -axis_stat_y22*axis_stat_psi_zz_box)*axis_stat_delta_z_box
    axis_stat_oracle_delta_r = sym(arena, "oracle_delta_R")
    axis_stat_oracle_delta_z = sym(arena, "oracle_delta_Z")
    axis_stat_oracle_psi_r0 = &
        axis_stat_psi_rr0*axis_stat_oracle_delta_r &
        +axis_stat_psi_rz0*axis_stat_oracle_delta_z
    axis_stat_oracle_psi_z0 = &
        axis_stat_psi_rz0*axis_stat_oracle_delta_r &
        +axis_stat_psi_zz0*axis_stat_oracle_delta_z
    axis_stat_newton_r_oracle = subs(subs(axis_stat_newton_r, &
        axis_stat_psi_r0, axis_stat_oracle_psi_r0), axis_stat_psi_z0, &
        axis_stat_oracle_psi_z0)
    axis_stat_newton_z_oracle = subs(subs(axis_stat_newton_z, &
        axis_stat_psi_r0, axis_stat_oracle_psi_r0), axis_stat_psi_z0, &
        axis_stat_oracle_psi_z0)
    axis_stat_krawczyk_r_oracle = axis_stat_krawczyk_r
    axis_stat_krawczyk_z_oracle = axis_stat_krawczyk_z
    axis_stat_krawczyk_r_oracle = subs(axis_stat_krawczyk_r_oracle, &
        axis_stat_psi_r0, axis_stat_oracle_psi_r0)
    axis_stat_krawczyk_r_oracle = subs(axis_stat_krawczyk_r_oracle, &
        axis_stat_psi_z0, axis_stat_oracle_psi_z0)
    axis_stat_krawczyk_r_oracle = subs(axis_stat_krawczyk_r_oracle, &
        axis_stat_psi_rr_box, axis_stat_psi_rr0)
    axis_stat_krawczyk_r_oracle = subs(axis_stat_krawczyk_r_oracle, &
        axis_stat_psi_rz_box, axis_stat_psi_rz0)
    axis_stat_krawczyk_r_oracle = subs(axis_stat_krawczyk_r_oracle, &
        axis_stat_psi_zz_box, axis_stat_psi_zz0)
    axis_stat_krawczyk_z_oracle = subs(axis_stat_krawczyk_z_oracle, &
        axis_stat_psi_r0, axis_stat_oracle_psi_r0)
    axis_stat_krawczyk_z_oracle = subs(axis_stat_krawczyk_z_oracle, &
        axis_stat_psi_z0, axis_stat_oracle_psi_z0)
    axis_stat_krawczyk_z_oracle = subs(axis_stat_krawczyk_z_oracle, &
        axis_stat_psi_rr_box, axis_stat_psi_rr0)
    axis_stat_krawczyk_z_oracle = subs(axis_stat_krawczyk_z_oracle, &
        axis_stat_psi_rz_box, axis_stat_psi_rz0)
    axis_stat_krawczyk_z_oracle = subs(axis_stat_krawczyk_z_oracle, &
        axis_stat_psi_zz_box, axis_stat_psi_zz0)

    ! ------------------------------------------------------------------
    ! Positive action, cyclotron frequency, and exact phase-space candidate.
    q_abs = abs(charge)
    q_phi_energy = charge*electrostatic_potential
    jk = mass*c_light*mu/q_abs
    omega_c = q_abs*bmod/(mass*c_light)
    jk_omega_c = jk*omega_c
    jk_candidate = (h - q_phi_energy)/omega_c
    jk_max = func("max", [zero, jk_candidate])
    d_jk_d_h = diff(jk_candidate, h)
    d_jk_d_phi = diff(jk_candidate, electrostatic_potential)
    omega_c_var = sym(arena, "omega_c_value")
    d_jk_d_omega = diff((h - q_phi_energy)/omega_c_var, omega_c_var)
    eq4_coeff = 2*(h - q_phi_energy) - jk_omega_c
    psi_star = c_light/charge*p_phi
    dpsi_dp_phi = diff(psi_star, p_phi)

    ! Frequency delta weights.  Derive F' from an x-dependent local model at
    ! the explicit root x_r.  The root side condition is g(x_r)=0; locally
    ! g(x)=g'(x_r)(x-x_r), tau(x)=tau_r+tau'(x_r)(x-x_r).  The emitted
    ! runtime value is the exact first Taylor coefficient at x_r, not a
    ! definition inserted after the fact.
    resonance_x = sym(arena, "x")
    resonance_xr = sym(arena, "x_r")
    ! Keep the emitted symbol identical to the public kernel argument.  A
    ! different local spelling leaves an undeclared symbol in generated code
    ! even though the Fortsym identity itself proves successfully.
    g_prime_x = g_prime
    tau_prime_x = sym(arena, "tau_prime_at_xr")
    resonance_delta = resonance_x - resonance_xr
    g_local = g_prime_x*resonance_delta
    tau_local = tau + tau_prime_x*resonance_delta
    f_local = n_mode*g_local/tau_local
    g_at_root = g_prime_x*(resonance_xr-resonance_xr)
    resonance_series = simplify_engine%series_coeff(f_local, resonance_x, &
        resonance_xr, 1)
    if (.not. resonance_series%ok) then
        error stop "fortsym could not derive resonance root derivative"
    end if
    f_prime_from_series = resonance_series%value
    f_prime = f_prime_from_series
    freq_weight = n_mode**2*tau/abs(f_prime)
    phase_weight = abs(n_mode)*tau**2/abs(g_prime)

    ! Eq. 17 force bracket and its single outer Phi_eff factor.
    force_bracket = a1 + (h - q_phi_energy)/temperature*a2
    e_ref = sym(arena, "Eref")
    n0 = sym(arena, "n0")
    residence = sym(arena, "residence")
    eq17_outer_unit_phi = -(pi_expr(arena)**rat(arena,3_int64,2_int64))/4* &
        e_ref*n0/(temperature/e_ref)**rat(arena,3_int64,2_int64)* &
        exp((q_phi_energy-h)/temperature)*residence
    eq17_outer_factor = phi_eff*eq17_outer_unit_phi
    gauge_constant = sym(arena, "gauge_C")
    eq17_shifted_energy = exact_expand( &
        charge*(electrostatic_potential+gauge_constant) - &
        (h+charge*gauge_constant), "Eq17 gauge-shifted energy")
    eq17_outer_shift = -(pi_expr(arena)**rat(arena,3_int64,2_int64))/4* &
        e_ref*phi_eff*n0/(temperature/e_ref)**rat(arena,3_int64,2_int64)* &
        exp(eq17_shifted_energy/temperature)*residence
    eq17_gauge_residual = eq17_outer_shift - eq17_outer_factor
    ! The complete outer-factor equality is proved above.  The runtime ledger
    ! emits its defining gauge-invariant exponent identity so it has no hidden
    ! Eref/n0/Phi_eff/residence inputs.
    gauge_outer_ledger_residual = &
        (charge*(electrostatic_potential+gauge_constant) - &
        (h+charge*gauge_constant)) - &
        (charge*electrostatic_potential-h)
    gauge_force_shift = a1 + ((h+charge*gauge_constant) - &
        charge*(electrostatic_potential+gauge_constant))/temperature*a2
    gauge_force_residual = gauge_force_shift-force_bracket

    ! Eq. 17--18 normalization and paired-domain Jacobians.  Consumers pass
    ! physical values and use these emitted conversions exactly once.
    norm_velocity = sym(arena, "reference_velocity")
    norm_jk = sym(arena, "J_K_physical")
    norm_psi_star = sym(arena, "psi_star_physical")
    norm_dpsi_star_dx = sym(arena, "dpsi_star_dx_physical")
    norm_tau_b = sym(arena, "tau_b_physical")
    norm_omega_b = sym(arena, "omega_b_physical")
    norm_omega_phi = sym(arena, "omega_phi_physical")
    norm_domega_b_dx = sym(arena, "domega_b_dx_physical")
    norm_domega_phi_dx = sym(arena, "domega_phi_dx_physical")
    norm_hm_real = sym(arena, "Hm_real_physical")
    norm_hm_imag = sym(arena, "Hm_imag_physical")
    norm_phi_eff = c_light*e_ref/(abs(charge)*norm_velocity)
    norm_jk_scale = mass*c_light*e_ref/abs(charge)
    norm_h_hat = h/e_ref
    norm_jk_hat = norm_jk/norm_jk_scale
    norm_psi_hat = norm_psi_star/norm_phi_eff
    norm_dpsi_hat_dx = norm_dpsi_star_dx/norm_phi_eff
    norm_tau_hat = norm_velocity*norm_tau_b
    norm_omega_b_hat = norm_omega_b/norm_velocity
    norm_omega_phi_hat = norm_omega_phi/norm_velocity
    norm_domega_b_hat_dx = norm_domega_b_dx/norm_velocity
    norm_domega_phi_hat_dx = norm_domega_phi_dx/norm_velocity
    norm_hm_real_hat = norm_hm_real/e_ref
    norm_hm_imag_hat = norm_hm_imag/e_ref

    quad_t = sym(arena, "energy_unit_node")
    quad_h_weight = sym(arena, "energy_unit_weight")
    quad_j_unit = sym(arena, "action_unit_node")
    quad_j_weight = sym(arena, "action_unit_weight")
    quad_qphi_min = sym(arena, "qPhi_min")
    quad_jk_max = sym(arena, "J_K_max_physical")
    quad_h_hat = quad_qphi_min/e_ref+quad_t/(1-quad_t)
    quad_h_phys = e_ref*quad_h_hat
    quad_dh_hat_dt = diff(quad_h_hat, quad_t)
    quad_weight_h = quad_h_weight*quad_dh_hat_dt
    quad_jk_hat_max = quad_jk_max/norm_jk_scale
    quad_jk_hat = quad_jk_hat_max*quad_j_unit
    quad_jk_phys = norm_jk_scale*quad_jk_hat
    quad_weight_j = quad_j_weight*quad_jk_hat_max
    quad_weight = quad_weight_h*quad_weight_j

    ! The Gauss-Legendre implementation supplies its canonical rule on
    ! [-1,1].  Derive the complete affine push-forward, including the
    ! measure Jacobian, rather than duplicating that normalization in the
    ! runtime consumer.
    gauss_standard_node = sym(arena, "standard_node")
    gauss_standard_weight = sym(arena, "standard_weight")
    gauss_lower_bound = sym(arena, "lower_bound")
    gauss_upper_bound = sym(arena, "upper_bound")
    gauss_half_span = (gauss_upper_bound-gauss_lower_bound)/2
    gauss_midpoint = (gauss_upper_bound+gauss_lower_bound)/2
    gauss_mapped_node = gauss_midpoint + &
        gauss_half_span*gauss_standard_node
    gauss_mapped_weight = gauss_half_span*gauss_standard_weight

    ! Conservative degree-five local-polynomial enclosure.  The inequality
    ! is the triangle bound supplied by the interval provider; Fortsym emits
    ! the radius and endpoints so signed c0 values are handled identically.
    poly_c0 = sym(arena, "coefficient_0")
    poly_c1 = sym(arena, "coefficient_1")
    poly_c2 = sym(arena, "coefficient_2")
    poly_c3 = sym(arena, "coefficient_3")
    poly_c4 = sym(arena, "coefficient_4")
    poly_c5 = sym(arena, "coefficient_5")
    poly_width = sym(arena, "cell_width")
    poly_tail_bound = abs(poly_c1)*poly_width + &
        abs(poly_c2)*poly_width**2 + abs(poly_c3)*poly_width**3 + &
        abs(poly_c4)*poly_width**4 + abs(poly_c5)*poly_width**5
    poly_lower = poly_c0-poly_tail_bound
    poly_upper = poly_c0+poly_tail_bound

    ! ------------------------------------------------------------------
    ! Conversion kernels at the three invariant boundaries.  These are kept
    ! separate so a runtime caller never evaluates an unrelated division.
    ! The wrappers own positivity/nonzero guards; Fortsym owns every formula.
    conversion_vperp = sym(arena, "v_perp")
    conversion_mu = sym(arena, "mu_phys")
    conversion_jk = sym(arena, "J_K")
    conversion_jpotato = sym(arena, "J_potato")
    conversion_v0 = sym(arena, "reference_velocity")
    mu_from_vperp = mass*conversion_vperp**2/(2*bmod)
    vperp_squared_from_mu = 2*conversion_mu*bmod/mass
    specific_energy_from_mu = conversion_mu*bmod/mass
    jk_from_mu = mass*c_light*conversion_mu/abs(charge)
    mu_from_jk = abs(charge)*conversion_jk/(mass*c_light)
    conversion_omega_c = abs(charge)*bmod/(mass*c_light)
    vperp_squared_from_jk = 2*conversion_jk*conversion_omega_c/mass
    specific_energy_from_jk = conversion_jk*conversion_omega_c/mass
    jpotato_from_vperp = conversion_vperp**2/(conversion_v0**2*bmod)
    vperp_squared_from_jpotato = conversion_jpotato*conversion_v0**2*bmod
    jpotato_from_mu = 2*conversion_mu/(mass*conversion_v0**2)
    mu_from_jpotato = mass*conversion_v0**2*conversion_jpotato/2

    ! Axisymmetric invariant and launch algebra.  The signed charge remains
    ! in q*Phi and psi_star; only J_K conversion above uses abs(q).
    invariant_h = p_parallel**2/(2*mass) + mu*bmod + q_phi_energy
    invariant_pphi = charge/c_light*psi + p_parallel*radius*bhat2
    invariant_psistar = c_light/charge*invariant_pphi
    invariant_vparallel_squared = 2*(h-mu*bmod-q_phi_energy)/mass
    invariant_ppar_from_pphi = (p_phi-charge*psi/c_light)/(radius*bhat2)
    invariant_ppar_squared = mass**2*invariant_vparallel_squared
    invariant_launch_residual = invariant_ppar_from_pphi**2 - &
        invariant_ppar_squared

    ! ------------------------------------------------------------------
    ! Central cylindrical Littlejohn kernel.  Components use (R,phi,Z)
    ! ordering; the public orbit derivative ordering is (R,Z,phi,p_parallel,mu).
    bstar1 = b1 + c_light/charge*p_parallel*curl1
    bstar2 = b2 + c_light/charge*p_parallel*curl2
    bstar3 = b3 + c_light/charge*p_parallel*curl3
    bparallel_star = bhat1*bstar1 + bhat2*bstar2 + bhat3*bstar3
    force1 = mu*gradb1 + charge*gradphi1
    force2 = mu*gradb2 + charge*gradphi2
    force3 = mu*gradb3 + charge*gradphi3
    cross1 = bhat2*force3 - bhat3*force2
    cross2 = bhat3*force1 - bhat1*force3
    cross3 = bhat1*force2 - bhat2*force1
    v1 = (p_parallel/mass*bstar1 + c_light/charge*cross1) &
        /bparallel_star
    v2 = (p_parallel/mass*bstar2 + c_light/charge*cross2) &
        /bparallel_star
    v3 = (p_parallel/mass*bstar3 + c_light/charge*cross3) &
        /bparallel_star
    dot_r = v1
    dot_z = v3
    dot_phi = v2/radius
    dot_p = -(bstar1*force1 + bstar2*force2 + bstar3*force3) &
        /bparallel_star
    dot_mu = zero
    cylindrical_measure = radius*bparallel_star
    hamiltonian_dot = p_parallel/mass*dot_p + force1*v1 + force2*v2 + force3*v3

    ! Genuine axisymmetric Noether identity along the generated Littlejohn
    ! dynamics.  Here B_R/B_phi/B_Z are full field components, while bhat_phi
    ! is the unit-vector component.  The flux convention is psi_R=R*B_Z and
    ! psi_Z=-R*B_R; curl_R=-d(bhat_phi)/dZ and
    ! curl_Z=d(bhat_phi)/dR+bhat_phi/R.
    axis_bmod = sym(arena, "Bmod_axis")
    axis_bhat_r = sym(arena, "bhat_R_axis")
    axis_bhat_phi = sym(arena, "bhat_phi_axis")
    axis_bhat_z = sym(arena, "bhat_Z_axis")
    axis_b_r = axis_bmod*axis_bhat_r
    axis_b_phi = axis_bmod*axis_bhat_phi
    axis_b_z = axis_bmod*axis_bhat_z
    axis_dbphi_d_r = sym(arena, "d_bhat_phi_d_R_axis")
    axis_dbphi_d_z = sym(arena, "d_bhat_phi_d_Z_axis")
    axis_curl_phi = sym(arena, "curl_bhat_phi_axis")
    axis_grad_b_r = sym(arena, "grad_b_R_axis")
    axis_grad_b_z = sym(arena, "grad_b_Z_axis")
    axis_phi_psi = sym(arena, "Phi_psi_axis")
    axis_bstar_r = axis_b_r - c_light/charge*p_parallel*axis_dbphi_d_z
    axis_bstar_phi = axis_b_phi + c_light/charge*p_parallel*axis_curl_phi
    axis_bstar_z = axis_b_z + c_light/charge*p_parallel* &
        (axis_dbphi_d_r + axis_bhat_phi/radius)
    axis_bparallel_star = axis_bhat_r*axis_bstar_r + &
        axis_bhat_phi*axis_bstar_phi + axis_bhat_z*axis_bstar_z
    axis_force_r = mu*axis_grad_b_r + charge*axis_phi_psi*radius*axis_b_z
    axis_force_z = mu*axis_grad_b_z - charge*axis_phi_psi*radius*axis_b_r
    axis_dot_r = (p_parallel/mass*axis_bstar_r + &
        c_light/charge*axis_bhat_phi*axis_force_z)/axis_bparallel_star
    axis_dot_z = (p_parallel/mass*axis_bstar_z - &
        c_light/charge*axis_bhat_phi*axis_force_r)/axis_bparallel_star
    axis_dot_p = -(axis_bstar_r*axis_force_r + &
        axis_bstar_z*axis_force_z)/axis_bparallel_star
    axis_psi = sym(arena, "psi_axis")
    axis_pphi = charge/c_light*axis_psi + &
        p_parallel*radius*axis_bhat_phi
    axis_pphi_dot = charge/c_light*(radius*axis_b_z*axis_dot_r - &
        radius*axis_b_r*axis_dot_z) + &
        p_parallel*(axis_bhat_phi + radius*axis_dbphi_d_r)*axis_dot_r + &
        p_parallel*radius*axis_dbphi_d_z*axis_dot_z + &
        radius*axis_bhat_phi*axis_dot_p
    axis_pphi_residual = axis_pphi_dot-zero

    ! On-shell perturbation limit C/E = 2 - eta B.
    energy_on_shell = p_parallel**2/(2*mass) + mu*bmod
    eta_on_shell = mu/energy_on_shell
    c_on_shell = 2*energy_on_shell - jk_omega_c
    perturbation_ratio = c_on_shell/energy_on_shell
    boozer_ratio = 2 - eta_on_shell*bmod

    ! Canonical temporal Fourier integrand for a single real-field toroidal
    ! amplitude A(R,Z), with delta B=Re[A exp(i*n*phi)].  The direct orbit
    ! state uses physical time and cylindrical phi; no Boozer angle or q
    ! reduction enters this phase.
    harmonic_delta_b_real = sym(arena, "delta_b_real")
    harmonic_delta_b_imag = sym(arena, "delta_b_imag")
    harmonic_phi = sym(arena, "orbit_phi")
    harmonic_phi_launch = sym(arena, "launch_phi")
    harmonic_time = sym(arena, "orbit_time")
    harmonic_omega_b = sym(arena, "omega_b")
    harmonic_omega_phi = sym(arena, "omega_phi")
    harmonic_phase_argument = n_mode*(harmonic_phi-harmonic_phi_launch) - &
        (m_mode*harmonic_omega_b+n_mode*harmonic_omega_phi)*harmonic_time
    harmonic_phase_real = cos(harmonic_phase_argument)
    harmonic_phase_imag = sin(harmonic_phase_argument)
    harmonic_shell_residual = eq4_coeff-c_on_shell
    harmonic_shell_identity = &
        2*((energy_on_shell+q_phi_energy)-q_phi_energy) - mu*bmod - &
        c_on_shell
    harmonic_factor = eq4_coeff/bmod
    harmonic_integrand_real = harmonic_factor* &
        (harmonic_delta_b_real*harmonic_phase_real - &
        harmonic_delta_b_imag*harmonic_phase_imag)
    harmonic_integrand_imag = harmonic_factor* &
        (harmonic_delta_b_real*harmonic_phase_imag + &
        harmonic_delta_b_imag*harmonic_phase_real)
    harmonic_complete_phase = n_mode*(harmonic_omega_phi*tau) - &
        (m_mode*(2*pi_expr(arena)/tau) + &
        n_mode*harmonic_omega_phi)*tau
    harmonic_complete_phase_residual = &
        harmonic_complete_phase + 2*pi_expr(arena)*m_mode
    harmonic_modulus_residual = harmonic_integrand_real**2 + &
        harmonic_integrand_imag**2 - harmonic_factor**2* &
        (harmonic_delta_b_real**2+harmonic_delta_b_imag**2)

    ! ------------------------------------------------------------------
    ! Real-field phase, sign, conjugation, and simultaneous (m,n) reversal.
    chi = m_mode*theta + n_mode*phi
    field_original = a_real*cos(chi) - a_imag*sin(chi)
    ! Parameterize the complete unit circle by homogeneous half-phase
    ! coordinates (u:v), excluding only the invalid pair (0:0).  This is the
    ! exact algebraic form of exp(i*beta) for
    ! (u,v) proportional to (cos(beta/2),sin(beta/2)); unlike a trigonometric
    ! rewrite rule, its unit norm and the SO(2) relabelling are both proved by
    ! Fortsym over the same expressions.
    phase_half_norm = phase_half_real**2+phase_half_imag**2
    phase_cos = (phase_half_real**2-phase_half_imag**2)/phase_half_norm
    phase_sin = 2*phase_half_real*phase_half_imag/phase_half_norm
    phase_unit_residual = exact_expand(phase_cos**2+phase_sin**2-one, &
        "projective phase unit norm")
    chi_cos = cos(chi)
    chi_sin = sin(chi)
    shifted_chi_cos = chi_cos*phase_cos + chi_sin*phase_sin
    shifted_chi_sin = chi_sin*phase_cos - chi_cos*phase_sin
    field_phase_shifted = &
        (a_real*phase_cos-a_imag*phase_sin)*shifted_chi_cos - &
        (a_real*phase_sin+a_imag*phase_cos)*shifted_chi_sin
    phase_relabel_residual = exact_expand( &
        field_phase_shifted-field_original, &
        "real-field phase relabelling residual")
    field_sign_sum = (-a_real)*cos(chi) - (-a_imag)*sin(chi) &
        + field_original
    field_reversal = a_real*cos(-chi) - (-a_imag)*sin(-chi)
    field_fixed_conjugate = a_real*cos(chi) - (-a_imag)*sin(chi)
    fixed_conjugation_difference = field_fixed_conjugate-field_original
    amp_square = a_real**2 + a_imag**2
    amp_square_sign = (-a_real)**2 + (-a_imag)**2
    amp_square_conj = a_real**2 + (-a_imag)**2
    amplitude_phase = sym(arena, "amplitude_phase")
    amp_rotated_real = a_real*cos(amplitude_phase) - &
        a_imag*sin(amplitude_phase)
    amp_rotated_imag = a_real*sin(amplitude_phase) + &
        a_imag*cos(amplitude_phase)
    amp_square_rotated = amp_rotated_real**2 + amp_rotated_imag**2
    sign_amplitude_residual = amp_square_sign-amp_square
    phase_amplitude_residual = amp_square_rotated-amp_square
    pair_field_residual = field_reversal-field_original
    conjugate_modulus_residual = amp_square_conj-amp_square

    ! Dimensionless refinement scales.  Each reference is supplied in the
    ! same physical unit as its observable and is required positive by the
    ! runtime caller; no mixed-unit absolute tolerance is encoded here.
    state_r = sym(arena, "state_R")
    state_z = sym(arena, "state_Z")
    state_phi = sym(arena, "state_phi")
    state_p_parallel = sym(arena, "state_p_parallel")
    state_mu = sym(arena, "state_mu")
    obs_tau = sym(arena, "tau_observable")
    obs_omega_b = sym(arena, "omega_b_observable")
    obs_omega_phi = sym(arena, "omega_phi_observable")
    hm_real = sym(arena, "Hm_real")
    hm_imag = sym(arena, "Hm_imag")
    root_derivative = sym(arena, "root_derivative")
    torque_density = sym(arena, "torque_density")
    torque = sym(arena, "torque")
    hm_squared = hm_real**2 + hm_imag**2
    r_ref = sym(arena, "R_ref")
    z_ref = sym(arena, "Z_ref")
    phi_ref = sym(arena, "phi_ref")
    p_parallel_ref = sym(arena, "p_parallel_ref")
    mu_ref = sym(arena, "mu_ref")
    tau_ref = sym(arena, "tau_ref")
    omega_b_ref = sym(arena, "omega_b_ref")
    omega_phi_ref = sym(arena, "omega_phi_ref")
    hm_real_ref = sym(arena, "Hm_real_ref")
    hm_imag_ref = sym(arena, "Hm_imag_ref")
    hm_squared_ref = sym(arena, "Hm_squared_ref")
    root_ref = sym(arena, "root_derivative_ref")
    torque_density_ref = sym(arena, "torque_density_ref")
    torque_ref = sym(arena, "torque_ref")
    scale_r = abs(state_r)/r_ref
    scale_z = abs(state_z)/z_ref
    scale_phi = abs(state_phi)/phi_ref
    scale_p_parallel = abs(state_p_parallel)/p_parallel_ref
    scale_mu = abs(state_mu)/mu_ref
    scale_tau = abs(obs_tau)/tau_ref
    scale_omega_b = abs(obs_omega_b)/omega_b_ref
    scale_omega_phi = abs(obs_omega_phi)/omega_phi_ref
    scale_hm_real = abs(hm_real)/hm_real_ref
    scale_hm_imag = abs(hm_imag)/hm_imag_ref
    scale_hm_squared = abs(hm_squared)/hm_squared_ref
    scale_root = abs(root_derivative)/root_ref
    scale_torque_density = abs(torque_density)/torque_density_ref
    scale_torque = abs(torque)/torque_ref

    ! Buchholz boundary forms.  Delta is a positive distance supplied after
    ! topology classification.  The residual outputs are generated class
    ! predicates: zero means the selected boundary class and coefficient agree.
    boundary_distance = sym(arena, "boundary_distance")
    boundary_reference_scale = sym(arena, "boundary_reference_scale")
    boundary_delta = boundary_distance/boundary_reference_scale
    regular_coefficient = sym(arena, "regular_coefficient")
    reflecting_coefficient = sym(arena, "reflecting_coefficient")
    separatrix_coefficient = sym(arena, "separatrix_coefficient")
    lambda_positive = sym(arena, "lambda_positive")
    regular_form = regular_coefficient
    reflecting_form = reflecting_coefficient*sqrt(boundary_delta)
    reflecting_derivative = diff(reflecting_form, boundary_distance)
    c_tau = 1/lambda_positive
    separatrix_coefficient = -c_tau
    separatrix_form = separatrix_coefficient*log(boundary_delta)
    separatrix_tau = c_tau*log(1/boundary_delta)
    xpoint_form = 2*separatrix_coefficient*log(boundary_delta)
    xpoint_tau = 2*c_tau*log(1/boundary_delta)
    regular_classification_residual = regular_form - regular_coefficient
    reflecting_classification_residual = reflecting_derivative - &
        reflecting_coefficient/(2*boundary_reference_scale*sqrt(boundary_delta))
    separatrix_classification_residual = separatrix_form - &
        separatrix_coefficient*log(boundary_delta)
    xpoint_pair_coefficient_residual = xpoint_form - 2*separatrix_form
    xpoint_pair_tau_residual = xpoint_tau - 2*separatrix_tau

    ! Buchholz cut and orientation in cylindrical components:
    ! C = (grad B cross grad psi) dot grad phi.
    grad_psi1 = sym(arena, "grad_psi1")
    grad_psi2 = sym(arena, "grad_psi2")
    grad_psi3 = sym(arena, "grad_psi3")
    grad_phi_coordinate1 = sym(arena, "grad_phi_coordinate1")
    grad_phi_coordinate2 = sym(arena, "grad_phi_coordinate2")
    grad_phi_coordinate3 = sym(arena, "grad_phi_coordinate3")
    cut_cross1 = gradb2*grad_psi3 - gradb3*grad_psi2
    cut_cross2 = gradb3*grad_psi1 - gradb1*grad_psi3
    cut_cross3 = gradb1*grad_psi2 - gradb2*grad_psi1
    cut_c = cut_cross1*grad_phi_coordinate1 + &
        cut_cross2*grad_phi_coordinate2 + cut_cross3*grad_phi_coordinate3
    d_cut_c_d_r = sym(arena, "d_cut_C_d_R")
    d_cut_c_d_arc_phi = sym(arena, "d_cut_C_d_arc_phi")
    d_cut_c_d_z = sym(arena, "d_cut_C_d_Z")
    dot_cut_r = sym(arena, "dot_cut_R")
    dot_cut_arc_phi = sym(arena, "dot_cut_arc_phi")
    dot_cut_z = sym(arena, "dot_cut_Z")
    cut_cdot = d_cut_c_d_r*dot_cut_r + d_cut_c_d_arc_phi*dot_cut_arc_phi + &
        d_cut_c_d_z*dot_cut_z
    cutdot_transversality = abs(cut_cdot)
    cutdot_orientation = cut_cdot

    dpsi_dx_section = sym(arena, "dpsi_star_dx")
    f_prime_section = sym(arena, "F_prime")
    cdot_section = sym(arena, "Cdot")
    signed_measure_section = sym(arena, "signed_R_Bparallel_star")
    period_section = sym(arena, "period")
    frequency_section = sym(arena, "frequency")
    dpsi_dx_reversed = -dpsi_dx_section
    f_prime_reversed = -f_prime_section
    cdot_reversed = -cdot_section
    root_weight_section = abs(dpsi_dx_section/f_prime_section)
    root_weight_reversed = abs(dpsi_dx_reversed/f_prime_reversed)
    positive_crossing_density = abs(signed_measure_section)*abs(cdot_section)
    positive_crossing_reversed = abs(signed_measure_section)*abs(cdot_reversed)
    period_reversed = period_section
    frequency_reversed = frequency_section

    crossing_signed_measure = sym(arena, "R_Bparallel_star")
    crossing_cdot = sym(arena, "Cdot")
    crossing_signed_density = crossing_signed_measure*crossing_cdot
    crossing_abs_measure = abs(crossing_signed_measure)
    crossing_abs_cdot = abs(crossing_cdot)
    crossing_positive_density = crossing_abs_measure*crossing_abs_cdot

    ! Exact sign/symmetry ledger.  Zero residuals are genuine relabellings;
    ! the explicitly nonzero differences are signs which retain physical
    ! content and must never be hidden by a post-hoc absolute value.
    section_root_residual = root_weight_reversed-root_weight_section
    crossing_orientation_residual = positive_crossing_reversed - &
        positive_crossing_density
    sign_omega_b = sym(arena, "sign_omega_b")
    sign_omega_phi = sym(arena, "sign_omega_phi")
    sign_torque_phi = sym(arena, "sign_torque_phi")
    sign_rotation_phi = sym(arena, "sign_rotation_phi")
    toroidal_resonance = m_mode*sign_omega_b+n_mode*sign_omega_phi
    toroidal_resonance_relabelled = m_mode*sign_omega_b + &
        (-n_mode)*(-sign_omega_phi)
    toroidal_resonance_residual = toroidal_resonance_relabelled - &
        toroidal_resonance
    toroidal_power = sign_torque_phi*sign_rotation_phi
    toroidal_power_relabelled = (-sign_torque_phi)*(-sign_rotation_phi)
    toroidal_power_residual = toroidal_power_relabelled-toroidal_power
    torque_component_difference = (-sign_torque_phi)-sign_torque_phi
    charge_reversal_energy_difference = (-charge)*electrostatic_potential - &
        charge*electrostatic_potential
    potential_reversal_energy_difference = charge*(-electrostatic_potential) - &
        charge*electrostatic_potential

    dpsi_star_dx_contribution = sym(arena, "dpsi_star_dx")
    abs_hm_squared = sym(arena, "abs_Hm_squared")
    tau_b = sym(arena, "tau_b")
    f_prime_contribution = n_mode*g_prime/tau_b
    frequency_contribution = abs(dpsi_star_dx_contribution)*abs_hm_squared* &
        tau_b/abs(f_prime_contribution)
    phase_contribution = abs(dpsi_star_dx_contribution)*abs_hm_squared* &
        abs(n_mode)*tau_b**2/abs(g_prime)
    n_squared_frequency_contribution = n_mode**2*frequency_contribution
    phase_contribution_identity = phase_contribution
    n_abs = sym(arena, "n_abs")
    tau_pos = sym(arena, "tau_pos")
    g_prime_abs = sym(arena, "g_prime_abs")
    f_prime_abs = n_abs*g_prime_abs/tau_pos
    dpsi_abs = sym(arena, "abs_dpsi_star_dx")
    hm_squared_positive = sym(arena, "abs_Hm_squared_positive")
    frequency_contribution_positive = dpsi_abs*hm_squared_positive*tau_pos/ &
        f_prime_abs
    phase_contribution_positive = dpsi_abs*hm_squared_positive* &
        n_abs*tau_pos**2/g_prime_abs
    ! Direct simple-root delta weight used by the class integral.  Absolute
    ! values belong only to the positive coordinate measure, complex harmonic
    ! modulus, and delta-function Jacobian.  The thermodynamic force remains
    ! signed and the transport layer, not this root kernel, owns n_mode**2.
    frequency_residual_derivative = &
        sym(arena, "frequency_residual_derivative")
    root_force_value = sym(arena, "force_value")
    root_hm_squared = hm_real**2 + hm_imag**2
    simple_root_weight = abs(dpsi_star_dx_contribution)*root_hm_squared* &
        tau_b/abs(frequency_residual_derivative)
    simple_root_force_contribution = simple_root_weight*root_force_value
    simple_root_phase_identity = simple_root_force_contribution - &
        simple_root_weight*root_force_value

    ! Cylindrical geometry kernel.  Both vector components and physical
    ! derivative directions are explicit: (R,phi,Z) and (R,arc_phi,Z).
    ! The arc_phi derivative is already physical (1/R*d/dphi); the only
    ! remaining cylindrical connection term is +bhat_phi/R in curl_Z.
    geom_radius = sym(arena, "radius")
    geom_arc_phi = sym(arena, "geometry_arc_phi")
    geom_z = sym(arena, "geometry_Z")
    geom_b_r = sym(arena, "b_R")
    geom_b_phi = sym(arena, "b_phi")
    geom_b_z = sym(arena, "b_Z")
    geom_dbr_dr = sym(arena, "d_b_R_d_R")
    geom_dbphi_dr = sym(arena, "d_b_phi_d_R")
    geom_dbz_dr = sym(arena, "d_b_Z_d_R")
    geom_dbr_darc_phi = sym(arena, "d_b_R_d_arc_phi")
    geom_dbphi_darc_phi = sym(arena, "d_b_phi_d_arc_phi")
    geom_dbz_darc_phi = sym(arena, "d_b_Z_d_arc_phi")
    geom_dbr_dz = sym(arena, "d_b_R_d_Z")
    geom_dbphi_dz = sym(arena, "d_b_phi_d_Z")
    geom_dbz_dz = sym(arena, "d_b_Z_d_Z")
    geom_b = [geom_b_r, geom_b_phi, geom_b_z]
    geom_d_b(:,1) = [geom_dbr_dr, geom_dbphi_dr, geom_dbz_dr]
    geom_d_b(:,2) = [geom_dbr_darc_phi, geom_dbphi_darc_phi, &
        geom_dbz_darc_phi]
    geom_d_b(:,3) = [geom_dbr_dz, geom_dbphi_dz, geom_dbz_dz]
    call build_cylindrical_geometry(geom_b, geom_d_b, geom_radius, geom_bmod, &
        geom_bhat, geom_grad_b, geom_dbhat, geom_curl)

    ! Manufactured analytic oracle: a unit-magnitude rotating field with
    ! alpha=R+arc_phi+Z.  Its derivatives are obtained by Fortsym diff,
    ! including the cylindrical connection term in curl_Z.
    oracle_coordinate = [geom_radius, geom_arc_phi, geom_z]
    oracle_b = [cos(geom_radius + geom_arc_phi + geom_z), &
        sin(geom_radius + geom_arc_phi + geom_z), zero]
    do k = 1, 3
        do i = 1, 3
            oracle_d_b(k,i) = diff(oracle_b(k), oracle_coordinate(i))
        end do
    end do
    call build_cylindrical_geometry(oracle_b, oracle_d_b, geom_radius, &
        oracle_bmod, oracle_bhat, oracle_grad_b, oracle_dbhat, oracle_curl)

    ! Complex bilinear interpolation on a physical R-Z cell.  The generated
    ! kernel owns both coordinate normalization and corner weights.
    interp_r = sym(arena, "coordinate_R")
    interp_z = sym(arena, "coordinate_Z")
    interp_r0 = sym(arena, "cell_R0")
    interp_r1 = sym(arena, "cell_R1")
    interp_z0 = sym(arena, "cell_Z0")
    interp_z1 = sym(arena, "cell_Z1")
    interp_amplitude_scale = sym(arena, "amplitude_scale")
    interp_u = (interp_r-interp_r0)/(interp_r1-interp_r0)
    interp_v = (interp_z-interp_z0)/(interp_z1-interp_z0)
    interp_u0_residual = (interp_r0-interp_r0)/(interp_r1-interp_r0)
    interp_u1_residual = (interp_r1-interp_r0)/(interp_r1-interp_r0)-1
    interp_v0_residual = (interp_z0-interp_z0)/(interp_z1-interp_z0)
    interp_v1_residual = (interp_z1-interp_z0)/(interp_z1-interp_z0)-1
    interp_w00 = (1-interp_u)*(1-interp_v)
    interp_w10 = interp_u*(1-interp_v)
    interp_w01 = (1-interp_u)*interp_v
    interp_w11 = interp_u*interp_v
    interp_r00 = sym(arena, "value_real_00")
    interp_r10 = sym(arena, "value_real_10")
    interp_r01 = sym(arena, "value_real_01")
    interp_r11 = sym(arena, "value_real_11")
    interp_i00 = sym(arena, "value_imag_00")
    interp_i10 = sym(arena, "value_imag_10")
    interp_i01 = sym(arena, "value_imag_01")
    interp_i11 = sym(arena, "value_imag_11")
    interp_value_real = interp_amplitude_scale*(interp_w00*interp_r00 + &
        interp_w10*interp_r10 + interp_w01*interp_r01 + interp_w11*interp_r11)
    interp_value_imag = interp_amplitude_scale*(interp_w00*interp_i00 + &
        interp_w10*interp_i10 + interp_w01*interp_i01 + interp_w11*interp_i11)
    interp_partition_residual = interp_w00 + interp_w10 + interp_w01 + &
        interp_w11 - 1
    interp_corner00 = (1-zero)*(1-zero)*interp_r00 + zero*(1-zero)* &
        interp_r10 + (1-zero)*zero*interp_r01 + zero*zero*interp_r11 - &
        interp_r00
    interp_corner10 = (1-one)*(1-zero)*interp_r00 + one*(1-zero)* &
        interp_r10 + (1-one)*zero*interp_r01 + one*zero*interp_r11 - &
        interp_r10
    interp_corner01 = (1-zero)*(1-one)*interp_r00 + zero*(1-one)* &
        interp_r10 + (1-zero)*one*interp_r01 + zero*one*interp_r11 - &
        interp_r01
    interp_corner11 = (1-one)*(1-one)*interp_r00 + one*(1-one)* &
        interp_r10 + (1-one)*one*interp_r01 + one*one*interp_r11 - &
        interp_r11

    ! Unique affine-in-s reconstruction from two monotone cell-center nodes.
    ! The runtime seam must reject s1==s0, non-monotone input, and invalid
    ! nonpositive references before calling this pure algebraic kernel.
    endpoint_s0 = sym(arena, "s0")
    endpoint_s1 = sym(arena, "s1")
    endpoint_f0 = sym(arena, "f0")
    endpoint_f1 = sym(arena, "f1")
    endpoint_slope = (endpoint_f1-endpoint_f0)/(endpoint_s1-endpoint_s0)
    endpoint_intercept = endpoint_f0-endpoint_slope*endpoint_s0
    endpoint_value_zero = endpoint_intercept
    endpoint_value_one = endpoint_intercept+endpoint_slope
    endpoint_derivative_zero = endpoint_slope
    endpoint_derivative_one = endpoint_slope
    endpoint_rho = sym(arena, "rho")
    endpoint_dfd_rho = 2*endpoint_rho*endpoint_derivative_zero
    endpoint_axis_residual = 2*zero*endpoint_derivative_zero
    affine_a = sym(arena, "affine_a")
    affine_b = sym(arena, "affine_b")
    affine_node0 = affine_a + affine_b*endpoint_s0
    affine_node1 = affine_a + affine_b*endpoint_s1
    affine_endpoint_slope = (affine_node1-affine_node0)/ &
        (endpoint_s1-endpoint_s0)
    affine_endpoint_zero = affine_node0-affine_endpoint_slope*endpoint_s0
    affine_endpoint_one = affine_endpoint_zero+affine_endpoint_slope

    ! Exact segment integral for an affine Omega_E(psi) profile.
    potential_psi0 = sym(arena, "psi0")
    potential_psi1 = sym(arena, "psi1")
    potential_omega0 = sym(arena, "Omega_E0")
    potential_omega1 = sym(arena, "Omega_E1")
    potential_c = sym(arena, "c_light")
    delta_phi_segment = (potential_psi1-potential_psi0)* &
        (potential_omega0+potential_omega1)/(2*potential_c)
    delta_phi_reversed = (potential_psi0-potential_psi1)* &
        (potential_omega1+potential_omega0)/(2*potential_c)
    omega_constant = sym(arena, "Omega_E_constant")
    delta_phi_constant = (potential_psi1-potential_psi0)*omega_constant/ &
        potential_c

    ! Exact potential and derivatives inside the same affine Omega_E segment.
    ! Phi is quadratic when Omega_E is affine; interpolating the integrated
    ! endpoint values linearly would replace the local electric field by its
    ! segment average and is therefore not the same physical profile.
    potential_psi = sym(arena, "psi")
    potential_phi0 = sym(arena, "Phi0")
    potential_delta_psi = potential_psi-potential_psi0
    potential_omega_slope = (potential_omega1-potential_omega0)/ &
        (potential_psi1-potential_psi0)
    potential_phi = potential_phi0 + &
        (potential_omega0*potential_delta_psi + &
        potential_omega_slope*potential_delta_psi**2/2)/potential_c
    potential_dphi_dpsi = diff(potential_phi, potential_psi)
    potential_d2phi_dpsi2 = diff(potential_dphi_dpsi, potential_psi)
    potential_phi1 = potential_phi0+delta_phi_segment
    potential_phi_reversed = potential_phi1 + &
        (potential_omega1*(potential_psi-potential_psi1) + &
        (potential_omega0-potential_omega1)/ &
        (potential_psi0-potential_psi1) * &
        (potential_psi-potential_psi1)**2/2)/potential_c
    potential_reversal_residual = potential_phi_reversed-potential_phi

    ! Explicit piecewise-affine map from native toroidal flux s_tor to the
    ! unscaled EQDSK normalized poloidal flux used by the cut atlas.  The
    ! factory stores field_scale*psi at profile nodes, hence field_scale is
    ! divided exactly once here and psi_hat is never identified with s_tor.
    flux_segment_s = sym(arena, "s_tor")
    flux_segment_s0 = sym(arena, "s_tor0")
    flux_segment_s1 = sym(arena, "s_tor1")
    flux_segment_psi0 = sym(arena, "scaled_psi0")
    flux_segment_psi1 = sym(arena, "scaled_psi1")
    flux_segment_field_scale = sym(arena, "field_scale")
    flux_segment_psi_sep = sym(arena, "psi_sep")
    flux_segment_target_psihat = sym(arena, "target_psihat")
    flux_segment_weight = (flux_segment_s-flux_segment_s0) &
        /(flux_segment_s1-flux_segment_s0)
    flux_segment_psihat = ((one-flux_segment_weight)*flux_segment_psi0 &
        +flux_segment_weight*flux_segment_psi1) &
        /(flux_segment_field_scale*flux_segment_psi_sep)
    flux_segment_dpsihat_dstor = (flux_segment_psi1-flux_segment_psi0) &
        /((flux_segment_s1-flux_segment_s0)*flux_segment_field_scale &
        *flux_segment_psi_sep)
    flux_segment_inverse_s = flux_segment_s0 &
        +(flux_segment_s1-flux_segment_s0) &
        *(flux_segment_target_psihat*flux_segment_field_scale &
        *flux_segment_psi_sep-flux_segment_psi0) &
        /(flux_segment_psi1-flux_segment_psi0)
    flux_segment_axis_residual = subs(flux_segment_psihat, flux_segment_s, &
        flux_segment_s0)-flux_segment_psi0 &
        /(flux_segment_field_scale*flux_segment_psi_sep)
    flux_segment_edge_residual = subs(flux_segment_psihat, flux_segment_s, &
        flux_segment_s1)-flux_segment_psi1 &
        /(flux_segment_field_scale*flux_segment_psi_sep)
    flux_segment_derivative_residual = diff(flux_segment_psihat, &
        flux_segment_s)-flux_segment_dpsihat_dstor
    flux_segment_forward_inverse_residual = subs(flux_segment_inverse_s, &
        flux_segment_target_psihat, flux_segment_psihat)-flux_segment_s
    flux_segment_inverse_forward_residual = subs(flux_segment_psihat, &
        flux_segment_s, flux_segment_inverse_s)-flux_segment_target_psihat
    scaled_flux_value = sym(arena, "scaled_psi")
    scaled_flux_psihat = scaled_flux_value &
        /(flux_segment_field_scale*flux_segment_psi_sep)
    scaled_flux_round_trip_residual = scaled_flux_psihat &
        *flux_segment_field_scale*flux_segment_psi_sep-scaled_flux_value
    physical_flux_normalized_value = sym(arena, "psi_physical")
    physical_flux_psihat = physical_flux_normalized_value &
        /(flux_segment_field_scale*flux_segment_psi_sep)
    physical_flux_normalization_residual = physical_flux_psihat &
        *flux_segment_field_scale*flux_segment_psi_sep &
        -physical_flux_normalized_value
    physical_flux_raw_psi = sym(arena, "psi")
    physical_flux_raw_psi_r = sym(arena, "psi_R")
    physical_flux_raw_psi_z = sym(arena, "psi_Z")
    physical_flux_dz_dr = sym(arena, "dZ_dR")
    physical_flux_raw_dpsi_dr = physical_flux_raw_psi_r + &
        physical_flux_raw_psi_z*physical_flux_dz_dr
    physical_flux_field_scale = sym(arena, "field_scale")
    physical_flux_psi = physical_flux_field_scale*physical_flux_raw_psi
    physical_flux_dpsi_dr = physical_flux_field_scale* &
        physical_flux_raw_dpsi_dr

    ! Physical B, allowed energy, and their derivatives on the regular R
    ! chart of the Eq. 13 cut.  Raw EQDSK psi derivatives and F enter here;
    ! field_scale is applied exactly once to B and physical psi.
    allowed_field_scale = sym(arena, "field_scale")
    allowed_dz_dr = sym(arena, "dZ_dR")
    allowed_d2z_dr2 = sym(arena, "d2Z_dR2")
    allowed_jk = sym(arena, "J_K")
    allowed_potential = sym(arena, "electrostatic_potential")
    allowed_dphi_dpsi = sym(arena, "dPhi_dpsi")
    allowed_d2phi_dpsi2 = sym(arena, "d2Phi_dpsi2")
    allowed_dpsi_dr = eqcut_psi_r+eqcut_psi_z*allowed_dz_dr
    allowed_equilibrium_psi_physical = allowed_field_scale*eqcut_psi0
    allowed_equilibrium_dpsi_physical_dr = &
        allowed_field_scale*allowed_dpsi_dr
    allowed_d2psi_dr2 = eqcut_psi_rr + &
        2*eqcut_psi_rz*allowed_dz_dr + &
        eqcut_psi_zz*allowed_dz_dr**2 + &
        eqcut_psi_z*allowed_d2z_dr2
    allowed_dpsi_r_dr = eqcut_psi_rr+eqcut_psi_rz*allowed_dz_dr
    allowed_dpsi_z_dr = eqcut_psi_rz+eqcut_psi_zz*allowed_dz_dr
    allowed_d2psi_r_dr2 = eqcut_psi_rrr + &
        2*eqcut_psi_rrz*allowed_dz_dr + &
        eqcut_psi_rzz*allowed_dz_dr**2 + &
        eqcut_psi_rz*allowed_d2z_dr2
    allowed_d2psi_z_dr2 = eqcut_psi_rrz + &
        2*eqcut_psi_rzz*allowed_dz_dr + &
        eqcut_psi_zzz*allowed_dz_dr**2 + &
        eqcut_psi_zz*allowed_d2z_dr2
    allowed_df_dr = eqcut_f_hat_first/eqcut_psi_sep*allowed_dpsi_dr
    allowed_d2f_dr2 = eqcut_f_hat_second/eqcut_psi_sep**2 * &
        allowed_dpsi_dr**2 + eqcut_f_hat_first/eqcut_psi_sep* &
        allowed_d2psi_dr2
    allowed_g = eqcut_psi_r**2+eqcut_psi_z**2+eqcut_f0**2
    allowed_sqrt_g = sqrt(allowed_g)
    allowed_dg_dr = 2*eqcut_psi_r*allowed_dpsi_r_dr + &
        2*eqcut_psi_z*allowed_dpsi_z_dr + &
        2*eqcut_f0*allowed_df_dr
    allowed_d2g_dr2 = 2*(allowed_dpsi_r_dr**2 + &
        eqcut_psi_r*allowed_d2psi_r_dr2 + &
        allowed_dpsi_z_dr**2 + &
        eqcut_psi_z*allowed_d2psi_z_dr2 + allowed_df_dr**2 + &
        eqcut_f0*allowed_d2f_dr2)
    allowed_dsqrt_g_dr = allowed_dg_dr/(2*allowed_sqrt_g)
    allowed_d2sqrt_g_dr2 = allowed_d2g_dr2/(2*allowed_sqrt_g) - &
        allowed_dg_dr**2/(4*allowed_sqrt_g**3)
    allowed_bmod = allowed_field_scale*allowed_sqrt_g/eqcut_radius
    allowed_dbmod_dr = allowed_field_scale* &
        (allowed_dsqrt_g_dr/eqcut_radius - &
        allowed_sqrt_g/eqcut_radius**2)
    allowed_d2bmod_dr2 = allowed_field_scale* &
        (allowed_d2sqrt_g_dr2/eqcut_radius - &
        2*allowed_dsqrt_g_dr/eqcut_radius**2 + &
        2*allowed_sqrt_g/eqcut_radius**3)
    allowed_omega_c = abs(charge)*allowed_bmod/(mass*c_light)
    allowed_domega_c_dr = abs(charge)*allowed_dbmod_dr/(mass*c_light)
    allowed_d2omega_c_dr2 = abs(charge)*allowed_d2bmod_dr2/ &
        (mass*c_light)
    allowed_energy_margin = h-allowed_jk*allowed_omega_c - &
        charge*allowed_potential
    allowed_denergy_dr = -allowed_jk*allowed_domega_c_dr - &
        charge*allowed_dphi_dpsi*allowed_field_scale*allowed_dpsi_dr
    allowed_d2energy_dr2 = -allowed_jk*allowed_d2omega_c_dr2 - charge* &
        (allowed_d2phi_dpsi2*(allowed_field_scale*allowed_dpsi_dr)**2 + &
        allowed_dphi_dpsi*allowed_field_scale*allowed_d2psi_dr2)
    allowed_vparallel_squared = 2*allowed_energy_margin/mass
    allowed_dvparallel_squared_dr = 2*allowed_denergy_dr/mass
    allowed_d2vparallel_squared_dr2 = 2*allowed_d2energy_dr2/mass
    allowed_bphi_covariant = eqcut_radius*eqcut_f0/allowed_sqrt_g
    allowed_dbphi_covariant_dr = &
        (eqcut_f0+eqcut_radius*allowed_df_dr)/allowed_sqrt_g - &
        eqcut_radius*eqcut_f0*allowed_dg_dr/(2*allowed_sqrt_g**3)

    ! Independent local path model used only by the exact proof suite.  Its
    ! first two derivatives at x=0 must reproduce the emitted R-chart jet.
    allowed_x = sym(arena, "allowed_path_x")
    allowed_radius_local = eqcut_radius+allowed_x
    allowed_psi_r_local = eqcut_psi_r + &
        allowed_dpsi_r_dr*allowed_x + &
        allowed_d2psi_r_dr2*allowed_x**2/2
    allowed_psi_z_local = eqcut_psi_z + &
        allowed_dpsi_z_dr*allowed_x + &
        allowed_d2psi_z_dr2*allowed_x**2/2
    allowed_f_local = eqcut_f0+allowed_df_dr*allowed_x + &
        allowed_d2f_dr2*allowed_x**2/2
    allowed_g_local = allowed_psi_r_local**2 + &
        allowed_psi_z_local**2+allowed_f_local**2
    allowed_bmod_local = allowed_field_scale*sqrt(allowed_g_local)/ &
        allowed_radius_local
    allowed_potential_local = allowed_potential + &
        allowed_dphi_dpsi*allowed_field_scale* &
        (allowed_dpsi_dr*allowed_x + &
        allowed_d2psi_dr2*allowed_x**2/2) + &
        allowed_d2phi_dpsi2*(allowed_field_scale* &
        allowed_dpsi_dr*allowed_x)**2/2
    allowed_energy_local = h-allowed_jk*abs(charge)*allowed_bmod_local/ &
        (mass*c_light)-charge*allowed_potential_local
    allowed_g_second_series_result = simplify_engine%series_coeff( &
        allowed_g_local, allowed_x, zero, 2)
    if (.not. allowed_g_second_series_result%ok) then
        error stop "fortsym could not derive Eq13 field-norm second jet"
    end if
    allowed_g_second_series = allowed_g_second_series_result%value
    allowed_potential_second_series_result = &
        simplify_engine%series_coeff(allowed_potential_local, allowed_x, &
        zero, 2)
    if (.not. allowed_potential_second_series_result%ok) then
        error stop "fortsym could not derive Eq13 potential second jet"
    end if
    allowed_potential_second_series = &
        allowed_potential_second_series_result%value

    ! The regular canonical chart is evaluated only where v_parallel^2>0.
    ! Its square-root singularity at a simple turning point is replaced by
    ! the generated y=sqrt(|x-x_t|) chart below.
    allowed_parallel_sign = sym(arena, "parallel_sign")
    allowed_psi_physical = sym(arena, "psi_physical")
    allowed_dpsi_physical_dr = sym(arena, "dpsi_physical_dR")
    allowed_vparallel_squared_input = sym(arena, "v_parallel_squared")
    allowed_dvparallel_squared_dr_input = &
        sym(arena, "dv_parallel_squared_dR")
    allowed_bphi_covariant_input = sym(arena, "bphi_covariant")
    allowed_dbphi_covariant_dr_input = &
        sym(arena, "dbphi_covariant_dR")
    allowed_vparallel = allowed_parallel_sign* &
        sqrt(allowed_vparallel_squared_input)
    allowed_dvparallel_dr = allowed_parallel_sign* &
        allowed_dvparallel_squared_dr_input/ &
        (2*sqrt(allowed_vparallel_squared_input))
    allowed_psi_star = allowed_psi_physical + mass*c_light/charge* &
        allowed_vparallel*allowed_bphi_covariant_input
    allowed_dpsi_star_dr = allowed_dpsi_physical_dr + &
        mass*c_light/charge*(allowed_dvparallel_dr* &
        allowed_bphi_covariant_input + &
        allowed_vparallel*allowed_dbphi_covariant_dr_input)

    turning_y = sym(arena, "root_coordinate_y")
    turning_direction = sym(arena, "allowed_side_sign")
    turning_dv2_dx = sym(arena, "dvparallel_squared_dx_root")
    turning_psi_root = sym(arena, "psi_star_flux_root")
    turning_dpsi_dx = sym(arena, "dpsi_flux_dx_root")
    turning_ratio_root = sym(arena, "bphi_covariant_root")
    turning_dratio_dx = sym(arena, "dbphi_covariant_dx_root")
    turning_delta_x = turning_direction*turning_y**2
    turning_vparallel = allowed_parallel_sign*turning_y* &
        sqrt(turning_direction*turning_dv2_dx)
    turning_psi_local = turning_psi_root + &
        turning_dpsi_dx*turning_delta_x
    turning_ratio_local = turning_ratio_root + &
        turning_dratio_dx*turning_delta_x
    turning_psi_star = turning_psi_local + mass*c_light/charge* &
        turning_vparallel*turning_ratio_local
    turning_dpsi_star_dy = diff(turning_psi_star, turning_y)
    turning_root_derivative = subs(turning_dpsi_star_dy, turning_y, zero)

    ! Generic canonical-momentum derivative numerator.  Let
    !
    !   F = v_parallel**2,  G = R*B_phi/B,  P = d psi/dR,
    !   kappa = mass*c_light*sigma/charge.
    !
    ! The local polynomials are only a symbolic jet device: Fortsym derives
    ! the value and the first two R derivatives of
    ! N = 2*sqrt(F)*P + kappa*(2*F*G' + F'*G).
    ! F and G require three jets and P requires two.  No finite-turning-point
    ! remainder is inferred here; the existing y chart above remains the
    ! only turning-point limit emitted by this generator.
    canonical_x = sym(arena, "canonical_R_offset")
    canonical_sigma = sym(arena, "sigma")
    canonical_kappa = mass*c_light*canonical_sigma/charge
    canonical_f0 = sym(arena, "F")
    canonical_f1 = sym(arena, "F_R")
    canonical_f2 = sym(arena, "F_RR")
    canonical_f3 = sym(arena, "F_RRR")
    canonical_g0 = sym(arena, "G")
    canonical_g1 = sym(arena, "G_R")
    canonical_g2 = sym(arena, "G_RR")
    canonical_g3 = sym(arena, "G_RRR")
    canonical_p0 = sym(arena, "P")
    canonical_p1 = sym(arena, "P_R")
    canonical_p2 = sym(arena, "P_RR")
    canonical_f_local = canonical_f0 + canonical_f1*canonical_x + &
        canonical_f2*canonical_x**2/2 + canonical_f3*canonical_x**3/6
    canonical_g_local = canonical_g0 + canonical_g1*canonical_x + &
        canonical_g2*canonical_x**2/2 + canonical_g3*canonical_x**3/6
    canonical_p_local = canonical_p0 + canonical_p1*canonical_x + &
        canonical_p2*canonical_x**2/2
    canonical_psi_local = psi + canonical_p0*canonical_x + &
        canonical_p1*canonical_x**2/2 + canonical_p2*canonical_x**3/6
    canonical_n_local = 2*sqrt(canonical_f_local)*canonical_p_local + &
        canonical_kappa*(2*canonical_f_local*exact_derivative( &
        canonical_g_local, canonical_x, "canonical G jet") + &
        exact_derivative(canonical_f_local, canonical_x, &
        "canonical F jet")*canonical_g_local)
    canonical_n = exact_simplify(subs(canonical_n_local, canonical_x, zero), &
        "canonical momentum numerator")
    canonical_n_r = exact_simplify(subs(exact_derivative(canonical_n_local, &
        canonical_x, "canonical numerator first jet"), canonical_x, zero), &
        "canonical momentum numerator first jet")
    canonical_n_rr = exact_simplify(subs(exact_derivative( &
        exact_derivative(canonical_n_local, canonical_x, &
        "canonical numerator first jet for second derivative"), &
        canonical_x, "canonical numerator second jet"), canonical_x, zero), &
        "canonical momentum numerator second jet")
    canonical_psi_star_local = canonical_psi_local + canonical_kappa* &
        sqrt(canonical_f_local)*canonical_g_local
    canonical_identity_residual = exact_simplify(canonical_n - &
        2*sqrt(canonical_f0)*subs(exact_derivative( &
        canonical_psi_star_local, canonical_x, &
        "canonical psi-star derivative"), canonical_x, zero), &
        "canonical psi-star derivative identity")

    ! Coordinate reversal changes only derivative orientation.  The scalar
    ! canonical flux is unchanged, so N, N_R, N_RR have alternating parity.
    canonical_n_coordinate = exact_simplify(-subs(canonical_n_local, &
        canonical_x, -canonical_x), "coordinate-reversed canonical N")
    canonical_n_r_coordinate = exact_simplify(subs( &
        exact_derivative(canonical_n_coordinate, canonical_x, &
        "coordinate-reversed canonical N first jet"), canonical_x, zero), &
        "coordinate-reversed canonical N first jet")
    canonical_n_rr_coordinate = exact_simplify(subs( &
        exact_derivative(exact_derivative(canonical_n_coordinate, &
        canonical_x, "coordinate-reversed canonical N first jet for second"), &
        canonical_x, "coordinate-reversed canonical N second jet"), &
        canonical_x, zero), "coordinate-reversed canonical N second jet")
    canonical_coordinate_residual = subs(canonical_n_coordinate, &
        canonical_x, zero)+canonical_n
    canonical_coordinate_r_residual = canonical_n_r_coordinate-canonical_n_r
    canonical_coordinate_rr_residual = &
        canonical_n_rr_coordinate+canonical_n_rr
    canonical_n_coordinate = subs(canonical_n_coordinate, canonical_x, zero)

    ! Simultaneous reversal of the orbit-time parameter and both frequency
    ! coordinates reverses the resonance scalar while preserving its zero
    ! set.  This is an algebraic reparametrization, not a physical
    ! time-reversal symmetry of fixed-B guiding-centre dynamics.
    canonical_time_phase = (m_mode*harmonic_omega_b + &
        n_mode*harmonic_omega_phi)*harmonic_time
    canonical_time_phase_reversed = subs(subs(subs( &
        canonical_time_phase, harmonic_omega_b, -harmonic_omega_b), &
        harmonic_omega_phi, -harmonic_omega_phi), harmonic_time, -harmonic_time)
    canonical_time_resonance = m_mode*harmonic_omega_b + &
        n_mode*harmonic_omega_phi
    canonical_time_resonance_reversed = subs(subs( &
        canonical_time_resonance, harmonic_omega_b, -harmonic_omega_b), &
        harmonic_omega_phi, -harmonic_omega_phi)
    canonical_time_phase_residual = cos(canonical_time_phase_reversed)- &
        cos(canonical_time_phase)
    canonical_time_resonance_residual = &
        canonical_time_resonance_reversed+canonical_time_resonance

    ! The Fourier-pair relabelling is a distinct algebraic identity: it reverses
    ! (m,n) and conjugates the complex amplitude at fixed orbit frequencies.
    ! The scalar label changes sign, while the real field and zero set remain
    ! the same.  Combining this relabelling with the orbit-time/frequency-
    ! coordinate reparametrization restores the original signed label.
    canonical_fourier_pair_phase = subs(subs( &
        canonical_time_phase, m_mode, -m_mode), n_mode, -n_mode)
    canonical_fourier_pair_resonance = subs(subs( &
        canonical_time_resonance, m_mode, -m_mode), n_mode, -n_mode)
    canonical_combined_resonance = subs(subs( &
        canonical_fourier_pair_resonance, harmonic_omega_b, &
        -harmonic_omega_b), harmonic_omega_phi, -harmonic_omega_phi)
    canonical_fourier_pair_phase_residual = &
        cos(canonical_fourier_pair_phase)-cos(canonical_time_phase)
    canonical_fourier_pair_resonance_residual = &
        canonical_fourier_pair_resonance+canonical_time_resonance
    canonical_combined_resonance_residual = &
        canonical_combined_resonance-canonical_time_resonance

    ! Toroidal-angle reversal relabels every toroidal covariant together:
    ! phi, n, omega_phi, psi_cov, and G.  The covariant one-form is the same
    ! physical object after relabelling; its coordinate component P_phi changes
    ! sign.
    canonical_toroidal_phase = n_mode*phi - &
        (m_mode*harmonic_omega_b+n_mode*harmonic_omega_phi)*harmonic_time
    canonical_toroidal_phase_relabelled = subs(subs(subs( &
        canonical_toroidal_phase, n_mode, -n_mode), phi, -phi), &
        harmonic_omega_phi, -harmonic_omega_phi)
    canonical_toroidal_phase_residual = &
        canonical_toroidal_phase_relabelled-canonical_toroidal_phase
    canonical_p_phi_covariant = charge/c_light*psi + &
        p_parallel*canonical_g0
    canonical_p_phi_covariant_relabelled = subs(subs( &
        canonical_p_phi_covariant, psi, -psi), canonical_g0, -canonical_g0)
    canonical_toroidal_p_phi_residual = &
        canonical_p_phi_covariant_relabelled+canonical_p_phi_covariant

    ! Coordinate reversal preserves the total variation, even though every
    ! signed endpoint difference reverses orientation.
    canonical_endpoint_delta_01 = sym(arena, "endpoint_delta_01")
    canonical_endpoint_delta_12 = sym(arena, "endpoint_delta_12")
    canonical_endpoint_variation = abs(canonical_endpoint_delta_01) + &
        abs(canonical_endpoint_delta_12)
    canonical_endpoint_variation_reversed = &
        abs(-canonical_endpoint_delta_01) + abs(-canonical_endpoint_delta_12)
    canonical_endpoint_variation_residual = &
        canonical_endpoint_variation_reversed-canonical_endpoint_variation

    ! Reversing sigma alone is not a symmetry: it reverses only the magnetic
    ! canonical offset while leaving the flux-gradient term unchanged.
    canonical_sigma_difference = subs(canonical_n, canonical_sigma, &
        -canonical_sigma)-canonical_n
    canonical_sigma_identity_residual = canonical_sigma_difference + &
        2*canonical_kappa*(2*canonical_f0*canonical_g1 + &
        canonical_f1*canonical_g0)

    ! The existing simple-turning chart supplies F=a*y**2, P=P_t, and
    ! G=G_t+G'_t*delta_R with delta_R=side*y**2.  It therefore proves the
    ! finite value N_t and the y derivative at the root, but not a finite
    ! R-derivative or a remainder.  Only those two exact limits are emitted.
    canonical_turning_f_local = turning_direction*turning_dv2_dx* &
        turning_y**2
    canonical_turning_g_local = turning_ratio_root + &
        turning_dratio_dx*turning_direction*turning_y**2
    canonical_turning_n_local = 2*turning_y*sqrt(turning_direction* &
        turning_dv2_dx)*turning_dpsi_dx + &
        canonical_kappa*(2*canonical_turning_f_local*turning_dratio_dx + &
        turning_dv2_dx*canonical_turning_g_local)
    canonical_turning_n_root = exact_simplify(subs( &
        canonical_turning_n_local, turning_y, zero), &
        "simple-turn canonical N limit")
    canonical_turning_n_y_root = exact_simplify(subs(exact_derivative( &
        canonical_turning_n_local, turning_y, &
        "simple-turn canonical N y derivative"), turning_y, zero), &
        "simple-turn canonical N y-derivative limit")
    canonical_turning_n_root_residual = canonical_turning_n_root - &
        canonical_kappa*turning_dv2_dx*turning_ratio_root
    canonical_turning_n_y_root_residual = canonical_turning_n_y_root - &
        2*sqrt(turning_direction*turning_dv2_dx)* &
        turning_dpsi_dx

    ! Axisymmetric phase-space one-form and Noether construction.  The
    ! explicit convention is A_phi_cov=psi and b_phi_cov=R*b_phi.
    b_phi_cov = sym(arena, "b_phi_cov")
    phi_coordinate = sym(arena, "phi_coordinate")
    phi_dot_symbol = sym(arena, "phi_dot")
    canonical_p_phi = charge/c_light*psi + p_parallel*b_phi_cov
    canonical_p_phi_cylindrical = charge/c_light*psi + &
        p_parallel*radius*bhat2
    phase_space_lagrangian = canonical_p_phi*phi_dot_symbol - h
    d_lagrangian_d_phi = diff(phase_space_lagrangian, phi_coordinate)
    d_lagrangian_d_phi_dot = diff(phase_space_lagrangian, phi_dot_symbol)
    canonical_residual = d_lagrangian_d_phi_dot - canonical_p_phi

    ! Manufactured axisymmetric flux convention in cylindrical coordinates.
    psi_r = sym(arena, "psi_r")
    psi_z = sym(arena, "psi_z")
    b_r_from_psi = -psi_z/radius
    b_z_from_psi = psi_r/radius
    flux_b_r_residual = b_r_from_psi - (-psi_z/radius)
    flux_b_z_residual = b_z_from_psi - psi_r/radius

    ! ------------------------------------------------------------------
    ! Strict proof suite.  The root identity is expressed in positive
    ! magnitudes because |n| and |g'| are not algebraically reducible without
    ! those domain facts.  This records the exact side conditions instead of
    ! asking an engine to guess them.
    n_abs = sym(arena, "n_abs")
    tau_pos = sym(arena, "tau_pos")
    g_prime_abs = sym(arena, "g_prime_abs")
    f_prime_abs = n_abs*g_prime_abs/tau_pos
    root_weight_abs = n_abs**2*tau_pos/f_prime_abs
    phase_weight_abs = n_abs*tau_pos**2/g_prime_abs
    h_axisymmetric = p_parallel**2/(2*mass) + mu*bmod + q_phi_energy

    call suite_begin(proofs, "NEO-RT direct full-FOW symbolic derivations")
    call check_identity(proofs, proof_engine, "J_K omega_c = mu B", &
        jk_omega_c - mu*bmod)
    call check_identity(proofs, proof_engine, &
        "J_K candidate times omega_c recovers H-q Phi", &
        jk_candidate*omega_c - (h - q_phi_energy))
    call check_identity(proofs, proof_engine, "resonance root side condition g(x_r)=0", &
        g_at_root)
    call check_identity(proofs, proof_engine, &
        "resonance root derivative from local expansion", &
        f_prime_from_series - n_mode*g_prime_x/tau)
    call check_identity(proofs, proof_engine, "Eq4 coefficient definition", &
        eq4_coeff - (2*(h - q_phi_energy) - jk*omega_c))
    call check_identity(proofs, proof_engine, &
        "canonical harmonic phase closes over one bounce", &
        harmonic_complete_phase_residual)
    call check_identity(proofs, proof_engine, &
        "canonical harmonic phase preserves amplitude modulus", &
        harmonic_modulus_residual)
    call check_identity(proofs, proof_engine, &
        "canonical harmonic Eq4 coefficient has the on-shell limit", &
        harmonic_shell_identity)
    call check_identity(proofs, proof_engine, "psi_star definition", &
        psi_star - c_light/charge*p_phi)
    call check_identity(proofs, proof_engine, &
        "canonical N is 2 sqrt(F) times d psi_star/dR", &
        canonical_identity_residual)
    call check_identity(proofs, proof_engine, &
        "simple-turn canonical N has the exact finite limit", &
        canonical_turning_n_root_residual)
    call check_identity(proofs, proof_engine, &
        "simple-turn canonical N y derivative has the exact limit", &
        canonical_turning_n_y_root_residual)
    call check_identity(proofs, proof_engine, &
        "coordinate reversal changes canonical N orientation only", &
        canonical_coordinate_residual)
    call check_identity(proofs, proof_engine, &
        "coordinate reversal preserves canonical N first derivative", &
        canonical_coordinate_r_residual)
    call check_identity(proofs, proof_engine, &
        "coordinate reversal changes canonical N second derivative orientation", &
        canonical_coordinate_rr_residual)
    call check_identity(proofs, proof_engine, &
        "orbit-time/frequency-coordinate reparametrization preserves phase", &
        canonical_time_phase_residual)
    call check_identity(proofs, proof_engine, &
        "orbit-time/frequency-coordinate reparametrization reverses resonance", &
        canonical_time_resonance_residual)
    call check_identity(proofs, proof_engine, &
        "Fourier-pair relabelling preserves the real phase", &
        canonical_fourier_pair_phase_residual)
    call check_identity(proofs, proof_engine, &
        "Fourier-pair relabelling reverses the resonance scalar", &
        canonical_fourier_pair_resonance_residual)
    call check_identity(proofs, proof_engine, &
        "orbit-time reparametrization plus Fourier-pair relabelling restores label", &
        canonical_combined_resonance_residual)
    call check_identity(proofs, proof_engine, &
        "toroidal reversal relabels the harmonic phase", &
        canonical_toroidal_phase_residual)
    call check_identity(proofs, proof_engine, &
        "toroidal reversal flips the covariant P_phi component", &
        canonical_toroidal_p_phi_residual)
    call check_identity(proofs, proof_engine, &
        "sigma reversal has the exact signed canonical difference", &
        canonical_sigma_identity_residual)
    call check_identity(proofs, proof_engine, &
        "coordinate reversal preserves endpoint total variation", &
        canonical_endpoint_variation_residual)
    call check_identity(proofs, proof_engine, &
        "toroidal reversal preserves torque power", &
        toroidal_power_residual)
    call check_identity(proofs, proof_engine, &
        "toroidal reversal flips the torque component", &
        torque_component_difference+2*sign_torque_phi)
    call check_identity(proofs, proof_engine, "d psi_star / d p_phi", &
        dpsi_dp_phi - c_light/charge)
    call check_identity(proofs, proof_engine, &
        "frequency-root and phase-root weights", &
        root_weight_abs - phase_weight_abs)
    call check_identity(proofs, proof_engine, &
        "frequency-root and phase-root contribution weights", &
        n_abs**2*frequency_contribution_positive - &
        phase_contribution_positive)
    call check_identity(proofs, proof_engine, &
        "simple-root harmonic modulus from complex components", &
        root_hm_squared-(hm_real**2+hm_imag**2))
    call check_identity(proofs, proof_engine, &
        "simple-root contribution preserves signed force", &
        simple_root_phase_identity)
    call check_identity(proofs, proof_engine, "Eq17 force bracket", &
        force_bracket - (a1 + (h - q_phi_energy)/temperature*a2))
    call check_identity(proofs, proof_engine, &
        "Eq17 outer factor owns exactly one Phi_eff", &
        eq17_outer_factor-phi_eff*eq17_outer_unit_phi)
    call check_identity(proofs, proof_engine, &
        "Eq17 outer factor gauge invariance", eq17_gauge_residual)
    call check_identity(proofs, proof_engine, &
        "Phi_eff normalization", &
        norm_phi_eff-c_light*e_ref/(abs(charge)*norm_velocity))
    call check_identity(proofs, proof_engine, &
        "J_K normalization scale", &
        norm_jk_scale-mass*c_light*e_ref/abs(charge))
    call check_identity(proofs, proof_engine, &
        "normalized psi_star round trip", norm_psi_hat*norm_phi_eff - &
        norm_psi_star)
    call check_identity(proofs, proof_engine, &
        "normalized bounce-time/frequency product", &
        norm_tau_hat*norm_omega_b_hat-norm_tau_b*norm_omega_b)
    call check_identity(proofs, proof_engine, &
        "semi-infinite energy-map Jacobian", &
        quad_dh_hat_dt-diff(quad_qphi_min/e_ref + &
        quad_t/(1-quad_t), quad_t))
    call check_identity(proofs, proof_engine, &
        "paired action map round trip", quad_jk_phys - &
        quad_jk_max*quad_j_unit)
    call check_identity(proofs, proof_engine, &
        "paired quadrature weight owns both normalized Jacobians", &
        quad_weight-quad_h_weight*quad_dh_hat_dt* &
        quad_j_weight*quad_jk_hat_max)
    call check_identity(proofs, proof_engine, &
        "Gauss interval affine coordinate", &
        2*gauss_mapped_node-gauss_lower_bound-gauss_upper_bound- &
        (gauss_upper_bound-gauss_lower_bound)*gauss_standard_node)
    call check_identity(proofs, proof_engine, &
        "Gauss interval measure Jacobian", &
        2*gauss_mapped_weight- &
        (gauss_upper_bound-gauss_lower_bound)*gauss_standard_weight)
    call check_identity(proofs, proof_engine, &
        "polynomial enclosure midpoint is signed c0", &
        (poly_lower+poly_upper)/2-poly_c0)
    call check_identity(proofs, proof_engine, &
        "polynomial enclosure width is twice the tail bound", &
        poly_upper-poly_lower-2*poly_tail_bound)
    call check_identity(proofs, proof_engine, &
        "physical magnetic-moment velocity round trip", &
        (2*mu_from_vperp*bmod/mass)-conversion_vperp**2)
    call check_identity(proofs, proof_engine, &
        "Buchholz action magnetic-moment round trip", &
        abs(charge)*jk_from_mu/(mass*c_light)-conversion_mu)
    call check_identity(proofs, proof_engine, &
        "Buchholz action owns perpendicular energy", &
        jk_from_mu*conversion_omega_c-conversion_mu*bmod)
    call check_identity(proofs, proof_engine, &
        "POTATO action magnetic-moment equivalence", &
        jpotato_from_mu-vperp_squared_from_mu/(conversion_v0**2*bmod))
    call check_identity(proofs, proof_engine, &
        "POTATO action magnetic-moment round trip", &
        mass*conversion_v0**2*jpotato_from_mu/2-conversion_mu)
    call check_identity(proofs, proof_engine, &
        "axisymmetric invariant Hamiltonian agrees with Littlejohn kernel", &
        invariant_h-(p_parallel**2/(2*mass)+mu*bmod+q_phi_energy))
    call check_identity(proofs, proof_engine, &
        "axisymmetric invariant canonical momentum agrees with Noether kernel", &
        invariant_pphi-canonical_p_phi_cylindrical)
    call check_identity(proofs, proof_engine, &
        "axisymmetric invariant canonical flux definition", &
        invariant_psistar-c_light/charge*invariant_pphi)
    call check_identity(proofs, proof_engine, &
        "fixed-H launch compares parallel momentum squared", &
        invariant_ppar_squared-2*mass*(h-mu*bmod-q_phi_energy))
    call check_identity(proofs, proof_engine, "force dot b cross force", &
        force1*cross1 + force2*cross2 + force3*cross3)
    call check_identity(proofs, proof_engine, "B-star cancellation in dH/dt", &
        p_parallel/mass*dot_p + &
        p_parallel/mass*(bstar1*force1 + bstar2*force2 + bstar3*force3) &
        /bparallel_star)
    call check_identity(proofs, proof_engine, "dH/dt = 0", hamiltonian_dot)
    call check_identity(proofs, proof_engine, "mu_dot = 0", dot_mu)
    call check_identity(proofs, proof_engine, "phi_dot = v_phi / R", &
        dot_phi - v2/radius)
    call check_identity(proofs, proof_engine, &
        "on-shell Hamiltonian perturbation ratio", &
        perturbation_ratio - boozer_ratio)
    call check_identity(proofs, proof_engine, "field sign reversal is odd", &
        field_sign_sum)
    call check_identity(proofs, proof_engine, &
        "projective phase parameterization has unit norm", &
        phase_unit_residual)
    call check_identity(proofs, proof_engine, &
        "phase shift is a coordinate relabelling", &
        phase_relabel_residual)
    call check_identity(proofs, proof_engine, &
        "conjugation plus (m,n) reversal preserves real field", &
        field_reversal - field_original)
    call check_identity(proofs, proof_engine, &
        "amplitude square ignores sign", amp_square_sign - amp_square)
    call check_identity(proofs, proof_engine, &
        "pointwise amplitude square ignores conjugation", &
        amp_square_conj - amp_square)
    call check_identity(proofs, proof_engine, &
        "global complex amplitude phase leaves squared amplitude", &
        amp_square_rotated-amp_square)
    call check_identity(proofs, proof_engine, &
        "fixed-mode conjugation changes the real field in general", &
        fixed_conjugation_difference-2*a_imag*sin(chi))
    call check_identity(proofs, proof_engine, &
        "gauge shift leaves the thermodynamic force bracket", &
        gauge_force_residual)
    call check_identity(proofs, proof_engine, &
        "toroidal coordinate reversal relabels the resonance", &
        toroidal_resonance_residual)
    call check_identity(proofs, proof_engine, &
        "toroidal coordinate reversal preserves torque power", &
        toroidal_power_residual)
    call check_identity(proofs, proof_engine, &
        "toroidal torque component changes sign under coordinate reversal", &
        torque_component_difference+2*sign_torque_phi)
    call check_identity(proofs, proof_engine, &
        "charge reversal alone changes potential energy", &
        charge_reversal_energy_difference+ &
        2*charge*electrostatic_potential)
    call check_identity(proofs, proof_engine, &
        "potential sign reversal is not a gauge shift", &
        potential_reversal_energy_difference+ &
        2*charge*electrostatic_potential)
    call check_identity(proofs, proof_engine, &
        "axisymmetric Hamiltonian has no phi derivative", &
        diff(h_axisymmetric, phi))
    call check_identity(proofs, proof_engine, &
        "axisymmetric Noether d p_phi / dt = dL/dphi", d_lagrangian_d_phi)
    call check_identity(proofs, proof_engine, &
        "phase-space one-form p_phi", canonical_residual)
    call check_identity(proofs, proof_engine, &
        "canonical p_phi is dL/dphi_dot", &
        d_lagrangian_d_phi_dot - canonical_p_phi)
    call check_identity(proofs, proof_engine, &
        "manufactured B_R flux convention", flux_b_r_residual)
    call check_identity(proofs, proof_engine, &
        "manufactured B_Z flux convention", flux_b_z_residual)
    call check_identity(proofs, proof_engine, &
        "cylindrical canonical p_phi uses unit-vector b_phi", &
        canonical_p_phi_cylindrical - (charge/c_light*psi + &
        p_parallel*radius*bhat2))
    call check_identity(proofs, proof_engine, &
        "axisymmetric P_phi conserved along Littlejohn dynamics", &
        axis_pphi_dot)
    call check_identity(proofs, proof_engine, &
        "axisymmetric P_phi residual is P_phi_dot minus zero", &
        axis_pphi_residual-axis_pphi_dot)
    call check_identity(proofs, proof_engine, &
        "axisymmetric grad Phi follows Phi_psi grad psi", &
        axis_force_r - (mu*axis_grad_b_r + charge*axis_phi_psi*radius*axis_b_z) + &
        axis_force_z - (mu*axis_grad_b_z - charge*axis_phi_psi*radius*axis_b_r))
    call check_identity(proofs, proof_engine, "cylindrical cross R component", &
        cut_cross1 - (gradb2*grad_psi3 - gradb3*grad_psi2))
    call check_identity(proofs, proof_engine, "cylindrical cross phi component", &
        cut_cross2 - (gradb3*grad_psi1 - gradb1*grad_psi3))
    call check_identity(proofs, proof_engine, "cylindrical cross Z component", &
        cut_cross3 - (gradb1*grad_psi2 - gradb2*grad_psi1))
    call check_identity(proofs, proof_engine, "Buchholz cut contraction", &
        cut_c - (cut_cross1*grad_phi_coordinate1 + &
        cut_cross2*grad_phi_coordinate2 + cut_cross3*grad_phi_coordinate3))
    call check_identity(proofs, proof_engine, "regular boundary classification", &
        regular_classification_residual)
    call check_identity(proofs, proof_engine, &
        "reflecting square-root boundary classification", &
        reflecting_classification_residual)
    call check_identity(proofs, proof_engine, &
        "separatrix logarithmic boundary classification", &
        separatrix_classification_residual)
    call check_identity(proofs, proof_engine, &
        "paired homoclinic X-point coefficient", &
        xpoint_pair_coefficient_residual)
    call check_identity(proofs, proof_engine, &
        "paired homoclinic X-point form", xpoint_pair_tau_residual)
    call check_identity(proofs, proof_engine, "cut directional derivative", &
        cut_cdot - (d_cut_c_d_r*dot_cut_r + d_cut_c_d_arc_phi*dot_cut_arc_phi + &
        d_cut_c_d_z*dot_cut_z))
    call check_identity(proofs, proof_engine, &
        "EQDSK cut jet directional derivative", &
        eqcut_cdot-(eqcut_dc_dr*eqcut_dot_r + &
        eqcut_dc_darc_phi*eqcut_dot_arc_phi + eqcut_dc_dz*eqcut_dot_z))
    call check_identity(proofs, proof_engine, &
        "EQDSK cut matches compact axisymmetric identity", &
        eqcut_compact_residual)
    call check_identity(proofs, proof_engine, &
        "emitted EQDSK cut matches compact axisymmetric identity", &
        eqcut_emitted_residual)
    call check_identity(proofs, proof_engine, &
        "EQDSK numerator matches compact N", eqcut_n_residual)
    call check_identity(proofs, proof_engine, &
        "EQDSK C numerator identity", eqcut_c_numerator_residual)
    call check_identity(proofs, proof_engine, &
        "EQDSK numerator value is independent of first F derivative", &
        eqcut_f1_n_residual)
    call check_identity(proofs, proof_engine, &
        "EQDSK numerator R derivative does not require second F derivative", &
        eqcut_f2_n_r_residual)
    call check_identity(proofs, proof_engine, &
        "EQDSK numerator Z derivative does not require second F derivative", &
        eqcut_f2_n_z_residual)
    call check_identity(proofs, proof_engine, &
        "simultaneous F jet reversal preserves EQDSK numerator", &
        eqcut_f_reversal_n_residual)
    call check_identity(proofs, proof_engine, &
        "simultaneous F jet reversal preserves EQDSK numerator RR", &
        eqcut_f_reversal_n_rr_residual)
    call check_identity(proofs, proof_engine, &
        "simultaneous F jet reversal preserves EQDSK numerator RZ", &
        eqcut_f_reversal_n_rz_residual)
    call check_identity(proofs, proof_engine, &
        "simultaneous F jet reversal preserves EQDSK numerator ZZ", &
        eqcut_f_reversal_n_zz_residual)
    call check_identity(proofs, proof_engine, &
        "EQDSK cut value is independent of first F derivative", &
        eqcut_f1_c_residual)
    call check_identity(proofs, proof_engine, &
        "EQDSK cut orientation reversal flips C", &
        eqcut_reversed+eqcut_c)
    call check_identity(proofs, proof_engine, &
        "EQDSK cut orientation reversal flips Cdot", &
        eqcut_cdot_reversed+eqcut_cdot)
    call check_identity(proofs, proof_engine, &
        "EQDSK cut positive transversality ignores orientation", &
        abs(eqcut_cdot_reversed)-eqcut_abs_cdot)
    call check_identity(proofs, proof_engine, &
        "EQDSK cut has exact field-scale squared normalization", &
        eqcut_scaled-eqcut_scale_multiplier**2*eqcut_c)
    call check_identity(proofs, proof_engine, &
        "EQDSK cut does not require second F derivative", &
        eqcut_f2_c_residual)
    call check_identity(proofs, proof_engine, &
        "EQDSK cut rate does not require second F derivative", &
        eqcut_f2_cdot_residual)
    call check_identity(proofs, proof_engine, &
        "up-down symmetric midplane is an Eq13 cut", eqcut_midplane)
    call check_identity(proofs, proof_engine, &
        "R-chart tangent is orthogonal to grad N", &
        eqcut_r_chart_tangent_residual)
    call check_identity(proofs, proof_engine, &
        "Z-chart tangent is orthogonal to grad N", &
        eqcut_z_chart_tangent_residual)
    call check_identity(proofs, proof_engine, &
        "R-chart normalized-flux derivative is the cut chain rule", &
        eqcut_r_chart_flux_residual)
    call check_identity(proofs, proof_engine, &
        "R-chart second derivative satisfies implicit cut equation", &
        eqcut_r_chart_second_tangent_residual)
    call check_identity(proofs, proof_engine, &
        "R-chart normalized-flux curvature is the second chain rule", &
        eqcut_r_chart_flux_second_residual)
    call check_identity(proofs, proof_engine, &
        "profile F second derivative contribution to cut Hessian", &
        eqcut_f2_hessian_contraction_residual)
    call check_identity(proofs, proof_engine, &
        "overall cut-numerator sign leaves R-chart slope unchanged", &
        eqcut_r_chart_slope_reversed-eqcut_r_chart_slope)
    call check_identity(proofs, proof_engine, &
        "overall cut-numerator sign leaves R-chart curvature unchanged", &
        eqcut_r_chart_second_reversed-eqcut_r_chart_second)
    call check_identity(proofs, proof_engine, &
        "EQDSK mean-value enclosure reproduces its center", &
        eqcut_mean_center_residual)
    call check_identity(proofs, proof_engine, &
        "EQDSK axis derivative has square-root limiting form", &
        eqcut_axis_derivative_limit_residual)
    call check_identity(proofs, proof_engine, &
        "EQDSK axis inverse has square-root limiting form", &
        eqcut_axis_inverse_limit_residual)
    call check_identity(proofs, proof_engine, &
        "rho_tor squared is s_tor", eqcut_s_tor-eqcut_rho**2)
    call check_identity(proofs, proof_engine, &
        "rho_tor map derivative is exact", eqcut_rho_map_derivative_residual)
    call check_identity(proofs, proof_engine, &
        "nonnegative s_tor square-root map is the inverse branch", &
        eqcut_s_tor_inverse_residual)
    call check_identity(proofs, proof_engine, &
        "rho_tor to normalized-flux derivative is the chain rule", &
        eqcut_flux_coordinate_chain_residual)
    call check_identity(proofs, proof_engine, &
        "interval cut-flux coordinate uses the exact chain rule", &
        eqcut_flux_coordinate_chain_residual)
    call check_identity(proofs, proof_engine, &
        "inverse cut-flux chart preserves normalized-flux derivative", &
        eqcut_flux_inverse_chain_residual)
    call check_identity(proofs, proof_engine, &
        "inverse cut-flux chart preserves second derivative", &
        eqcut_flux_inverse_second_chain_residual)
    call check_identity(proofs, proof_engine, &
        "allowed energy rho derivative is the inverse-chart chain", &
        eqcut_allowed_first_chain_residual)
    call check_identity(proofs, proof_engine, &
        "allowed energy rho curvature is the inverse-chart chain", &
        eqcut_allowed_second_chain_residual)
    call check_identity(proofs, proof_engine, &
        "canonical rho derivative is the inverse-chart chain", &
        eqcut_allowed_psistar_chain_residual)
    call check_identity(proofs, proof_engine, &
        "axis rho_tor derivative has finite square-root limit", &
        eqcut_axis_rho_limit_residual)
    call check_identity(proofs, proof_engine, &
        "axis allowed-energy rho derivative uses the finite chart", &
        eqcut_allowed_axis_g_chain_residual)
    call check_identity(proofs, proof_engine, &
        "axis canonical rho derivative uses the finite chart", &
        eqcut_allowed_axis_psistar_chain_residual)
    call check_identity(proofs, proof_engine, &
        "endpoint inverse Jacobian row one column one", &
        endpoint_y11*endpoint_nr0+endpoint_y12*endpoint_gr0-one)
    call check_identity(proofs, proof_engine, &
        "endpoint inverse Jacobian row one column two", &
        endpoint_y11*endpoint_nz0+endpoint_y12*endpoint_gz0)
    call check_identity(proofs, proof_engine, &
        "endpoint inverse Jacobian row two column one", &
        endpoint_y21*endpoint_nr0+endpoint_y22*endpoint_gr0)
    call check_identity(proofs, proof_engine, &
        "endpoint inverse Jacobian row two column two", &
        endpoint_y21*endpoint_nz0+endpoint_y22*endpoint_gz0-one)
    call check_identity(proofs, proof_engine, &
        "endpoint Newton map solves affine R system", &
        endpoint_newton_r_oracle-(endpoint_r0-endpoint_oracle_delta_r))
    call check_identity(proofs, proof_engine, &
        "endpoint Newton map solves affine Z system", &
        endpoint_newton_z_oracle-(endpoint_z0-endpoint_oracle_delta_z))
    call check_identity(proofs, proof_engine, &
        "endpoint Krawczyk map contracts affine R system", &
        endpoint_krawczyk_r_oracle-(endpoint_r0-endpoint_oracle_delta_r))
    call check_identity(proofs, proof_engine, &
        "endpoint Krawczyk map contracts affine Z system", &
        endpoint_krawczyk_z_oracle-(endpoint_z0-endpoint_oracle_delta_z))
    call check_identity(proofs, proof_engine, &
        "axis inverse Hessian row one column one", &
        axis_stat_y11*axis_stat_psi_rr0 &
        +axis_stat_y12*axis_stat_psi_rz0-one)
    call check_identity(proofs, proof_engine, &
        "axis inverse Hessian row one column two", &
        axis_stat_y11*axis_stat_psi_rz0 &
        +axis_stat_y12*axis_stat_psi_zz0)
    call check_identity(proofs, proof_engine, &
        "axis inverse Hessian row two column one", &
        axis_stat_y21*axis_stat_psi_rr0 &
        +axis_stat_y22*axis_stat_psi_rz0)
    call check_identity(proofs, proof_engine, &
        "axis inverse Hessian row two column two", &
        axis_stat_y21*axis_stat_psi_rz0 &
        +axis_stat_y22*axis_stat_psi_zz0-one)
    call check_identity(proofs, proof_engine, &
        "axis Newton map solves affine R stationarity", &
        axis_stat_newton_r_oracle &
        -(axis_stat_r0-axis_stat_oracle_delta_r))
    call check_identity(proofs, proof_engine, &
        "axis Newton map solves affine Z stationarity", &
        axis_stat_newton_z_oracle &
        -(axis_stat_z0-axis_stat_oracle_delta_z))
    call check_identity(proofs, proof_engine, &
        "axis Krawczyk map contracts affine R stationarity", &
        axis_stat_krawczyk_r_oracle &
        -(axis_stat_r0-axis_stat_oracle_delta_r))
    call check_identity(proofs, proof_engine, &
        "axis Krawczyk map contracts affine Z stationarity", &
        axis_stat_krawczyk_z_oracle &
        -(axis_stat_z0-axis_stat_oracle_delta_z))
    call check_identity(proofs, proof_engine, &
        "section reversal flips dpsi_star/dx", &
        dpsi_dx_reversed + dpsi_dx_section)
    call check_identity(proofs, proof_engine, &
        "section reversal flips F prime", f_prime_reversed + f_prime_section)
    call check_identity(proofs, proof_engine, &
        "section reversal flips Cdot", cdot_reversed + cdot_section)
    call check_identity(proofs, proof_engine, &
        "absolute root weight is section-orientation invariant", &
        root_weight_reversed - root_weight_section)
    call check_identity(proofs, proof_engine, &
        "positive crossing density is section-orientation invariant", &
        positive_crossing_reversed - positive_crossing_density)
    call check_identity(proofs, proof_engine, &
        "root delta weight is section-orientation independent", &
        section_root_residual)
    call check_identity(proofs, proof_engine, &
        "positive crossing density is cut-orientation independent", &
        crossing_orientation_residual)
    call check_identity(proofs, proof_engine, &
        "section period is orientation invariant", period_reversed-period_section)
    call check_identity(proofs, proof_engine, &
        "section frequency is orientation invariant", &
        frequency_reversed-frequency_section)
    call check_identity(proofs, proof_engine, "signed crossing density", &
        crossing_signed_density - crossing_signed_measure*crossing_cdot)
    call check_identity(proofs, proof_engine, &
        "positive crossing density uses absolute factors", &
        crossing_positive_density - crossing_abs_measure*crossing_abs_cdot)
    call check_identity(proofs, proof_engine, "bilinear partition of unity", &
        interp_partition_residual)
    call check_identity(proofs, proof_engine, "bilinear physical R lower corner", &
        interp_u0_residual)
    call check_identity(proofs, proof_engine, "bilinear physical R upper corner", &
        interp_u1_residual)
    call check_identity(proofs, proof_engine, "bilinear physical Z lower corner", &
        interp_v0_residual)
    call check_identity(proofs, proof_engine, "bilinear physical Z upper corner", &
        interp_v1_residual)
    call check_identity(proofs, proof_engine, "bilinear corner (0,0) reproduction", &
        interp_corner00)
    call check_identity(proofs, proof_engine, "bilinear corner (1,0) reproduction", &
        interp_corner10)
    call check_identity(proofs, proof_engine, "bilinear corner (0,1) reproduction", &
        interp_corner01)
    call check_identity(proofs, proof_engine, "bilinear corner (1,1) reproduction", &
        interp_corner11)
    call check_identity(proofs, proof_engine, "affine endpoint s=0 exactness", &
        affine_endpoint_zero - affine_a)
    call check_identity(proofs, proof_engine, "affine endpoint s=1 exactness", &
        affine_endpoint_one - (affine_a + affine_b))
    call check_identity(proofs, proof_engine, "affine endpoint derivative exactness", &
        affine_endpoint_slope - affine_b)
    call check_identity(proofs, proof_engine, "axis rho derivative regularity", &
        diff(affine_a + affine_b*endpoint_rho**2, endpoint_rho) - &
        2*endpoint_rho*affine_b)
    call check_identity(proofs, proof_engine, "axis derivative vanishes at rho=0", &
        endpoint_axis_residual)
    call check_identity(proofs, proof_engine, &
        "manufactured cylindrical b magnitude", oracle_bmod**2 - &
        (oracle_b(1)**2 + oracle_b(2)**2 + oracle_b(3)**2))
    call check_identity(proofs, proof_engine, &
        "manufactured cylindrical grad magnitude", oracle_grad_b(1) - &
        diff(oracle_bmod, oracle_coordinate(1)))
    call check_identity(proofs, proof_engine, &
        "manufactured cylindrical unit-vector derivative", oracle_dbhat(1,1) - &
        diff(oracle_bhat(1), oracle_coordinate(1)))
    call check_identity(proofs, proof_engine, &
        "manufactured cylindrical curl R", oracle_curl(1) - &
        (diff(oracle_bhat(3), oracle_coordinate(2)) - &
        diff(oracle_bhat(2), oracle_coordinate(3))))
    call check_identity(proofs, proof_engine, &
        "affine Omega_E segment integral", delta_phi_segment - &
        (potential_psi1-potential_psi0)* &
        (potential_omega0+potential_omega1)/(2*potential_c))
    call check_identity(proofs, proof_engine, &
        "profile potential endpoint reversal antisymmetry", &
        delta_phi_reversed + delta_phi_segment)
    call check_identity(proofs, proof_engine, &
        "constant Omega_E segment limit", &
        ((potential_psi1-potential_psi0)* &
        (omega_constant+omega_constant)/(2*potential_c)) - &
        delta_phi_constant)
    call check_identity(proofs, proof_engine, &
        "quadratic profile potential reproduces left endpoint", &
        subs(potential_phi, potential_psi, potential_psi0)-potential_phi0)
    call check_identity(proofs, proof_engine, &
        "quadratic profile potential reproduces integrated right endpoint", &
        subs(potential_phi, potential_psi, potential_psi1)-potential_phi1)
    call check_identity(proofs, proof_engine, &
        "profile potential derivative is local affine Omega_E over c", &
        potential_dphi_dpsi - &
        (potential_omega0+potential_omega_slope* &
        potential_delta_psi)/potential_c)
    call check_identity(proofs, proof_engine, &
        "profile potential curvature is affine Omega_E slope over c", &
        potential_d2phi_dpsi2-potential_omega_slope/potential_c)
    call check_identity(proofs, proof_engine, &
        "profile potential segment reversal leaves physical Phi unchanged", &
        potential_reversal_residual)
    call check_identity(proofs, proof_engine, &
        "flux segment reproduces axis-side endpoint", &
        flux_segment_axis_residual)
    call check_identity(proofs, proof_engine, &
        "flux segment reproduces edge-side endpoint", &
        flux_segment_edge_residual)
    call check_identity(proofs, proof_engine, &
        "flux segment normalized derivative is exact", &
        flux_segment_derivative_residual)
    call check_identity(proofs, proof_engine, &
        "flux segment forward map inverts exactly", &
        flux_segment_forward_inverse_residual)
    call check_identity(proofs, proof_engine, &
        "flux segment inverse map reproduces target", &
        flux_segment_inverse_forward_residual)
    call check_identity(proofs, proof_engine, &
        "scaled equilibrium flux normalization round trip", &
        scaled_flux_round_trip_residual)
    call check_identity(proofs, proof_engine, &
        "physical equilibrium flux normalization round trip", &
        physical_flux_normalization_residual)
    call check_identity(proofs, proof_engine, &
        "physical equilibrium flux applies field scale exactly once", &
        physical_flux_psi/physical_flux_field_scale-physical_flux_raw_psi)
    call check_identity(proofs, proof_engine, &
        "physical equilibrium flux derivative applies the same scale", &
        physical_flux_dpsi_dr/physical_flux_field_scale - &
        physical_flux_raw_dpsi_dr)
    call check_identity(proofs, proof_engine, &
        "Eq13 path field magnitude value", &
        subs(allowed_bmod_local, allowed_x, zero)-allowed_bmod)
    call check_identity(proofs, proof_engine, &
        "Eq13 path field magnitude first derivative", &
        subs(diff(allowed_bmod_local, allowed_x), allowed_x, zero) - &
        allowed_dbmod_dr)
    call check_identity(proofs, proof_engine, &
        "Eq13 path field norm second derivative", &
        2*allowed_g_second_series-allowed_d2g_dr2)
    call check_identity(proofs, proof_engine, &
        "Eq13 square-root field norm second chain rule", &
        2*(allowed_dsqrt_g_dr**2 + &
        allowed_sqrt_g*allowed_d2sqrt_g_dr2)-allowed_d2g_dr2)
    call check_identity(proofs, proof_engine, &
        "Eq13 field magnitude second product rule", &
        eqcut_radius*allowed_d2bmod_dr2 + 2*allowed_dbmod_dr - &
        allowed_field_scale*allowed_d2sqrt_g_dr2)
    call check_identity(proofs, proof_engine, &
        "Eq13 allowed energy value", &
        subs(allowed_energy_local, allowed_x, zero)-allowed_energy_margin)
    call check_identity(proofs, proof_engine, &
        "Eq13 allowed energy first derivative", &
        subs(diff(allowed_energy_local, allowed_x), allowed_x, zero) - &
        allowed_denergy_dr)
    call check_identity(proofs, proof_engine, &
        "Eq13 potential second path chain rule", &
        2*allowed_potential_second_series - &
        (allowed_d2phi_dpsi2*(allowed_field_scale*allowed_dpsi_dr)**2 + &
        allowed_dphi_dpsi*allowed_field_scale*allowed_d2psi_dr2))
    call check_identity(proofs, proof_engine, &
        "Eq13 allowed energy second derivative composition", &
        allowed_d2energy_dr2 + allowed_jk*allowed_d2omega_c_dr2 + &
        charge*(allowed_d2phi_dpsi2* &
        (allowed_field_scale*allowed_dpsi_dr)**2 + &
        allowed_dphi_dpsi*allowed_field_scale*allowed_d2psi_dr2))
    call check_identity(proofs, proof_engine, &
        "simple-turn root coordinate has quadratic physical displacement", &
        turning_delta_x-turning_direction*turning_y**2)
    call check_identity(proofs, proof_engine, &
        "simple-turn parallel speed squares to local vparallel squared", &
        turning_vparallel**2-allowed_parallel_sign**2* &
        turning_direction*turning_dv2_dx*turning_y**2)
    call check_identity(proofs, proof_engine, &
        "simple-turn canonical coordinate reaches flux root", &
        subs(turning_psi_star, turning_y, zero)-turning_psi_root)
    call check_identity(proofs, proof_engine, &
        "simple-turn canonical root derivative is finite", &
        turning_root_derivative-allowed_parallel_sign*mass*c_light/charge* &
        turning_ratio_root*sqrt(turning_direction*turning_dv2_dx))
    if (proofs%failed /= 0) error stop "full-FOW symbolic proof failed"
    call suite_end(proofs)

    roots = [jk, omega_c, jk_omega_c, jk_candidate, jk_max, d_jk_d_h, &
        d_jk_d_phi, d_jk_d_omega, eq4_coeff, psi_star, dpsi_dp_phi, &
        q_phi_energy, energy_on_shell, c_on_shell, eta_on_shell, &
        perturbation_ratio, boozer_ratio, f_prime, freq_weight, phase_weight, &
        force_bracket, eq17_outer_factor, bstar1, bstar2, bstar3, &
        bparallel_star, force1, force2, force3, cross1, cross2, cross3, &
        v1, v2, v3, dot_r, dot_z, dot_phi, dot_p, dot_mu, &
        cylindrical_measure, hamiltonian_dot, h_axisymmetric, &
        canonical_p_phi_cylindrical, &
        canonical_p_phi, d_lagrangian_d_phi, d_lagrangian_d_phi_dot, &
        canonical_residual]
    refinement_roots = [scale_r, scale_z, scale_phi, scale_p_parallel, scale_mu, &
        scale_tau, scale_omega_b, scale_omega_phi, scale_hm_real, scale_hm_imag, &
        scale_hm_squared, scale_root, scale_torque_density, scale_torque]
    boundary_roots = [c_tau, regular_form, reflecting_form, reflecting_derivative, &
        separatrix_form, separatrix_tau, xpoint_form, xpoint_tau, &
        regular_classification_residual, reflecting_classification_residual, &
        separatrix_classification_residual, xpoint_pair_coefficient_residual, &
        xpoint_pair_tau_residual]
    cut_roots = [cut_cross1, cut_cross2, cut_cross3, cut_c, &
        cut_c - (cut_cross1*grad_phi_coordinate1 + &
        cut_cross2*grad_phi_coordinate2 + cut_cross3*grad_phi_coordinate3)]
    cutdot_roots = [cut_cdot, cutdot_transversality, cutdot_orientation]
    section_roots = [dpsi_dx_reversed, f_prime_reversed, cdot_reversed, &
        root_weight_section, root_weight_reversed, positive_crossing_density, &
        positive_crossing_reversed, period_reversed, frequency_reversed]
    crossing_roots = [crossing_signed_measure, crossing_signed_density, &
        crossing_abs_measure, crossing_abs_cdot, crossing_positive_density]
    geometry_roots(1) = geom_bmod
    geometry_roots(2:4) = geom_bhat
    geometry_roots(5:7) = geom_grad_b
    geometry_roots(8) = geom_dbhat(1,1)
    geometry_roots(9) = geom_dbhat(2,1)
    geometry_roots(10) = geom_dbhat(3,1)
    geometry_roots(11) = geom_dbhat(1,2)
    geometry_roots(12) = geom_dbhat(2,2)
    geometry_roots(13) = geom_dbhat(3,2)
    geometry_roots(14) = geom_dbhat(1,3)
    geometry_roots(15) = geom_dbhat(2,3)
    geometry_roots(16) = geom_dbhat(3,3)
    geometry_roots(17:19) = geom_curl
    interpolation_roots = [interp_u, interp_v, interp_w00, interp_w10, &
        interp_w01, interp_w11, interp_value_real, interp_value_imag, &
        interp_partition_residual]
    endpoint_roots = [endpoint_value_zero, endpoint_value_one, &
        endpoint_derivative_zero, endpoint_derivative_one, endpoint_slope, &
        endpoint_intercept, endpoint_dfd_rho, endpoint_axis_residual]
    profile_potential_roots = [delta_phi_segment, delta_phi_reversed, &
        delta_phi_constant]
    profile_potential_map_roots = [potential_phi, potential_dphi_dpsi, &
        potential_d2phi_dpsi2]
    eqdsk_flux_profile_segment_roots = [flux_segment_psihat, &
        flux_segment_dpsihat_dstor, flux_segment_inverse_s]
    eqdsk_scaled_flux_normalization_roots = [scaled_flux_psihat]
    eqdsk_physical_flux_normalization_roots = [physical_flux_psihat]
    eqdsk_physical_flux_map_roots = [physical_flux_psi, &
        physical_flux_dpsi_dr]
    eqdsk_allowed_energy_roots = [allowed_g, &
        allowed_equilibrium_psi_physical, &
        allowed_equilibrium_dpsi_physical_dr, allowed_bmod, &
        allowed_dbmod_dr, allowed_d2bmod_dr2, allowed_omega_c, &
        allowed_domega_c_dr, allowed_d2omega_c_dr2, &
        allowed_energy_margin, allowed_denergy_dr, &
        allowed_d2energy_dr2, allowed_vparallel_squared, &
        allowed_dvparallel_squared_dr, allowed_d2vparallel_squared_dr2, &
        allowed_bphi_covariant, allowed_dbphi_covariant_dr]
    eqdsk_canonical_cut_roots = [allowed_vparallel, &
        allowed_dvparallel_dr, allowed_psi_star, allowed_dpsi_star_dr]
    eqdsk_turning_chart_roots = [turning_delta_x, turning_vparallel, &
        turning_psi_star, turning_dpsi_star_dy, turning_root_derivative]
    eqdsk_cell_jet_roots = [cell_psi, cell_psi_r, cell_psi_z, cell_psi_rr, &
        cell_psi_rz, cell_psi_zz, cell_psi_rrr, cell_psi_rrz, &
        cell_psi_rzz, cell_psi_zzz]
    eqdsk_cell_fourth_jet_roots = [cell_psi_rrrr, cell_psi_rrrz, &
        cell_psi_rrzz, cell_psi_rzzz, cell_psi_zzzz]
    eqdsk_profile_jet_roots = [profile_value, profile_first, profile_second, &
        profile_vacuum_f]
    eqdsk_cut_jet_roots = [eqcut_c, eqcut_dc_dr, eqcut_dc_darc_phi, &
        eqcut_dc_dz, eqcut_cdot, eqcut_abs_cdot, &
        eqcut_orientation_scalar]
    eqdsk_cut_numerator_roots = [eqcut_n, eqcut_n_r, eqcut_n_z]
    eqdsk_cut_interval_roots = [eqcut_n, eqcut_n_r, eqcut_n_z, eqcut_g]
    canonical_numerator_roots = [canonical_n, canonical_n_r, canonical_n_rr]
    canonical_turning_roots = [canonical_turning_n_root, &
        canonical_turning_n_y_root]
    canonical_symmetry_roots = [canonical_n_coordinate, &
        canonical_n_r_coordinate, canonical_n_rr_coordinate, &
        canonical_time_phase_reversed, canonical_time_resonance_reversed, &
        canonical_fourier_pair_resonance, canonical_combined_resonance, &
        canonical_toroidal_phase_relabelled, &
        canonical_p_phi_covariant_relabelled, canonical_sigma_difference, &
        canonical_time_phase_residual, canonical_time_resonance_residual, &
        canonical_fourier_pair_phase_residual, &
        canonical_fourier_pair_resonance_residual, &
        canonical_combined_resonance_residual, &
        canonical_endpoint_variation_residual, &
        pair_field_residual, &
        canonical_toroidal_phase_residual, canonical_toroidal_p_phi_residual, &
        toroidal_power_residual, torque_component_difference]
    eqdsk_cut_r_chart_roots = [eqcut_r_chart_slope, eqcut_r_chart_ds]
    eqdsk_cut_z_chart_roots = [eqcut_z_chart_slope, eqcut_z_chart_ds]
    eqdsk_cut_r_flux_chart_roots = [eqcut_r_chart_slope, &
        eqcut_r_chart_dpsihat]
    eqdsk_cut_mean_value_roots = [eqcut_mean_enclosure]
    eqdsk_cut_numerator_hessian_roots = [eqcut_n_rr, eqcut_n_rz, eqcut_n_zz]
    eqdsk_cut_r_flux_curvature_roots = [eqcut_r_chart_slope, &
        eqcut_r_chart_second, eqcut_r_chart_dpsihat, &
        eqcut_r_chart_d2psihat]
    eqdsk_cut_axis_curvature_roots = [eqcut_axis_curvature, &
        eqcut_axis_hessian_determinant]
    eqdsk_cut_axis_limit_roots = [eqcut_axis_curvature, &
        eqcut_axis_abs_dpsihat, eqcut_axis_inverse_delta_r]
    eqdsk_rho_tor_map_roots = [eqcut_s_tor, eqcut_rho_dstor_drho]
    eqdsk_s_tor_to_rho_roots = [eqcut_rho_from_s]
    eqdsk_flux_profile_rho_chain_roots = [eqcut_dpsihat_drho]
    eqdsk_cut_flux_coordinate_roots = [eqcut_dpsihat_drho, eqcut_dR_drho, &
        eqcut_dZ_drho]
    eqdsk_cut_axis_rho_limit_roots = [eqcut_axis_curvature, &
        eqcut_axis_rho_dR, eqcut_axis_rho_dZ]
    eqdsk_allowed_rho_chain_roots = [eqcut_dpsihat_drho, &
        eqcut_d2psihat_drho2, eqcut_dR_drho, eqcut_d2R_drho2, &
        eqcut_allowed_g_rho, eqcut_allowed_g_rhorho, &
        eqcut_allowed_psistar_rho]
    eqdsk_allowed_axis_rho_chain_roots = [eqcut_allowed_axis_g_rho, &
        eqcut_allowed_axis_psistar_rho]
    eqdsk_cut_endpoint_system_roots = [endpoint_n, endpoint_flux_residual, &
        endpoint_n_r, endpoint_n_z, endpoint_flux_r, endpoint_flux_z]
    eqdsk_cut_endpoint_newton_roots = [endpoint_determinant, &
        endpoint_newton_r, endpoint_newton_z, endpoint_y11, endpoint_y12, &
        endpoint_y21, endpoint_y22]
    eqdsk_cut_endpoint_krawczyk_roots = [endpoint_krawczyk_r, &
        endpoint_krawczyk_z]
    eqdsk_axis_stationarity_system_roots = [axis_stat_psi_r, &
        axis_stat_psi_z, axis_stat_psi_rr, axis_stat_psi_rz, &
        axis_stat_psi_rz, axis_stat_psi_zz]
    eqdsk_axis_stationarity_newton_roots = [axis_stat_determinant, &
        axis_stat_newton_r, axis_stat_newton_z, axis_stat_y11, axis_stat_y12, &
        axis_stat_y21, axis_stat_y22]
    eqdsk_axis_stationarity_krawczyk_roots = [axis_stat_krawczyk_r, &
        axis_stat_krawczyk_z]
    eq17_outer_roots = [eq17_outer_factor]
    frequency_contribution_roots = [frequency_contribution, phase_contribution]
    frequency_identity_roots = [n_squared_frequency_contribution, &
        phase_contribution_identity]
    simple_root_force_roots = [root_hm_squared, simple_root_weight, &
        simple_root_force_contribution]
    axisymmetric_pphi_roots = [axis_pphi, axis_pphi_dot, axis_pphi_residual]
    sign_symmetry_roots = [sign_amplitude_residual, phase_amplitude_residual, &
        pair_field_residual, fixed_conjugation_difference, &
        conjugate_modulus_residual, section_root_residual, &
        crossing_orientation_residual, gauge_force_residual, &
        gauge_outer_ledger_residual, toroidal_resonance_residual, &
        toroidal_power_residual, torque_component_difference, &
        charge_reversal_energy_difference, &
        potential_reversal_energy_difference]
    harmonic_integrand_roots = [eq4_coeff, c_on_shell, &
        harmonic_shell_residual, harmonic_phase_argument, &
        harmonic_phase_real, harmonic_phase_imag, harmonic_integrand_real, &
        harmonic_integrand_imag]
    normalization_roots = [norm_phi_eff, norm_jk_scale, norm_h_hat, &
        norm_jk_hat, norm_psi_hat, norm_dpsi_hat_dx, norm_tau_hat, &
        norm_omega_b_hat, norm_omega_phi_hat, norm_domega_b_hat_dx, &
        norm_domega_phi_hat_dx, norm_hm_real_hat, norm_hm_imag_hat]
    quadrature_map_roots = [quad_h_hat, quad_h_phys, quad_weight_h, &
        quad_jk_hat_max, quad_jk_hat, quad_jk_phys, quad_weight_j, &
        quad_weight]
    gauss_interval_map_roots = [gauss_mapped_node, gauss_mapped_weight]
    polynomial_enclosure_roots = [poly_tail_bound, poly_lower, poly_upper]
    cylindrical_bstar_roots = [bstar1, bstar2, bstar3, &
        bparallel_star, cylindrical_measure]
    physical_mu_roots = [mu_from_vperp, vperp_squared_from_mu, &
        specific_energy_from_mu]
    buchholz_action_roots = [jk_from_mu, mu_from_jk]
    buchholz_energy_roots = [conversion_omega_c, &
        vperp_squared_from_jk, specific_energy_from_jk]
    buchholz_specific_energy_roots = [conversion_jk*omega_c_var/mass]
    potato_velocity_roots = [jpotato_from_vperp, &
        vperp_squared_from_jpotato]
    potato_mu_roots = [jpotato_from_mu, mu_from_jpotato]
    cylindrical_hamiltonian_roots = [invariant_h]
    cylindrical_canonical_roots = [ &
        subs(invariant_pphi, bhat2, cylindrical_bhat_phi), &
        subs(invariant_psistar, bhat2, cylindrical_bhat_phi)]
    cylindrical_vparallel_roots = [invariant_vparallel_squared]
    cylindrical_launch_roots = [ &
        subs(invariant_vparallel_squared, bhat2, cylindrical_bhat_phi), &
        subs(invariant_ppar_from_pphi, bhat2, cylindrical_bhat_phi), &
        subs(invariant_launch_residual, bhat2, cylindrical_bhat_phi)]

    do k = 1, size(roots)
        simplified = simplify_engine%simplify(roots(k))
        if (.not. simplified%ok) then
            write (output_unit, "(a,i0)") "simplification failed for output ", k
            error stop 1
        end if
        roots(k) = simplified%value
    end do
    call simplify_array(refinement_roots)
    call simplify_array(boundary_roots)
    call simplify_array(cut_roots)
    call simplify_array(cutdot_roots)
    call simplify_array(section_roots)
    call simplify_array(crossing_roots)
    call simplify_array(geometry_roots)
    call simplify_array(interpolation_roots)
    call simplify_array(endpoint_roots)
    call simplify_array(profile_potential_roots)
    call simplify_array(profile_potential_map_roots)
    call simplify_array(eqdsk_flux_profile_segment_roots)
    call simplify_array(eqdsk_scaled_flux_normalization_roots)
    call simplify_array(eqdsk_physical_flux_normalization_roots)
    call simplify_array(eqdsk_physical_flux_map_roots)
    call simplify_array(eqdsk_allowed_energy_roots)
    call simplify_array(eqdsk_canonical_cut_roots)
    call simplify_array(eqdsk_turning_chart_roots)
    call simplify_array(eqdsk_cell_jet_roots)
    call simplify_array(eqdsk_cell_fourth_jet_roots)
    call simplify_array(eqdsk_profile_jet_roots)
    call simplify_array(eqdsk_cut_jet_roots)
    call simplify_array(eqdsk_cut_numerator_roots)
    call simplify_array(eqdsk_cut_interval_roots)
    call simplify_array(canonical_numerator_roots)
    call simplify_array(canonical_turning_roots)
    call simplify_array(canonical_symmetry_roots)
    call simplify_array(eqdsk_cut_r_chart_roots)
    call simplify_array(eqdsk_cut_z_chart_roots)
    call simplify_array(eqdsk_cut_r_flux_chart_roots)
    call simplify_array(eqdsk_cut_mean_value_roots)
    call simplify_array(eqdsk_cut_numerator_hessian_roots)
    call simplify_array(eqdsk_cut_r_flux_curvature_roots)
    call simplify_array(eqdsk_cut_axis_curvature_roots)
    call simplify_array(eqdsk_cut_axis_limit_roots)
    call simplify_array(eqdsk_rho_tor_map_roots)
    call simplify_array(eqdsk_s_tor_to_rho_roots)
    call simplify_array(eqdsk_flux_profile_rho_chain_roots)
    call simplify_array(eqdsk_cut_flux_coordinate_roots)
    call simplify_array(eqdsk_cut_axis_rho_limit_roots)
    call simplify_array(eqdsk_allowed_rho_chain_roots)
    call simplify_array(eqdsk_allowed_axis_rho_chain_roots)
    call simplify_array(eqdsk_cut_endpoint_system_roots)
    call simplify_array(eqdsk_cut_endpoint_newton_roots)
    call simplify_array(eqdsk_cut_endpoint_krawczyk_roots)
    call simplify_array(eqdsk_axis_stationarity_system_roots)
    call simplify_array(eqdsk_axis_stationarity_newton_roots)
    call simplify_array(eqdsk_axis_stationarity_krawczyk_roots)
    call simplify_array(eq17_outer_roots)
    call simplify_array(frequency_contribution_roots)
    call simplify_array(frequency_identity_roots)
    call simplify_array(simple_root_force_roots)
    call simplify_array(axisymmetric_pphi_roots)
    call simplify_array(sign_symmetry_roots)
    call simplify_array(harmonic_integrand_roots)
    call simplify_array(normalization_roots)
    call simplify_array(quadrature_map_roots)
    call simplify_array(gauss_interval_map_roots)
    call simplify_array(polynomial_enclosure_roots)
    call simplify_array(cylindrical_bstar_roots)
    call simplify_array(physical_mu_roots)
    call simplify_array(buchholz_action_roots)
    call simplify_array(buchholz_energy_roots)
    call simplify_array(buchholz_specific_energy_roots)
    call simplify_array(potato_velocity_roots)
    call simplify_array(potato_mu_roots)
    call simplify_array(cylindrical_hamiltonian_roots)
    call simplify_array(cylindrical_canonical_roots)
    call simplify_array(cylindrical_vparallel_roots)
    call simplify_array(cylindrical_launch_roots)

    action_roots = roots(1:11)
    perturbation_roots = roots(12:17)
    resonance_roots = roots(18:20)
    eq17_roots = [roots(12), roots(21)]
    cylindrical_roots = roots(23:44)
    noether_roots = roots(45:48)

    call emit_kernel_file(trim(output_path)// &
        "/neort_full_fow_action_symbolic.f90", &
        "neort_full_fow_action_symbolic", "evaluate_neort_action_normalization", &
        [character(len=64) :: "mass", "charge", "c_light", "mu", "bmod", "h", &
        "electrostatic_potential", "p_phi", "omega_c_value"], action_roots, &
        [character(len=64) :: "j_k", "omega_c", "j_k_omega_c", "j_k_candidate", "j_k_max", &
        "d_j_k_candidate_d_h", "d_j_k_candidate_d_potential", &
        "d_j_k_candidate_d_omega_c", "eq4_coefficient", "psi_star", &
        "d_psi_star_d_p_phi"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_full_fow_perturbation_symbolic.f90", &
        "neort_full_fow_perturbation_symbolic", &
        "evaluate_neort_perturbation_coefficient", &
        [character(len=64) :: "mass", "charge", "mu", "bmod", "h", &
        "electrostatic_potential", "p_parallel", "temperature"], &
        perturbation_roots, [character(len=64) :: "q_phi_energy", "e_on_shell", &
        "c_on_shell", "eta_on_shell", "perturbation_ratio", "boozer_ratio"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_full_fow_harmonic_symbolic.f90", &
        "neort_full_fow_harmonic_symbolic", &
        "evaluate_neort_full_fow_harmonic_integrand", &
        [character(len=64) :: "mass", "charge", "mu", "bmod", "h", &
        "electrostatic_potential", "p_parallel", &
        "delta_b_real", "delta_b_imag", "m_mode", "n_mode", &
        "orbit_phi", "launch_phi", "orbit_time", "omega_b", "omega_phi"], &
        harmonic_integrand_roots, [character(len=64) :: "eq4_coefficient", &
        "c_on_shell", "on_shell_residual", "phase_argument", &
        "phase_real", "phase_imag", &
        "integrand_real", "integrand_imag"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_full_fow_resonance_symbolic.f90", &
        "neort_full_fow_resonance_symbolic", &
        "evaluate_neort_resonance_weights", [character(len=64) :: "n_mode", "tau", "g_prime"], &
        resonance_roots, [character(len=64) :: "frequency_phase_derivative", &
        "frequency_root_weight", "phase_root_weight"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_full_fow_frequency_contribution_symbolic.f90", &
        "neort_full_fow_frequency_contribution_symbolic", &
        "evaluate_neort_frequency_root_contribution", &
        [character(len=64) :: "dpsi_star_dx", "abs_Hm_squared", "tau_b", "n_mode", "g_prime"], &
        frequency_contribution_roots, [character(len=64) :: "frequency_root_contribution", &
        "phase_root_contribution"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_full_fow_frequency_identity_symbolic.f90", &
        "neort_full_fow_frequency_identity_symbolic", &
        "evaluate_neort_frequency_phase_identity", [character(len=64) :: "dpsi_star_dx", &
        "abs_Hm_squared", "tau_b", "n_mode", "g_prime"], &
        frequency_identity_roots, [character(len=64) :: "n_squared_frequency_contribution", &
        "phase_contribution_identity"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_full_fow_simple_root_symbolic.f90", &
        "neort_full_fow_simple_root_symbolic", &
        "evaluate_neort_full_fow_simple_root_force", &
        [character(len=64) :: "dpsi_star_dx", "Hm_real", "Hm_imag", &
        "tau_b", "frequency_residual_derivative", "force_value"], &
        simple_root_force_roots, [character(len=64) :: "abs_Hm_squared", &
        "simple_root_weight", "force_contribution"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_full_fow_sign_symmetry_symbolic.f90", &
        "neort_full_fow_sign_symmetry_symbolic", &
        "evaluate_neort_full_fow_sign_symmetry_ledger", &
        [character(len=64) :: "a_real", "a_imag", &
        "m_mode", "n_mode", "theta", "phi", "amplitude_phase", &
        "dpsi_star_dx", "F_prime", "Cdot", "signed_R_Bparallel_star", &
        "h", "charge", "electrostatic_potential", "temperature", "a1", &
        "a2", "gauge_C", "sign_omega_b", "sign_omega_phi", &
        "sign_torque_phi", "sign_rotation_phi"], sign_symmetry_roots, &
        [character(len=64) :: "amplitude_sign_squared_residual", &
        "amplitude_phase_squared_residual", "real_field_pair_residual", &
        "fixed_mode_conjugation_field_difference", &
        "fixed_mode_conjugation_modulus_residual", &
        "section_root_weight_residual", &
        "cut_crossing_density_residual", "gauge_force_residual", &
        "gauge_outer_residual", "toroidal_resonance_relabel_residual", &
        "toroidal_power_relabel_residual", &
        "toroidal_torque_component_difference", &
        "charge_reversal_potential_energy_difference", &
        "potential_sign_reversal_energy_difference"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_full_fow_normalization_symbolic.f90", &
        "neort_full_fow_normalization_symbolic", &
        "evaluate_neort_full_fow_normalization", [character(len=64) :: "mass", "charge", &
        "c_light", "Eref", "reference_velocity", "h", "J_K_physical", &
        "psi_star_physical", "dpsi_star_dx_physical", "tau_b_physical", &
        "omega_b_physical", "omega_phi_physical", &
        "domega_b_dx_physical", "domega_phi_dx_physical", &
        "Hm_real_physical", "Hm_imag_physical"], normalization_roots, &
        [character(len=64) :: "Phi_eff", "J_K_scale", "H_hat", "J_K_hat", "psi_star_hat", &
        "dpsi_star_hat_dx", "tau_b_hat", "omega_b_hat", "omega_phi_hat", &
        "domega_b_hat_dx", "domega_phi_hat_dx", "Hm_real_hat", &
        "Hm_imag_hat"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_full_fow_quadrature_map_symbolic.f90", &
        "neort_full_fow_quadrature_map_symbolic", &
        "evaluate_neort_full_fow_quadrature_map", [character(len=64) :: "mass", "charge", &
        "c_light", "Eref", "energy_unit_node", "energy_unit_weight", &
        "action_unit_node", "action_unit_weight", "qPhi_min", &
        "J_K_max_physical"], quadrature_map_roots, [character(len=64) :: "H_hat", "H_physical", &
        "energy_normalized_weight", "J_K_hat_max", "J_K_hat", &
        "J_K_physical", "action_normalized_weight", "paired_weight"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_gauss_legendre_interval_map_symbolic.f90", &
        "neort_gauss_legendre_interval_map_symbolic", &
        "evaluate_neort_gauss_legendre_interval_map", &
        [character(len=64) :: "standard_node", "standard_weight", &
        "lower_bound", "upper_bound"], gauss_interval_map_roots, &
        [character(len=64) :: "mapped_node", "mapped_weight"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_polynomial_cell_enclosure_symbolic.f90", &
        "neort_polynomial_cell_enclosure_symbolic", &
        "evaluate_neort_polynomial_cell_enclosure", [character(len=64) :: "coefficient_0", &
        "coefficient_1", "coefficient_2", "coefficient_3", &
        "coefficient_4", "coefficient_5", "cell_width"], &
        polynomial_enclosure_roots, [character(len=64) :: "tail_absolute_bound", &
        "polynomial_lower_bound", "polynomial_upper_bound"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_cylindrical_bstar_symbolic.f90", &
        "neort_cylindrical_bstar_symbolic", &
        "evaluate_neort_cylindrical_bstar", [character(len=64) :: "charge", "c_light", &
        "radius", "p_parallel", "b1", "b2", "b3", "bhat1", "bhat2", &
        "bhat3", "curl_bhat1", "curl_bhat2", "curl_bhat3"], &
        cylindrical_bstar_roots, [character(len=64) :: "b_star_1", "b_star_2", "b_star_3", &
        "b_parallel_star", "cylindrical_measure"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_physical_mu_symbolic.f90", "neort_physical_mu_symbolic", &
        "evaluate_neort_physical_mu", [character(len=64) :: "v_perp", "mass", "bmod", &
        "mu_phys"], physical_mu_roots, [character(len=64) :: "mu_phys_from_v_perp", &
        "v_perp_squared_from_mu_phys", "specific_energy_from_mu_phys"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_buchholz_action_symbolic.f90", &
        "neort_buchholz_action_symbolic", "evaluate_neort_buchholz_action", &
        [character(len=64) :: "mu_phys", "J_K", "mass", "charge", "c_light"], &
        buchholz_action_roots, [character(len=64) :: "J_K_from_mu_phys", "mu_phys_from_J_K"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_buchholz_energy_symbolic.f90", &
        "neort_buchholz_energy_symbolic", "evaluate_neort_buchholz_energy", &
        [character(len=64) :: "J_K", "mass", "charge", "c_light", "bmod"], &
        buchholz_energy_roots, [character(len=64) :: "omega_c", "v_perp_squared_from_J_K", &
        "specific_energy_from_J_K"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_buchholz_specific_energy_symbolic.f90", &
        "neort_buchholz_specific_energy_symbolic", &
        "evaluate_neort_buchholz_specific_energy", [character(len=64) :: "J_K", &
        "omega_c_value", "mass"], buchholz_specific_energy_roots, &
        [character(len=64) :: "specific_energy_from_J_K"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_potato_velocity_symbolic.f90", &
        "neort_potato_velocity_symbolic", "evaluate_neort_potato_velocity", &
        [character(len=64) :: "v_perp", "J_potato", "reference_velocity", "bmod"], &
        potato_velocity_roots, [character(len=64) :: "J_potato_from_v_perp", &
        "v_perp_squared_from_J_potato"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_potato_mu_symbolic.f90", "neort_potato_mu_symbolic", &
        "evaluate_neort_potato_mu", [character(len=64) :: "mu_phys", "J_potato", "mass", &
        "reference_velocity"], potato_mu_roots, &
        [character(len=64) :: "J_potato_from_mu_phys", "mu_phys_from_J_potato"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_cylindrical_hamiltonian_symbolic.f90", &
        "neort_cylindrical_hamiltonian_symbolic", &
        "evaluate_neort_cylindrical_hamiltonian", [character(len=64) :: "mass", "charge", &
        "mu", "bmod", "electrostatic_potential", "p_parallel"], &
        cylindrical_hamiltonian_roots, [character(len=64) :: "hamiltonian"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_cylindrical_canonical_symbolic.f90", &
        "neort_cylindrical_canonical_symbolic", &
        "evaluate_neort_cylindrical_canonical", [character(len=64) :: "charge", "c_light", &
        "radius", "p_parallel", "psi", "bhat_phi"], &
        cylindrical_canonical_roots, [character(len=64) :: "canonical_p_phi", "psi_star"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_cylindrical_vparallel_symbolic.f90", &
        "neort_cylindrical_vparallel_symbolic", &
        "evaluate_neort_cylindrical_vparallel", &
        [character(len=64) :: "h", "mu", "bmod", &
        "electrostatic_potential", "mass", "charge"], &
        cylindrical_vparallel_roots, [character(len=64) :: &
        "v_parallel_squared"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_cylindrical_launch_symbolic.f90", &
        "neort_cylindrical_launch_symbolic", &
        "evaluate_neort_cylindrical_launch", [character(len=64) :: "mass", "charge", &
        "c_light", "mu", "bmod", "electrostatic_potential", "radius", &
        "psi", "bhat_phi", "h", "p_phi"], cylindrical_launch_roots, &
        [character(len=64) :: "v_parallel_squared", "p_parallel_from_p_phi", &
        "launch_energy_residual"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_full_fow_eq17_symbolic.f90", "neort_full_fow_eq17_symbolic", &
        "evaluate_neort_eq17_force", [character(len=64) :: "h", "charge", &
        "electrostatic_potential", "temperature", "a1", "a2"], &
        eq17_roots, &
        [character(len=64) :: "q_phi_energy", "force_bracket"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_full_fow_eq17_outer_symbolic.f90", &
        "neort_full_fow_eq17_outer_symbolic", &
        "evaluate_neort_eq17_outer_factor", [character(len=64) :: "Eref", "Phi_eff", "n0", &
        "temperature", "charge", "electrostatic_potential", "h", &
        "residence"], eq17_outer_roots, [character(len=64) :: "eq17_outer_factor"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_cylindrical_littlejohn_symbolic.f90", &
        "neort_cylindrical_littlejohn_symbolic", &
        "evaluate_neort_cylindrical_littlejohn", &
        [character(len=64) :: "mass", "charge", "c_light", "mu", "bmod", &
        "electrostatic_potential", "radius", "p_parallel", "psi", "b1", &
        "b2", "b3", "bhat1", "bhat2", "bhat3", "curl_bhat1", &
        "curl_bhat2", "curl_bhat3", "grad_b1", "grad_b2", "grad_b3", &
        "grad_phi1", "grad_phi2", "grad_phi3"], cylindrical_roots, &
        [character(len=64) :: "b_star_1", "b_star_2", "b_star_3", "b_parallel_star", "force_1", &
        "force_2", "force_3", "cross_1", "cross_2", "cross_3", "v_r", &
        "v_phi", "v_z", "dot_r", "dot_z", "dot_phi", "dot_p_parallel", &
        "dot_mu", "cylindrical_measure", "d_hamiltonian_dt", &
        "hamiltonian", "canonical_p_phi"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_axisymmetric_pphi_symbolic.f90", &
        "neort_axisymmetric_pphi_symbolic", &
        "evaluate_neort_axisymmetric_p_phi", [character(len=64) :: "mass", "charge", "c_light", &
        "mu", "radius", "p_parallel", "psi_axis", "Bmod_axis", &
        "bhat_R_axis", "bhat_phi_axis", "bhat_Z_axis", &
        "curl_bhat_phi_axis", "d_bhat_phi_d_R_axis", &
        "d_bhat_phi_d_Z_axis", "grad_b_R_axis", "grad_b_Z_axis", &
        "Phi_psi_axis"], &
        axisymmetric_pphi_roots, [character(len=64) :: "canonical_p_phi", &
        "d_canonical_p_phi_dt", "p_phi_dot_minus_zero"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_axisymmetric_noether_symbolic.f90", &
        "neort_axisymmetric_noether_symbolic", &
        "evaluate_neort_axisymmetric_noether", [character(len=64) :: "charge", "c_light", "psi", &
        "b_phi_cov", "p_parallel", "h", "phi_dot"], noether_roots, &
        [character(len=64) :: "canonical_p_phi", "d_canonical_p_phi_dt_from_noether", &
        "d_lagrangian_d_phi_dot", "canonical_p_phi_residual"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_full_fow_refinement_symbolic.f90", &
        "neort_full_fow_refinement_symbolic", &
        "evaluate_neort_full_fow_refinement_scales", &
        [character(len=64) :: "state_R", "state_Z", "state_phi", "state_p_parallel", "state_mu", &
        "tau_observable", "omega_b_observable", "omega_phi_observable", &
        "Hm_real", "Hm_imag", "root_derivative", "torque_density", "torque", &
        "R_ref", "Z_ref", "phi_ref", "p_parallel_ref", "mu_ref", "tau_ref", &
        "omega_b_ref", "omega_phi_ref", "Hm_real_ref", "Hm_imag_ref", &
        "Hm_squared_ref", "root_derivative_ref", "torque_density_ref", &
        "torque_ref"], refinement_roots, &
        [character(len=64) :: "scale_R", "scale_Z", "scale_phi", "scale_p_parallel", "scale_mu", &
        "scale_tau", "scale_omega_b", "scale_omega_phi", "scale_Hm_real", &
        "scale_Hm_imag", "scale_Hm_squared", "scale_root_derivative", &
        "scale_torque_density", "scale_torque"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_buchholz_boundary_symbolic.f90", &
        "neort_buchholz_boundary_symbolic", &
        "evaluate_neort_buchholz_boundary_limits", [character(len=64) :: "boundary_distance", &
        "boundary_reference_scale", "regular_coefficient", "reflecting_coefficient", &
        "lambda_positive"], boundary_roots, &
        [character(len=64) :: "C_tau", "regular_form", "reflecting_sqrt_form", &
        "reflecting_derivative_form", "separatrix_log_form", &
        "separatrix_positive_tau_form", "xpoint_log_form", &
        "xpoint_positive_tau_form", &
        "regular_classification_residual", &
        "reflecting_classification_residual", &
        "separatrix_classification_residual", &
        "xpoint_pair_coefficient_residual", "xpoint_pair_tau_residual"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_buchholz_cut_symbolic.f90", "neort_buchholz_cut_symbolic", &
        "evaluate_neort_buchholz_cut", [character(len=64) :: "grad_b1", "grad_b2", "grad_b3", &
        "grad_psi1", "grad_psi2", "grad_psi3", "grad_phi_coordinate1", &
        "grad_phi_coordinate2", "grad_phi_coordinate3"], cut_roots, &
        [character(len=64) :: "cut_cross_1", "cut_cross_2", "cut_cross_3", "cut_C", &
        "cut_identity_residual"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_buchholz_cdot_symbolic.f90", "neort_buchholz_cdot_symbolic", &
        "evaluate_neort_buchholz_cdot", [character(len=64) :: "d_cut_C_d_R", &
        "d_cut_C_d_arc_phi", "d_cut_C_d_Z", "dot_cut_R", &
        "dot_cut_arc_phi", "dot_cut_Z"], cutdot_roots, &
        [character(len=64) :: "cut_Cdot", "transversality_scalar", "orientation_scalar"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_cylindrical_crossing_symbolic.f90", &
        "neort_cylindrical_crossing_symbolic", &
        "evaluate_neort_cylindrical_crossing_density", &
        [character(len=64) :: "R_Bparallel_star", "Cdot"], crossing_roots, &
        [character(len=64) :: "signed_R_Bparallel_star", "signed_crossing_density", &
        "absolute_R_Bparallel_star", "absolute_Cdot", &
        "positive_phase_space_crossing_density"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_full_fow_section_reversal_symbolic.f90", &
        "neort_full_fow_section_reversal_symbolic", &
        "evaluate_neort_section_reversal", [character(len=64) :: "dpsi_star_dx", "F_prime", &
        "Cdot", "signed_R_Bparallel_star", "period", "frequency"], &
        section_roots, [character(len=64) :: "dpsi_star_dx_reversed", "F_prime_reversed", &
        "Cdot_reversed", "absolute_root_weight", &
        "absolute_root_weight_reversed", "positive_crossing_density", &
        "positive_crossing_density_reversed", "period_reversed", &
        "frequency_reversed"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_cylindrical_geometry_symbolic.f90", &
        "neort_cylindrical_geometry_symbolic", &
        "evaluate_neort_cylindrical_geometry", [character(len=64) :: "radius", "b_R", &
        "b_phi", "b_Z", "d_b_R_d_R", "d_b_phi_d_R", "d_b_Z_d_R", &
        "d_b_R_d_arc_phi", "d_b_phi_d_arc_phi", "d_b_Z_d_arc_phi", &
        "d_b_R_d_Z", "d_b_phi_d_Z", "d_b_Z_d_Z"], geometry_roots, &
        [character(len=64) :: "bmod", "bhat_R", "bhat_phi", "bhat_Z", "grad_b_R", &
        "grad_b_arc_phi", "grad_b_Z", "dbhat_R_dR", "dbhat_phi_dR", &
        "dbhat_Z_dR", "dbhat_R_darc_phi", "dbhat_phi_darc_phi", &
        "dbhat_Z_darc_phi", "dbhat_R_dZ", "dbhat_phi_dZ", "dbhat_Z_dZ", &
        "curl_bhat_R", "curl_bhat_phi", "curl_bhat_Z"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_cylindrical_bilinear_complex_symbolic.f90", &
        "neort_cylindrical_bilinear_complex_symbolic", &
        "evaluate_neort_cylindrical_bilinear_complex", &
        [character(len=64) :: "coordinate_R", "coordinate_Z", "cell_R0", &
        "cell_R1", "cell_Z0", "cell_Z1", "amplitude_scale", "value_real_00", &
        "value_imag_00", "value_real_10", "value_imag_10", &
        "value_real_01", "value_imag_01", "value_real_11", "value_imag_11"], &
        interpolation_roots, [character(len=64) :: "unit_R", "unit_Z", &
        "weight_00", "weight_10", "weight_01", "weight_11", &
        "value_real", "value_imag", &
        "partition_of_unity_residual"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_profile_endpoint_symbolic.f90", &
        "neort_profile_endpoint_symbolic", &
        "evaluate_neort_profile_endpoints", [character(len=64) :: "s0", "s1", "f0", "f1", &
        "rho"], endpoint_roots, [character(len=64) :: "endpoint_value_s0", "endpoint_value_s1", &
        "endpoint_derivative_s0", "endpoint_derivative_s1", "affine_slope", &
        "affine_intercept", "df_drho", "axis_regular_residual"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_profile_potential_segment_symbolic.f90", &
        "neort_profile_potential_segment_symbolic", &
        "evaluate_neort_profile_potential_segment", [character(len=64) :: "psi0", "psi1", &
        "Omega_E0", "Omega_E1", "Omega_E_constant", "c_light"], &
        profile_potential_roots, &
        [character(len=64) :: "delta_Phi", "delta_Phi_reversed", "delta_Phi_constant_limit"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_profile_potential_map_symbolic.f90", &
        "neort_profile_potential_map_symbolic", &
        "evaluate_neort_profile_potential_map", &
        [character(len=64) :: "psi", "psi0", "psi1", "Phi0", &
        "Omega_E0", "Omega_E1", "c_light"], &
        profile_potential_map_roots, &
        [character(len=64) :: "Phi", "dPhi_dpsi", "d2Phi_dpsi2"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_profile_potential_map_interval_symbolic.f90", &
        "neort_profile_potential_map_interval_symbolic", &
        "evaluate_neort_profile_potential_map_interval", &
        [character(len=64) :: "psi", "psi0", "psi1", "Phi0", &
        "Omega_E0", "Omega_E1", "c_light"], &
        profile_potential_map_roots, &
        [character(len=64) :: "Phi", "dPhi_dpsi", "d2Phi_dpsi2"], &
        interval_kernel=.true.)
    call emit_kernel_file(trim(output_path)// &
        "/neort_eqdsk_allowed_energy_symbolic.f90", &
        "neort_eqdsk_allowed_energy_symbolic", &
        "evaluate_neort_eqdsk_allowed_energy", &
        [character(len=64) :: "radius", "field_scale", "psi", &
        "psi_R", "psi_Z", &
        "psi_RR", "psi_RZ", "psi_ZZ", "psi_RRR", "psi_RRZ", &
        "psi_RZZ", "psi_ZZZ", "F", "dF_dpsihat", &
        "d2F_dpsihat2", "psi_sep", "dZ_dR", "d2Z_dR2", "h", &
        "J_K", "mass", "charge", "c_light", &
        "electrostatic_potential", "dPhi_dpsi", "d2Phi_dpsi2"], &
        eqdsk_allowed_energy_roots, &
        [character(len=64) :: "field_norm_squared", "psi_physical", &
        "dpsi_physical_dR", "bmod", "dbmod_dR", &
        "d2bmod_dR2", "omega_c", "domega_c_dR", "d2omega_c_dR2", &
        "energy_margin", "denergy_margin_dR", "d2energy_margin_dR2", &
        "v_parallel_squared", "dv_parallel_squared_dR", &
        "d2v_parallel_squared_dR2", "bphi_covariant", &
        "dbphi_covariant_dR"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_eqdsk_allowed_energy_interval_symbolic.f90", &
        "neort_eqdsk_allowed_energy_interval_symbolic", &
        "evaluate_neort_eqdsk_allowed_energy_interval", &
        [character(len=64) :: "radius", "field_scale", "psi", &
        "psi_R", "psi_Z", &
        "psi_RR", "psi_RZ", "psi_ZZ", "psi_RRR", "psi_RRZ", &
        "psi_RZZ", "psi_ZZZ", "F", "dF_dpsihat", &
        "d2F_dpsihat2", "psi_sep", "dZ_dR", "d2Z_dR2", "h", &
        "J_K", "mass", "charge", "c_light", &
        "electrostatic_potential", "dPhi_dpsi", "d2Phi_dpsi2"], &
        eqdsk_allowed_energy_roots, &
        [character(len=64) :: "field_norm_squared", "psi_physical", &
        "dpsi_physical_dR", "bmod", "dbmod_dR", &
        "d2bmod_dR2", "omega_c", "domega_c_dR", "d2omega_c_dR2", &
        "energy_margin", "denergy_margin_dR", "d2energy_margin_dR2", &
        "v_parallel_squared", "dv_parallel_squared_dR", &
        "d2v_parallel_squared_dR2", "bphi_covariant", &
        "dbphi_covariant_dR"], interval_kernel=.true.)
    call emit_kernel_file(trim(output_path)// &
        "/neort_eqdsk_canonical_cut_symbolic.f90", &
        "neort_eqdsk_canonical_cut_symbolic", &
        "evaluate_neort_eqdsk_canonical_cut", &
        [character(len=64) :: "v_parallel_squared", &
        "dv_parallel_squared_dR", "mass", "charge", "c_light", &
        "parallel_sign", "psi_physical", "dpsi_physical_dR", &
        "bphi_covariant", "dbphi_covariant_dR"], &
        eqdsk_canonical_cut_roots, &
        [character(len=64) :: "v_parallel", "dv_parallel_dR", &
        "psi_star", "dpsi_star_dR"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_eqdsk_canonical_cut_interval_symbolic.f90", &
        "neort_eqdsk_canonical_cut_interval_symbolic", &
        "evaluate_neort_eqdsk_canonical_cut_interval", &
        [character(len=64) :: "v_parallel_squared", &
        "dv_parallel_squared_dR", "mass", "charge", "c_light", &
        "parallel_sign", "psi_physical", "dpsi_physical_dR", &
        "bphi_covariant", "dbphi_covariant_dR"], &
        eqdsk_canonical_cut_roots, &
        [character(len=64) :: "v_parallel", "dv_parallel_dR", &
        "psi_star", "dpsi_star_dR"], interval_kernel=.true.)
    call emit_kernel_file(trim(output_path)// &
        "/neort_eqdsk_turning_chart_symbolic.f90", &
        "neort_eqdsk_turning_chart_symbolic", &
        "evaluate_neort_eqdsk_turning_chart", &
        [character(len=64) :: "root_coordinate_y", "allowed_side_sign", &
        "dvparallel_squared_dx_root", "psi_star_flux_root", &
        "dpsi_flux_dx_root", "bphi_covariant_root", &
        "dbphi_covariant_dx_root", "mass", "charge", "c_light", &
        "parallel_sign"], eqdsk_turning_chart_roots, &
        [character(len=64) :: "delta_x", "v_parallel", "psi_star", &
        "dpsi_star_dy", "dpsi_star_dy_root"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_full_fow_canonical_numerator_symbolic.f90", &
        "neort_full_fow_canonical_numerator_symbolic", &
        "evaluate_neort_full_fow_canonical_numerator", &
        [character(len=64) :: "F", "F_R", "F_RR", "F_RRR", "G", "G_R", &
        "G_RR", "G_RRR", "P", "P_R", "P_RR", "mass", "c_light", &
        "sigma", "charge"], canonical_numerator_roots, &
        [character(len=64) :: "N", "N_R", "N_RR"], line_limit=88)
    call emit_kernel_file(trim(output_path)// &
        "/neort_full_fow_canonical_numerator_interval_symbolic.f90", &
        "neort_full_fow_canonical_numerator_interval_symbolic", &
        "evaluate_neort_full_fow_canonical_numerator_interval", &
        [character(len=64) :: "F", "F_R", "F_RR", "F_RRR", "G", "G_R", &
        "G_RR", "G_RRR", "P", "P_R", "P_RR", "mass", "c_light", &
        "sigma", "charge"], canonical_numerator_roots, &
        [character(len=64) :: "N", "N_R", "N_RR"], &
        interval_kernel=.true., line_limit=88)
    call emit_kernel_file(trim(output_path)// &
        "/neort_full_fow_canonical_turning_symbolic.f90", &
        "neort_full_fow_canonical_turning_symbolic", &
        "evaluate_neort_full_fow_canonical_turning", &
        [character(len=64) :: "allowed_side_sign", &
        "dvparallel_squared_dx_root", "bphi_covariant_root", &
        "dbphi_covariant_dx_root", "dpsi_flux_dx_root", "mass", &
        "c_light", "sigma", "charge"], &
        canonical_turning_roots, [character(len=64) :: "N_turning", &
        "dN_dy_turning"], line_limit=88)
    call emit_kernel_file(trim(output_path)// &
        "/neort_full_fow_canonical_turning_interval_symbolic.f90", &
        "neort_full_fow_canonical_turning_interval_symbolic", &
        "evaluate_neort_full_fow_canonical_turning_interval", &
        [character(len=64) :: "allowed_side_sign", &
        "dvparallel_squared_dx_root", "bphi_covariant_root", &
        "dbphi_covariant_dx_root", "dpsi_flux_dx_root", "mass", &
        "c_light", "sigma", "charge"], &
        canonical_turning_roots, [character(len=64) :: "N_turning", &
        "dN_dy_turning"], interval_kernel=.true., line_limit=88)
    call emit_kernel_file(trim(output_path)// &
        "/neort_full_fow_canonical_symmetry_symbolic.f90", &
        "neort_full_fow_canonical_symmetry_symbolic", &
        "evaluate_neort_full_fow_canonical_symmetry", &
        [character(len=64) :: "F", "F_R", "F_RR", "F_RRR", "G", "G_R", &
        "G_RR", "G_RRR", "P", "P_R", "P_RR", "mass", "c_light", &
        "sigma", "charge", "m_mode", "n_mode", "omega_b", "omega_phi", &
        "orbit_time", "phi", "psi", "p_parallel", "theta", "a_real", &
        "a_imag", "endpoint_delta_01", "endpoint_delta_12", &
        "sign_torque_phi", "sign_rotation_phi"], &
        canonical_symmetry_roots, [character(len=64) :: "N_coordinate", &
        "N_R_coordinate", "N_RR_coordinate", "time_phase_reversed", &
        "time_resonance_reversed", "fourier_pair_resonance", &
        "combined_resonance", "toroidal_phase_relabelled", &
        "P_phi_relabelled", "sigma_reversal_difference", &
        "time_phase_residual", "time_resonance_residual", &
        "fourier_pair_phase_residual", "fourier_pair_resonance_residual", &
        "combined_resonance_residual", "endpoint_variation_residual", &
        "fourier_pair_field_residual", "toroidal_phase_residual", &
        "toroidal_P_phi_residual", "toroidal_power_residual", &
        "torque_component_difference"], line_limit=88)
    call emit_kernel_file(trim(output_path)// &
        "/neort_eqdsk_flux_profile_segment_symbolic.f90", &
        "neort_eqdsk_flux_profile_segment_symbolic", &
        "evaluate_neort_eqdsk_flux_profile_segment", &
        [character(len=64) :: "s_tor", "target_psihat", "s_tor0", "s_tor1", &
        "scaled_psi0", "scaled_psi1", "field_scale", "psi_sep"], &
        eqdsk_flux_profile_segment_roots, &
        [character(len=64) :: "psihat", "dpsihat_dstor", &
        "s_tor_from_psihat"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_eqdsk_flux_profile_segment_interval_symbolic.f90", &
        "neort_eqdsk_flux_profile_segment_interval_symbolic", &
        "evaluate_neort_eqdsk_flux_profile_segment_interval", &
        [character(len=64) :: "s_tor", "target_psihat", "s_tor0", &
        "s_tor1", "scaled_psi0", "scaled_psi1", "field_scale", &
        "psi_sep"], eqdsk_flux_profile_segment_roots, &
        [character(len=64) :: "psihat", "dpsihat_dstor", &
        "s_tor_from_psihat"], interval_kernel=.true.)
    call emit_kernel_file(trim(output_path)// &
        "/neort_eqdsk_scaled_flux_normalization_symbolic.f90", &
        "neort_eqdsk_scaled_flux_normalization_symbolic", &
        "evaluate_neort_eqdsk_scaled_flux_normalization", &
        [character(len=64) :: "scaled_psi", "field_scale", "psi_sep"], &
        eqdsk_scaled_flux_normalization_roots, &
        [character(len=64) :: "psihat"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_eqdsk_physical_flux_normalization_symbolic.f90", &
        "neort_eqdsk_physical_flux_normalization_symbolic", &
        "evaluate_neort_eqdsk_physical_flux_normalization", &
        [character(len=64) :: "psi_physical", "field_scale", "psi_sep"], &
        eqdsk_physical_flux_normalization_roots, &
        [character(len=64) :: "psihat"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_eqdsk_physical_flux_normalization_interval_symbolic.f90", &
        "neort_eqdsk_physical_flux_normalization_interval_symbolic", &
        "evaluate_neort_eqdsk_physical_flux_normalization_interval", &
        [character(len=64) :: "psi_physical", "field_scale", "psi_sep"], &
        eqdsk_physical_flux_normalization_roots, &
        [character(len=64) :: "psihat"], interval_kernel=.true.)
    call emit_kernel_file(trim(output_path)// &
        "/neort_eqdsk_physical_flux_map_symbolic.f90", &
        "neort_eqdsk_physical_flux_map_symbolic", &
        "evaluate_neort_eqdsk_physical_flux_map", &
        [character(len=64) :: "psi", "psi_R", "psi_Z", "dZ_dR", &
        "field_scale"], &
        eqdsk_physical_flux_map_roots, &
        [character(len=64) :: "psi_physical", "dpsi_physical_dR"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_eqdsk_physical_flux_map_interval_symbolic.f90", &
        "neort_eqdsk_physical_flux_map_interval_symbolic", &
        "evaluate_neort_eqdsk_physical_flux_map_interval", &
        [character(len=64) :: "psi", "psi_R", "psi_Z", "dZ_dR", &
        "field_scale"], &
        eqdsk_physical_flux_map_roots, &
        [character(len=64) :: "psi_physical", "dpsi_physical_dR"], &
        interval_kernel=.true.)
    call emit_kernel_file(trim(output_path)// &
        "/neort_eqdsk_quintic_cell_jet_symbolic.f90", &
        "neort_eqdsk_quintic_cell_jet_symbolic", &
        "evaluate_neort_eqdsk_quintic_cell_jet", eqdsk_cell_arg_names, &
        eqdsk_cell_jet_roots, [character(len=64) :: "psi", "psi_R", &
        "psi_Z", "psi_RR", "psi_RZ", "psi_ZZ", "psi_RRR", &
        "psi_RRZ", "psi_RZZ", "psi_ZZZ"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_eqdsk_quintic_cell_jet_interval_symbolic.f90", &
        "neort_eqdsk_quintic_cell_jet_interval_symbolic", &
        "evaluate_neort_eqdsk_quintic_cell_jet_interval", &
        eqdsk_cell_arg_names, eqdsk_cell_jet_roots, &
        [character(len=64) :: "psi", "psi_R", "psi_Z", "psi_RR", &
        "psi_RZ", "psi_ZZ", "psi_RRR", "psi_RRZ", "psi_RZZ", &
        "psi_ZZZ"], interval_kernel=.true.)
    call emit_kernel_file(trim(output_path)// &
        "/neort_eqdsk_quintic_cell_fourth_jet_symbolic.f90", &
        "neort_eqdsk_quintic_cell_fourth_jet_symbolic", &
        "evaluate_neort_eqdsk_quintic_cell_fourth_jet", &
        eqdsk_cell_arg_names, eqdsk_cell_fourth_jet_roots, &
        [character(len=64) :: "psi_RRRR", "psi_RRRZ", "psi_RRZZ", &
        "psi_RZZZ", "psi_ZZZZ"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_eqdsk_quintic_cell_fourth_jet_interval_symbolic.f90", &
        "neort_eqdsk_quintic_cell_fourth_jet_interval_symbolic", &
        "evaluate_neort_eqdsk_quintic_cell_fourth_jet_interval", &
        eqdsk_cell_arg_names, eqdsk_cell_fourth_jet_roots, &
        [character(len=64) :: "psi_RRRR", "psi_RRRZ", "psi_RRZZ", &
        "psi_RZZZ", "psi_ZZZZ"], interval_kernel=.true.)
    call emit_kernel_file(trim(output_path)// &
        "/neort_eqdsk_quintic_profile_jet_symbolic.f90", &
        "neort_eqdsk_quintic_profile_jet_symbolic", &
        "evaluate_neort_eqdsk_quintic_profile_jet", &
        eqdsk_profile_arg_names, eqdsk_profile_jet_roots, &
        [character(len=64) :: "profile_value", "profile_first", &
        "profile_second", "vacuum_F"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_eqdsk_quintic_profile_jet_interval_symbolic.f90", &
        "neort_eqdsk_quintic_profile_jet_interval_symbolic", &
        "evaluate_neort_eqdsk_quintic_profile_jet_interval", &
        eqdsk_profile_arg_names, eqdsk_profile_jet_roots, &
        [character(len=64) :: "profile_value", "profile_first", &
        "profile_second", "vacuum_F"], interval_kernel=.true.)
    call emit_kernel_file(trim(output_path)// &
        "/neort_eqdsk_cut_jet_symbolic.f90", &
        "neort_eqdsk_cut_jet_symbolic", &
        "evaluate_neort_eqdsk_cut_jet", &
        [character(len=64) :: "radius", "field_scale", "cut_orientation", &
        "psi", "psi_R", "psi_Z", "psi_RR", "psi_RZ", "psi_ZZ", &
        "psi_RRR", "psi_RRZ", "psi_RZZ", "psi_ZZZ", "F", &
        "dF_dpsihat", "psi_sep", "dot_R", &
        "dot_phi", "dot_Z"], eqdsk_cut_jet_roots, &
        [character(len=64) :: "cut_C", "d_cut_C_d_R", &
        "d_cut_C_d_arc_phi", "d_cut_C_d_Z", "cut_Cdot", &
        "absolute_cut_Cdot", "orientation_scalar"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_eqdsk_cut_numerator_symbolic.f90", &
        "neort_eqdsk_cut_numerator_symbolic", &
        "evaluate_neort_eqdsk_cut_numerator", &
        [character(len=64) :: "radius", "psi_R", "psi_Z", &
        "psi_RR", "psi_RZ", "psi_ZZ", "psi_RRR", "psi_RRZ", &
        "psi_RZZ", "psi_ZZZ", "F", "dF_dpsihat", "psi_sep"], &
        eqdsk_cut_numerator_roots, [character(len=64) :: "N", "N_R", "N_Z"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_eqdsk_cut_numerator_interval_symbolic.f90", &
        "neort_eqdsk_cut_numerator_interval_symbolic", &
        "evaluate_neort_eqdsk_cut_numerator_interval", &
        [character(len=64) :: "radius", "psi_R", "psi_Z", &
        "psi_RR", "psi_RZ", "psi_ZZ", "psi_RRR", "psi_RRZ", &
        "psi_RZZ", "psi_ZZZ", "F", "dF_dpsihat", "psi_sep"], &
        eqdsk_cut_interval_roots, [character(len=64) :: "N", "N_R", &
        "N_Z", "Q"], interval_kernel=.true.)
    call emit_kernel_file(trim(output_path)// &
        "/neort_eqdsk_cut_numerator_hessian_symbolic.f90", &
        "neort_eqdsk_cut_numerator_hessian_symbolic", &
        "evaluate_neort_eqdsk_cut_numerator_hessian", &
        [character(len=64) :: "radius", "psi_R", "psi_Z", "psi_RR", &
        "psi_RZ", "psi_ZZ", "psi_RRR", "psi_RRZ", "psi_RZZ", &
        "psi_ZZZ", "psi_RRRR", "psi_RRRZ", "psi_RRZZ", "psi_RZZZ", &
        "psi_ZZZZ", "F", "dF_dpsihat", "d2F_dpsihat2", "psi_sep"], &
        eqdsk_cut_numerator_hessian_roots, &
        [character(len=64) :: "N_RR", "N_RZ", "N_ZZ"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_eqdsk_cut_numerator_hessian_interval_symbolic.f90", &
        "neort_eqdsk_cut_numerator_hessian_interval_symbolic", &
        "evaluate_neort_eqdsk_cut_numerator_hessian_interval", &
        [character(len=64) :: "radius", "psi_R", "psi_Z", "psi_RR", &
        "psi_RZ", "psi_ZZ", "psi_RRR", "psi_RRZ", "psi_RZZ", &
        "psi_ZZZ", "psi_RRRR", "psi_RRRZ", "psi_RRZZ", "psi_RZZZ", &
        "psi_ZZZZ", "F", "dF_dpsihat", "d2F_dpsihat2", "psi_sep"], &
        eqdsk_cut_numerator_hessian_roots, &
        [character(len=64) :: "N_RR", "N_RZ", "N_ZZ"], &
        interval_kernel=.true.)
    call emit_kernel_file(trim(output_path)// &
        "/neort_eqdsk_cut_r_chart_symbolic.f90", &
        "neort_eqdsk_cut_r_chart_symbolic", &
        "evaluate_neort_eqdsk_cut_r_chart", &
        [character(len=64) :: "N_R", "N_Z"], eqdsk_cut_r_chart_roots, &
        [character(len=64) :: "dZ_dR", "ds_dR"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_eqdsk_cut_z_chart_symbolic.f90", &
        "neort_eqdsk_cut_z_chart_symbolic", &
        "evaluate_neort_eqdsk_cut_z_chart", &
        [character(len=64) :: "N_R", "N_Z"], eqdsk_cut_z_chart_roots, &
        [character(len=64) :: "dR_dZ", "ds_dZ"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_eqdsk_cut_r_flux_chart_symbolic.f90", &
        "neort_eqdsk_cut_r_flux_chart_symbolic", &
        "evaluate_neort_eqdsk_cut_r_flux_chart", &
        [character(len=64) :: "N_R", "N_Z", "psi_R", "psi_Z", &
        "psi_sep"], eqdsk_cut_r_flux_chart_roots, &
        [character(len=64) :: "dZ_dR", "dpsihat_dR"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_eqdsk_cut_r_flux_chart_interval_symbolic.f90", &
        "neort_eqdsk_cut_r_flux_chart_interval_symbolic", &
        "evaluate_neort_eqdsk_cut_r_flux_chart_interval", &
        [character(len=64) :: "N_R", "N_Z", "psi_R", "psi_Z", &
        "psi_sep"], eqdsk_cut_r_flux_chart_roots, &
        [character(len=64) :: "dZ_dR", "dpsihat_dR"], &
        interval_kernel=.true.)
    call emit_kernel_file(trim(output_path)// &
        "/neort_eqdsk_cut_r_flux_curvature_symbolic.f90", &
        "neort_eqdsk_cut_r_flux_curvature_symbolic", &
        "evaluate_neort_eqdsk_cut_r_flux_curvature", &
        [character(len=64) :: "N_R", "N_Z", "N_RR", "N_RZ", "N_ZZ", &
        "psi_R", "psi_Z", "psi_RR", "psi_RZ", "psi_ZZ", "psi_sep"], &
        eqdsk_cut_r_flux_curvature_roots, &
        [character(len=64) :: "dZ_dR", "d2Z_dR2", "dpsihat_dR", &
        "d2psihat_dR2"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_eqdsk_cut_r_flux_curvature_interval_symbolic.f90", &
        "neort_eqdsk_cut_r_flux_curvature_interval_symbolic", &
        "evaluate_neort_eqdsk_cut_r_flux_curvature_interval", &
        [character(len=64) :: "N_R", "N_Z", "N_RR", "N_RZ", "N_ZZ", &
        "psi_R", "psi_Z", "psi_RR", "psi_RZ", "psi_ZZ", "psi_sep"], &
        eqdsk_cut_r_flux_curvature_roots, &
        [character(len=64) :: "dZ_dR", "d2Z_dR2", "dpsihat_dR", &
        "d2psihat_dR2"], interval_kernel=.true.)
    call emit_kernel_file(trim(output_path)// &
        "/neort_eqdsk_cut_mean_value_interval_symbolic.f90", &
        "neort_eqdsk_cut_mean_value_interval_symbolic", &
        "evaluate_neort_eqdsk_cut_mean_value_interval", &
        [character(len=64) :: "N_mid", "N_R_box", "N_Z_box", &
        "delta_R", "delta_Z"], eqdsk_cut_mean_value_roots, &
        [character(len=64) :: "N_mean_value"], interval_kernel=.true.)
    call emit_kernel_file(trim(output_path)// &
        "/neort_eqdsk_cut_axis_curvature_interval_symbolic.f90", &
        "neort_eqdsk_cut_axis_curvature_interval_symbolic", &
        "evaluate_neort_eqdsk_cut_axis_curvature_interval", &
        [character(len=64) :: "dZ_dR", "psi_RR", "psi_RZ", "psi_ZZ", &
        "psi_sep"], eqdsk_cut_axis_curvature_roots, &
        [character(len=64) :: "axis_flux_curvature", &
        "axis_hessian_determinant"], interval_kernel=.true.)
    call emit_kernel_file(trim(output_path)// &
        "/neort_eqdsk_cut_axis_limit_symbolic.f90", &
        "neort_eqdsk_cut_axis_limit_symbolic", &
        "evaluate_neort_eqdsk_cut_axis_limit", &
        [character(len=64) :: "dZ_dR", "psi_RR", "psi_RZ", "psi_ZZ", &
        "psi_sep", "delta_psihat", "branch_sign"], &
        eqdsk_cut_axis_limit_roots, &
        [character(len=64) :: "axis_flux_curvature", &
        "absolute_dpsihat_dR_limit", "delta_R_limit"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_eqdsk_rho_tor_map_symbolic.f90", &
        "neort_eqdsk_rho_tor_map_symbolic", &
        "evaluate_neort_eqdsk_rho_tor_map", &
        [character(len=64) :: "rho_tor"], eqdsk_rho_tor_map_roots, &
        [character(len=64) :: "s_tor", "dstor_drho_tor"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_eqdsk_rho_tor_map_interval_symbolic.f90", &
        "neort_eqdsk_rho_tor_map_interval_symbolic", &
        "evaluate_neort_eqdsk_rho_tor_map_interval", &
        [character(len=64) :: "rho_tor"], eqdsk_rho_tor_map_roots, &
        [character(len=64) :: "s_tor", "dstor_drho_tor"], &
        interval_kernel=.true.)
    call emit_kernel_file(trim(output_path)// &
        "/neort_eqdsk_s_tor_to_rho_symbolic.f90", &
        "neort_eqdsk_s_tor_to_rho_symbolic", &
        "evaluate_neort_eqdsk_s_tor_to_rho", &
        [character(len=64) :: "s_tor"], eqdsk_s_tor_to_rho_roots, &
        [character(len=64) :: "rho_tor"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_eqdsk_s_tor_to_rho_interval_symbolic.f90", &
        "neort_eqdsk_s_tor_to_rho_interval_symbolic", &
        "evaluate_neort_eqdsk_s_tor_to_rho_interval", &
        [character(len=64) :: "s_tor"], eqdsk_s_tor_to_rho_roots, &
        [character(len=64) :: "rho_tor"], interval_kernel=.true.)
    call emit_kernel_file(trim(output_path)// &
        "/neort_eqdsk_flux_profile_rho_chain_symbolic.f90", &
        "neort_eqdsk_flux_profile_rho_chain_symbolic", &
        "evaluate_neort_eqdsk_flux_profile_rho_chain", &
        [character(len=64) :: "dstor_drho_tor", "dpsihat_dstor"], &
        eqdsk_flux_profile_rho_chain_roots, &
        [character(len=64) :: "dpsihat_drho_tor"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_eqdsk_cut_flux_coordinate_symbolic.f90", &
        "neort_eqdsk_cut_flux_coordinate_symbolic", &
        "evaluate_neort_eqdsk_cut_flux_coordinate", &
        [character(len=64) :: "dstor_drho_tor", "dpsihat_dstor", &
        "dpsihat_dR", "dZ_dR"], eqdsk_cut_flux_coordinate_roots, &
        [character(len=64) :: "dpsihat_drho_tor", "dR_drho_tor", &
        "dZ_drho_tor"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_eqdsk_cut_flux_coordinate_interval_symbolic.f90", &
        "neort_eqdsk_cut_flux_coordinate_interval_symbolic", &
        "evaluate_neort_eqdsk_cut_flux_coordinate_interval", &
        [character(len=64) :: "dstor_drho_tor", "dpsihat_dstor", &
        "dpsihat_dR", "dZ_dR"], eqdsk_cut_flux_coordinate_roots, &
        [character(len=64) :: "dpsihat_drho_tor", "dR_drho_tor", &
        "dZ_drho_tor"], interval_kernel=.true., line_limit=88)
    call emit_kernel_file(trim(output_path)// &
        "/neort_eqdsk_cut_axis_rho_limit_symbolic.f90", &
        "neort_eqdsk_cut_axis_rho_limit_symbolic", &
        "evaluate_neort_eqdsk_cut_axis_rho_limit", &
        [character(len=64) :: "dZ_dR", "psi_RR", "psi_RZ", "psi_ZZ", &
        "psi_sep", "dpsihat_dstor", "branch_sign"], &
        eqdsk_cut_axis_rho_limit_roots, &
        [character(len=64) :: "axis_flux_curvature", &
        "dR_drho_tor_limit", "dZ_drho_tor_limit"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_eqdsk_cut_axis_rho_limit_interval_symbolic.f90", &
        "neort_eqdsk_cut_axis_rho_limit_interval_symbolic", &
        "evaluate_neort_eqdsk_cut_axis_rho_limit_interval", &
        [character(len=64) :: "dZ_dR", "psi_RR", "psi_RZ", "psi_ZZ", &
        "psi_sep", "dpsihat_dstor", "branch_sign"], &
        eqdsk_cut_axis_rho_limit_roots, &
        [character(len=64) :: "axis_flux_curvature", &
        "dR_drho_tor_limit", "dZ_drho_tor_limit"], &
        interval_kernel=.true.)
    call emit_kernel_file(trim(output_path)// &
        "/neort_eqdsk_allowed_rho_chain_symbolic.f90", &
        "neort_eqdsk_allowed_rho_chain_symbolic", &
        "evaluate_neort_eqdsk_allowed_rho_chain", &
        [character(len=64) :: "dstor_drho_tor", "dpsihat_dstor", &
        "d2psihat_dstor2", "dpsihat_dR", &
        "d2psihat_dR2", "dvparallel_squared_dR", &
        "d2vparallel_squared_dR2", "dpsi_star_dR"], &
        eqdsk_allowed_rho_chain_roots, &
        [character(len=64) :: "dpsihat_drho_tor", &
        "d2psihat_drho_tor2", "dR_drho_tor", "d2R_drho_tor2", &
        "dvparallel_squared_drho_tor", &
        "d2vparallel_squared_drho_tor2", "dpsi_star_drho_tor"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_eqdsk_allowed_rho_chain_interval_symbolic.f90", &
        "neort_eqdsk_allowed_rho_chain_interval_symbolic", &
        "evaluate_neort_eqdsk_allowed_rho_chain_interval", &
        [character(len=64) :: "dstor_drho_tor", "dpsihat_dstor", &
        "d2psihat_dstor2", "dpsihat_dR", &
        "d2psihat_dR2", "dvparallel_squared_dR", &
        "d2vparallel_squared_dR2", "dpsi_star_dR"], &
        eqdsk_allowed_rho_chain_roots, &
        [character(len=64) :: "dpsihat_drho_tor", &
        "d2psihat_drho_tor2", "dR_drho_tor", "d2R_drho_tor2", &
        "dvparallel_squared_drho_tor", &
        "d2vparallel_squared_drho_tor2", "dpsi_star_drho_tor"], &
        interval_kernel=.true.)
    call emit_kernel_file(trim(output_path)// &
        "/neort_eqdsk_allowed_axis_rho_chain_symbolic.f90", &
        "neort_eqdsk_allowed_axis_rho_chain_symbolic", &
        "evaluate_neort_eqdsk_allowed_axis_rho_chain", &
        [character(len=64) :: "dR_drho_tor_limit", &
        "dvparallel_squared_dR", "dpsi_star_dR"], &
        eqdsk_allowed_axis_rho_chain_roots, &
        [character(len=64) :: "dvparallel_squared_drho_tor_axis", &
        "dpsi_star_drho_tor_axis"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_eqdsk_allowed_axis_rho_chain_interval_symbolic.f90", &
        "neort_eqdsk_allowed_axis_rho_chain_interval_symbolic", &
        "evaluate_neort_eqdsk_allowed_axis_rho_chain_interval", &
        [character(len=64) :: "dR_drho_tor_limit", &
        "dvparallel_squared_dR", "dpsi_star_dR"], &
        eqdsk_allowed_axis_rho_chain_roots, &
        [character(len=64) :: "dvparallel_squared_drho_tor_axis", &
        "dpsi_star_drho_tor_axis"], interval_kernel=.true.)
    call emit_kernel_file(trim(output_path)// &
        "/neort_eqdsk_cut_endpoint_system_symbolic.f90", &
        "neort_eqdsk_cut_endpoint_system_symbolic", &
        "evaluate_neort_eqdsk_cut_endpoint_system", &
        [character(len=64) :: "cut_numerator", "d_cut_numerator_d_R", &
        "d_cut_numerator_d_Z", "psi", "psi_R", "psi_Z", "psi_sep", &
        "target_psihat"], eqdsk_cut_endpoint_system_roots, &
        [character(len=64) :: "cut_residual", "flux_residual", &
        "cut_jacobian_R", "cut_jacobian_Z", "flux_jacobian_R", &
        "flux_jacobian_Z"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_eqdsk_cut_endpoint_system_interval_symbolic.f90", &
        "neort_eqdsk_cut_endpoint_system_interval_symbolic", &
        "evaluate_neort_eqdsk_cut_endpoint_system_interval", &
        [character(len=64) :: "cut_numerator", "d_cut_numerator_d_R", &
        "d_cut_numerator_d_Z", "psi", "psi_R", "psi_Z", "psi_sep", &
        "target_psihat"], eqdsk_cut_endpoint_system_roots, &
        [character(len=64) :: "cut_residual", "flux_residual", &
        "cut_jacobian_R", "cut_jacobian_Z", "flux_jacobian_R", &
        "flux_jacobian_Z"], interval_kernel=.true.)
    call emit_kernel_file(trim(output_path)// &
        "/neort_eqdsk_cut_endpoint_newton_symbolic.f90", &
        "neort_eqdsk_cut_endpoint_newton_symbolic", &
        "evaluate_neort_eqdsk_cut_endpoint_newton", &
        [character(len=64) :: "R0", "Z0", "cut_numerator0", &
        "flux_residual0", "d_cut_numerator_d_R0", &
        "d_cut_numerator_d_Z0", "d_flux_residual_d_R0", &
        "d_flux_residual_d_Z0"], eqdsk_cut_endpoint_newton_roots, &
        [character(len=64) :: "jacobian_determinant", "newton_R", &
        "newton_Z", "inverse_11", "inverse_12", "inverse_21", &
        "inverse_22"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_eqdsk_cut_endpoint_krawczyk_interval_symbolic.f90", &
        "neort_eqdsk_cut_endpoint_krawczyk_interval_symbolic", &
        "evaluate_neort_eqdsk_cut_endpoint_krawczyk_interval", &
        [character(len=64) :: "R0", "Z0", "cut_numerator0", &
        "flux_residual0", "d_cut_numerator_d_R0", &
        "d_cut_numerator_d_Z0", "d_flux_residual_d_R0", &
        "d_flux_residual_d_Z0", "d_cut_numerator_d_R_box", &
        "d_cut_numerator_d_Z_box", "d_flux_residual_d_R_box", &
        "d_flux_residual_d_Z_box", "delta_R_box", "delta_Z_box"], &
        eqdsk_cut_endpoint_krawczyk_roots, &
        [character(len=64) :: "krawczyk_R", "krawczyk_Z"], &
        interval_kernel=.true.)
    call emit_kernel_file(trim(output_path)// &
        "/neort_eqdsk_axis_stationarity_system_symbolic.f90", &
        "neort_eqdsk_axis_stationarity_system_symbolic", &
        "evaluate_neort_eqdsk_axis_stationarity_system", &
        [character(len=64) :: "psi_R", "psi_Z", "psi_RR", "psi_RZ", &
        "psi_ZZ"], eqdsk_axis_stationarity_system_roots, &
        [character(len=64) :: "residual_R", "residual_Z", "jacobian_11", &
        "jacobian_12", "jacobian_21", "jacobian_22"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_eqdsk_axis_stationarity_system_interval_symbolic.f90", &
        "neort_eqdsk_axis_stationarity_system_interval_symbolic", &
        "evaluate_neort_eqdsk_axis_stationarity_system_interval", &
        [character(len=64) :: "psi_R", "psi_Z", "psi_RR", "psi_RZ", &
        "psi_ZZ"], eqdsk_axis_stationarity_system_roots, &
        [character(len=64) :: "residual_R", "residual_Z", "jacobian_11", &
        "jacobian_12", "jacobian_21", "jacobian_22"], &
        interval_kernel=.true.)
    call emit_kernel_file(trim(output_path)// &
        "/neort_eqdsk_axis_stationarity_newton_symbolic.f90", &
        "neort_eqdsk_axis_stationarity_newton_symbolic", &
        "evaluate_neort_eqdsk_axis_stationarity_newton", &
        [character(len=64) :: "R0", "Z0", "psi_R0", "psi_Z0", &
        "psi_RR0", "psi_RZ0", "psi_ZZ0"], &
        eqdsk_axis_stationarity_newton_roots, &
        [character(len=64) :: "hessian_determinant", "newton_R", &
        "newton_Z", "inverse_11", "inverse_12", "inverse_21", &
        "inverse_22"])
    call emit_kernel_file(trim(output_path)// &
        "/neort_eqdsk_axis_stationarity_krawczyk_interval_symbolic.f90", &
        "neort_eqdsk_axis_stationarity_krawczyk_interval_symbolic", &
        "evaluate_neort_eqdsk_axis_stationarity_krawczyk_interval", &
        [character(len=64) :: "R0", "Z0", "psi_R0", "psi_Z0", &
        "psi_RR0", "psi_RZ0", "psi_ZZ0", "psi_RR_box", "psi_RZ_box", &
        "psi_ZZ_box", "delta_R_box", "delta_Z_box"], &
        eqdsk_axis_stationarity_krawczyk_roots, &
        [character(len=64) :: "krawczyk_R", "krawczyk_Z"], &
        interval_kernel=.true.)
    call emit_certificate_registry(trim(output_path)// &
        "/neort_generated_certificate_registry.f90")

contains

    function exact_simplify(expression, context) result(value)
        type(expr_t), intent(in) :: expression
        character(*), intent(in) :: context
        type(expr_t) :: value
        type(engine_result_t) :: result

        result = simplify_engine%simplify(expression)
        if (.not. result%ok) then
            write (output_unit, "(a,1x,a)") &
                "fortsym exact simplification failed:", trim(context)
            error stop 1
        end if
        value = result%value
    end function exact_simplify

    function exact_derivative(expression, variable, context) result(value)
        type(expr_t), intent(in) :: expression, variable
        character(*), intent(in) :: context
        type(expr_t) :: value
        type(engine_result_t) :: result

        result = simplify_engine%diff(expression, variable)
        if (.not. result%ok) then
            write (output_unit, "(a,1x,a)") &
                "fortsym exact differentiation failed:", trim(context)
            error stop 1
        end if
        value = result%value
    end function exact_derivative

    function exact_expand(expression, context) result(value)
        type(expr_t), intent(in) :: expression
        character(*), intent(in) :: context
        type(expr_t) :: value
        type(engine_result_t) :: result

        result = simplify_engine%expand(expression)
        if (.not. result%ok) then
            write (output_unit, "(a,1x,a)") &
                "fortsym exact expansion failed:", trim(context)
            error stop 1
        end if
        value = result%value
    end function exact_expand

    subroutine build_cylindrical_geometry(b, d_b, radius_in, bmod_out, &
            bhat_out, grad_b_out, dbhat_out, curl_out)
        type(expr_t), intent(in) :: b(3), d_b(3,3), radius_in
        type(expr_t), intent(out) :: bmod_out, bhat_out(3), grad_b_out(3)
        type(expr_t), intent(out) :: dbhat_out(3,3), curl_out(3)
        integer :: component, direction

        bmod_out = sqrt(b(1)**2 + b(2)**2 + b(3)**2)
        do component = 1, 3
            bhat_out(component) = b(component)/bmod_out
        end do
        do direction = 1, 3
            grad_b_out(direction) = bhat_out(1)*d_b(1,direction) + &
                bhat_out(2)*d_b(2,direction) + &
                bhat_out(3)*d_b(3,direction)
            do component = 1, 3
                dbhat_out(component,direction) = &
                    (d_b(component,direction) - &
                    bhat_out(component)*grad_b_out(direction))/bmod_out
            end do
        end do
        curl_out(1) = dbhat_out(3,2) - dbhat_out(2,3)
        curl_out(2) = dbhat_out(1,3) - dbhat_out(3,1)
        curl_out(3) = dbhat_out(2,1) + bhat_out(2)/radius_in - &
            dbhat_out(1,2)
    end subroutine build_cylindrical_geometry

    subroutine simplify_array(expressions)
        type(expr_t), intent(inout) :: expressions(:)
        integer :: i

        do i = 1, size(expressions)
            simplified = simplify_engine%simplify(expressions(i))
            if (.not. simplified%ok) then
                error stop "fortsym simplification failed"
            end if
            expressions(i) = simplified%value
        end do
    end subroutine simplify_array

    subroutine emit_kernel_file(path, module_name, procedure_name, arg_names, &
            kernel_roots, output_names, interval_kernel, line_limit)
        character(*), intent(in) :: path, module_name, procedure_name
        character(*), intent(in) :: arg_names(:), output_names(:)
        type(expr_t), intent(in) :: kernel_roots(:)
        logical, intent(in), optional :: interval_kernel
        integer, intent(in), optional :: line_limit
        type(kernel_spec_t) :: kernel_spec
        type(str_t) :: emitted
        character(:), allocatable :: emitted_text
        integer :: unit, ios, i, j, maximum_allowed_length
        logical :: ok, make_interval

        if (size(kernel_roots) /= size(output_names)) then
            error stop "fortsym generator arity mismatch"
        end if
        if (size(arg_names) == 0 .or. size(output_names) == 0) then
            error stop "fortsym generator emitted empty interface"
        end if
        do i = 1, size(arg_names)
            if (len_trim(arg_names(i)) == 0) error stop "empty generator argument"
            do j = 1, i - 1
                if (trim(arg_names(i)) == trim(arg_names(j))) then
                    error stop "duplicate generator argument name"
                end if
            end do
        end do
        do i = 1, size(output_names)
            if (len_trim(output_names(i)) == 0) error stop "empty generator output"
            do j = 1, i - 1
                if (trim(output_names(i)) == trim(output_names(j))) then
                    error stop "duplicate generator output name"
                end if
            end do
        end do

        kernel_spec%name = str(procedure_name)
        kernel_spec%module_name = str(module_name)
        kernel_spec%mode = KERNEL_SUBROUTINE
        kernel_spec%pure_procedure = .true.
        kernel_spec%generator = str(&
            "tools/gc_symbolics/app/gen_full_fow_physics.f90")
        kernel_spec%generator_revision = str(FORTSYM_REVISION)
        kernel_spec%regenerate_command = str(&
            "cd tools/gc_symbolics && fo exec gen_full_fow_physics "// &
            "../../src/generated")
        allocate (kernel_spec%args(size(arg_names)), &
            kernel_spec%outputs(size(output_names)))
        do i = 1, size(arg_names)
            kernel_spec%args(i) = str(trim(arg_names(i)))
        end do
        do i = 1, size(output_names)
            kernel_spec%outputs(i) = str(trim(output_names(i)))
        end do

        open (newunit=unit, file=trim(path), status="replace", action="write", &
            iostat=ios)
        if (ios /= 0) error stop "cannot open generated kernel output"
        emitted = emit_kernel(kernel_roots, kernel_spec, ok)
        if (.not. ok .or. len(chars(emitted)) == 0) then
            close (unit)
            error stop "fortsym refused generated kernel"
        end if
        emitted_text = chars(emitted)
        make_interval = .false.
        if (present(interval_kernel)) make_interval = interval_kernel
        if (make_interval) emitted_text = intervalize_kernel(emitted_text)
        if (present(line_limit)) emitted_text = &
            wrap_generated_lines(emitted_text, line_limit)
        maximum_allowed_length = 132
        if (present(line_limit)) maximum_allowed_length = line_limit
        if (maximum_line_length(emitted_text) > maximum_allowed_length) then
            close (unit)
            error stop "fortsym emitted a nonconforming overlong line"
        end if
        write (unit, "(a)") emitted_text
        close (unit)
        write (output_unit, "(a)") "wrote "//trim(path)
    end subroutine emit_kernel_file

    function intervalize_kernel(source) result(interval_source)
        character(*), intent(in) :: source
        character(:), allocatable :: interval_source
        character(*), parameter :: interval_use = &
            "        use neort_gc_outward_interval, only: "// &
            "gc_outward_interval_t, &"//new_line('a')// &
            "            operator(+), operator(-), operator(*), "// &
            "operator(/), operator(**), abs, sqrt"//new_line('a')// &
            "        implicit none"

        interval_source = replace_all(source, "        implicit none", &
            interval_use)
        interval_source = replace_all(interval_source, "real(dp)", &
            "type(gc_outward_interval_t)")
    end function intervalize_kernel

    function replace_all(source, target, replacement) result(replaced)
        character(*), intent(in) :: source, target, replacement
        character(:), allocatable :: replaced
        integer :: first, relative

        if (len(target) == 0) error stop "empty generated-text replacement"
        replaced = ''
        first = 1
        do
            relative = index(source(first:), target)
            if (relative == 0) then
                replaced = replaced//source(first:)
                exit
            end if
            if (relative > 1) then
                replaced = replaced//source(first:first+relative-2)
            end if
            replaced = replaced//replacement
            first = first+relative-1+len(target)
            if (first > len(source)) exit
        end do
    end function replace_all

    function wrap_generated_lines(source, limit) result(wrapped)
        character(*), intent(in) :: source
        integer, intent(in) :: limit
        character(:), allocatable :: wrapped
        character(:), allocatable :: line
        integer :: start, finish

        wrapped = ''
        start = 1
        do while (start <= len(source))
            finish = index(source(start:), new_line('a'))
            if (finish == 0) then
                line = source(start:)
                start = len(source)+1
            else
                if (finish == 1) then
                    line = ''
                else
                    line = source(start:start+finish-2)
                end if
                start = start+finish
            end if
            wrapped = wrapped//wrap_generated_line(line, limit)// &
                new_line('a')
        end do
    end function wrap_generated_lines

    function wrap_generated_line(line, limit) result(wrapped)
        character(*), intent(in) :: line
        integer, intent(in) :: limit
        character(:), allocatable :: wrapped, rest
        character(:), allocatable :: prefix
        integer :: cut, available
        logical :: comment_line

        wrapped = ''
        rest = line
        comment_line = .false.
        prefix = adjustl(rest)
        if (len_trim(prefix) > 0) comment_line = prefix(1:1) == '!'
        do while (len(rest) > limit)
            available = limit-2
            cut = generated_line_break(rest, available)
            if (cut < 1) cut = available
            if (comment_line) then
                wrapped = wrapped//rest(:cut)//new_line('a')// &
                    '! '
            else
                wrapped = wrapped//rest(:cut)//' &'//new_line('a')// &
                    '            '
            end if
            rest = adjustl(rest(cut+1:))
        end do
        wrapped = wrapped//rest
    end function wrap_generated_line

    pure integer function generated_line_break(line, available) result(cut)
        character(*), intent(in) :: line
        integer, intent(in) :: available
        integer :: i

        cut = 0
        do i = min(available, len(line)), 1, -1
            if (line(i:i) == ' ' .or. line(i:i) == ',' .or. &
                    line(i:i) == '*') then
                cut = i
                exit
            end if
        end do
    end function generated_line_break

    pure integer function maximum_line_length(text) result(maximum_length)
        character(*), intent(in) :: text
        integer :: line_start, i

        maximum_length = 0
        line_start = 1
        do i = 1, len(text)
            if (text(i:i) /= new_line('a')) cycle
            maximum_length = max(maximum_length, i-line_start)
            line_start = i+1
        end do
        maximum_length = max(maximum_length, len(text)-line_start+1)
    end function maximum_line_length

    subroutine emit_certificate_registry(path)
        character(*), intent(in) :: path
        integer :: unit, ios

        open (newunit=unit, file=trim(path), status="replace", action="write", &
            iostat=ios)
        if (ios /= 0) error stop "cannot open certificate registry output"
        write (unit, "(a)") "module neort_generated_certificate_registry"
        write (unit, "(a)") "    implicit none"
        write (unit, "(a)") "    character(*), parameter :: fortsym_revision = &"
        write (unit, "(a)") "        'fortsym@545788453a204d58705f735b519c3863c2f734c8'"
        write (unit, "(a)") "    character(*), parameter :: regenerate_command = &"
        write (unit, "(a)") "        'cd tools/gc_symbolics && fo exec gen_full_fow_physics ../../src/generated'"
        write (unit, "(a)") "    integer, parameter :: certificate_count = 49"
        write (unit, "(a)") "    character(len=32), parameter :: certificate_id(certificate_count) = &"
        write (unit, "(a)") "        [character(len=32) :: 'geometry', 'littlejohn', 'eq13_cdot', 'boundary_limits', &"
        write (unit, "(a)") "        'root_enclosures', 'interpolation', 'profile_endpoints', &"
        write (unit, "(a)") "        'refinement', 'harmonic_integrand', 'simple_root_force', &"
        write (unit, "(a)") "        'eqdsk_cell_jet', 'eqdsk_cell_fourth_jet', &"
        write (unit, "(a)") "        'eqdsk_profile_jet', 'eqdsk_cut_jet', &"
        write (unit, "(a)") "        'eqdsk_cut_numerator_jet', 'eqdsk_cut_r_chart', &"
        write (unit, "(a)") "        'eqdsk_cut_z_chart', 'eqdsk_cut_r_flux_chart', &"
        write (unit, "(a)") "        'eqdsk_cut_numerator_hessian', &"
        write (unit, "(a)") "        'eqdsk_cut_r_flux_curvature', 'eqdsk_cut_mean_value', &"
        write (unit, "(a)") "        'eqdsk_cut_axis_curvature', &"
        write (unit, "(a)") "        'eqdsk_cut_axis_limit', &"
        write (unit, "(a)") "        'eqdsk_rho_tor_map', &"
        write (unit, "(a)") "        'eqdsk_s_tor_to_rho', &"
        write (unit, "(a)") "        'eqdsk_flux_profile_rho_chain', &"
        write (unit, "(a)") "        'eqdsk_cut_flux_coordinate', &"
        write (unit, "(a)") "        'eqdsk_flux_coordinate_interval', &"
        write (unit, "(a)") "        'eqdsk_cut_axis_rho_limit', &"
        write (unit, "(a)") "        'eqdsk_allowed_rho_chain', &"
        write (unit, "(a)") "        'eqdsk_allowed_axis_rho_chain', &"
        write (unit, "(a)") "        'eqdsk_cut_endpoint_system', &"
        write (unit, "(a)") "        'eqdsk_cut_endpoint_newton', &"
        write (unit, "(a)") "        'eqdsk_cut_endpoint_krawczyk', &"
        write (unit, "(a)") "        'eqdsk_axis_stationarity_system', &"
        write (unit, "(a)") "        'eqdsk_axis_stationarity_newton', &"
        write (unit, "(a)") "        'eqdsk_axis_stationarity_krawczyk', &"
        write (unit, "(a)") "        'eqdsk_flux_profile_segment', &"
        write (unit, "(a)") "        'eqdsk_scaled_flux_normalization', &"
        write (unit, "(a)") "        'eqdsk_physical_flux_norm', &"
        write (unit, "(a)") "        'gauss_interval_map', &"
        write (unit, "(a)") "        'profile_potential_map', &"
        write (unit, "(a)") "        'eqdsk_allowed_energy', &"
        write (unit, "(a)") "        'eqdsk_canonical_cut', &"
        write (unit, "(a)") "        'eqdsk_turning_chart', &"
        write (unit, "(a)") "        'eqdsk_physical_flux_map', &"
        write (unit, "(a)") "        'full_fow_canonical_numerator', &"
        write (unit, "(a)") "        'full_fow_canonical_turning', &"
        write (unit, "(a)") "        'full_fow_canonical_symmetry' ]"
        write (unit, "(a)") "    character(len=64), parameter :: certificate_fingerprint(certificate_count) = &"
        write (unit, "(a)") "        [character(len=64) :: 'neort-cert-v1:geometry:19:fortsym-5457884', &"
        write (unit, "(a)") "        'neort-cert-v1:littlejohn:22:fortsym-5457884', &"
        write (unit, "(a)") "        'neort-cert-v1:eq13_cdot:3:fortsym-5457884', &"
        write (unit, "(a)") "        'neort-cert-v1:boundary_limits:13:fortsym-5457884', &"
        write (unit, "(a)") "        'neort-cert-v1:root_enclosures:3:fortsym-5457884', &"
        write (unit, "(a)") "        'neort-cert-v1:interpolation:9:fortsym-5457884', &"
        write (unit, "(a)") "        'neort-cert-v1:profile_endpoints:8:fortsym-5457884', &"
        write (unit, "(a)") "        'neort-cert-v1:refinement:14:fortsym-5457884', &"
        write (unit, "(a)") "        'neort-cert-v1:harmonic_integrand:8:fortsym-5457884', &"
        write (unit, "(a)") "        'neort-cert-v1:simple_root_force:3:fortsym-5457884', &"
        write (unit, "(a)") "        'neort-cert-v1:eqdsk_cell_jet:10:fortsym-5457884', &"
        write (unit, "(a)") "        'neort-cert-v1:eqdsk_cell_fourth_jet:5:fortsym-5457884', &"
        write (unit, "(a)") "        'neort-cert-v1:eqdsk_profile_jet:4:fortsym-5457884', &"
        write (unit, "(a)") "        'neort-cert-v1:eqdsk_cut_jet:7:fortsym-5457884', &"
        write (unit, "(a)") "        'neort-cert-v1:eqdsk_cut_numerator_jet:3:fortsym-5457884', &"
        write (unit, "(a)") "        'neort-cert-v1:eqdsk_cut_r_chart:2:fortsym-5457884', &"
        write (unit, "(a)") "        'neort-cert-v1:eqdsk_cut_z_chart:2:fortsym-5457884', &"
        write (unit, "(a)") "        'neort-cert-v1:eqdsk_cut_r_flux_chart:2:fortsym-5457884', &"
        write (unit, "(a)") "        'neort-cert-v1:eqdsk_cut_numerator_hessian:3:fortsym-5457884', &"
        write (unit, "(a)") "        'neort-cert-v1:eqdsk_cut_r_flux_curvature:4:fortsym-5457884', &"
        write (unit, "(a)") "        'neort-cert-v1:eqdsk_cut_mean_value:1:fortsym-5457884', &"
        write (unit, "(a)") "        'neort-cert-v1:eqdsk_cut_axis_curvature:2:fortsym-5457884', &"
        write (unit, "(a)") "        'neort-cert-v1:eqdsk_cut_axis_limit:3:fortsym-5457884', &"
        write (unit, "(a)") "        'neort-cert-v1:eqdsk_rho_tor_map:2:fortsym-5457884', &"
        write (unit, "(a)") &
            "        'neort-cert-v1:eqdsk_s_tor_to_rho:1:fortsym-5457884', &"
        write (unit, "(a)") "        'neort-cert-v1:eqdsk_flux_profile_rho_chain:1:fortsym-5457884', &"
        write (unit, "(a)") "        'neort-cert-v1:eqdsk_cut_flux_coordinate:3:fortsym-5457884', &"
        write (unit, "(a)") &
            "        'neort-cert-v1:eqdsk_flux_coordinate_interval:3:fortsym-5457884', &"
        write (unit, "(a)") "        'neort-cert-v1:eqdsk_cut_axis_rho_limit:3:fortsym-5457884', &"
        write (unit, "(a)") &
            "        'neort-cert-v1:eqdsk_allowed_rho_chain:7:fortsym-5457884', &"
        write (unit, "(a)") &
            "        'neort-cert-v1:eqdsk_allowed_axis_rho_chain:2:fortsym-5457884', &"
        write (unit, "(a)") "        'neort-cert-v1:eqdsk_cut_endpoint_system:6:fortsym-5457884', &"
        write (unit, "(a)") "        'neort-cert-v1:eqdsk_cut_endpoint_newton:7:fortsym-5457884', &"
        write (unit, "(a)") "        'neort-cert-v1:eqdsk_cut_endpoint_krawczyk:2:fortsym-5457884', &"
        write (unit, "(a)") "        'neort-cert-v1:eqdsk_axis_stationarity_system:6:fortsym-5457884', &"
        write (unit, "(a)") "        'neort-cert-v1:eqdsk_axis_stationarity_newton:7:fortsym-5457884', &"
        write (unit, "(a)") "        'neort-cert-v1:eqdsk_axis_stationarity_krawczyk:2:fortsym-5457884', &"
        write (unit, "(a)") "        'neort-cert-v1:eqdsk_flux_profile_segment:3:fortsym-5457884', &"
        write (unit, "(a)") &
            "        'neort-cert-v1:eqdsk_scaled_flux_normalization:1:fortsym-5457884', &"
        write (unit, "(a)") &
            "        'neort-cert-v1:eqdsk_physical_flux_norm:1:fortsym-5457884', &"
        write (unit, "(a)") &
            "        'neort-cert-v1:gauss_interval_map:2:fortsym-5457884', &"
        write (unit, "(a)") &
            "        'neort-cert-v1:profile_potential_map:3:fortsym-5457884', &"
        write (unit, "(a)") &
            "        'neort-cert-v1:eqdsk_allowed_energy:17:fortsym-5457884', &"
        write (unit, "(a)") &
            "        'neort-cert-v1:eqdsk_canonical_cut:4:fortsym-5457884', &"
        write (unit, "(a)") &
            "        'neort-cert-v1:eqdsk_turning_chart:5:fortsym-5457884', &"
        write (unit, "(a)") &
            "        'neort-cert-v1:eqdsk_physical_flux_map:2:fortsym-5457884', &"
        write (unit, "(a)") &
            "        'neort-cert-v1:full_fow_canonical_numerator:3:fortsym-5457884', &"
        write (unit, "(a)") &
            "        'neort-cert-v1:full_fow_canonical_turning:2:fortsym-5457884', &"
        write (unit, "(a)") &
            "        'neort-cert-v1:full_fow_canonical_symmetry:21:fortsym-5457884' ]"
        write (unit, "(a)") "    ! Fingerprints are provenance/arity manifests, not algebraic proofs."
        write (unit, "(a)") "    ! Root multiplicity and crossing counts require interval/theorem gates."
        write (unit, "(a)") "contains"
        write (unit, "(a)") "    pure integer function certificate_index(id) result(index)"
        write (unit, "(a)") "        character(*), intent(in) :: id"
        write (unit, "(a)") "        integer :: k"
        write (unit, "(a)") "        index = 0"
        write (unit, "(a)") "        do k = 1, certificate_count"
        write (unit, "(a)") "            if (id == certificate_id(k)) index = k"
        write (unit, "(a)") "        end do"
        write (unit, "(a)") "    end function certificate_index"
        write (unit, "(a)") "    pure logical function certificate_matches(id, fingerprint)"
        write (unit, "(a)") "        character(*), intent(in) :: id, fingerprint"
        write (unit, "(a)") "        integer :: k"
        write (unit, "(a)") "        k = certificate_index(id)"
        write (unit, "(a)") "        certificate_matches = .false."
        write (unit, "(a)") "        if (k > 0) then"
        write (unit, "(a)") "            certificate_matches = fingerprint == &"
        write (unit, "(a)") "                certificate_fingerprint(k)"
        write (unit, "(a)") "        end if"
        write (unit, "(a)") "    end function certificate_matches"
        write (unit, "(a)") "end module neort_generated_certificate_registry"
        close (unit)
        write (output_unit, "(a)") "wrote "//trim(path)
    end subroutine emit_certificate_registry

end program gen_full_fow_physics
