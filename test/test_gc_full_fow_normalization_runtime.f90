program test_gc_full_fow_normalization_runtime
    use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_gc_full_fow_normalization_runtime, only: &
        GC_FULL_FOW_NORMALIZATION_CERTIFICATE_REQUIRED, &
        GC_FULL_FOW_NORMALIZATION_DOMAIN_ERROR, &
        GC_FULL_FOW_NORMALIZATION_NONFINITE_INPUT, &
        GC_FULL_FOW_NORMALIZATION_NOT_INITIALIZED, &
        GC_FULL_FOW_NORMALIZATION_REFERENCE_MISMATCH, &
        GC_FULL_FOW_NORMALIZATION_SUCCESS, &
        evaluate_gc_full_fow_canonical_flux, &
        evaluate_gc_full_fow_degree5_cell_enclosure, &
        evaluate_gc_full_fow_eq17, evaluate_gc_full_fow_jk_envelope, &
        evaluate_gc_full_fow_scaled_magnitude, &
        gc_full_fow_degree5_enclosure_t, gc_full_fow_eq17_result_t, &
        gc_full_fow_energy_quadrature_t, gc_full_fow_jk_envelope_t, &
        gc_full_fow_normalized_sample_t, gc_full_fow_physical_sample_t, &
        gc_full_fow_phase_space_bound_certificate_t, &
        gc_full_fow_paired_quadrature_t, gc_full_fow_reference_scales_t, &
        initialize_gc_full_fow_reference, map_gc_full_fow_energy_quadrature, &
        map_gc_full_fow_paired_quadrature, &
        normalize_gc_full_fow_sample
    use util_for_test, only: pass_test
    implicit none

    type(gc_full_fow_reference_scales_t) :: reference, other_reference
    type(gc_full_fow_reference_scales_t) :: uninitialized_reference
    type(gc_full_fow_physical_sample_t) :: physical, bad_physical
    type(gc_full_fow_normalized_sample_t) :: normalized
    type(gc_full_fow_phase_space_bound_certificate_t) :: certificate
    type(gc_full_fow_phase_space_bound_certificate_t) :: bad_certificate
    type(gc_full_fow_phase_space_bound_certificate_t) :: zero_field_certificate
    type(gc_full_fow_jk_envelope_t) :: envelope, bad_envelope
    type(gc_full_fow_energy_quadrature_t) :: energy_mapped
    type(gc_full_fow_paired_quadrature_t) :: mapped
    type(gc_full_fow_eq17_result_t) :: eq17
    type(gc_full_fow_degree5_enclosure_t) :: enclosure
    real(dp) :: coefficients(6), nan_value
    integer :: status
    character(len=256) :: message
    real(dp) :: expected_outer, psi_star, dpsi_star_dp_phi, scaled_magnitude

    call initialize_gc_full_fow_reference(2.0_dp, -3.0_dp, 5.0_dp, 2.0_dp, &
        4.0_dp, reference, status, message)
    call require_status('reference initialization', status, &
        GC_FULL_FOW_NORMALIZATION_SUCCESS, message)
    call require_close('generated Phi_eff', reference%phi_eff, 5.0_dp/6.0_dp)
    call require_close('generated J_K scale', reference%jk_scale, 20.0_dp/3.0_dp)
    call evaluate_gc_full_fow_canonical_flux(reference, -2.0_dp, psi_star, &
        dpsi_star_dp_phi, status, message)
    call require_status('canonical flux', status, &
        GC_FULL_FOW_NORMALIZATION_SUCCESS, message)
    call require_close('canonical psi_star', psi_star, 10.0_dp/3.0_dp)
    call require_close('canonical derivative', dpsi_star_dp_phi, -5.0_dp/3.0_dp)
    call evaluate_gc_full_fow_scaled_magnitude(-3.0_dp, 2.0_dp, &
        scaled_magnitude, status, message)
    call require_status('refinement scaling', status, &
        GC_FULL_FOW_NORMALIZATION_SUCCESS, message)
    call require_close('refinement scaled magnitude', scaled_magnitude, 1.5_dp)

    physical%h_physical = 9.0_dp
    physical%jk_physical = 8.0_dp
    physical%psi_star_physical = 2.0_dp
    physical%dpsi_star_dx_physical = -1.0_dp
    physical%tau_b_physical = 3.0_dp
    physical%omega_b_physical = 6.0_dp
    physical%omega_phi_physical = -2.0_dp
    physical%domega_b_dx_physical = 4.0_dp
    physical%domega_phi_dx_physical = -6.0_dp
    physical%hm_physical = cmplx(4.0_dp, -6.0_dp, kind=dp)
    call normalize_gc_full_fow_sample(reference, physical, normalized, status, &
        message)
    call require_status('sample normalization', status, &
        GC_FULL_FOW_NORMALIZATION_SUCCESS, message)
    call require_close('H hat', normalized%h_hat, 4.5_dp)
    call require_close('J_K hat', normalized%jk_hat, 1.2_dp)
    call require_close('psi_star hat', normalized%psi_star_hat, 2.4_dp)
    call require_close('d psi_star hat/dx', normalized%dpsi_star_hat_dx, -1.2_dp)
    call require_close('tau_b hat', normalized%tau_b_hat, 12.0_dp)
    call require_close('omega_b hat', normalized%omega_b_hat, 1.5_dp)
    call require_close('omega_phi hat', normalized%omega_phi_hat, -0.5_dp)
    call require_close('d omega_b hat/dx', normalized%domega_b_hat_dx, 1.0_dp)
    call require_close('d omega_phi hat/dx', normalized%domega_phi_hat_dx, -1.5_dp)
    call require_close('Hm real hat', real(normalized%hm_hat, dp), 2.0_dp)
    call require_close('Hm imag hat', aimag(normalized%hm_hat), -3.0_dp)

    certificate%qphi_min_certified = .true.
    certificate%bmod_min_certified = .true.
    certificate%qphi_min = 3.0_dp
    certificate%bmod_min = 4.0_dp
    certificate%method = 'independent manufactured phase-space certificate'
    certificate%certificate_id = 'normalization-adapter-test'
    call evaluate_gc_full_fow_jk_envelope(reference, 9.0_dp, certificate, &
        envelope, status, message)
    call require_status('J_K envelope', status, &
        GC_FULL_FOW_NORMALIZATION_SUCCESS, message)
    ! Independent oracle: omega_c,min=|q|B_min/(mc)=1.2 and the positive
    ! action candidate is (H-qPhi_min)/omega_c,min=(9-3)/1.2=5.
    call require_close('omega_c minimum', envelope%omega_c_min, 1.2_dp)
    call require_close('J_K candidate', envelope%jk_candidate_physical, 5.0_dp)
    call require_close('J_K maximum', envelope%jk_max_physical, 5.0_dp)

    call map_gc_full_fow_energy_quadrature(reference, certificate%qphi_min, &
        0.25_dp, 0.4_dp, energy_mapped, status, message)
    call require_status('energy quadrature map', status, &
        GC_FULL_FOW_NORMALIZATION_SUCCESS, message)
    ! Independent oracle for qPhi_min/Eref=3/2 and t=1/4.
    call require_close('energy-only H hat', energy_mapped%h_hat, 11.0_dp/6.0_dp)
    call require_close('energy-only H physical', energy_mapped%h_physical, &
        11.0_dp/3.0_dp)
    call require_close('energy-only map weight', energy_mapped%normalized_weight, &
        32.0_dp/45.0_dp)

    call map_gc_full_fow_paired_quadrature(reference, envelope, 0.25_dp, &
        0.4_dp, 0.3_dp, 0.6_dp, mapped, status, message)
    call require_status('paired quadrature map', status, &
        GC_FULL_FOW_NORMALIZATION_SUCCESS, message)
    ! Independent oracle for t=0.25, u=0.3 and J_K,max=5.
    call require_close('mapped H hat', mapped%h_hat, 11.0_dp/6.0_dp)
    call require_close('mapped H physical', mapped%h_physical, 11.0_dp/3.0_dp)
    call require_close('energy map weight', mapped%energy_normalized_weight, &
        32.0_dp/45.0_dp)
    call require_close('mapped J_K hat maximum', mapped%jk_hat_max, 0.75_dp)
    call require_close('mapped J_K hat', mapped%jk_hat, 0.225_dp)
    call require_close('mapped J_K physical', mapped%jk_physical, 1.5_dp)
    call require_close('action map weight', mapped%action_normalized_weight, &
        0.45_dp)
    call require_close('paired map weight', mapped%paired_weight, 0.32_dp)

    call evaluate_gc_full_fow_eq17(reference, 9.0_dp, 0.2_dp, 1.5_dp, 0.2_dp, &
        -0.4_dp, 3.0_dp, 0.75_dp, eq17, status, message)
    call require_status('Eq. 17 adapter', status, &
        GC_FULL_FOW_NORMALIZATION_SUCCESS, message)
    call require_close('Eq. 17 q Phi', eq17%qphi_energy, -0.6_dp)
    call require_close('Eq. 17 force bracket', eq17%force_bracket, -2.36_dp)
    expected_outer = -(acos(-1.0_dp)**1.5_dp/4.0_dp)*2.0_dp*(5.0_dp/6.0_dp)* &
        3.0_dp/(1.5_dp/2.0_dp)**1.5_dp*exp((-0.6_dp-9.0_dp)/1.5_dp)*0.75_dp
    call require_close('Eq. 17 outer factor', eq17%outer_factor, expected_outer)

    coefficients = [-10.0_dp, 1.0_dp, -2.0_dp, 0.5_dp, 0.0_dp, 1.0_dp]
    call evaluate_gc_full_fow_degree5_cell_enclosure(coefficients, 0.5_dp, &
        enclosure, status, message)
    call require_status('degree-five enclosure', status, &
        GC_FULL_FOW_NORMALIZATION_SUCCESS, message)
    call require_close('degree-five tail bound', enclosure%tail_absolute_bound, &
        1.09375_dp)
    call require_close('degree-five lower bound', enclosure%lower_bound, &
        -11.09375_dp)
    call require_close('degree-five upper bound', enclosure%upper_bound, &
        -8.90625_dp)

    call initialize_gc_full_fow_reference(2.0_dp, 0.0_dp, 5.0_dp, 2.0_dp, &
        4.0_dp, other_reference, status, message)
    call require_status('zero charge rejection', status, &
        GC_FULL_FOW_NORMALIZATION_DOMAIN_ERROR, message)
    nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
    call initialize_gc_full_fow_reference(2.0_dp, -3.0_dp, 5.0_dp, nan_value, &
        4.0_dp, other_reference, status, message)
    call require_status('nonfinite reference rejection', status, &
        GC_FULL_FOW_NORMALIZATION_NONFINITE_INPUT, message)

    call normalize_gc_full_fow_sample(uninitialized_reference, physical, &
        normalized, status, message)
    call require_status('uninitialized sample rejection', status, &
        GC_FULL_FOW_NORMALIZATION_NOT_INITIALIZED, message)
    bad_physical = physical
    bad_physical%jk_physical = -1.0_dp
    call normalize_gc_full_fow_sample(reference, bad_physical, normalized, &
        status, message)
    call require_status('negative action rejection', status, &
        GC_FULL_FOW_NORMALIZATION_DOMAIN_ERROR, message)

    call evaluate_gc_full_fow_jk_envelope(reference, 9.0_dp, bad_certificate, &
        bad_envelope, status, message)
    call require_status('uncertified bound rejection', status, &
        GC_FULL_FOW_NORMALIZATION_CERTIFICATE_REQUIRED, message)
    zero_field_certificate = certificate
    zero_field_certificate%bmod_min = 0.0_dp
    call evaluate_gc_full_fow_jk_envelope(reference, 9.0_dp, &
        zero_field_certificate, bad_envelope, status, message)
    call require_status('zero B_min rejection', status, &
        GC_FULL_FOW_NORMALIZATION_DOMAIN_ERROR, message)

    call map_gc_full_fow_paired_quadrature(reference, bad_envelope, 0.25_dp, &
        0.4_dp, 0.3_dp, 0.6_dp, mapped, status, message)
    call require_status('uncertified map rejection', status, &
        GC_FULL_FOW_NORMALIZATION_CERTIFICATE_REQUIRED, message)
    call map_gc_full_fow_paired_quadrature(reference, envelope, 1.0_dp, 0.4_dp, &
        0.3_dp, 0.6_dp, mapped, status, message)
    call require_status('singular energy-node rejection', status, &
        GC_FULL_FOW_NORMALIZATION_DOMAIN_ERROR, message)
    call initialize_gc_full_fow_reference(3.0_dp, -3.0_dp, 5.0_dp, 2.0_dp, &
        4.0_dp, other_reference, status, message)
    call require_status('second reference initialization', status, &
        GC_FULL_FOW_NORMALIZATION_SUCCESS, message)
    call map_gc_full_fow_paired_quadrature(other_reference, envelope, 0.25_dp, &
        0.4_dp, 0.3_dp, 0.6_dp, mapped, status, message)
    call require_status('reference mismatch rejection', status, &
        GC_FULL_FOW_NORMALIZATION_REFERENCE_MISMATCH, message)

    call evaluate_gc_full_fow_eq17(reference, 9.0_dp, 0.2_dp, 0.0_dp, 0.2_dp, &
        -0.4_dp, 3.0_dp, 0.75_dp, eq17, status, message)
    call require_status('nonpositive temperature rejection', status, &
        GC_FULL_FOW_NORMALIZATION_DOMAIN_ERROR, message)
    nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
    coefficients(3) = nan_value
    call evaluate_gc_full_fow_degree5_cell_enclosure(coefficients, 0.5_dp, &
        enclosure, status, message)
    call require_status('nonfinite coefficient rejection', status, &
        GC_FULL_FOW_NORMALIZATION_NONFINITE_INPUT, message)
    coefficients(3) = -2.0_dp
    call evaluate_gc_full_fow_degree5_cell_enclosure(coefficients, -0.5_dp, &
        enclosure, status, message)
    call require_status('negative cell width rejection', status, &
        GC_FULL_FOW_NORMALIZATION_DOMAIN_ERROR, message)

    write (*, '(A)') 'test_gc_full_fow_normalization_runtime OK'
    call pass_test

contains

    subroutine require_status(label, actual, expected, message)
        character(*), intent(in) :: label, message
        integer, intent(in) :: actual, expected

        if (actual /= expected) then
            write (*, '(A,2(1X,I0),1X,A)') trim(label), actual, expected, &
                trim(message)
            error stop 'normalization runtime status oracle failed'
        end if
        if (len_trim(message) == 0) then
            error stop 'normalization runtime did not report a status message'
        end if
    end subroutine require_status

    subroutine require_close(label, actual, expected)
        character(*), intent(in) :: label
        real(dp), intent(in) :: actual, expected
        real(dp) :: scale

        scale = max(1.0_dp, max(abs(actual), abs(expected)))
        if (abs(actual-expected) > 2.0e-12_dp*scale) then
            write (*, '(A,2(1X,ES24.16))') trim(label), actual, expected
            error stop 'normalization runtime behavioral oracle failed'
        end if
    end subroutine require_close

end program test_gc_full_fow_normalization_runtime
