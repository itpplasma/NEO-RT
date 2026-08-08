program test_gc_eqdsk_operational_point
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_eqdsk_allowed_energy_symbolic, only: &
        evaluate_neort_eqdsk_allowed_energy
    use neort_eqdsk_canonical_cut_symbolic, only: &
        evaluate_neort_eqdsk_canonical_cut
    use neort_eqdsk_cut_r_flux_curvature_symbolic, only: &
        evaluate_neort_eqdsk_cut_r_flux_curvature
    use neort_gc_eqdsk_cut_jet, only: eqdsk_cut_jet_t
    use neort_gc_eqdsk_operational_point, only: &
        EQDSK_OPERATIONAL_POINT_INVALID_CANONICAL_SQRT, &
        EQDSK_OPERATIONAL_POINT_INVALID_CUT_DENOMINATOR, &
        EQDSK_OPERATIONAL_POINT_INVALID_CUT_JET, &
        EQDSK_OPERATIONAL_POINT_INVALID_INPUT, &
        EQDSK_OPERATIONAL_POINT_SUCCESS, eqdsk_operational_point_t, &
        evaluate_eqdsk_operational_point
    use util_for_test, only: pass_test
    implicit none

    type(eqdsk_cut_jet_t) :: cut, invalid_cut
    type(eqdsk_operational_point_t) :: actual
    real(dp) :: expected(32), observed(32)
    real(dp) :: dZ_dR, d2Z_dR2, dpsihat_dR, d2psihat_dR2
    real(dp) :: field_norm_squared, psi_physical, dpsi_physical_dR
    real(dp) :: d2psi_physical_dR2, bmod, dbmod_dR, d2bmod_dR2
    real(dp) :: omega_c, domega_c_dR, d2omega_c_dR2
    real(dp) :: energy_margin, denergy_margin_dR, d2energy_margin_dR2
    real(dp) :: v_parallel_squared, dv_parallel_squared_dR
    real(dp) :: d2v_parallel_squared_dR2, bphi_covariant
    real(dp) :: dbphi_covariant_dR, d2bphi_covariant_dR2
    real(dp) :: v_parallel, dv_parallel_dR, d2v_parallel_dR2
    real(dp) :: psi_star, dpsi_star_dR, d2psi_star_dR2
    integer :: status

    call make_cut(cut)
    call evaluate_eqdsk_operational_point(cut, 2.3_dp, 1.7_dp, 5.7_dp, &
        20.0_dp, 0.8_dp, 2.5_dp, 1.3_dp, 4.0_dp, 0.2_dp, 0.05_dp, &
        0.01_dp, 1, actual, status)
    call require(status == EQDSK_OPERATIONAL_POINT_SUCCESS, &
        'valid operational point failed')
    call require(actual%status == status, 'result status is not stable')

    call evaluate_neort_eqdsk_cut_r_flux_curvature( &
        cut%d_cut_numerator_d_R, cut%d_cut_numerator_d_Z, &
        cut%d2_cut_numerator_d_R2, cut%d2_cut_numerator_d_RdZ, &
        cut%d2_cut_numerator_d_Z2, cut%psi_jet(2), cut%psi_jet(3), &
        cut%psi_jet(4), cut%psi_jet(5), cut%psi_jet(6), 5.7_dp, dZ_dR, &
        d2Z_dR2, dpsihat_dR, d2psihat_dR2)
    call evaluate_neort_eqdsk_allowed_energy(2.3_dp, 1.7_dp, &
        cut%psi_jet(1), cut%psi_jet(2), cut%psi_jet(3), cut%psi_jet(4), &
        cut%psi_jet(5), cut%psi_jet(6), cut%psi_jet(7), cut%psi_jet(8), &
        cut%psi_jet(9), cut%psi_jet(10), cut%f_jet(1), cut%f_jet(2), &
        cut%f_jet(3), 5.7_dp, dZ_dR, d2Z_dR2, 20.0_dp, 0.8_dp, 2.5_dp, &
        1.3_dp, 4.0_dp, 0.2_dp, 0.05_dp, 0.01_dp, field_norm_squared, &
        psi_physical, dpsi_physical_dR, bmod, dbmod_dR, d2bmod_dR2, &
        omega_c, domega_c_dR, d2omega_c_dR2, energy_margin, &
        denergy_margin_dR, d2energy_margin_dR2, v_parallel_squared, &
        dv_parallel_squared_dR, d2v_parallel_squared_dR2, bphi_covariant, &
        dbphi_covariant_dR, d2psi_physical_dR2, d2bphi_covariant_dR2)
    call evaluate_neort_eqdsk_canonical_cut(v_parallel_squared, &
        dv_parallel_squared_dR, 2.5_dp, 1.3_dp, 4.0_dp, 1.0_dp, &
        psi_physical, dpsi_physical_dR, bphi_covariant, &
        dbphi_covariant_dR, d2v_parallel_squared_dR2, &
        d2psi_physical_dR2, d2bphi_covariant_dR2, v_parallel, &
        dv_parallel_dR, psi_star, dpsi_star_dR, d2v_parallel_dR2, &
        d2psi_star_dR2)

    expected = [dZ_dR, d2Z_dR2, dpsihat_dR, d2psihat_dR2, &
        field_norm_squared, psi_physical, dpsi_physical_dR, &
        d2psi_physical_dR2, bmod, dbmod_dR, d2bmod_dR2, omega_c, &
        domega_c_dR, d2omega_c_dR2, energy_margin, denergy_margin_dR, &
        d2energy_margin_dR2, v_parallel_squared, dv_parallel_squared_dR, &
        d2v_parallel_squared_dR2, bphi_covariant, dbphi_covariant_dR, &
        d2bphi_covariant_dR2, v_parallel, dv_parallel_dR, &
        d2v_parallel_dR2, psi_star, dpsi_star_dR, d2psi_star_dR2, &
        cut%psi_jet(1), cut%f_jet(1), cut%cut_numerator]
    observed = [actual%dZ_dR, actual%d2Z_dR2, actual%dpsihat_dR, &
        actual%d2psihat_dR2, actual%field_norm_squared, actual%psi_physical, &
        actual%dpsi_physical_dR, actual%d2psi_physical_dR2, actual%bmod, &
        actual%dbmod_dR, actual%d2bmod_dR2, actual%omega_c, &
        actual%domega_c_dR, actual%d2omega_c_dR2, actual%energy_margin, &
        actual%denergy_margin_dR, actual%d2energy_margin_dR2, &
        actual%v_parallel_squared, actual%dv_parallel_squared_dR, &
        actual%d2v_parallel_squared_dR2, actual%bphi_covariant, &
        actual%dbphi_covariant_dR, actual%d2bphi_covariant_dR2, &
        actual%v_parallel, actual%dv_parallel_dR, actual%d2v_parallel_dR2, &
        actual%psi_star, actual%dpsi_star_dR, actual%d2psi_star_dR2, &
        cut%psi_jet(1), cut%f_jet(1), cut%cut_numerator]
    call require(maxval(abs(observed-expected)) <= &
        2.0e-12_dp*max(1.0_dp, maxval(abs(expected))), &
        'orchestration differs from direct generated kernels')

    invalid_cut = cut
    invalid_cut%d_cut_numerator_d_Z = 0.0_dp
    call evaluate_eqdsk_operational_point(invalid_cut, 2.3_dp, 1.7_dp, 5.7_dp, &
        20.0_dp, 0.8_dp, 2.5_dp, 1.3_dp, 4.0_dp, 0.2_dp, 0.05_dp, &
        0.01_dp, 1, actual, status)
    call require(status == EQDSK_OPERATIONAL_POINT_INVALID_CUT_DENOMINATOR, &
        'zero cut denominator was not rejected')

    invalid_cut = cut
    invalid_cut%psi_jet(3) = ieee_value(0.0_dp, ieee_quiet_nan)
    call evaluate_eqdsk_operational_point(invalid_cut, 2.3_dp, 1.7_dp, 5.7_dp, &
        20.0_dp, 0.8_dp, 2.5_dp, 1.3_dp, 4.0_dp, 0.2_dp, 0.05_dp, &
        0.01_dp, 1, actual, status)
    call require(status == EQDSK_OPERATIONAL_POINT_INVALID_CUT_JET, &
        'nonfinite cut jet was not rejected')

    call evaluate_eqdsk_operational_point(cut, 2.3_dp, 1.7_dp, 5.7_dp, &
        -1.0_dp, 0.8_dp, 2.5_dp, 1.3_dp, 4.0_dp, 0.2_dp, 0.05_dp, &
        0.01_dp, 1, actual, status)
    call require(status == EQDSK_OPERATIONAL_POINT_INVALID_CANONICAL_SQRT, &
        'negative allowed energy was not rejected before canonical kernel')

    call evaluate_eqdsk_operational_point(cut, 0.0_dp, 1.7_dp, 5.7_dp, &
        20.0_dp, 0.8_dp, 2.5_dp, 1.3_dp, 4.0_dp, 0.2_dp, 0.05_dp, &
        0.01_dp, 1, actual, status)
    call require(status == EQDSK_OPERATIONAL_POINT_INVALID_INPUT, &
        'invalid radius was not rejected')

    call evaluate_eqdsk_operational_point(cut, 2.3_dp, 1.7_dp, 5.7_dp, &
        20.0_dp, 0.8_dp, 2.5_dp, 1.3_dp, 4.0_dp, 0.2_dp, 0.05_dp, &
        0.01_dp, 0, actual, status)
    call require(status == EQDSK_OPERATIONAL_POINT_INVALID_INPUT, &
        'invalid parallel branch was not rejected')

    write (*, '(a)') 'test_gc_eqdsk_operational_point: PASS'
    call pass_test

contains

    subroutine make_cut(value)
        type(eqdsk_cut_jet_t), intent(out) :: value

        value = eqdsk_cut_jet_t()
        value%psi_jet = [1.2_dp, 0.7_dp, 0.4_dp, 0.2_dp, -0.1_dp, 0.3_dp, &
            0.05_dp, -0.02_dp, 0.04_dp, -0.03_dp, 0.01_dp, -0.02_dp, &
            0.03_dp, 0.02_dp, -0.01_dp]
        value%f_jet = [1.1_dp, 0.08_dp, -0.02_dp]
        value%cut_value = 0.4_dp
        value%d_cut_d_R = 0.9_dp
        value%d_cut_d_arc_phi = -0.2_dp
        value%d_cut_d_Z = 1.6_dp
        value%cut_rate = 0.3_dp
        value%absolute_cut_rate = 0.3_dp
        value%orientation_scalar = 1.0_dp
        value%cut_numerator = 0.6_dp
        value%d_cut_numerator_d_R = 0.9_dp
        value%d_cut_numerator_d_Z = 1.6_dp
        value%d2_cut_numerator_d_R2 = 0.3_dp
        value%d2_cut_numerator_d_RdZ = -0.2_dp
        value%d2_cut_numerator_d_Z2 = 0.4_dp
    end subroutine make_cut

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message

        if (.not. condition) error stop message
    end subroutine require

end program test_gc_eqdsk_operational_point
