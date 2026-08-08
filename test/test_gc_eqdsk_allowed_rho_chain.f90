program test_gc_eqdsk_allowed_rho_chain
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_eqdsk_allowed_axis_rho_chain_interval_symbolic, only: &
        evaluate_neort_eqdsk_allowed_axis_rho_chain_interval
    use neort_eqdsk_allowed_axis_rho_chain_symbolic, only: &
        evaluate_neort_eqdsk_allowed_axis_rho_chain
    use neort_eqdsk_allowed_rho_chain_interval_symbolic, only: &
        evaluate_neort_eqdsk_allowed_rho_chain_interval
    use neort_eqdsk_allowed_rho_chain_symbolic, only: &
        evaluate_neort_eqdsk_allowed_rho_chain
    use neort_eqdsk_flux_profile_segment_interval_symbolic, only: &
        evaluate_neort_eqdsk_flux_profile_segment_interval
    use neort_eqdsk_rho_tor_map_interval_symbolic, only: &
        evaluate_neort_eqdsk_rho_tor_map_interval
    use neort_gc_outward_interval, only: gc_outward_interval, &
        gc_outward_interval_t
    implicit none

    real(dp) :: p_first, p_second, r_first, r_second
    real(dp) :: g_first, g_second, canonical_first
    real(dp) :: axis_g_first, axis_canonical_first
    type(gc_outward_interval_t) :: interval_output(7)
    type(gc_outward_interval_t) :: axis_output(2)
    type(gc_outward_interval_t) :: s_tor, ds_drho
    type(gc_outward_interval_t) :: psihat, dpsihat_ds, inverse_s
    real(dp), parameter :: expected_p_first = 2.4_dp
    real(dp), parameter :: expected_p_second = 6.32_dp
    real(dp), parameter :: expected_r_first = 1.2_dp
    real(dp), parameter :: expected_r_second = 3.34_dp
    real(dp), parameter :: expected_g_first = 6.0_dp
    real(dp), parameter :: expected_g_second = 26.78_dp
    real(dp), parameter :: expected_canonical_first = -4.8_dp

    call evaluate_neort_eqdsk_allowed_rho_chain(0.8_dp, 3.0_dp, 0.5_dp, &
        2.0_dp, -0.25_dp, 5.0_dp, 7.0_dp, -4.0_dp, p_first, p_second, &
        r_first, r_second, g_first, g_second, canonical_first)
    call require_close(p_first, expected_p_first, 'flux first chain')
    call require_close(p_second, expected_p_second, 'flux second chain')
    call require_close(r_first, expected_r_first, 'inverse first chain')
    call require_close(r_second, expected_r_second, 'inverse second chain')
    call require_close(g_first, expected_g_first, 'energy first chain')
    call require_close(g_second, expected_g_second, 'energy second chain')
    call require_close(canonical_first, expected_canonical_first, &
        'canonical first chain')

    call evaluate_neort_eqdsk_allowed_rho_chain_interval(point(0.8_dp), &
        point(3.0_dp), point(0.5_dp), point(2.0_dp), point(-0.25_dp), &
        point(5.0_dp), point(7.0_dp), point(-4.0_dp), interval_output(1), &
        interval_output(2), interval_output(3), interval_output(4), &
        interval_output(5), interval_output(6), interval_output(7))
    call require_contains(interval_output(1), expected_p_first, &
        'interval flux first chain')
    call require_contains(interval_output(2), expected_p_second, &
        'interval flux second chain')
    call require_contains(interval_output(3), expected_r_first, &
        'interval inverse first chain')
    call require_contains(interval_output(4), expected_r_second, &
        'interval inverse second chain')
    call require_contains(interval_output(5), expected_g_first, &
        'interval energy first chain')
    call require_contains(interval_output(6), expected_g_second, &
        'interval energy second chain')
    call require_contains(interval_output(7), expected_canonical_first, &
        'interval canonical first chain')

    call evaluate_neort_eqdsk_allowed_axis_rho_chain(1.5_dp, 5.0_dp, &
        -4.0_dp, axis_g_first, axis_canonical_first)
    call require_close(axis_g_first, 7.5_dp, 'axis energy chain')
    call require_close(axis_canonical_first, -6.0_dp, &
        'axis canonical chain')
    call evaluate_neort_eqdsk_allowed_axis_rho_chain_interval( &
        point(1.5_dp), point(5.0_dp), point(-4.0_dp), axis_output(1), &
        axis_output(2))
    call require_contains(axis_output(1), 7.5_dp, &
        'interval axis energy chain')
    call require_contains(axis_output(2), -6.0_dp, &
        'interval axis canonical chain')

    call evaluate_neort_eqdsk_rho_tor_map_interval( &
        gc_outward_interval(0.3_dp, 0.4_dp), s_tor, ds_drho)
    call require_contains(s_tor, 0.09_dp, 'rho map lower endpoint')
    call require_contains(s_tor, 0.16_dp, 'rho map upper endpoint')
    call require_contains(ds_drho, 0.6_dp, 'rho derivative lower endpoint')
    call require_contains(ds_drho, 0.8_dp, 'rho derivative upper endpoint')

    call evaluate_neort_eqdsk_flux_profile_segment_interval( &
        gc_outward_interval(0.2_dp, 0.3_dp), point(0.06_dp), &
        point(0.0_dp), point(1.0_dp), point(0.0_dp), point(10.0_dp), &
        point(2.0_dp), point(20.0_dp), psihat, dpsihat_ds, inverse_s)
    call require_contains(psihat, 0.05_dp, 'flux profile lower endpoint')
    call require_contains(psihat, 0.075_dp, 'flux profile upper endpoint')
    call require_contains(dpsihat_ds, 0.25_dp, 'flux profile slope')
    call require_contains(inverse_s, 0.24_dp, 'flux profile inverse')

    write (*, '(A)') 'test_gc_eqdsk_allowed_rho_chain: PASS'

contains

    pure function point(value) result(interval)
        real(dp), intent(in) :: value
        type(gc_outward_interval_t) :: interval

        interval = gc_outward_interval(value, value)
    end function point

    subroutine require_close(actual, expected, label)
        real(dp), intent(in) :: actual, expected
        character(len=*), intent(in) :: label

        call require(abs(actual-expected) <= &
            1.0e-13_dp*max(1.0_dp, abs(expected)), &
            trim(label)//' does not match the independent oracle')
    end subroutine require_close

    subroutine require_contains(interval, expected, label)
        type(gc_outward_interval_t), intent(in) :: interval
        real(dp), intent(in) :: expected
        character(len=*), intent(in) :: label

        call require(interval%lo <= expected .and. interval%hi >= expected, &
            trim(label)//' is not enclosed')
    end subroutine require_contains

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) error stop message
    end subroutine require

end program test_gc_eqdsk_allowed_rho_chain
