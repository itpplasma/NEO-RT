program test_gc_eqdsk_allowed_region_signs
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_gc_eqdsk_allowed_region_interval, only: &
        EQDSK_ALLOWED_INTERVAL_SUCCESS, eqdsk_allowed_interval_result_t, &
        evaluate_eqdsk_allowed_region_interval
    use neort_gc_eqdsk_cut_interval, only: eqdsk_cut_interval_result_t
    use neort_gc_outward_interval, only: gc_outward_interval, &
        gc_outward_interval_t
    implicit none

    type(eqdsk_cut_interval_result_t) :: cut(1)
    type(eqdsk_allowed_interval_result_t) :: sigma_plus, sigma_minus
    type(eqdsk_allowed_interval_result_t) :: charge_flip
    type(eqdsk_allowed_interval_result_t) :: zero_charge_plus, zero_charge_minus
    real(dp), parameter :: radius = 2.0_dp
    real(dp), parameter :: field_scale = 2.0_dp
    real(dp), parameter :: separatrix_psi = 4.0_dp
    real(dp), parameter :: h0 = 20.0_dp
    real(dp), parameter :: j_k = 1.2_dp
    real(dp), parameter :: mass = 1.5_dp
    real(dp), parameter :: charge = 2.0_dp
    real(dp), parameter :: c_light = 4.0_dp
    real(dp), parameter :: psi_nodes(2) = [0.0_dp, 4.0_dp]
    real(dp), parameter :: phi_nodes(2) = [3.0_dp, 3.0_dp]
    real(dp), parameter :: zero_phi_nodes(2) = [0.0_dp, 0.0_dp]
    real(dp), parameter :: omega_nodes(2) = [0.0_dp, 0.0_dp]
    ! Hand-evaluated values for the manufactured point: the 3-4-0 field
    ! triangle has norm 5, the constant potential is 3, and the resulting
    ! positive/negative-charge energy margins are 12 and 24.
    real(dp), parameter :: field_norm = 5.0_dp
    real(dp), parameter :: bmod = 5.0_dp
    real(dp), parameter :: bphi = 1.6_dp
    real(dp), parameter :: psi_physical = 2.5_dp
    real(dp), parameter :: expected_energy = 12.0_dp
    real(dp), parameter :: expected_flip_energy = 24.0_dp
    real(dp), parameter :: expected_sigma_plus = 21.7_dp
    real(dp), parameter :: expected_sigma_minus = -16.7_dp
    real(dp), parameter :: expected_zero_canonical_offset = 23.5153013442625_dp
    integer :: status

    call make_manufactured_cut(cut(1))
    call evaluate_case(cut, phi_nodes, 1, charge, sigma_plus, status)
    call require(status == EQDSK_ALLOWED_INTERVAL_SUCCESS, &
        'positive sigma evaluation failed')
    call evaluate_case(cut, phi_nodes, -1, charge, sigma_minus, status)
    call require(status == EQDSK_ALLOWED_INTERVAL_SUCCESS, &
        'negative sigma evaluation failed')
    call require(sigma_plus%canonical_chart_certified .and. &
        sigma_minus%canonical_chart_certified, &
        'sigma cases did not certify the canonical chart')

    call require_contains(sigma_plus%field_norm_squared, field_norm**2, &
        'manufactured field norm')
    call require_contains(sigma_plus%bmod, bmod, 'manufactured field magnitude')
    call require_contains(sigma_plus%bphi_covariant, bphi, &
        'manufactured covariant field')
    call require_contains(sigma_plus%energy_margin, expected_energy, &
        'independent positive-charge energy')
    call require_contains(sigma_plus%psi_star, expected_sigma_plus, &
        'independent positive-sigma canonical flux')
    call require_contains(sigma_minus%psi_star, expected_sigma_minus, &
        'independent negative-sigma canonical flux')
    call require_close(midpoint(sigma_plus%energy_margin), &
        midpoint(sigma_minus%energy_margin), 1.0e-12_dp, &
        'sigma reversal changed allowed energy')
    call require_close(midpoint(sigma_plus%psi_star)+ &
        midpoint(sigma_minus%psi_star), 2.0_dp*psi_physical, 1.0e-11_dp, &
        'sigma reversal did not reverse canonical offset')
    call require(abs(midpoint(sigma_plus%psi_star)- &
        midpoint(sigma_minus%psi_star)) > 1.0_dp, &
        'sigma reversal left canonical flux unchanged')

    call evaluate_case(cut, phi_nodes, 1, -charge, charge_flip, status)
    call require(status == EQDSK_ALLOWED_INTERVAL_SUCCESS, &
        'charge-reversed evaluation failed')
    call require_contains(charge_flip%energy_margin, expected_flip_energy, &
        'independent charge-reversed energy')
    call require(abs(expected_flip_energy-expected_energy) > 1.0_dp, &
        'manufactured potential did not distinguish charge reversal')
    call require(abs(midpoint(charge_flip%energy_margin)- &
        midpoint(sigma_plus%energy_margin)) > 1.0_dp, &
        'charge reversal was incorrectly treated as a symmetry')
    call require_close(midpoint(charge_flip%v_parallel_squared), &
        32.0_dp, 1.0e-11_dp, &
        'charge-reversed speed is inconsistent with its energy')

    call evaluate_case(cut, zero_phi_nodes, 1, charge, zero_charge_plus, status)
    call require(status == EQDSK_ALLOWED_INTERVAL_SUCCESS, &
        'positive charge zero-potential evaluation failed')
    call evaluate_case(cut, zero_phi_nodes, 1, -charge, zero_charge_minus, status)
    call require(status == EQDSK_ALLOWED_INTERVAL_SUCCESS, &
        'negative charge zero-potential evaluation failed')
    call require_close(midpoint(zero_charge_plus%energy_margin), &
        midpoint(zero_charge_minus%energy_margin), 1.0e-12_dp, &
        'zero potential did not isolate charge parity')
    call require_close(midpoint(zero_charge_plus%psi_star)-psi_physical, &
        expected_zero_canonical_offset, 1.0e-11_dp, &
        'positive-charge zero-potential canonical offset is wrong')
    call require_close(midpoint(zero_charge_minus%psi_star)-psi_physical, &
        -expected_zero_canonical_offset, 1.0e-11_dp, &
        'negative-charge zero-potential canonical offset is wrong')

    write (*, '(A)') 'test_gc_eqdsk_allowed_region_signs: PASS'

contains

    subroutine evaluate_case(cut, phi, sigma, particle_charge, result, status)
        type(eqdsk_cut_interval_result_t), intent(in) :: cut(:)
        real(dp), intent(in) :: phi(:)
        integer, intent(in) :: sigma
        real(dp), intent(in) :: particle_charge
        type(eqdsk_allowed_interval_result_t), intent(out) :: result
        integer, intent(out) :: status

        call evaluate_eqdsk_allowed_region_interval( &
            gc_outward_interval(radius, radius), field_scale, separatrix_psi, &
            cut, psi_nodes, phi, omega_nodes, h0, j_k, mass, particle_charge, &
            c_light, sigma, result, status)
    end subroutine evaluate_case

    subroutine make_manufactured_cut(value)
        type(eqdsk_cut_interval_result_t), intent(out) :: value

        value = eqdsk_cut_interval_result_t()
        value%psi = point(1.25_dp)
        value%psi_R = point(3.0_dp)
        value%psi_Z = point(0.0_dp)
        value%psi_RR = point(0.0_dp)
        value%psi_RZ = point(0.0_dp)
        value%psi_ZZ = point(0.0_dp)
        value%psi_RRR = point(0.0_dp)
        value%psi_RRZ = point(0.0_dp)
        value%psi_RZZ = point(0.0_dp)
        value%psi_ZZZ = point(0.0_dp)
        value%F = point(4.0_dp)
        value%dF_dpsi_hat = point(0.0_dp)
        value%d2F_dpsi_hat2 = point(0.0_dp)
        value%dZ_dR = point(0.0_dp)
        value%d2Z_dR2 = point(0.0_dp)
        value%denominator_positive_certified = .true.
        value%r_chart_certified = .true.
    end subroutine make_manufactured_cut

    pure function point(value) result(interval)
        real(dp), intent(in) :: value
        type(gc_outward_interval_t) :: interval

        interval = gc_outward_interval(value, value)
    end function point

    pure function midpoint(interval) result(value)
        type(gc_outward_interval_t), intent(in) :: interval
        real(dp) :: value

        value = 0.5_dp*(interval%lo+interval%hi)
    end function midpoint

    subroutine require_contains(interval, expected, label)
        type(gc_outward_interval_t), intent(in) :: interval
        real(dp), intent(in) :: expected
        character(len=*), intent(in) :: label

        call require(interval%lo <= expected .and. interval%hi >= expected, &
            trim(label)//' is not enclosed')
    end subroutine require_contains

    subroutine require_close(actual, expected, tolerance, message)
        real(dp), intent(in) :: actual, expected, tolerance
        character(len=*), intent(in) :: message

        call require(abs(actual-expected) <= tolerance*max(1.0_dp, abs(expected)), &
            message)
    end subroutine require_close

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) error stop message
    end subroutine require

end program test_gc_eqdsk_allowed_region_signs
