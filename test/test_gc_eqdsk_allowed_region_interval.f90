program test_gc_eqdsk_allowed_region_interval
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_gc_eqdsk_allowed_region_interval, only: &
        EQDSK_ALLOWED_INTERVAL_PROFILE_GAP, &
        EQDSK_ALLOWED_INTERVAL_NONFINITE, &
        EQDSK_ALLOWED_INTERVAL_SUCCESS, &
        EQDSK_ALLOWED_INTERVAL_UNCERTIFIED_CUT, &
        eqdsk_allowed_interval_result_t, &
        evaluate_eqdsk_allowed_region_interval, &
        evaluate_eqdsk_profile_potential_interval
    use neort_gc_eqdsk_cut_interval, only: eqdsk_cut_interval_result_t
    use neort_gc_outward_interval, only: gc_outward_interval, &
        gc_outward_interval_t
    implicit none

    type(eqdsk_cut_interval_result_t) :: cut(1)
    type(eqdsk_cut_interval_result_t) :: partial_cut(2)
    type(eqdsk_allowed_interval_result_t) :: result, second_result
    type(gc_outward_interval_t) :: potential, first, second
    real(dp), parameter :: sqrt_ten = sqrt(10.0_dp)
    real(dp), parameter :: psi_nodes(2) = [0.0_dp, 4.0_dp]
    real(dp), parameter :: phi_nodes(2) = [0.0_dp, 8.0_dp]
    real(dp), parameter :: omega_nodes(2) = [2.0_dp, 6.0_dp]
    real(dp) :: energy, energy_first, energy_second, speed
    real(dp) :: canonical, canonical_first
    integer :: segments, status

    call make_manufactured_cut(cut(1))
    call evaluate_eqdsk_allowed_region_interval(point(2.0_dp), 2.0_dp, &
        10.0_dp, cut, psi_nodes, phi_nodes, omega_nodes, 10.0_dp, 1.0_dp, &
        2.0_dp, -1.0_dp, 2.0_dp, 1, result, status)
    call require(status == EQDSK_ALLOWED_INTERVAL_SUCCESS, &
        'manufactured allowed jet did not certify')
    call require(result%canonical_chart_certified, &
        'positive manufactured energy did not certify canonical chart')
    call require(result%cut_enclosures_covered == 1, &
        'cut enclosure count is wrong')
    call require(result%profile_segments_covered == 1, &
        'profile segment count is wrong')

    energy = 13.0_dp-sqrt_ten/4.0_dp
    energy_first = 4.0_dp+sqrt_ten/8.0_dp
    energy_second = 2.0_dp-sqrt_ten/8.0_dp
    speed = sqrt(energy)
    canonical = 2.0_dp-24.0_dp*speed/sqrt_ten
    canonical_first = 2.0_dp-4.0_dp*( &
        energy_first/(2.0_dp*speed)*6.0_dp/sqrt_ten + &
        speed*3.0_dp/sqrt_ten)

    call require_contains(result%field_norm_squared, 10.0_dp, &
        'manufactured field norm')
    call require_contains(result%psi_physical, 2.0_dp, 'physical psi')
    call require_contains(result%dpsi_physical_dR, 2.0_dp, &
        'physical psi derivative')
    call require_contains(result%bmod, sqrt_ten, 'field magnitude')
    call require_contains(result%dbmod_dR, -sqrt_ten/2.0_dp, &
        'field first derivative')
    call require_contains(result%d2bmod_dR2, sqrt_ten/2.0_dp, &
        'field second derivative')
    call require_contains(result%energy_margin, energy, 'allowed energy')
    call require_contains(result%denergy_margin_dR, energy_first, &
        'allowed energy first derivative')
    call require_contains(result%d2energy_margin_dR2, energy_second, &
        'allowed energy second derivative')
    call require_contains(result%v_parallel_squared, energy, &
        'parallel speed squared')
    call require_contains(result%bphi_covariant, 6.0_dp/sqrt_ten, &
        'covariant toroidal field ratio')
    call require_contains(result%dbphi_covariant_dR, 3.0_dp/sqrt_ten, &
        'covariant ratio derivative')
    call require_contains(result%v_parallel, speed, 'parallel speed')
    call require_contains(result%psi_star, canonical, 'canonical flux')
    call require_contains(result%dpsi_star_dR, canonical_first, &
        'canonical flux derivative')

    call evaluate_eqdsk_allowed_region_interval(point(2.0_dp), 2.0_dp, &
        10.0_dp, cut, psi_nodes, phi_nodes, omega_nodes, 10.0_dp, 2.0_dp, &
        2.0_dp, -1.0_dp, 2.0_dp, 1, second_result, status)
    call require(status == EQDSK_ALLOWED_INTERVAL_SUCCESS, &
        'second J_K evaluation failed')
    call require_contains(second_result%energy_margin, &
        13.0_dp-sqrt_ten/2.0_dp, 'J_K energy dependence')

    call evaluate_eqdsk_allowed_region_interval(point(2.0_dp), 2.0_dp, &
        10.0_dp, cut, psi_nodes, phi_nodes, omega_nodes, -10.0_dp, 1.0_dp, &
        2.0_dp, -1.0_dp, 2.0_dp, 1, second_result, status)
    call require(status == EQDSK_ALLOWED_INTERVAL_SUCCESS, &
        'forbidden-region enclosure failed')
    call require(.not. second_result%canonical_chart_certified, &
        'forbidden region entered regular square-root chart')

    call test_profile_knot(potential, first, second, segments, status)

    partial_cut = [cut(1), cut(1)]
    partial_cut(2)%psi = gc_outward_interval(1.0_dp, 0.0_dp)
    call evaluate_eqdsk_allowed_region_interval(point(2.0_dp), 2.0_dp, &
        10.0_dp, partial_cut, psi_nodes, phi_nodes, omega_nodes, 10.0_dp, &
        1.0_dp, 2.0_dp, -1.0_dp, 2.0_dp, 1, second_result, status)
    call require(status == EQDSK_ALLOWED_INTERVAL_NONFINITE, &
        'invalid cut interval was accepted')
    call require(second_result%cut_enclosures_covered == 0 .and. &
        .not. second_result%canonical_chart_certified, &
        'failed multi-cut evaluation leaked a partial result')

    cut(1)%r_chart_certified = .false.
    call evaluate_eqdsk_allowed_region_interval(point(2.0_dp), 2.0_dp, &
        10.0_dp, cut, psi_nodes, phi_nodes, omega_nodes, 10.0_dp, 1.0_dp, &
        2.0_dp, -1.0_dp, 2.0_dp, 1, second_result, status)
    call require(status == EQDSK_ALLOWED_INTERVAL_UNCERTIFIED_CUT, &
        'uncertified graph enclosure was accepted')

    write (*, '(A)') 'test_gc_eqdsk_allowed_region_interval: PASS'

contains

    subroutine make_manufactured_cut(value)
        type(eqdsk_cut_interval_result_t), intent(out) :: value

        value = eqdsk_cut_interval_result_t()
        value%psi = point(1.0_dp)
        value%psi_R = point(1.0_dp)
        value%psi_Z = point(0.0_dp)
        value%psi_RR = point(0.0_dp)
        value%psi_RZ = point(0.0_dp)
        value%psi_ZZ = point(0.0_dp)
        value%psi_RRR = point(0.0_dp)
        value%psi_RRZ = point(0.0_dp)
        value%psi_RZZ = point(0.0_dp)
        value%psi_ZZZ = point(0.0_dp)
        value%F = point(3.0_dp)
        value%dF_dpsi_hat = point(0.0_dp)
        value%d2F_dpsi_hat2 = point(0.0_dp)
        value%dZ_dR = point(0.0_dp)
        value%d2Z_dR2 = point(0.0_dp)
        value%denominator_positive_certified = .true.
        value%r_chart_certified = .true.
    end subroutine make_manufactured_cut

    subroutine test_profile_knot(potential, first, second, segments, status)
        type(gc_outward_interval_t), intent(out) :: potential, first, second
        integer, intent(out) :: segments, status
        real(dp), parameter :: knot_psi(3) = [0.0_dp, 2.0_dp, 4.0_dp]
        real(dp), parameter :: knot_phi(3) = [0.0_dp, 4.0_dp, 8.0_dp]
        real(dp), parameter :: knot_omega(3) = [2.0_dp, 6.0_dp, 2.0_dp]

        call evaluate_eqdsk_profile_potential_interval(point(2.0_dp), &
            knot_psi, knot_phi, knot_omega, 2.0_dp, potential, first, second, &
            segments, status)
        call require(status == EQDSK_ALLOWED_INTERVAL_SUCCESS, &
            'profile knot enclosure failed')
        call require(segments == 2, 'profile knot did not hull both segments')
        call require_contains(potential, 4.0_dp, 'profile knot potential')
        call require_contains(first, 3.0_dp, 'profile knot first derivative')
        call require_contains(second, -1.0_dp, &
            'profile knot right curvature')
        call require_contains(second, 1.0_dp, 'profile knot left curvature')

        call evaluate_eqdsk_profile_potential_interval( &
            gc_outward_interval(-0.1_dp, 0.0_dp), knot_psi, knot_phi, &
            knot_omega, 2.0_dp, potential, first, second, segments, status)
        call require(status == EQDSK_ALLOWED_INTERVAL_PROFILE_GAP, &
            'out-of-profile interval was accepted')
    end subroutine test_profile_knot

    pure function point(value) result(interval)
        real(dp), intent(in) :: value
        type(gc_outward_interval_t) :: interval

        interval = gc_outward_interval(value, value)
    end function point

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

end program test_gc_eqdsk_allowed_region_interval
