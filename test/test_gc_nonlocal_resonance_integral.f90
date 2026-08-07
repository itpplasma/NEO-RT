program test_gc_nonlocal_resonance_integral
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_gc_nonlocal_resonance_integral, only: &
        integrate_gc_nonlocal_resonance
    use neort_gc_nonlocal_resonance_types, only: &
        GC_NONLOCAL_COMPONENT_IDENTITY, GC_NONLOCAL_INVALID_INPUT, &
        GC_NONLOCAL_PARTIAL, GC_NONLOCAL_SAMPLE_VALID, &
        GC_NONLOCAL_SAMPLE_WALL, GC_NONLOCAL_SUCCESS, &
        gc_nonlocal_component_t, gc_nonlocal_orbit_sample_t, &
        gc_nonlocal_resonance_options_t, gc_nonlocal_resonance_result_t
    implicit none

    type(gc_nonlocal_component_t) :: components(2)
    type(gc_nonlocal_component_t) :: coincident_components(2)
    type(gc_nonlocal_resonance_options_t) :: options
    type(gc_nonlocal_resonance_result_t) :: result
    integer :: status

    components(1) = gc_nonlocal_component_t(101, 1, -1.0_dp, 1.0_dp)
    components(2) = gc_nonlocal_component_t(202, -1, 2.0_dp, 4.0_dp)
    options = gc_nonlocal_resonance_options_t()
    options%scan_intervals = 8
    options%max_roots = 8
    options%force_count = 2
    options%residual_tolerance = 1.0e-12_dp
    options%x_tolerance = 1.0e-12_dp

    coincident_components(1) = gc_nonlocal_component_t(303, 1, -1.0_dp, 1.0_dp)
    coincident_components(2) = gc_nonlocal_component_t(303, -1, -1.0_dp, 1.0_dp)

    call integrate_gc_nonlocal_resonance(manufactured_orbit, 7.0_dp, 11.0_dp, &
        2, 1, components, options, result, status)
    call require(status == GC_NONLOCAL_SUCCESS, &
        'manufactured nonlocal integral did not certify')
    call require(result%certified, 'successful integral was not certified')
    call require(result%nroots == 2, 'both disconnected resonance roots missing')
    call require(result%root_component_id(1) == 101 .and. &
        result%root_component_id(2) == 202, &
        'root component identity was not retained')
    call require(result%root_sigma(1) == 1 .and. result%root_sigma(2) == -1, &
        'native sigma labels were not retained')
    call require_close(result%contribution(1), 11.25_dp, &
        'native signed first force', 2.0e-10_dp)
    call require_close(result%contribution(2), -51.25_dp, &
        'native signed second force', 2.0e-10_dp)
    call require_close(result%component_contribution(1, 1), -3.75_dp, &
        'first component first force', 2.0e-10_dp)
    call require_close(result%component_contribution(1, 2), 15.0_dp, &
        'second component first force', 2.0e-10_dp)

    call integrate_gc_nonlocal_resonance(coincident_orbit, 7.0_dp, 11.0_dp, &
        1, 0, coincident_components, options, result, status)
    call require(status == GC_NONLOCAL_SUCCESS, &
        'coincident sigma branches were rejected')
    call require(result%nroots == 2, &
        'coincident sigma branches were merged')
    call require(result%root_component_id(1) == 303 .and. &
        result%root_component_id(2) == 303, &
        'tuple component identity lost the shared component id')
    call require(result%root_sigma(1) == 1 .and. result%root_sigma(2) == -1, &
        'coincident sigma branch identity was lost')
    call require_close(result%contribution(1), -4.0_dp, &
        'coincident branch signed force', 2.0e-10_dp)

    call integrate_gc_nonlocal_resonance(wall_orbit, 7.0_dp, 11.0_dp, 2, 1, &
        components, options, result, status)
    call require(status == GC_NONLOCAL_PARTIAL, &
        'unresolved wall sample was silently accepted')
    call require(.not. result%certified, &
        'partial wall scan was incorrectly certified')

    call check_fail_closed_contract(components, options)

    write (*, '(a)') 'test_gc_nonlocal_resonance_integral OK'

contains

    subroutine manufactured_orbit(h0, jperp, x, sigma, component_id, sample, &
            callback_status)
        real(dp), intent(in) :: h0, jperp, x
        integer, intent(in) :: sigma, component_id
        type(gc_nonlocal_orbit_sample_t), intent(out) :: sample
        integer, intent(out) :: callback_status

        sample = gc_nonlocal_orbit_sample_t()
        if (.not. all(ieee_is_finite([h0, jperp]))) then
            callback_status = 1
            return
        end if
        sample%status = GC_NONLOCAL_SAMPLE_VALID
        sample%sigma = sigma
        sample%component_id = component_id
        sample%psi_star = x + real(component_id, dp)
        sample%tau_b = 0.5_dp
        sample%h_m = cmplx(1.0_dp, 2.0_dp, kind=dp)
        sample%nforce = 2
        sample%thermodynamic_force(1:2) = [-3.0_dp, 7.0_dp]
        if (component_id == 202) then
            sample%omega_b = 3.63_dp - x
            sample%domega_b_dx = -1.0_dp
            sample%dpsi_star_dx = -3.0_dp
            sample%tau_b = 2.0_dp
            sample%thermodynamic_force(1:2) = [0.5_dp, -2.0_dp]
        else
            sample%omega_b = x + 0.5_dp
            sample%domega_b_dx = 1.0_dp
            sample%dpsi_star_dx = 2.0_dp
        end if
        sample%omega_phi = -1.0_dp
        sample%domega_phi_dx = 0.0_dp
        sample%derivatives_available = .true.
        callback_status = GC_NONLOCAL_SUCCESS
    end subroutine manufactured_orbit

    subroutine wall_orbit(h0, jperp, x, sigma, component_id, sample, &
            callback_status)
        real(dp), intent(in) :: h0, jperp, x
        integer, intent(in) :: sigma, component_id
        type(gc_nonlocal_orbit_sample_t), intent(out) :: sample
        integer, intent(out) :: callback_status

        call manufactured_orbit(h0, jperp, x, sigma, component_id, sample, &
            callback_status)
        if (component_id == 202) then
            if (x > 3.5_dp) sample%status = GC_NONLOCAL_SAMPLE_WALL
        end if
    end subroutine wall_orbit

    subroutine coincident_orbit(h0, jperp, x, sigma, component_id, sample, &
            callback_status)
        real(dp), intent(in) :: h0, jperp, x
        integer, intent(in) :: sigma, component_id
        type(gc_nonlocal_orbit_sample_t), intent(out) :: sample
        integer, intent(out) :: callback_status

        sample = gc_nonlocal_orbit_sample_t()
        if (.not. all(ieee_is_finite([h0, jperp]))) then
            callback_status = 1
            return
        end if
        sample%status = GC_NONLOCAL_SAMPLE_VALID
        sample%sigma = sigma
        sample%component_id = component_id
        sample%psi_star = x
        sample%dpsi_star_dx = 1.0_dp
        sample%tau_b = 1.0_dp
        sample%h_m = cmplx(1.0_dp, 0.0_dp, kind=dp)
        sample%omega_phi = 0.0_dp
        sample%domega_phi_dx = 0.0_dp
        sample%nforce = 2
        if (sigma == 1) then
            sample%omega_b = x + 0.25_dp
            sample%thermodynamic_force(1:2) = [2.0_dp, 0.0_dp]
        else
            sample%omega_b = x - 0.25_dp
            sample%dpsi_star_dx = 2.0_dp
            sample%thermodynamic_force(1:2) = [-3.0_dp, 0.0_dp]
        end if
        sample%domega_b_dx = 1.0_dp
        sample%derivatives_available = .true.
        callback_status = GC_NONLOCAL_SUCCESS
    end subroutine coincident_orbit

    subroutine check_fail_closed_contract(valid_components, valid_options)
        type(gc_nonlocal_component_t), intent(in) :: valid_components(:)
        type(gc_nonlocal_resonance_options_t), intent(in) :: valid_options

        type(gc_nonlocal_component_t) :: bad_components(2)
        type(gc_nonlocal_resonance_options_t) :: bad_options

        bad_components = valid_components
        bad_components(2)%component_id = bad_components(1)%component_id
        bad_components(2)%sigma = bad_components(1)%sigma
        call integrate_gc_nonlocal_resonance(manufactured_orbit, 7.0_dp, 11.0_dp, &
            2, 1, bad_components, valid_options, result, status)
        call require(status == GC_NONLOCAL_COMPONENT_IDENTITY, &
            'duplicate component identity was accepted')
        bad_options = valid_options
        bad_options%scan_intervals = 32769
        call integrate_gc_nonlocal_resonance(manufactured_orbit, 7.0_dp, 11.0_dp, &
            2, 1, valid_components, bad_options, result, status)
        call require(status == GC_NONLOCAL_INVALID_INPUT, &
            'unbounded scan option was accepted')
    end subroutine check_fail_closed_contract

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message

        if (.not. condition) then
            write (*, '(a)') trim(message)
            error stop 1
        end if
    end subroutine require

    subroutine require_close(value, reference, label, tolerance)
        real(dp), intent(in) :: value, reference, tolerance
        character(*), intent(in) :: label

        if (abs(value - reference) > tolerance) then
            write (*, '(a,2(1x,es24.16))') trim(label), value, reference
            error stop 'manufactured oracle mismatch'
        end if
    end subroutine require_close

end program test_gc_nonlocal_resonance_integral
