module test_physical_return_models
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_gc_cylindrical_model, only: GC_CYL_SUCCESS, &
        gc_cylindrical_field_sample_t, gc_cylindrical_field_t, &
        gc_cylindrical_wall_t, make_gc_cylindrical_field_sample

    implicit none

    type, extends(gc_cylindrical_field_t) :: circular_helical_field_t
    contains
        procedure :: evaluate => evaluate_circular_helical_field
    end type circular_helical_field_t

    type, extends(gc_cylindrical_wall_t) :: z_limit_wall_t
        real(dp) :: limit = 0.25_dp
    contains
        procedure :: evaluate => evaluate_z_limit_wall
    end type z_limit_wall_t

contains

    subroutine evaluate_circular_helical_field(self, position, sample, status)
        class(circular_helical_field_t), intent(in) :: self
        real(dp), intent(in) :: position(3)
        type(gc_cylindrical_field_sample_t), intent(out) :: sample
        integer, intent(out) :: status

        real(dp), parameter :: toroidal_amplitude = 6.0_dp
        real(dp) :: b(3), db_dq(3, 3)

        associate (unused_self => self)
        end associate
        b = [0.0_dp, toroidal_amplitude/position(1), 0.0_dp]
        db_dq = 0.0_dp
        db_dq(2, 1) = -toroidal_amplitude/position(1)**2
        call make_gc_cylindrical_field_sample(position(1), b, db_dq, 0.0_dp, &
            [0.0_dp, 0.0_dp, 0.0_dp], sample, status)
    end subroutine evaluate_circular_helical_field

    subroutine evaluate_z_limit_wall(self, position, distance, status)
        class(z_limit_wall_t), intent(in) :: self
        real(dp), intent(in) :: position(3)
        real(dp), intent(out) :: distance
        integer, intent(out) :: status

        associate (unused_phi => position(3))
        end associate
        distance = self%limit - position(2)
        status = GC_CYL_SUCCESS
    end subroutine evaluate_z_limit_wall

end module test_physical_return_models

program test_gc_cylindrical_physical_return
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_gc_cylindrical_model, only: GC_CYL_EQUILIBRIUM_DOMAIN, &
        GC_CYL_NO_RETURN, GC_CYL_SUCCESS, &
        GC_CYL_WALL_LOSS, gc_cylindrical_field_sample_t, &
        gc_cylindrical_invariants_t, gc_cylindrical_state_t, &
        gc_cylindrical_zero_potential_t, invariants_from_cylindrical_state
    use neort_gc_cylindrical_physical_return, only: &
        GC_CYL_PHYSICAL_EVENT_RADIAL_DOMAIN, GC_CYL_PHYSICAL_EVENT_RETURN, &
        GC_CYL_PHYSICAL_EVENT_WALL, &
        gc_cylindrical_physical_return_options_t, &
        gc_cylindrical_physical_return_t, compute_gc_cylindrical_physical_return
    use test_physical_return_models, only: circular_helical_field_t, &
        z_limit_wall_t
    use util_for_test, only: pass_test

    implicit none

    real(dp), parameter :: c_light = 1.0_dp
    real(dp), parameter :: mass = 1.0_dp
    real(dp), parameter :: charge = 1.0_dp
    real(dp), parameter :: radius = 2.0_dp
    real(dp), parameter :: expected_circular_period = 4.0_dp*acos(-1.0_dp)
    real(dp), parameter :: expected_helical_period = 12.0_dp*acos(-1.0_dp)

    type(circular_helical_field_t) :: field
    type(gc_cylindrical_zero_potential_t) :: potential
    type(z_limit_wall_t) :: wall
    type(gc_cylindrical_state_t) :: initial_state
    type(gc_cylindrical_field_sample_t) :: field_sample
    type(gc_cylindrical_invariants_t) :: invariants
    type(gc_cylindrical_physical_return_options_t) :: options
    type(gc_cylindrical_physical_return_t) :: result, coarse, fine
    integer :: status

    initial_state = gc_cylindrical_state_t()
    initial_state%R = radius
    initial_state%phi = 0.0_dp
    initial_state%p_parallel = 1.0_dp
    call field%evaluate([initial_state%R, initial_state%Z, initial_state%phi], &
        field_sample, status)
    if (status /= GC_CYL_SUCCESS) error stop 'manufactured field failed'
    call invariants_from_cylindrical_state(initial_state, field_sample, 0.0_dp, &
        mass, charge, c_light, invariants, status)
    if (status /= GC_CYL_SUCCESS) error stop 'manufactured invariants failed'

    call make_options(options)
    call check_circular_return(initial_state, invariants, options)
    call check_disarm_restart(initial_state, invariants, options)
    call check_helical_return(initial_state, invariants, options)
    call check_orientation_reversal(initial_state, field, potential, options)
    call check_event_time_convergence(initial_state, invariants, options, coarse, fine)
    call check_wall_before_return(initial_state, invariants, options)
    call check_radial_domain_before_return(initial_state, invariants, options)
    call check_no_return(initial_state, invariants, options)
    call pass_test

contains

    subroutine make_options(local_options)
        type(gc_cylindrical_physical_return_options_t), intent(out) :: local_options

        local_options = gc_cylindrical_physical_return_options_t()
        local_options%relative_tolerance = 2.0e-10_dp
        local_options%absolute_tolerance = 1.0e-12_dp
        local_options%invariant_relative_tolerance = 2.0e-8_dp
        local_options%invariant_absolute_tolerance = 1.0e-9_dp
        local_options%event_time_tolerance = 1.0e-11_dp
        local_options%event_value_tolerance = 1.0e-10_dp
        local_options%minimum_return_time = 0.1_dp
        local_options%maximum_time = 50.0_dp
        local_options%maximum_step = 0.1_dp
    end subroutine make_options

    subroutine check_circular_return(initial, local_invariants, local_options)
        type(gc_cylindrical_state_t), intent(in) :: initial
        type(gc_cylindrical_invariants_t), intent(in) :: local_invariants
        type(gc_cylindrical_physical_return_options_t), intent(in) :: local_options

        type(gc_cylindrical_physical_return_t) :: local_result

        call compute_gc_cylindrical_physical_return(field, potential, initial, &
            local_invariants, mass, charge, c_light, circular_event, &
            local_options, local_result)
        if (local_result%status /= GC_CYL_SUCCESS) then
            error stop 'circular physical return failed'
        end if
        if (local_result%event_kind /= GC_CYL_PHYSICAL_EVENT_RETURN) then
            error stop 'circular event kind failed'
        end if
        if (.not. local_result%physical_return_found) then
            error stop 'circular return was not certified'
        end if
        call require_close('circular period', local_result%period, &
            expected_circular_period, 3.0e-8_dp)
        call require_close('circular delta phi', local_result%delta_phi, &
            2.0_dp*acos(-1.0_dp), 3.0e-9_dp)
        if (local_result%return_orientation /= 1) then
            error stop 'circular orientation was not rising'
        end if
        if (local_result%maximum_invariant_scaled_drift > 1.0_dp) then
            error stop 'circular invariant drift was accepted'
        end if
    end subroutine check_circular_return

    subroutine check_disarm_restart(initial, local_invariants, local_options)
        type(gc_cylindrical_state_t), intent(in) :: initial
        type(gc_cylindrical_invariants_t), intent(in) :: local_invariants
        type(gc_cylindrical_physical_return_options_t), intent(in) :: local_options

        type(gc_cylindrical_physical_return_options_t) :: disarm_options
        type(gc_cylindrical_physical_return_t) :: local_result

        disarm_options = local_options
        disarm_options%require_opposite_intersection = .true.
        disarm_options%require_transverse_intersection = .true.
        ! This is manufactured complete-atlas evidence for the generic event
        ! driver only.  The production EQDSK adapter cannot set this flag
        ! without the committed Fortsym cut-atlas theorem certificate.
        disarm_options%complete_atlas_multiplicity_certified = .true.
        call compute_gc_cylindrical_physical_return(field, potential, initial, &
            local_invariants, mass, charge, c_light, circular_event, &
            disarm_options, local_result, return_event_rate=circular_event_rate)
        if (local_result%status /= GC_CYL_SUCCESS) then
            error stop 'two-stage disarm return failed'
        end if
        if (.not. local_result%physical_return_found) then
            error stop 'two-stage disarm return was not certified'
        end if
        if (local_result%intersection_count /= 2) then
            error stop 'two-stage disarm did not record two crossings'
        end if
        if (local_result%intersection_orientations(1) /= -1 .or. &
                local_result%intersection_orientations(2) /= 1) then
            error stop 'two-stage crossing orientations are wrong'
        end if
        if (local_result%intersection_times(2) <= &
                local_result%intersection_times(1) + &
                disarm_options%event_time_tolerance) then
            error stop 'same-oriented event rediscovered at t=start'
        end if
        call require_close('disarmed full period', local_result%period, &
            expected_circular_period, 3.0e-8_dp)
        call require_close('true same-oriented event time', &
            local_result%intersection_times(2), expected_circular_period, &
            3.0e-8_dp)
    end subroutine check_disarm_restart

    subroutine check_helical_return(initial, local_invariants, local_options)
        type(gc_cylindrical_state_t), intent(in) :: initial
        type(gc_cylindrical_invariants_t), intent(in) :: local_invariants
        type(gc_cylindrical_physical_return_options_t), intent(in) :: local_options

        type(gc_cylindrical_physical_return_t) :: local_result

        call compute_gc_cylindrical_physical_return(field, potential, initial, &
            local_invariants, mass, charge, c_light, helical_event, &
            local_options, local_result)
        if (local_result%status /= GC_CYL_SUCCESS) then
            error stop 'helical physical return failed'
        end if
        call require_close('helical period', local_result%period, &
            expected_helical_period, 2.0e-7_dp)
        call require_close('helical delta phi', local_result%delta_phi, &
            6.0_dp*acos(-1.0_dp), 3.0e-8_dp)
        if (local_result%return_orientation /= 1) then
            error stop 'helical orientation was not rising'
        end if
    end subroutine check_helical_return

    subroutine check_orientation_reversal(initial, local_field, local_potential, &
            local_options)
        type(gc_cylindrical_state_t), intent(in) :: initial
        type(circular_helical_field_t), intent(in) :: local_field
        type(gc_cylindrical_zero_potential_t), intent(in) :: local_potential
        type(gc_cylindrical_physical_return_options_t), intent(in) :: local_options

        type(gc_cylindrical_state_t) :: reversed_state
        type(gc_cylindrical_field_sample_t) :: reversed_field
        type(gc_cylindrical_invariants_t) :: reversed_invariants
        type(gc_cylindrical_physical_return_t) :: local_result
        integer :: local_status

        reversed_state = initial
        reversed_state%p_parallel = -initial%p_parallel
        call local_field%evaluate([reversed_state%R, reversed_state%Z, &
            reversed_state%phi], reversed_field, local_status)
        if (local_status /= GC_CYL_SUCCESS) error stop 'reversed field failed'
        call invariants_from_cylindrical_state(reversed_state, reversed_field, 0.0_dp, &
            mass, charge, c_light, reversed_invariants, local_status)
        if (local_status /= GC_CYL_SUCCESS) then
            error stop 'reversed invariants failed'
        end if
        call compute_gc_cylindrical_physical_return(local_field, local_potential, &
            reversed_state, reversed_invariants, mass, charge, c_light, &
            circular_event, local_options, local_result)
        if (local_result%status /= GC_CYL_SUCCESS) then
            error stop 'reversed physical return failed'
        end if
        if (local_result%return_orientation /= -1) then
            error stop 'orientation reversal did not reverse event direction'
        end if
        call require_close('reversed period', local_result%period, &
            expected_circular_period, 3.0e-8_dp)
        call require_close('reversed delta phi', local_result%delta_phi, &
            -2.0_dp*acos(-1.0_dp), 3.0e-9_dp)
    end subroutine check_orientation_reversal

    subroutine check_event_time_convergence(initial, local_invariants, local_options, &
            coarse_result, fine_result)
        type(gc_cylindrical_state_t), intent(in) :: initial
        type(gc_cylindrical_invariants_t), intent(in) :: local_invariants
        type(gc_cylindrical_physical_return_options_t), intent(in) :: local_options
        type(gc_cylindrical_physical_return_t), intent(out) :: coarse_result, fine_result

        type(gc_cylindrical_physical_return_options_t) :: coarse_options, fine_options

        coarse_options = local_options
        fine_options = local_options
        coarse_options%maximum_step = 0.4_dp
        fine_options%maximum_step = 0.1_dp
        call compute_gc_cylindrical_physical_return(field, potential, initial, &
            local_invariants, mass, charge, c_light, circular_event, &
            coarse_options, coarse_result)
        call compute_gc_cylindrical_physical_return(field, potential, initial, &
            local_invariants, mass, charge, c_light, circular_event, &
            fine_options, fine_result)
        if (coarse_result%status /= GC_CYL_SUCCESS .or. &
                fine_result%status /= GC_CYL_SUCCESS) then
            error stop 'event convergence run failed'
        end if
        if (abs(fine_result%period-expected_circular_period) > &
                abs(coarse_result%period-expected_circular_period)) then
            error stop 'event-time refinement did not improve'
        end if
        if (abs(fine_result%period-expected_circular_period) > 3.0e-8_dp) then
            error stop 'event-time refinement is not converged'
        end if
    end subroutine check_event_time_convergence

    subroutine check_wall_before_return(initial, local_invariants, local_options)
        type(gc_cylindrical_state_t), intent(in) :: initial
        type(gc_cylindrical_invariants_t), intent(in) :: local_invariants
        type(gc_cylindrical_physical_return_options_t), intent(in) :: local_options

        type(gc_cylindrical_physical_return_t) :: local_result
        type(z_limit_wall_t) :: local_wall

        local_wall%limit = 0.25_dp
        call compute_gc_cylindrical_physical_return(field, potential, initial, &
            local_invariants, mass, charge, c_light, circular_event, &
            local_options, local_result, wall_model=local_wall)
        if (local_result%status /= GC_CYL_WALL_LOSS) then
            error stop 'wall-before-return status failed'
        end if
        if (local_result%event_kind /= GC_CYL_PHYSICAL_EVENT_WALL) then
            error stop 'wall-before-return event kind failed'
        end if
        if (.not. local_result%wall_hit) error stop 'wall hit was not recorded'
        call require_close('wall event time', local_result%event_time, 1.5_dp, &
            2.0e-8_dp)
    end subroutine check_wall_before_return

    subroutine check_radial_domain_before_return(initial, local_invariants, &
            local_options)
        type(gc_cylindrical_state_t), intent(in) :: initial
        type(gc_cylindrical_invariants_t), intent(in) :: local_invariants
        type(gc_cylindrical_physical_return_options_t), intent(in) :: local_options

        type(gc_cylindrical_physical_return_t) :: local_result

        call compute_gc_cylindrical_physical_return(field, potential, initial, &
            local_invariants, mass, charge, c_light, circular_event, &
            local_options, local_result, radial_domain=radial_domain_event)
        if (local_result%status /= GC_CYL_EQUILIBRIUM_DOMAIN) then
            error stop 'radial-domain status failed'
        end if
        if (local_result%event_kind /= GC_CYL_PHYSICAL_EVENT_RADIAL_DOMAIN) then
            error stop 'radial-domain event kind failed'
        end if
        if (.not. local_result%radial_domain_hit) then
            error stop 'radial-domain hit was not recorded'
        end if
        call require_close('radial event time', local_result%event_time, 1.5_dp, &
            2.0e-8_dp)
    end subroutine check_radial_domain_before_return

    subroutine check_no_return(initial, local_invariants, local_options)
        type(gc_cylindrical_state_t), intent(in) :: initial
        type(gc_cylindrical_invariants_t), intent(in) :: local_invariants
        type(gc_cylindrical_physical_return_options_t), intent(in) :: local_options

        type(gc_cylindrical_physical_return_options_t) :: short_options
        type(gc_cylindrical_physical_return_t) :: local_result

        short_options = local_options
        short_options%maximum_time = 1.0_dp
        short_options%minimum_return_time = 0.1_dp
        call compute_gc_cylindrical_physical_return(field, potential, initial, &
            local_invariants, mass, charge, c_light, circular_event, &
            short_options, local_result)
        if (local_result%status /= GC_CYL_NO_RETURN) then
            error stop 'no-return status failed'
        end if
        if (local_result%physical_return_found) then
            error stop 'no-return fabricated a physical cycle'
        end if
    end subroutine check_no_return

    subroutine circular_event(position, state, sample, user_data, value, status)
        real(dp), intent(in) :: position(3)
        type(gc_cylindrical_state_t), intent(in) :: state
        type(gc_cylindrical_field_sample_t), intent(in) :: sample
        class(*), pointer, intent(inout) :: user_data
        real(dp), intent(out) :: value
        integer, intent(out) :: status

        associate (unused_state => state, unused_sample => sample, &
            unused_user_data => user_data)
        end associate
        value = sin(position(3))
        status = GC_CYL_SUCCESS
    end subroutine circular_event

    subroutine circular_event_rate(position, state, sample, user_data, rate, status)
        real(dp), intent(in) :: position(3)
        type(gc_cylindrical_state_t), intent(in) :: state
        type(gc_cylindrical_field_sample_t), intent(in) :: sample
        class(*), pointer, intent(inout) :: user_data
        real(dp), intent(out) :: rate
        integer, intent(out) :: status

        associate (unused_sample => sample, unused_user_data => user_data)
        end associate
        ! Independent manufactured oracle: for the pure toroidal test field,
        ! dphi/dt=p_parallel/(mass*R), and C=sin(phi).  Production Cdot is
        ! supplied by the generated geometry kernel, never by this test.
        rate = cos(position(3))*state%p_parallel/(mass*position(1))
        status = GC_CYL_SUCCESS
    end subroutine circular_event_rate

    subroutine helical_event(position, state, sample, user_data, value, status)
        real(dp), intent(in) :: position(3)
        type(gc_cylindrical_state_t), intent(in) :: state
        type(gc_cylindrical_field_sample_t), intent(in) :: sample
        class(*), pointer, intent(inout) :: user_data
        real(dp), intent(out) :: value
        integer, intent(out) :: status

        associate (unused_state => state, unused_sample => sample, &
            unused_user_data => user_data)
        end associate
        value = sin(position(3)-2.0_dp*position(2))
        status = GC_CYL_SUCCESS
    end subroutine helical_event

    subroutine radial_domain_event(position, state, sample, user_data, margin, &
            status)
        real(dp), intent(in) :: position(3)
        type(gc_cylindrical_state_t), intent(in) :: state
        type(gc_cylindrical_field_sample_t), intent(in) :: sample
        class(*), pointer, intent(inout) :: user_data
        real(dp), intent(out) :: margin
        integer, intent(out) :: status

        associate (unused_phi => position(3), unused_state => state, &
            unused_sample => sample, unused_user_data => user_data)
        end associate
        margin = 0.25_dp-position(2)
        status = GC_CYL_SUCCESS
    end subroutine radial_domain_event

    subroutine require_close(label, actual, expected, tolerance)
        character(*), intent(in) :: label
        real(dp), intent(in) :: actual, expected, tolerance

        if (abs(actual-expected) > tolerance) then
            write (*, '(a,2(1x,es24.16),1x,es12.4)') trim(label), actual, &
                expected, tolerance
            error stop 'physical return oracle failed'
        end if
    end subroutine require_close

end program test_gc_cylindrical_physical_return
