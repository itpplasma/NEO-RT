module test_gc_cylindrical_orbit_fields
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_gc_cylindrical_model, only: GC_CYL_SUCCESS, &
        gc_cylindrical_field_sample_t, gc_cylindrical_field_t, &
        make_gc_cylindrical_field_sample

    implicit none

    type, extends(gc_cylindrical_field_t) :: toroidal_field_t
    contains
        procedure :: evaluate => evaluate_toroidal_field
    end type toroidal_field_t

    type, extends(gc_cylindrical_field_t) :: axis_crossing_field_t
    contains
        procedure :: evaluate => evaluate_axis_crossing_field
    end type axis_crossing_field_t

contains

    subroutine evaluate_toroidal_field(self, position, sample, status)
        class(toroidal_field_t), intent(in) :: self
        real(dp), intent(in) :: position(3)
        type(gc_cylindrical_field_sample_t), intent(out) :: sample
        integer, intent(out) :: status
        real(dp) :: b(3), db(3, 3)

        associate (unused_self => self)
        end associate
        b = [0.0_dp, 6.0_dp/position(1), 0.0_dp]
        db = 0.0_dp
        db(2, 1) = -6.0_dp/position(1)**2
        call make_gc_cylindrical_field_sample(position(1), b, db, 0.0_dp, &
            [0.0_dp, 0.0_dp, 0.0_dp], sample, status)
    end subroutine evaluate_toroidal_field

    subroutine evaluate_axis_crossing_field(self, position, sample, status)
        class(axis_crossing_field_t), intent(in) :: self
        real(dp), intent(in) :: position(3)
        type(gc_cylindrical_field_sample_t), intent(out) :: sample
        integer, intent(out) :: status
        real(dp), parameter :: radius_axis = 10.0_dp, poloidal_scale = 1.0_dp
        real(dp), parameter :: toroidal_field = 10.0_dp
        real(dp) :: radius, height, b(3), db(3, 3), grad_psi(3), psi

        associate (unused_self => self)
        end associate
        radius = position(1)
        height = position(2)
        b = [-poloidal_scale*height/radius, toroidal_field, &
            poloidal_scale*(radius - radius_axis)/radius]
        db = 0.0_dp
        db(1, 1) = poloidal_scale*height/radius**2
        db(1, 3) = -poloidal_scale/radius
        db(3, 1) = poloidal_scale*radius_axis/radius**2
        psi = 0.5_dp*poloidal_scale*((radius - radius_axis)**2 + height**2)
        grad_psi = [poloidal_scale*(radius - radius_axis), 0.0_dp, &
            poloidal_scale*height]
        call make_gc_cylindrical_field_sample(radius, b, db, psi, grad_psi, &
            sample, status)
    end subroutine evaluate_axis_crossing_field

end module test_gc_cylindrical_orbit_fields

program test_gc_cylindrical_orbit
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_gc_cylindrical_model, only: GC_CYL_SECTION_PLANE, GC_CYL_SUCCESS, &
        GC_CYL_SECTION_PHI, gc_cylindrical_field_sample_t, &
        gc_cylindrical_field_t, gc_cylindrical_invariants_t, &
        gc_cylindrical_potential_t, gc_cylindrical_state_t, &
        gc_cylindrical_wall_t, gc_cylindrical_section_t, &
        gc_cylindrical_polygon_wall_t, &
        make_gc_cylindrical_field_sample, gc_cylindrical_zero_potential_t, &
        invariants_from_cylindrical_state
    use neort_gc_cylindrical_dynamics, only: gc_cylindrical_rhs
    use neort_gc_cylindrical_orbit, only: GC_CYL_WALL_LOSS, &
        GC_CYL_NO_RETURN, gc_cylindrical_orbit_options_t, &
        gc_cylindrical_phase_result_t, gc_cylindrical_return_t, &
        compute_gc_cylindrical_phase_average, compute_gc_cylindrical_return, &
        gc_cylindrical_status_is_numerical, gc_cylindrical_status_is_physical_loss
    use test_gc_cylindrical_orbit_fields, only: axis_crossing_field_t, &
        toroidal_field_t
    use util_for_test, only: pass_test

    implicit none

    real(dp), parameter :: c_light = 1.0_dp
    type(toroidal_field_t) :: toroidal_field
    type(axis_crossing_field_t) :: axis_field
    type(gc_cylindrical_zero_potential_t) :: zero_potential
    type(gc_cylindrical_polygon_wall_t) :: wall
    type(gc_cylindrical_invariants_t) :: invariants
    type(gc_cylindrical_orbit_options_t) :: options
    type(gc_cylindrical_return_t) :: result
    type(gc_cylindrical_phase_result_t) :: phase_result
    type(gc_cylindrical_state_t) :: initial
    type(gc_cylindrical_field_sample_t) :: field_sample
    type(gc_cylindrical_section_t) :: section
    integer :: status

    call check_toroidal_return_and_phase()
    call check_axis_crossing()
    call check_wall_status_seam()
    call check_no_return_status_seam()
    call pass_test

contains

    subroutine check_toroidal_return_and_phase()
        real(dp), parameter :: radius = 2.0_dp, field_strength = 3.0_dp
        real(dp), parameter :: expected_period = 4.0_dp*acos(-1.0_dp)

        initial = gc_cylindrical_state_t()
        initial%R = radius
        initial%p_parallel = 1.0_dp
        call toroidal_field%evaluate([initial%R, initial%Z, initial%phi], &
            field_sample, status)
        call invariants_from_cylindrical_state(initial, field_sample, 0.0_dp, &
            1.0_dp, 1.0_dp, c_light, invariants, status)
        options = gc_cylindrical_orbit_options_t()
        section = gc_cylindrical_section_t()
        section%kind = GC_CYL_SECTION_PHI
        section%point = [radius, 0.0_dp, 0.0_dp]
        section%winding = 1
        section%frequency_winding = 1
        section%fieldline_q = 1.0_dp
        options%maximum_time = 20.0_dp
        call compute_gc_cylindrical_return(toroidal_field, zero_potential, &
            invariants, [radius, 0.0_dp, 0.0_dp], 1, 1.0_dp, 1.0_dp, c_light, &
            section, options, result)
        if (result%status /= GC_CYL_SUCCESS) error stop 'toroidal return failed'
        call require_close('toroidal return period', result%period, &
            expected_period, 2.0e-9_dp)
        call require_close('toroidal return delta phi', result%delta_phi, &
            2.0_dp*acos(-1.0_dp), 2.0e-10_dp)
        call require_close('toroidal precession frequency', result%omega_phi, &
            0.5_dp, 2.0e-10_dp)
        call require_close('passing bounce frequency', result%omega_b, 0.5_dp, &
            2.0e-10_dp)
        call require_close('passing precession frequency', result%omega_prec, &
            0.0_dp, 2.0e-10_dp)
        if (result%delta_phi <= 1.5_dp*acos(-1.0_dp)) then
            error stop 'section phi was wrapped instead of unwrapped'
        end if

        call compute_gc_cylindrical_phase_average(toroidal_field, zero_potential, &
            invariants, [radius, 0.0_dp, 0.0_dp], 1, 1.0_dp, 1.0_dp, c_light, &
            0.0_dp, result%omega_phi, 0, 1, constant_perturbation, options, &
            section, phase_result)
        if (phase_result%status /= GC_CYL_SUCCESS) error stop 'phase return failed'
        call require_close('temporal phase average real', &
            real(phase_result%perturbation_average, dp), 1.0_dp, 2.0e-8_dp)
        call require_close('temporal phase average imaginary', &
            aimag(phase_result%perturbation_average), 0.0_dp, 2.0e-8_dp)
        associate (unused_field_strength => field_strength)
        end associate
    end subroutine check_toroidal_return_and_phase

    subroutine check_axis_crossing()
        real(dp), parameter :: radius_axis = 10.0_dp

        initial = gc_cylindrical_state_t()
        initial%R = radius_axis + 0.001_dp
        initial%Z = 0.05_dp
        initial%p_parallel = 1.0_dp
        call axis_field%evaluate([initial%R, initial%Z, initial%phi], &
            field_sample, status)
        call invariants_from_cylindrical_state(initial, field_sample, 0.0_dp, &
            1.0_dp, 1.0_dp, c_light, invariants, status)
        options = gc_cylindrical_orbit_options_t()
        section = gc_cylindrical_section_t()
        section%kind = GC_CYL_SECTION_PLANE
        section%point = [radius_axis, 0.0_dp, 0.0_dp]
        section%normal = [1.0_dp, 0.0_dp]
        section%direction = -1
        options%maximum_time = 4.0_dp
        call compute_gc_cylindrical_return(axis_field, zero_potential, invariants, &
            [initial%R, initial%Z, initial%phi], 1, 1.0_dp, 1.0_dp, c_light, &
            section, options, result)
        if (result%status /= GC_CYL_SUCCESS) error stop 'axis-crossing return failed'
        call require_close('axis-crossing section', result%state_at_return%R, &
            radius_axis, 1.0e-7_dp)
    end subroutine check_axis_crossing

    subroutine check_wall_status_seam()
        initial = gc_cylindrical_state_t()
        initial%R = 12.0_dp
        initial%p_parallel = 1.0_dp
        call toroidal_field%evaluate([initial%R, initial%Z, initial%phi], &
            field_sample, status)
        call invariants_from_cylindrical_state(initial, field_sample, 0.0_dp, &
            1.0_dp, 1.0_dp, c_light, invariants, status)
        call wall%set_vertices(reshape([9.0_dp, -1.0_dp, 11.0_dp, -1.0_dp, &
            11.0_dp, 1.0_dp, 9.0_dp, 1.0_dp], [2, 4]), status)
        if (status /= GC_CYL_SUCCESS) error stop 'polygon wall setup failed'
        options = gc_cylindrical_orbit_options_t()
        section = gc_cylindrical_section_t()
        section%kind = GC_CYL_SECTION_PHI
        section%point = [initial%R, initial%Z, initial%phi]
        section%winding = 1
        call compute_gc_cylindrical_return(toroidal_field, zero_potential, &
            invariants, [initial%R, initial%Z, initial%phi], 1, 1.0_dp, 1.0_dp, &
            c_light, section, options, result, wall)
        if (result%status /= GC_CYL_WALL_LOSS) error stop 'wall loss was not physical'
        if (.not. gc_cylindrical_status_is_physical_loss(result%status)) then
            error stop 'wall loss crossed the numerical seam'
        end if
        if (gc_cylindrical_status_is_numerical(result%status)) then
            error stop 'wall loss was classified as numerical'
        end if
    end subroutine check_wall_status_seam

    subroutine check_no_return_status_seam()
        initial = gc_cylindrical_state_t()
        initial%R = 2.0_dp
        initial%p_parallel = 1.0_dp
        call toroidal_field%evaluate([initial%R, initial%Z, initial%phi], &
            field_sample, status)
        call invariants_from_cylindrical_state(initial, field_sample, 0.0_dp, &
            1.0_dp, 1.0_dp, c_light, invariants, status)
        options = gc_cylindrical_orbit_options_t()
        section = gc_cylindrical_section_t()
        section%kind = GC_CYL_SECTION_PHI
        section%point = [initial%R, initial%Z, initial%phi]
        section%winding = 1
        options%maximum_time = 0.01_dp
        call compute_gc_cylindrical_return(toroidal_field, zero_potential, &
            invariants, [initial%R, initial%Z, initial%phi], 1, 1.0_dp, 1.0_dp, &
            c_light, section, options, result)
        if (result%status /= GC_CYL_NO_RETURN) error stop 'no-return status changed'
        if (.not. gc_cylindrical_status_is_numerical(result%status)) then
            error stop 'no-return status was not numerical'
        end if
    end subroutine check_no_return_status_seam

    subroutine constant_perturbation(position, sample, amplitude, status)
        real(dp), intent(in) :: position(3)
        type(gc_cylindrical_field_sample_t), intent(in) :: sample
        complex(dp), intent(out) :: amplitude
        integer, intent(out) :: status

        associate (unused_position => position, unused_sample => sample)
        end associate
        amplitude = cmplx(1.0_dp, 0.0_dp, dp)
        status = GC_CYL_SUCCESS
    end subroutine constant_perturbation

    subroutine require_close(label, actual, reference, tolerance)
        character(*), intent(in) :: label
        real(dp), intent(in) :: actual, reference, tolerance

        if (abs(actual - reference) > tolerance) then
            write(*, '(a,2(1x,es24.16),1x,es12.4)') trim(label), actual, reference, tolerance
            error stop 'cylindrical orbit oracle failed'
        end if
    end subroutine require_close

end program test_gc_cylindrical_orbit
