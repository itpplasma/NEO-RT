module test_uniform_gc_field
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_gc_dynamics, only: gc_field_sample_t
    use neort_gc_models, only: GC_MODEL_SUCCESS, gc_field_t

    implicit none

    type, extends(gc_field_t), public :: uniform_field_t
    contains
        procedure :: evaluate => evaluate_uniform_field
    end type uniform_field_t

contains

    subroutine evaluate_uniform_field(self, position, sample, status)
        class(uniform_field_t), intent(in) :: self
        real(dp), intent(in) :: position(3)
        type(gc_field_sample_t), intent(out) :: sample
        integer, intent(out) :: status

        associate (unused_self => self)
        end associate
        sample = gc_field_sample_t()
        if (position(1) <= 0.0_dp .or. position(1) >= 1.0_dp) then
            status = 1
            return
        end if
        sample%bmod = 1.0_dp
        sample%sqrtg = 1.0_dp
        sample%hcon = [0.0_dp, 0.5_dp, 1.0_dp]
        sample%hcov = sample%hcon
        sample%psi = position(1)
        sample%grad_psi = [1.0_dp, 0.0_dp, 0.0_dp]
        status = GC_MODEL_SUCCESS
    end subroutine evaluate_uniform_field

end module test_uniform_gc_field

program test_gc_transport_orbit
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use test_uniform_gc_field, only: uniform_field_t
    use neort_gc_dynamics, only: gc_field_sample_t
    use neort_gc_models, only: GC_MODEL_SUCCESS, &
        gc_invariants_t, gc_zero_potential_t, invariants_from_state
    use neort_gc_orbit_integrator, only: GC_ORBIT_PASSING, &
        GC_ORBIT_SUCCESS, GC_ORBIT_NO_RETURN, gc_orbit_average_t, &
        gc_orbit_options_t, compute_gc_orbit_average, &
        compute_gc_full_orbit_average, normalized_full_hamiltonian_factor

    implicit none

    type(uniform_field_t) :: field
    type(gc_zero_potential_t) :: potential
    type(gc_invariants_t) :: invariants
    type(gc_invariants_t) :: full_invariants
    type(gc_orbit_average_t) :: result
    type(gc_orbit_options_t) :: options
    type(gc_field_sample_t) :: sample_from_state
    real(dp) :: reference_position(3)
    real(dp) :: expected_hamiltonian, expected_period
    complex(dp) :: expected_full_average
    integer :: status

    reference_position = [0.5_dp, 0.0_dp, 0.0_dp]
    call field%evaluate(reference_position, sample_from_state, status)
    call invariants_from_state(sample_from_state, 0.0_dp, 1.0_dp, 0.0_dp, &
        1.0_dp, 0.8_dp, invariants, status)
    if (status /= GC_MODEL_SUCCESS) error stop "uniform invariant setup failed"
    call invariants_from_state(sample_from_state, 0.0_dp, 1.0_dp, 1.0_dp, &
        1.0_dp, 0.8_dp, full_invariants, status)
    if (status /= GC_MODEL_SUCCESS) error stop "full invariant setup failed"
    full_invariants%energy = 5.0_dp
    full_invariants%magnetic_moment = 0.75_dp
    if (abs(normalized_full_hamiltonian_factor(full_invariants, 1.0_dp, &
            2.0_dp, 4.0_dp) - 1.625_dp) > 1.0e-14_dp) then
        error stop "full Hamiltonian potential normalization mismatch"
    end if
    call invariants_from_state(sample_from_state, 0.0_dp, 1.0_dp, 1.0_dp, &
        1.0_dp, 0.8_dp, full_invariants, status)

    expected_period = 2.0_dp*acos(-1.0_dp)/0.8_dp
    call compute_gc_orbit_average(field, potential, invariants, &
        reference_position, 1, 1.0_dp, 1.0_dp, 0.36_dp, GC_ORBIT_PASSING, 1, &
        expected_period, 0.8_dp, 0.4_dp, 0.5_dp, 1, 1, &
        unit_poloidal_perturbation, options, result)
    if (result%status /= GC_ORBIT_SUCCESS) error stop "uniform orbit failed"
    expected_hamiltonian = 2.0_dp - 0.36_dp
    if (abs(result%period - expected_period) > 1.0e-8_dp) then
        error stop "uniform orbit period mismatch"
    end if
    if (abs(result%perturbation_average - &
        cmplx(expected_hamiltonian, 0.0_dp, dp)) > 1.0e-8_dp) then
        error stop "uniform orbit phase average mismatch"
    end if
    if (abs(result%inverse_b_average - 1.0_dp) > 1.0e-10_dp) then
        error stop "uniform inverse-B average mismatch"
    end if
    if (abs(result%b_average - 1.0_dp) > 1.0e-10_dp) then
        error stop "uniform B average mismatch"
    end if

    call compute_gc_full_orbit_average(field, potential, full_invariants, &
        reference_position, 1, 1.0_dp, 1.0_dp, 0.36_dp, GC_ORBIT_PASSING, 1, &
        expected_period, 0.8_dp, 0.4_dp, 0, 1, constant_perturbation, options, result)
    if (result%status /= GC_ORBIT_SUCCESS) error stop "full uniform orbit failed"
    expected_full_average = cmplx(2.0_dp - 0.36_dp, 0.0_dp, dp)
    if (abs(result%period - expected_period) > 1.0e-8_dp) then
        error stop "full uniform orbit period mismatch"
    end if
    if (abs(result%perturbation_average - expected_full_average) > 1.0e-8_dp) then
        error stop "full canonical phase mismatch"
    end if
    if (abs(result%inverse_b_average - 1.0_dp) > 1.0e-10_dp) then
        error stop "full inverse-B average mismatch"
    end if
    if (abs(result%b_average - 1.0_dp) > 1.0e-10_dp) then
        error stop "full B average mismatch"
    end if

    options%max_periods = 1.0e-6_dp
    call compute_gc_full_orbit_average(field, potential, full_invariants, &
        reference_position, 1, 1.0_dp, 1.0_dp, 0.36_dp, GC_ORBIT_PASSING, 1, &
        expected_period, 0.8_dp, 0.4_dp, 0, 1, constant_perturbation, options, result)
    if (result%status /= GC_ORBIT_NO_RETURN) then
        error stop "full no-return status was not preserved"
    end if
    write (*, '(A)') 'test_gc_transport_orbit OK'

contains

    subroutine unit_poloidal_perturbation(position, bmod, amplitude, local_status)
        real(dp), intent(in) :: position(3), bmod
        complex(dp), intent(out) :: amplitude
        integer, intent(out) :: local_status

        amplitude = cmplx(cos(position(3)), sin(position(3)), dp)
        local_status = 0
        if (bmod <= 0.0_dp) local_status = 1
    end subroutine unit_poloidal_perturbation

    subroutine constant_perturbation(position, bmod, amplitude, local_status)
        real(dp), intent(in) :: position(3), bmod
        complex(dp), intent(out) :: amplitude
        integer, intent(out) :: local_status

        associate (unused_position => position)
        end associate
        amplitude = cmplx(1.0_dp, 0.0_dp, dp)
        local_status = 0
        if (bmod <= 0.0_dp) local_status = 1
    end subroutine constant_perturbation

end program test_gc_transport_orbit
