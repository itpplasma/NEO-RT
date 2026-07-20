program test_potato_invariants
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use potato_invariants_mod, only: POTATO_INVARIANT_INVALID_INPUT, &
        POTATO_INVARIANT_NONZERO_POTENTIAL, POTATO_INVARIANT_SUCCESS, &
        potato_invariant_t, potato_invariant_from_local_state
    implicit none

    type(potato_invariant_t) :: invariant

    call potato_invariant_from_local_state(2.0_dp, -0.5_dp, 4.0_dp, 3.0_dp, &
        7.0_dp, 0.0_dp, 0.25_dp, 1.0_dp, 9.0_dp, 2.5e7_dp, invariant)
    call require(invariant%status == POTATO_INVARIANT_SUCCESS, 'valid status')
    call require_close(invariant%h0, 4.0_dp, 1.0e-14_dp, 'H0 normalization')
    call require_close(invariant%j_perp, 0.75_dp, 1.0e-14_dp, &
        'J_perp normalization')
    call require_close(invariant%psi_star, 6.25_dp, 1.0e-14_dp, &
        'psi_star normalization')
    call require(invariant%sigma == -1, 'parallel sign')
    call require_close(invariant%psi_axis, 1.0_dp, 1.0e-14_dp, 'axis flux')
    call require_close(invariant%psi_edge, 9.0_dp, 1.0e-14_dp, 'edge flux')
    call require_close(invariant%v0, 2.5e7_dp, 1.0e-7_dp, 'reference velocity')

    call potato_invariant_from_local_state(1.0_dp, 0.3_dp, 2.0_dp, 1.0_dp, &
        0.5_dp, 1.0e-3_dp, 0.1_dp, 0.0_dp, 1.0_dp, 1.0_dp, invariant)
    call require(invariant%status == POTATO_INVARIANT_NONZERO_POTENTIAL, &
        'electrostatic potential rejection')
    call require_close(invariant%h0, 1.001_dp, 1.0e-14_dp, &
        'POTATO total energy includes potential')

    call potato_invariant_from_local_state(1.0_dp, 1.1_dp, 2.0_dp, 1.0_dp, &
        0.5_dp, 0.0_dp, 0.1_dp, 0.0_dp, 1.0_dp, 1.0_dp, invariant)
    call require(invariant%status == POTATO_INVARIANT_INVALID_INPUT, &
        'invalid pitch rejection')

contains

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) then
            write (*, '(A)') 'FAILED: '//message
            error stop 1
        end if
    end subroutine require

    subroutine require_close(actual, expected, tolerance, message)
        real(dp), intent(in) :: actual, expected, tolerance
        character(len=*), intent(in) :: message

        call require(abs(actual - expected) <= tolerance, message)
    end subroutine require_close
end program test_potato_invariants
