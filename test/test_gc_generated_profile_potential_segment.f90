program test_gc_generated_profile_potential_segment
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_profile_potential_segment_symbolic, only: &
        evaluate_neort_profile_potential_segment
    use neort_eqdsk_flux_profile_segment_symbolic, only: &
        evaluate_neort_eqdsk_flux_profile_segment
    use util_for_test, only: pass_test
    implicit none
    real(dp) :: result(3)
    real(dp), parameter :: psi0 = 1.0_dp, psi1 = 5.0_dp
    real(dp), parameter :: omega0 = 2.0_dp, omega1 = 6.0_dp
    real(dp), parameter :: omega_constant = 4.0_dp, c_light = 10.0_dp
    real(dp) :: psihat, dpsihat_dstor

    call evaluate_neort_profile_potential_segment(psi0, psi1, omega0, omega1, &
        omega_constant, c_light, result(1), result(2), result(3))
    call require_close('affine profile segment', result(1), 1.6_dp)
    call require_close('endpoint reversal', result(2), -1.6_dp)
    call require_close('constant profile limit', result(3), 1.6_dp)
    call evaluate_neort_eqdsk_flux_profile_segment(0.25_dp, 0.0_dp, 0.5_dp, &
        0.0_dp, 6.0_dp, 2.0_dp, 3.0_dp, psihat, dpsihat_dstor)
    call require_close('explicit s_tor to psihat map', psihat, 0.5_dp)
    call require_close('explicit psihat derivative', dpsihat_dstor, 2.0_dp)
    write (*, '(A)') 'test_gc_generated_profile_potential_segment OK'
    call pass_test

contains
    subroutine require_close(label, actual, expected)
        character(*), intent(in) :: label
        real(dp), intent(in) :: actual, expected
        if (abs(actual-expected) > 2.0e-12_dp*max(1.0_dp, abs(expected))) then
            write (*, '(A,2(1X,ES24.16))') trim(label), actual, expected
            error stop 'profile potential segment oracle failed'
        end if
    end subroutine require_close
end program test_gc_generated_profile_potential_segment
