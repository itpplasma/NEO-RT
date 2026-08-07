program test_gc_generated_resonance_eq17_seams
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_full_fow_frequency_contribution_symbolic, only: &
        evaluate_neort_frequency_root_contribution
    use neort_full_fow_eq17_outer_symbolic, only: &
        evaluate_neort_eq17_outer_factor
    use util_for_test, only: pass_test
    implicit none
    real(dp) :: contribution(2), outer_factor
    real(dp), parameter :: dpsi = 2.0_dp, hm2 = 3.0_dp, tau = 4.0_dp
    real(dp), parameter :: n = 2.0_dp, gprime = 5.0_dp
    real(dp), parameter :: e_ref = 2.0_dp, phi_eff = 0.25_dp, n0 = 3.0_dp
    real(dp), parameter :: temp = 4.0_dp, charge = 2.0_dp, phi = 1.5_dp
    real(dp), parameter :: h = 5.0_dp, residence = 0.75_dp
    real(dp) :: expected, gauge_shifted

    call evaluate_neort_frequency_root_contribution(dpsi, hm2, tau, n, gprime, &
        contribution(1), contribution(2))
    expected = abs(dpsi)*hm2*tau/abs(n*gprime/tau)
    call require_close('frequency contribution', contribution(1), expected)
    call require_close('n-squared full phase contribution', contribution(2), &
        n**2*expected)

    call evaluate_neort_eq17_outer_factor(e_ref, phi_eff, n0, temp, charge, phi, h, &
        residence, outer_factor)
    expected = -(acos(-1.0_dp)**1.5_dp/4.0_dp)*e_ref*phi_eff*n0/ &
        (temp/e_ref)**1.5_dp*exp((charge*phi-h)/temp)*residence
    call require_close('signed Eq17 outer factor', outer_factor, expected)
    gauge_shifted = -(acos(-1.0_dp)**1.5_dp/4.0_dp)*e_ref*phi_eff*n0/ &
        (temp/e_ref)**1.5_dp*exp((charge*(phi+0.7_dp)-(h+charge*0.7_dp))/temp)*residence
    call require_close('Eq17 gauge invariance', gauge_shifted, expected)

    write (*, '(A)') 'test_gc_generated_resonance_eq17_seams OK'
    call pass_test

contains
    subroutine require_close(label, actual, expected)
        character(*), intent(in) :: label
        real(dp), intent(in) :: actual, expected
        if (abs(actual-expected) > 2.0e-12_dp*max(1.0_dp, abs(expected))) then
            write (*, '(A,2(1X,ES24.16))') trim(label), actual, expected
            error stop 'generated resonance/Eq17 seam oracle failed'
        end if
    end subroutine require_close
end program test_gc_generated_resonance_eq17_seams
