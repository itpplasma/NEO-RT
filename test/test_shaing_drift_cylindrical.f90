program test_shaing_drift_cylindrical
    ! Independent cylindrical oracle in the deeply trapped, large-aspect-ratio
    ! limit: B=B_axis*R_axis/R, so v_DZ=sigma_B*c*mu/(qi*R) and
    ! alpha_dot=-q*v_DZ/r. Neither the oracle nor its units use Boozer fluxes.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use do_magfie_mod, only: a, eps, iota, psi_pr, q, R0, s
    use driftorbit, only: B0
    use neort_profiles, only: Om_tE
    use shaing, only: omph_shaing
    use util, only: c, mi, qe, qi

    implicit none

    integer :: sigma_B, sigma_q

    do sigma_B = -1, 1, 2
        do sigma_q = -1, 1, 2
            call check_cylindrical_limit(sigma_B, sigma_q, 1.0_dp, 1000.0_dp, 2.0_dp)
            call check_cylindrical_limit(sigma_B, sigma_q, 2.0_dp, 2000.0_dp, 4.0_dp)
        end do
    end do

contains

    subroutine check_cylindrical_limit(sigma_B, sigma_q, radius_a, radius_R, q_abs)
        integer, intent(in) :: sigma_B, sigma_q
        real(dp), intent(in) :: radius_a, radius_R, q_abs

        real(dp), parameter :: speed = 1.0e8_dp
        real(dp) :: radius, eta, magnetic_moment, expected, actual, relative_error

        s = 0.25_dp
        a = radius_a
        R0 = radius_R
        radius = a*sqrt(s)
        eps = radius/R0
        q = real(sigma_q, dp)*q_abs
        iota = 1.0_dp/q
        B0 = 2.0e4_dp
        qi = qe
        psi_pr = real(sigma_B, dp)*B0*a**2/2.0_dp
        Om_tE = 0.0_dp
        eta = 1.0_dp/(B0*(1.0_dp - eps))
        magnetic_moment = 0.5_dp*mi*eta*speed**2
        expected = -q*real(sigma_B, dp)*c*magnetic_moment &
            /(qi*(R0 + radius)*radius)
        actual = omph_shaing(speed, eta)
        if (.not. ieee_is_finite(actual)) error stop "nonfinite Shaing drift"
        relative_error = abs(actual - expected)/abs(expected)
        print *, "cylindrical Shaing oracle:", sigma_B, sigma_q, expected, actual
        if (relative_error > 2.0_dp*eps) error stop "incorrect Shaing drift normalization"
    end subroutine check_cylindrical_limit

end program test_shaing_drift_cylindrical
