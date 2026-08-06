program test_trapping_parameter
    !! Oracle for the trapping parameter eps.
    !!
    !! eps parametrises B = B0*(1 - eps*cos(theta)), the large-aspect-ratio
    !! form every trapped-particle formula downstream is derived from.  It is
    !! therefore a property of the field strength on the surface, not of the
    !! angle used to label positions on it.  Reading it off as the m=1 cosine
    !! coefficient made it chart dependent: the same TC24 surface reported
    !! 0.0839 in Boozer coordinates and 0.0916 in the direct GEQDSK chart, a
    !! 9% error that entered only the trapped branch.
    !!
    !! The test reparametrises the poloidal angle and requires the reported
    !! eps to be unchanged, which the harmonic definition cannot satisfy.
    use, intrinsic :: iso_fortran_env, only: dp => real64, output_unit
    use math_constants, only: pi

    implicit none

    real(dp), parameter :: B0_ref = 5.0_dp, eps_ref = 0.12_dp
    real(dp), parameter :: tol = 1.0e-6_dp
    integer, parameter :: nth = 4096

    real(dp) :: eps_uniform, eps_skewed, eps_harmonic_uniform, eps_harmonic_skewed

    call sample(0.0_dp, eps_uniform, eps_harmonic_uniform)
    ! A monotonic, 2*pi-periodic relabelling of the same surface: the field
    ! strength at each physical point is untouched, only its angular label
    ! moves, exactly as it does between Boozer and geometric angles.
    call sample(0.45_dp, eps_skewed, eps_harmonic_skewed)

    write (output_unit, '(a,f10.6,a,f10.6)') &
        'extrema-based eps: uniform ', eps_uniform, '  relabelled ', eps_skewed
    write (output_unit, '(a,f10.6,a,f10.6)') &
        'harmonic eps     : uniform ', eps_harmonic_uniform, &
        '  relabelled ', eps_harmonic_skewed

    if (abs(eps_uniform - eps_ref) > tol) then
        write (output_unit, '(a)') 'FAIL: extrema-based eps does not recover the model value'
        error stop
    end if
    if (abs(eps_skewed - eps_uniform) > tol) then
        write (output_unit, '(a,es12.4)') &
            'FAIL: eps changed under a pure relabelling of the poloidal angle: ', &
            abs(eps_skewed - eps_uniform)
        error stop
    end if
    ! Guard the regression: the old definition must actually have been broken,
    ! otherwise this test proves nothing.
    if (abs(eps_harmonic_skewed - eps_harmonic_uniform) < 1.0e-3_dp) then
        write (output_unit, '(a)') &
            'FAIL: the relabelling is too weak to distinguish the definitions'
        error stop
    end if

contains

    subroutine sample(skew, eps_extrema, eps_harmonic)
        !! Walk a surface whose field strength is B0*(1 - eps_ref*cos(theta_true))
        !! while labelling points by theta = theta_true - skew*sin(theta_true).
        real(dp), intent(in) :: skew
        real(dp), intent(out) :: eps_extrema, eps_harmonic

        real(dp) :: theta, theta_true, bmod, dth, bmin, bmax, b_mean, b_cos
        integer :: k

        dth = 2.0_dp*pi/real(nth, dp)
        bmin = huge(1.0_dp)
        bmax = -huge(1.0_dp)
        b_mean = 0.0_dp
        b_cos = 0.0_dp
        do k = 1, nth
            theta = real(k - 1, dp)*dth
            theta_true = invert_label(theta, skew)
            bmod = B0_ref*(1.0_dp - eps_ref*cos(theta_true))
            bmin = min(bmin, bmod)
            bmax = max(bmax, bmod)
            b_mean = b_mean + bmod*dth
            b_cos = b_cos - cos(theta)*bmod*dth
        end do
        b_mean = b_mean/(2.0_dp*pi)
        eps_extrema = (bmax - bmin)/(bmax + bmin)
        eps_harmonic = b_cos/(b_mean*pi)
    end subroutine sample

    real(dp) function invert_label(theta, skew)
        !! Solve theta = theta_true - skew*sin(theta_true) for theta_true.
        real(dp), intent(in) :: theta, skew
        real(dp) :: guess, residual
        integer :: iter

        guess = theta
        do iter = 1, 60
            residual = guess - skew*sin(guess) - theta
            guess = guess - residual/(1.0_dp - skew*cos(guess))
            if (abs(residual) < 1.0e-14_dp) exit
        end do
        invert_label = guess
    end function invert_label

end program test_trapping_parameter
