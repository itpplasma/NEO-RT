program test_krook
    use iso_fortran_env, only: dp => real64
    use neort_krook, only: krook_rate, integrate_krook, lorentzian
    use fortnum_quadrature, only: gauss_legendre
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp) :: answer(1), expected, nu, a, b, slope, offset, amplitude(2), root
    real(dp) :: cfold, shift
    integer :: model, harmonic, status, trial, direction, placement, calls
    logical :: quadratic, cubic, layer, noisy, wiggle

    do model = 0, 2
        do harmonic = -4, 4
            expected = 3.0_dp
            if (model > 0) expected = expected/8.0_dp
            if (model == 1) expected = expected*(1.0_dp + harmonic**2/4.0_dp)
            call check(krook_rate(3.0_dp, 2.0_dp, real(harmonic, dp), model, 1.0_dp), &
                       expected, 2.0e-15_dp, 'MARS Krook rate')
        end do
    end do
    call check(krook_rate(2.0_dp, 0.5_dp, -2.25_dp, 1, 3.0_dp), &
               6.0_dp*8.0_dp*(1.0_dp + 2.25_dp**2/4.0_dp), &
               2.0e-15_dp, 'canonical harmonic and geometry factor')
    call check(lorentzian(0.0_dp, 2.0_dp), 1.0_dp/(2.0_dp*pi), &
               2.0e-15_dp, 'Lorentzian normalization')

    quadratic = .false.
    cubic = .false.
    layer = .false.
    noisy = .false.
    wiggle = .false.
    calls = 0
    amplitude = [2.0_dp, 0.2_dp]
    do direction = -1, 1, 2
        slope = 2.0_dp*direction
        do placement = 1, 5
            a = -0.7_dp
            b = 1.3_dp
            root = -2.0_dp + real(placement - 1, dp)
            if (placement == 2) root = a
            if (placement == 4) root = b
            offset = -slope*root
            do trial = 0, 5
                nu = 10.0_dp**(3 - 3*trial)
                expected = affine_oracle()
                call integrate_krook(residual, weight, a, b, nu, 12, &
                                     1.0e-8_dp, 0.0_dp, answer, status)
                if (status /= 0) error stop 'affine quadrature did not converge'
         call check(answer(1), expected, 3.0e-8_dp, 'affine amplitude and finite tails')
            end do
        end do
    end do
    ! At finite width, a domain with no resonance still has nonzero tails.
    a = 1.0_dp
    b = 2.0_dp
    nu = 1.0e-12_dp
    slope = 1.0_dp
    offset = 0.0_dp
    amplitude = [1.0_dp, 0.0_dp]
    expected = atan(nu*(b - a)/(nu**2 + a*b))/pi
    call integrate_krook(residual, weight, a, b, nu, 8, &
                         1.0e-8_dp, 1.0e-30_dp, answer, status)
    if (status /= 0) error stop 'tail quadrature did not converge'
    call check(answer(1), expected, 2.0e-8_dp, 'nonresonant finite-width tails')
    quadratic = .true.
    a = -2.0_dp
    b = 3.0_dp
    nu = 1.0e-6_dp
    shift = 0.0_dp
    cfold = 1.0_dp
    call check_fold(0.999999976796914125621098927688_dp)
    cfold = 0.0_dp
    call check_fold(707.106781169354860795238760524_dp)
    cfold = -1.0_dp
    call check_fold(4.86408460164014427167751549034e-7_dp)
    a = 2.0_dp
    b = 3.0_dp
    cfold = 1.0_dp
    nu = 1.0e-12_dp
    call check_fold(1.41543036317840934490040509795e-14_dp)
    a = -2.0_dp
    b = 3.0_dp
    cfold = 1.0e-6_dp
    nu = 1.0e-9_dp
    shift = 0.371_dp
    call check_fold(999.999625000259637633569136827237_dp)
    quadratic = .false.
    slope = 0.0_dp
    offset = 2.0_dp
    nu = 0.3_dp
    expected = (b - a)*nu/(pi*(offset**2 + nu**2))
    call integrate_krook(residual, weight, a, b, nu, 2, &
                         1.0e-8_dp, 0.0_dp, answer, status)
    if (status /= 0) error stop 'constant residual quadrature did not converge'
    call check(answer(1), expected, 2.0e-8_dp, 'constant detuning tails')
    cubic = .true.
    a = 0.0_dp
    b = 1.0_dp
    nu = 1.0e-8_dp
    call check_fold(6666.59252747415932782976684048705691_dp)
    nu = 1.0e-10_dp
    call check_fold(6666.66665859435780367541546164747952_dp)
    cubic = .false.
    call check_boundary_layer()
    call check_negligible_noise()
    call check_hidden_extrema()
    print *, 'test_krook PASSED'

contains

    subroutine residual(eta, res, dres)
        real(dp), intent(in) :: eta
        real(dp), intent(out) :: res, dres
        real(dp) :: bump
        if (wiggle) then
            ! Decreasing detuning with a narrow bump near the upper endpoint,
            ! like the spline wiggle of R next to the trapped-passing
            ! boundary: R' changes sign twice inside one seed cell, so the
            ! seed scan cannot see these extrema.
            res = offset + slope*eta
            dres = slope
            bump = (eta - 0.999_dp)/0.004_dp
            ! Below exp(-26**2) the bump is beneath rounding (no underflow).
            if (abs(bump) < 26.0_dp) then
                res = res + 0.05_dp*exp(-bump**2)
                dres = dres - 0.05_dp*2.0_dp*bump/0.004_dp*exp(-bump**2)
            end if
        else if (cubic) then
            res = (eta - 0.1_dp)**3 - 0.0003_dp*(eta - 0.1_dp)
            dres = 3.0_dp*(eta - 0.1_dp)**2 - 0.0003_dp
        else if (quadratic) then
            res = (eta - shift)**2 - cfold
            dres = 2.0_dp*(eta - shift)
        else
            res = slope*eta + offset
            dres = slope
        end if
    end subroutine residual

    subroutine weight(eta, values)
        real(dp), intent(in) :: eta
        real(dp), intent(out) :: values(:)
        calls = calls + 1
        if (layer) then
            ! Integrable endpoint layer of width 1e-8, like the bounce weight
            ! at the trapped-passing boundary eta_tp(1 -+ 1e-8).
            values = 1.0_dp/sqrt(b + 1.0e-8_dp*(b - a) - eta)
        else if (noisy) then
            ! Decaying weight with 10% deterministic jitter where it is below
            ! 1e-21 of its peak, a stand-in for ODE noise in a region that
            ! cannot affect the integral.
            values = exp(-eta/0.01_dp)
            if (eta > 0.5_dp .and. modulo(aint(eta*2.0_dp**40), 2.0_dp) == 1.0_dp) &
                values = 1.1_dp*values
        else
            values = amplitude(1) + amplitude(2)*eta
        end if
    end subroutine weight

    subroutine check_boundary_layer()
        ! Constant detuning inside the width selects the direct pitch rule.
        real(dp) :: delta, reference
        layer = .true.
        slope = 0.0_dp
        offset = 0.5_dp
        nu = 2.0_dp
        a = 0.0_dp
        b = 1.0_dp
        delta = 1.0e-8_dp*(b - a)
        reference = 2.0_dp*(sqrt(b - a + delta) - sqrt(delta)) &
                    *nu/(pi*(offset**2 + nu**2))
        call integrate_krook(residual, weight, a, b, nu, 4, &
                             1.0e-10_dp, 0.0_dp, answer, status)
        if (status /= 0) error stop 'endpoint layer quadrature did not converge'
        call check(answer(1), reference, 1.0e-9_dp, 'endpoint layer oracle')
        layer = .false.
    end subroutine check_boundary_layer

    subroutine check_negligible_noise()
        real(dp) :: reference
        noisy = .true.
        slope = 0.0_dp
        offset = 0.5_dp
        nu = 2.0_dp
        a = 0.0_dp
        b = 1.0_dp
        ! The jitter changes the integral by less than 1e-21 relative.
        reference = 0.01_dp*(1.0_dp - exp(-100.0_dp))*nu/(pi*(offset**2 + nu**2))
        calls = 0
        call integrate_krook(residual, weight, a, b, nu, 4, &
                             1.0e-10_dp, 0.0_dp, answer, status)
        if (status /= 0) error stop 'noisy negligible tail did not converge'
        call check(answer(1), reference, 1.0e-9_dp, 'noisy negligible tail oracle')
        ! Refinement is driven by the error of the whole integral, so the
        ! negligible noisy tail must not be bisected to the depth limit.
        if (calls > 20000) then
            print *, 'weight evaluations:', calls
            error stop 'negligible noisy tail drove the quadrature to its depth limit'
        end if
        noisy = .false.
    end subroutine check_negligible_noise

    subroutine check_hidden_extrema()
        ! A stationary point missed by the seed scan used to stop the
        ! inverse residual map with status 3 (wrong-sign derivative).
        real(dp) :: xg(16), wg(16), reference, eta, res, dres, h
        integer :: panel, k
        integer, parameter :: panels = 4000
        wiggle = .true.
        slope = -5.0_dp
        offset = 10.0_dp
        nu = 1.0e-3_dp
        a = 0.0_dp
        b = 1.0_dp
        amplitude = [1.0_dp, 0.0_dp]
        call gauss_legendre(16, xg, wg)
        reference = 0.0_dp
        h = (b - a)/panels
        do panel = 1, panels
            do k = 1, 16
                eta = a + h*(panel - 0.5_dp + 0.5_dp*xg(k))
                call residual(eta, res, dres)
                reference = reference + 0.5_dp*h*wg(k)*lorentzian(res, nu)
            end do
        end do
        call integrate_krook(residual, weight, a, b, nu, 4, &
                             1.0e-10_dp, 0.0_dp, answer, status)
        if (status /= 0) then
            print *, 'status', status
            error stop 'hidden extrema quadrature did not converge'
        end if
        call check(answer(1), reference, 1.0e-9_dp, 'hidden extrema oracle')
        wiggle = .false.
    end subroutine check_hidden_extrema

    subroutine check(actual, reference, tolerance, message)
        real(dp), intent(in) :: actual, reference, tolerance
        character(len=*), intent(in) :: message
        if (abs(actual - reference) > tolerance*max(abs(reference), tiny(1.0_dp))) then
            print *, message, actual, reference
            error stop 'Krook oracle mismatch'
        end if
    end subroutine check

    function affine_oracle() result(value)
        real(dp) :: value, ra, rb, angle, constant, linear
        ra = slope*a + offset
        rb = slope*b + offset
        angle = atan2(nu*(rb - ra), nu**2 + ra*rb)
        constant = amplitude(1) - amplitude(2)*offset/slope
        linear = amplitude(2)/slope
        value = (constant*angle + linear*nu/2.0_dp* &
                 log((rb**2 + nu**2)/(ra**2 + nu**2)))/(pi*slope)
    end function affine_oracle

    subroutine check_fold(reference)
        real(dp), intent(in) :: reference
        call integrate_krook(residual, weight, a, b, nu, 2, &
                             1.0e-8_dp, 0.0_dp, answer, status)
        if (status /= 0) error stop 'fold quadrature did not converge'
        call check(answer(1), reference, 1.0e-7_dp, 'finite-width fold oracle')
    end subroutine check_fold
end program test_krook
