program test_resonance_cutout
    ! Regression checks for the optional resonance cutout (resonance_cutout).
    ! Without an argument: layer formulas, all three variants and the
    ! transport dispatch on a manufactured fold. With an argument: one invalid
    ! configuration that must stop.
    use iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use neort_resonance_cutout, only: CUTOUT_OFF, CUTOUT_ZERO, CUTOUT_FILL, &
                  resonance_cutout, resonance_cutout_width, validate_resonance_cutout, &
                     reset_resonance_cutout, bounce_pitch_diffusion, fold_layer_width, &
                                      evaluate_cutout_weight => cutout_weight
    implicit none
    character(len=32) :: mode
    ! Manufactured fold: R = fold_curvature (eta - eta_f)**2 / 2 - ..., with
    ! roots at eta_f +- fold_half. Shared with the resonance_override seam.
    real(dp) :: fold_half = 0.0_dp, fold_curvature = 0.0_dp

    call get_command_argument(1, mode)
    select case (trim(mode))
    case ('')
        call check_pitch_diffusion()
        call check_layer_variants()
        call check_singular_weights()
        call check_transport_dispatch()
        print *, 'test_resonance_cutout PASSED'
    case ('invalid-mode')
        resonance_cutout = 3
        call validate_resonance_cutout()
    case ('narrow-width')
        resonance_cutout = CUTOUT_ZERO
        resonance_cutout_width = 0.5_dp
        call validate_resonance_cutout()
    case ('nan-width')
        resonance_cutout = CUTOUT_ZERO
        resonance_cutout_width = ieee_value(0.0_dp, ieee_quiet_nan)
        call validate_resonance_cutout()
    case ('krook-combination')
        call reject_combination(.true.)
    case ('layer-combination')
        call reject_combination(.false.)
    case default
        error stop 'unknown mode'
    end select

contains

    subroutine check_pitch_diffusion()
        ! Independent oracle: the orbit average <(1 - eta B)/B> by direct
        ! theta quadrature with weight 1/|v_par| in B = B0 (1 - eps cos theta).
        real(dp), parameter :: B0 = 2.0e4_dp, eps = 0.1_dp, nu_d = 3.0e2_dp
        real(dp) :: etatp, etadt, eta(12), expected, got
        integer :: k

        etatp = 1.0_dp / (B0 * (1.0_dp + eps))
        etadt = 1.0_dp / (B0 * (1.0_dp - eps))
        do k = 1, 6
            eta(k) = etatp * (0.2_dp + 0.15_dp * (k - 1)) ! passing
            eta(6 + k) = etatp + (etadt - etatp) * (0.1_dp + 0.15_dp * (k - 1)) ! trapped
        end do
        do k = 1, size(eta)
            expected = 2.0_dp * nu_d * eta(k) * orbit_average(eta(k), B0, eps, etatp)
            got = bounce_pitch_diffusion(nu_d, eta(k), etatp, etadt)
            ! The closed form shares the large-aspect-ratio O(eps) modulus
            ! approximation of layer_delta_eta; the bound is 2*eps.
            if (abs(got - expected) > 2.0_dp * eps * abs(expected)) then
                print *, 'eta =', eta(k), ' closed form', got, ' quadrature', expected
                error stop 'bounce_pitch_diffusion disagrees with quadrature'
            end if
        end do
    end subroutine check_pitch_diffusion

    function orbit_average(eta, B0, eps, etatp) result(avg)
        ! Substitution theta = 2 asin(sqrt(kk) sin(x)) for trapped orbits
        ! removes the turning-point singularity; passing orbits use a plain
        ! midpoint rule over one period.
        real(dp), intent(in) :: eta, B0, eps, etatp
        real(dp) :: avg, num, den, th, Bl, w, kk, x
        integer :: j
        integer, parameter :: n = 20000
        real(dp), parameter :: pi = acos(-1.0_dp)

        num = 0.0_dp
        den = 0.0_dp
        kk = (1.0_dp / eta - B0 * (1.0_dp - eps)) / (2.0_dp * B0 * eps)
        do j = 1, n
            if (eta < etatp) then
                th = 2.0_dp * pi * (j - 0.5_dp) / n
                w = 1.0_dp / sqrt(1.0_dp - eta * B0 * (1.0_dp - eps * cos(th)))
            else
                x = (pi / 2.0_dp) * (j - 0.5_dp) / n
                th = 2.0_dp * asin(sqrt(kk) * sin(x))
                ! d theta / sqrt(1 - eta B) is proportional to dx/cos(theta/2)
                w = 1.0_dp / cos(th / 2.0_dp)
            end if
            Bl = B0 * (1.0_dp - eps * cos(th))
            num = num + w * (1.0_dp - eta * Bl) / Bl
            den = den + w
        end do
        avg = num / den
    end function orbit_average

    subroutine check_layer_variants()
        ! Manufactured fold R = A (eta - eta_f)**2 - A c**2: roots at
        ! d = c from the fold with |R_eta| = 2 A c and |R_eta_eta| = 2 A.
        real(dp), parameter :: A = 3.0e9_dp, D = 4.0e-12_dp
        real(dp) :: layer, c, f0, f1, f2, f_wide

        layer = fold_layer_width(D, 2.0_dp * A)
        if (abs(layer - (D / (2.0_dp * A))**0.25_dp) > 1.0e-14_dp * layer) &
            error stop 'fold layer width is not (D/|R_eta_eta|)**(1/4)'
        if (fold_layer_width(0.0_dp, 2.0_dp*A) /= 0.0_dp) &
            error stop 'D = 0 must give no layer'
        if (fold_layer_width(D, 0.0_dp) /= 0.0_dp) &
            error stop 'flat R must give no layer'

        ! Inside: c = layer/2.
        c = 0.5_dp * layer
        f0 = cutout_weight(CUTOUT_OFF, 1.0_dp, 2.0_dp * A * c, 2.0_dp * A, D)
        f1 = cutout_weight(CUTOUT_ZERO, 1.0_dp, 2.0_dp * A * c, 2.0_dp * A, D)
        f2 = cutout_weight(CUTOUT_FILL, 1.0_dp, 2.0_dp * A * c, 2.0_dp * A, D)
        if (f0 /= 1.0_dp / (2.0_dp * A * c)) &
            error stop 'off must leave the weight unchanged'
        if (f1 /= 0.0_dp) error stop 'zero variant must remove a root inside the layer'
        ! Fill is the finite boundary weight 1/(|R_eta_eta| layer).
        if (abs(f2 * (2.0_dp * A * layer) - 1.0_dp) > 1.0e-14_dp) &
            error stop 'fill weight is not the boundary weight'
        ! Outside: c = 1.01 layer, all variants unchanged.
        c = 1.01_dp * layer
        f0 = 1.0_dp / (2.0_dp * A * c)
        f1 = cutout_weight(CUTOUT_ZERO, 1.0_dp, 2.0_dp * A * c, 2.0_dp * A, D)
        f2 = cutout_weight(CUTOUT_FILL, 1.0_dp, 2.0_dp * A * c, 2.0_dp * A, D)
        if (f1 /= f0 .or. f2 /= f0) &
            error stop 'root outside the layer must be unchanged'
        ! Continuity of the fill weight at the boundary.
        c = (1.0_dp - 1.0e-12_dp) * layer
        f2 = cutout_weight(CUTOUT_FILL, 1.0_dp, 2.0_dp * A * c, 2.0_dp * A, D)
        if (abs(f2 * (2.0_dp * A * c) - 1.0_dp) > 1.0e-11_dp) &
            error stop 'fill weight is discontinuous at the layer edge'
        ! Wider cutoff: the root at 1.01 layer is removed for w = 2, and the
        ! fill factor scales as 1/w.
        c = 1.01_dp * layer
        f_wide = cutout_weight(CUTOUT_ZERO, 2.0_dp, 2.0_dp * A * c, 2.0_dp * A, D)
        if (f_wide /= 0.0_dp) &
            error stop 'wider cutoff must remove roots out to w * layer'
        f_wide = cutout_weight(CUTOUT_FILL, 2.0_dp, 2.0_dp * A * c, 2.0_dp * A, D)
        if (abs(f_wide * (2.0_dp * A * 2.0_dp * layer) - 1.0_dp) > 1.0e-14_dp) &
            error stop 'wide fill weight is not the widened boundary weight'
    end subroutine check_layer_variants

    function cutout_weight(mode, width, slope, curvature, diffusion) result(weight)
        integer, intent(in) :: mode
        real(dp), intent(in) :: width, slope, curvature, diffusion
        real(dp) :: weight
        logical :: inside_layer

        call evaluate_cutout_weight(mode, width, slope, curvature, diffusion, &
                                    weight, inside_layer)
        if (.not. inside_layer) weight = 1.0_dp / slope
    end function cutout_weight

    subroutine check_singular_weights()
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        real(dp), parameter :: curvature = 6.0e9_dp, diffusion = 4.0e-12_dp
        real(dp) :: slope, weight, expected
        integer :: k

        expected = 1.0_dp / (curvature * (diffusion / curvature)**0.25_dp)
        do k = 1, 2
            slope = 0.0_dp
            if (k == 2) slope = tiny(1.0_dp) * epsilon(1.0_dp)
            weight = cutout_weight(CUTOUT_ZERO, 1.0_dp, slope, curvature, diffusion)
            if (weight /= 0.0_dp) error stop 'zero cutout must be finite at a fold'
            weight = cutout_weight(CUTOUT_FILL, 1.0_dp, slope, curvature, diffusion)
            if (.not. ieee_is_finite(weight)) &
                error stop 'fill weight diverged at a fold'
            if (abs(weight - expected) > 1.0e-14_dp * expected) &
                error stop 'fill weight must equal the finite boundary weight at a fold'
        end do
    end subroutine check_singular_weights

    subroutine check_transport_dispatch()
        ! Real geometry (examples/base), manufactured R with a fold in the
        ! trapped interval: verifies that the variants reach the transport sum
        ! and that the default stays bit-identical.
        use neort_lib, only: neort_init, neort_prepare_splines, neort_setup_at_s
        use neort_resonance, only: resonance_override
        use driftorbit, only: mth, comptorque, nopassing, etamin, etamax
        real(dp) :: base(3), off(3), zero(3), fill(3), wide(3), far(3), far_base(3)

        call neort_init('driftorbit.in', 'in_file')
        nopassing = .true.
        call neort_prepare_splines('plasma.in', 'profile.in')
        call neort_setup_at_s(0.5_dp)
        mth = 0
        comptorque = .true.
        resonance_override => manufactured_fold
        ! Roots 1e-9 of the interval apart on a gentle fold: deep inside the
        ! collisional fold layer for any positive collision frequency here.
        fold_half = 1.0e-9_dp * (etamax - etamin)
        fold_curvature = 2.0e3_dp / (etamax - etamin)**2

        call reset_resonance_cutout()
        call run(base)
        call run(off)
        if (any(off /= base)) error stop 'repeat of the default changed transport'
        if (base(1) <= 0.0_dp) error stop 'manufactured fold produced no transport'
        resonance_cutout = CUTOUT_ZERO
        call run(zero)
        if (any(zero /= 0.0_dp)) &
            error stop 'zero variant left a contribution inside the layer'
        resonance_cutout = CUTOUT_FILL
        call run(fill)
        if (.not. (fill(1) > 0.0_dp .and. fill(1) < base(1))) &
            error stop 'fill variant must give 0 < D11 < unregularized D11'
        resonance_cutout_width = 4.0_dp
        call run(wide)
        if (abs(wide(1) - fill(1) / 4.0_dp) > 1.0e-10_dp * fill(1)) &
            error stop 'fill weight must scale as 1/w inside the layer'
        ! Roots far from a sharp fold (layer ~ 1e-3 of the interval times
        ! the near-case layer): every variant equals the default.
        fold_half = 0.25_dp * (etamax - etamin)
        fold_curvature = 2.0e15_dp / (etamax - etamin)**2
        call reset_resonance_cutout()
        call run(far_base)
        resonance_cutout = CUTOUT_ZERO
        resonance_cutout_width = 4.0_dp
        call run(far)
        if (any(far /= far_base)) error stop 'cutout changed roots far from any fold'
        resonance_cutout_width = 1.0_dp
        resonance_cutout = CUTOUT_FILL
        call run(far)
        if (any(far /= far_base)) error stop 'fill changed roots far from any fold'
        resonance_override => null()
        call reset_resonance_cutout()
        print *, 'cutout D11 base/fill/wide-fill:', base(1), fill(1), wide(1)
    end subroutine check_transport_dispatch

    subroutine run(result)
        ! One thermal-speed node of the collisionless transport integral.
        use neort_transport, only: compute_transport_integral
        use neort_profiles, only: vth
        real(dp), intent(out) :: result(3)

        call compute_transport_integral(0.9_dp * vth, 1.1_dp * vth, 1, -1.0_dp, &
                                        result(1:2), result(3))
    end subroutine run

    subroutine manufactured_fold(v, eta, res, dresdeta)
        use driftorbit, only: etamin, etamax
        real(dp), intent(in) :: v, eta
        real(dp), intent(out) :: res, dresdeta
        real(dp) :: eta_f

        if (.not. v > 0.0_dp) error stop 'manufactured fold needs positive speed'
        eta_f = 0.5_dp * (etamin + etamax)
        res = 0.5_dp * fold_curvature * ((eta - eta_f)**2 - fold_half**2)
        dresdeta = fold_curvature * (eta - eta_f)
    end subroutine manufactured_fold

    subroutine reject_combination(with_krook)
        use neort_config, only: config_t, set_config
        logical, intent(in) :: with_krook
        type(config_t) :: config

        config%s = 0.5_dp
        config%qs = 1.0_dp
        config%ms = 2.0_dp
        config%vth = 4.0e7_dp
        config%mph = 3
        config%krook_collisions = with_krook
        if (with_krook) config%krook_nu = 1.0e3_dp
        config%collisional_layer = .not. with_krook
        config%resonance_cutout = CUTOUT_ZERO
        call set_config(config)
    end subroutine reject_combination

end program test_resonance_cutout
