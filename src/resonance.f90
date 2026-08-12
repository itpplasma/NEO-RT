module neort_resonance
    use iso_fortran_env, only: dp => real64
    use ieee_arithmetic, only: ieee_is_finite
    use neort_freq, only: Om_th, Om_ph, Om_ph_passing_from_omth
    use driftorbit, only: mth, mph, nlev, vth, sign_vpar, etatp
    implicit none

contains

    subroutine resonance_value(v, eta, res, dresdeta)
        ! R(v, eta) = mth*Omega_theta + mph*Omega_phi.
        ! Om_th and Om_ph already provide the analytic eta derivatives, so
        ! one callback supplies both the Newton residual and its Jacobian.
        real(dp), intent(in) :: v, eta
        real(dp), intent(out) :: res, dresdeta
        real(dp) :: Omph, dOmphdv, dOmphdeta
        real(dp) :: Omth, dOmthdv, dOmthdeta

        call Om_th(v, eta, Omth, dOmthdv, dOmthdeta)
        if (eta <= etatp) then
            ! Om_ph repeats this Om_th call in the passing branch.  Expand
            ! that identity here so each residual evaluation computes the
            ! expensive poloidal frequency only once. The shared completion
            ! routine retains the passing magnetic-drift path exactly.
            call Om_ph_passing_from_omth(v, eta, Omth, dOmthdv, dOmthdeta, &
                Omph, dOmphdv, dOmphdeta)
        else
            call Om_ph(v, eta, Omph, dOmphdv, dOmphdeta)
        end if
        ! Keep the historical/public ordering of the cancellation.  The
        ! resonance oracle evaluates the public frequency path in this order,
        ! and changing it under -ffast-math can move a poorly conditioned root
        ! by several ulps even though the algebra is identical.
        res = mph * Omph + mth * Omth
        dresdeta = mph * dOmphdeta + mth * dOmthdeta
    end subroutine resonance_value

    subroutine driftorbit_coarse(v, eta_min, eta_max, roots, nroots)
        real(dp), intent(in) :: v, eta_min, eta_max
        real(dp), intent(out) :: roots(:, :)
        integer, intent(out) :: nroots
        real(dp) :: deta
        real(dp) :: res, dresdeta
        real(dp) :: resold, dresdetaold
        real(dp) :: eta
        integer :: k, ninterv

        ninterv = size(roots, 1)

        deta = (eta_max - eta_min) / ninterv
        nroots = 0
        roots = 0.0_dp

        do k = 0, ninterv
            eta = eta_min + k * deta
            call resonance_value(v, eta, res, dresdeta)
            if (k > 0) then
                ! Keep exact grid-point roots as zero-width brackets.  Do not
                ! let an exact zero be lost because SIGN(1, 0) is positive.
                if (res == 0.0_dp) then
                    if (resold /= 0.0_dp .and. nroots < ninterv) then
                        nroots = nroots + 1
                        roots(nroots, 1) = eta
                        roots(nroots, 2) = eta
                        roots(nroots, 3) = dresdeta
                    end if
                else if (resold /= 0.0_dp) then
                    if ((res > 0.0_dp .and. resold < 0.0_dp) .or. &
                        (res < 0.0_dp .and. resold > 0.0_dp)) then
                        nroots = nroots + 1
                        roots(nroots, 1) = eta - deta
                        roots(nroots, 2) = eta
                        roots(nroots, 3) = dresdetaold
                    end if
                end if
            else if (res == 0.0_dp) then
                nroots = 1
                roots(1, 1) = eta
                roots(1, 2) = eta
                roots(1, 3) = dresdeta
            end if
            resold = res
            dresdetaold = dresdeta
        end do
    end subroutine driftorbit_coarse

    function driftorbit_nroot(v, eta_min, eta_max)
        integer :: driftorbit_nroot
        real(dp), intent(in) :: v
        real(dp), intent(in) :: eta_min, eta_max
        real(dp) :: roots(nlev, 3)

        call driftorbit_coarse(v, eta_min, eta_max, roots, driftorbit_nroot)
    end function driftorbit_nroot

    function driftorbit_root(v, tol, eta_min, eta_max)
        use logger, only: debug, warning

        ! Bisection can exhaust the representable eta values before a
        ! zero-frequency lane reaches an absolute residual of exactly zero.
        ! Keep the fallback well below the transport quadrature accuracy while
        ! acknowledging the conditioning of the separatrix frequency splines.
        real(dp), parameter :: relative_tolerance_floor = 1.0e-9_dp
        real(dp), parameter :: eta_relative_tolerance = 1.0e-12_dp
        integer, parameter :: maxit = 100
        real(dp) :: driftorbit_root(2)
        real(dp), intent(in) :: v, tol, eta_min, eta_max
        real(dp) :: res, dresdeta, dfa_deta, dfb_deta, resmin, resmax, tol_eff
        real(dp) :: a, b, fa, fb, eta, eta_next, denom
        real(dp) :: eta_scale, eta_tol, width_prev
        integer :: k
        logical :: converged

        character(len=1024) :: msg
        character, parameter :: TAB = char(9)
        character, parameter :: LF = char(10)

        ! The caller supplies one interval returned by driftorbit_coarse.
        ! Re-scanning all nlev intervals here was both redundant and the
        ! dominant cost of resonance refinement.  A sign check is sufficient
        ! for the documented bracket contract; the safeguarded iteration below
        ! never leaves this interval.
        a = min(eta_min, eta_max)
        b = max(eta_min, eta_max)
        call resonance_value(v, a, fa, dfa_deta)
        call resonance_value(v, b, fb, dfb_deta)
        resmin = fa
        resmax = fb
        driftorbit_root = [-1.0_dp, 0.0_dp]

        if (.not. ieee_is_finite(fa) .or. .not. ieee_is_finite(fb)) then
            write (msg, "(a,g0,a,g0,a,g0)") &
                "driftorbit_root couldn't bracket 0 for v/vth = ", v / vth, LF // &
                TAB // "etamin = ", a, ", etamax = ", b
            call warning(msg)
            return
        end if

        tol_eff = max(abs(tol), relative_tolerance_floor * &
            max(1.0_dp, abs(fa), abs(fb)))

        if (fa == 0.0_dp) then
            if (.not. ieee_is_finite(dfa_deta) .or. dfa_deta == 0.0_dp) then
                call warning("driftorbit_root: zero/invalid Jacobian at exact bracket endpoint")
                return
            end if
            driftorbit_root = [a, dfa_deta]
            return
        end if
        if (fb == 0.0_dp) then
            if (.not. ieee_is_finite(dfb_deta) .or. dfb_deta == 0.0_dp) then
                call warning("driftorbit_root: zero/invalid Jacobian at exact bracket endpoint")
                return
            end if
            driftorbit_root = [b, dfb_deta]
            return
        end if

        if ((fa > 0.0_dp .and. fb > 0.0_dp) .or. &
            (fa < 0.0_dp .and. fb < 0.0_dp)) then
            write (msg, "(a,g0,a,g0,a,g0)") &
                "driftorbit_root couldn't bracket 0 for v/vth = ", v / vth, LF // &
                TAB // "etamin = ", a, ", etamax = ", b
            call warning(msg)
            return
        end if

        eta_scale = max(abs(a), abs(b))
        eta_tol = 4.0_dp * epsilon(1.0_dp) * eta_scale

        ! Start with a secant estimate, then use Newton steps whenever the
        ! analytic step remains strictly inside the active bracket.  The
        ! midpoint fallback makes convergence unconditional for continuous
        ! residuals with a valid sign-changing bracket.
        denom = fb - fa
        if (ieee_is_finite(denom) .and. denom /= 0.0_dp) then
            eta = a - fa * (b - a) / denom
        else
            eta = 0.5_dp * (a + b)
        end if
        if (.not. ieee_is_finite(eta) .or. eta <= a .or. eta >= b) then
            eta = 0.5_dp * (a + b)
        end if

        converged = .false.
        do k = 1, maxit
            width_prev = abs(b - a)
            call resonance_value(v, eta, res, dresdeta)
            if (.not. ieee_is_finite(res)) exit
            driftorbit_root(1) = eta
            driftorbit_root(2) = dresdeta

            if (abs(res) <= tol_eff .and. abs(b - a) <= eta_relative_tolerance * eta_scale .and. &
                ieee_is_finite(dresdeta)) then
                converged = .true.
                exit
            end if

            if ((res > 0.0_dp .and. fa > 0.0_dp) .or. &
                (res < 0.0_dp .and. fa < 0.0_dp)) then
                a = eta
                fa = res
            else
                b = eta
                fb = res
            end if

            if (abs(b - a) <= eta_tol) then
                eta = 0.5_dp * (a + b)
                call resonance_value(v, eta, res, dresdeta)
                driftorbit_root = [eta, dresdeta]
                converged = ieee_is_finite(res) .and. ieee_is_finite(dresdeta) .and. dresdeta /= 0.0_dp
                exit
            end if

            ! False-position steps can retain one endpoint indefinitely.  Force
            ! a true bracket contraction before trying another interpolation;
            ! this is the safeguard that makes the hybrid globally convergent.
            if (abs(b - a) > 0.5_dp * width_prev) then
                eta = 0.5_dp * (a + b)
                cycle
            end if

            eta_next = eta
            if (ieee_is_finite(dresdeta) .and. dresdeta /= 0.0_dp) then
                eta_next = eta - res / dresdeta
            end if
            if (.not. ieee_is_finite(eta_next) .or. eta_next <= a .or. eta_next >= b .or. &
                eta_next == eta) then
                denom = fb - fa
                if (ieee_is_finite(denom) .and. denom /= 0.0_dp) then
                    eta_next = a - fa * (b - a) / denom
                else
                    eta_next = 0.5_dp * (a + b)
                end if
            end if
            if (.not. ieee_is_finite(eta_next) .or. eta_next <= a .or. eta_next >= b) then
                eta_next = 0.5_dp * (a + b)
            end if
            eta = eta_next
        end do

        if (.not. converged) then
            driftorbit_root(1) = -2.0_dp
            driftorbit_root(2) = 0.0_dp
            write (msg, "(*(g0))") &
                "driftorbit_root: did not converge within ", maxit, " iterations" // LF // &
                TAB // "v/vth = ", v / vth, ", mth = ", mth, ", sign_vpar = ", sign_vpar, LF // &
                TAB // "etamin = ", a, ", etamax = ", b, ", eta = ", eta, LF // &
                TAB // "resmin = ", resmin, ", resmax = ", resmax, ", res = ", res, LF // &
                TAB // "requested tol = ", tol, ", effective tol = ", tol_eff
            call warning(msg)
        end if
    end function driftorbit_root
end module neort_resonance
