module neort_resonance
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_next_after
    use iso_fortran_env, only: dp => real64
    use neort_freq, only: Om_th, Om_ph
    use driftorbit, only: mth, mph, nlev, vth, sign_vpar, etatp
    implicit none

    abstract interface
        subroutine resonance_evaluator(v, eta, residual, derivative)
            import :: dp
            real(dp), intent(in) :: v, eta
            real(dp), intent(out) :: residual, derivative
        end subroutine resonance_evaluator
    end interface

    public :: solve_bracketed_root

contains

    subroutine driftorbit_coarse(v, eta_min, eta_max, roots, nroots)
        real(dp), intent(in) :: v, eta_min, eta_max
        real(dp), intent(out) :: roots(:, :)
        integer, intent(out) :: nroots
        real(dp) :: residual, residual_old, dresdv, dresdeta
        real(dp) :: eta, eta_old, distance, distance_old
        real(dp) :: distance_min, distance_max
        integer :: k, nscan, side
        logical :: valid

        nscan = size(roots, 1)
        roots = 0.0_dp
        nroots = 0
        if (nscan < 2) return
        if (eta_max <= eta_min) return

        ! Resolve the caller's complete physical bracket in the logarithmic
        ! distance from the trapped-passing boundary.  This changes only the
        ! parametrisation of the scan; it does not discard a resonance lane.
        call class_distance_bounds(eta_min, eta_max, distance_min, &
            distance_max, side, valid)
        if (.not. valid) return
        if (distance_max <= distance_min) return

        distance_old = distance_min
        eta_old = eta_from_distance(distance_old, side)
        call resonance_value(v, eta_old, residual_old, dresdv, dresdeta)

        do k = 1, nscan - 1
            distance = exp(log(distance_min) + real(k, dp) * &
                (log(distance_max) - log(distance_min))/real(nscan - 1, dp))
            eta = eta_from_distance(distance, side)
            call resonance_value(v, eta, residual, dresdv, dresdeta)
            if (residual == 0.0_dp .or. residual*residual_old < 0.0_dp) then
                if (nroots < size(roots, 1)) then
                    nroots = nroots + 1
                    roots(nroots, 1) = min(eta_old, eta)
                    roots(nroots, 2) = max(eta_old, eta)
                end if
            end if
            distance_old = distance
            eta_old = eta
            residual_old = residual
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
        use logger, only: warning

        real(dp) :: driftorbit_root(2)
        real(dp), intent(in) :: v, tol, eta_min, eta_max
        real(dp) :: eta_root, derivative_root
        integer :: status
        character(len=1024) :: msg
        character, parameter :: TAB = char(9)
        character, parameter :: LF = char(10)

        call solve_bracketed_root(v, tol, eta_min, eta_max, resonance_eta, &
            eta_root, derivative_root, status)
        driftorbit_root = [eta_root, derivative_root]
        if (status == -1) then
            write (msg, "(*(g0))") &
                "driftorbit_root couldn't bracket 0 for v/vth = ", v/vth, &
                LF//TAB//"etamin = ", eta_min, ", etamax = ", eta_max
            call warning(msg)
        else if (status == -2) then
            write (msg, "(*(g0))") &
                "driftorbit_root: did not converge within the safeguarded iteration budget", LF// &
                TAB//"v/vth = ", v/vth, ", mth = ", mth, ", sign_vpar = ", sign_vpar, LF// &
                TAB//"etamin = ", eta_min, ", etamax = ", eta_max, ", requested tol = ", tol
            call warning(msg)
        end if
    end function driftorbit_root

    subroutine solve_bracketed_root(v, tol, eta_min, eta_max, evaluate, eta_root, &
            derivative_root, status)
        real(dp), intent(in) :: v, tol, eta_min, eta_max
        procedure(resonance_evaluator) :: evaluate
        real(dp), intent(out) :: eta_root, derivative_root
        integer, intent(out) :: status
        real(dp) :: residual, residual_left, residual_right
        real(dp) :: derivative, eta, eta_new, eta_left, eta_right
        real(dp) :: derivative_left, derivative_right
        real(dp) :: tol_eff, residual_scale
        real(dp) :: eta_scale, eta_tol, newton_candidate, derivative_scale
        real(dp) :: best_eta, best_residual, best_derivative
        integer, parameter :: maxit = 100
        integer :: k
        logical :: bracketed

        status = -2
        eta_root = -2.0_dp
        derivative_root = 0.0_dp
        eta_left = eta_min
        eta_right = eta_max

        call evaluate(v, eta_left, residual_left, derivative_left)
        call evaluate(v, eta_right, residual_right, derivative_right)
        if (.not. ieee_is_finite(residual_left) .or. .not. ieee_is_finite(residual_right)) return
        residual_scale = max(1.0_dp, abs(residual_left), abs(residual_right))
        tol_eff = max(abs(tol), 64.0_dp*epsilon(1.0_dp)*residual_scale)
        eta_scale = max(abs(eta_left), abs(eta_right), tiny(1.0_dp))
        derivative_scale = 0.0_dp
        if (ieee_is_finite(derivative_left)) derivative_scale = max(derivative_scale, abs(derivative_left))
        if (ieee_is_finite(derivative_right)) derivative_scale = max(derivative_scale, abs(derivative_right))
        tol_eff = max(tol_eff, 8.0_dp * derivative_scale * spacing(eta_scale))
        bracketed = residual_left == 0.0_dp .or. residual_right == 0.0_dp .or. &
            residual_left*residual_right < 0.0_dp
        if (.not. bracketed) then
            status = -1
            eta_root = -1.0_dp
            return
        end if

        if (abs(residual_left) <= tol_eff) then
            eta_root = eta_left
            derivative_root = derivative_left
            status = 1
            return
        end if
        if (abs(residual_right) <= tol_eff) then
            eta_root = eta_right
            derivative_root = derivative_right
            status = 1
            return
        end if

        eta = 0.5_dp*(eta_left + eta_right)
        do k = 1, maxit
            call evaluate(v, eta, residual, derivative)
            if (.not. ieee_is_finite(residual)) return
            eta_root = eta
            derivative_root = derivative
            if (abs(residual) <= tol_eff) then
                status = 1
                return
            end if

            if (residual_left*residual <= 0.0_dp) then
                eta_right = eta
                residual_right = residual
                derivative_right = derivative
            else
                eta_left = eta
                residual_left = residual
                derivative_left = derivative
            end if

            eta_scale = max(abs(eta_left), abs(eta_right), tiny(1.0_dp))
            eta_tol = 2.0_dp*spacing(eta_scale)
            if (eta_right - eta_left <= eta_tol) then
                eta = 0.5_dp*(eta_left + eta_right)
                call evaluate(v, eta, residual, derivative)
                best_eta = eta
                best_residual = residual
                best_derivative = derivative
                if (abs(residual_left) < abs(best_residual)) then
                    best_eta = eta_left
                    best_residual = residual_left
                    best_derivative = derivative_left
                end if
                if (abs(residual_right) < abs(best_residual)) then
                    best_eta = eta_right
                    best_residual = residual_right
                    best_derivative = derivative_right
                end if
                eta_root = best_eta
                derivative_root = best_derivative
                status = 1
                return
            end if

            ! Newton is used only when it remains inside the updated bracket;
            ! bisection is the deterministic fallback for a flat or noisy
            ! derivative.
            eta_new = 0.5_dp*(eta_left + eta_right)
            newton_candidate = eta
            if (ieee_is_finite(derivative)) then
                if (derivative /= 0.0_dp) newton_candidate = eta - residual/derivative
            end if
            if (ieee_is_finite(newton_candidate) .and. &
                    newton_candidate > eta_left .and. newton_candidate < eta_right) then
                eta_new = newton_candidate
            end if
            eta = eta_new
        end do
        eta_root = -2.0_dp
        derivative_root = 0.0_dp
    end subroutine solve_bracketed_root

    subroutine resonance_value(v, eta, residual, dresdv, dresdeta)
        real(dp), intent(in) :: v, eta
        real(dp), intent(out) :: residual, dresdv, dresdeta
        real(dp) :: Omph, dOmphdv, dOmphdeta
        real(dp) :: Omth, dOmthdv, dOmthdeta

        call Om_ph(v, eta, Omph, dOmphdv, dOmphdeta)
        call Om_th(v, eta, Omth, dOmthdv, dOmthdeta)
        residual = mth*Omth + mph*Omph
        dresdv = mth*dOmthdv + mph*dOmphdv
        dresdeta = mth*dOmthdeta + mph*dOmphdeta
    end subroutine resonance_value

    subroutine resonance_eta(v, eta, residual, derivative)
        real(dp), intent(in) :: v, eta
        real(dp), intent(out) :: residual, derivative
        real(dp) :: dresdv

        call resonance_value(v, eta, residual, dresdv, derivative)
    end subroutine resonance_eta

    subroutine class_distance_bounds(eta_min, eta_max, distance_min, &
            distance_max, side, valid)
        real(dp), intent(in) :: eta_min, eta_max
        real(dp), intent(out) :: distance_min, distance_max
        integer, intent(out) :: side
        logical, intent(out) :: valid
        real(dp) :: x_min, x_max

        distance_min = 0.0_dp
        distance_max = 0.0_dp
        side = 0
        valid = .false.
        if (etatp <= 0.0_dp .or. eta_max <= 0.0_dp) return

        x_min = (eta_min - etatp)/etatp
        x_max = (eta_max - etatp)/etatp
        if (x_min == 0.0_dp .and. x_max > 0.0_dp) then
            x_min = (ieee_next_after(etatp, huge(etatp)) - etatp)/etatp
        else if (x_max == 0.0_dp .and. x_min < 0.0_dp) then
            x_max = (ieee_next_after(etatp, 0.0_dp) - etatp)/etatp
        end if
        if (x_min > 0.0_dp .and. x_max > 0.0_dp) then
            side = 1
        else if (x_min < 0.0_dp .and. x_max < 0.0_dp) then
            side = -1
            x_min = -x_min
            x_max = -x_max
        else
            return
        end if
        distance_min = min(x_min, x_max)
        distance_max = max(x_min, x_max)
        valid = ieee_is_finite(distance_min) .and. ieee_is_finite(distance_max) &
            .and. distance_min > 0.0_dp .and. distance_max > distance_min
    end subroutine class_distance_bounds

    function eta_from_distance(distance, side) result(eta)
        real(dp), intent(in) :: distance
        integer, intent(in) :: side
        real(dp) :: eta

        if (side > 0) then
            eta = etatp + distance*etatp
            if (eta <= etatp) eta = ieee_next_after(etatp, huge(etatp))
        else
            eta = etatp - distance*etatp
            if (eta >= etatp) eta = ieee_next_after(etatp, 0.0_dp)
        end if
    end function eta_from_distance

end module neort_resonance
