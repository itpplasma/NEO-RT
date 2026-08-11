module neort_resonance
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use iso_fortran_env, only: dp => real64
    use neort_freq, only: Om_th, Om_ph
    use driftorbit, only: mth, mph, nlev, vth, sign_vpar, etatp
    implicit none

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
        real(dp) :: residual, residual_left, residual_right, residual_old
        real(dp) :: derivative, eta, eta_new, eta_left, eta_right
        real(dp) :: derivative_left, derivative_right, dresdv_tmp
        real(dp) :: tol_eff, residual_scale
        integer :: maxit, k, state
        logical :: bracketed
        character(len=1024) :: msg
        character, parameter :: TAB = char(9)
        character, parameter :: LF = char(10)

        maxit = 100
        state = -2
        residual_old = 0.0_dp
        driftorbit_root = [ -2.0_dp, 0.0_dp ]
        eta_left = eta_min
        eta_right = eta_max

        call resonance_value(v, eta_left, residual_left, dresdv_tmp, derivative_left)
        call resonance_value(v, eta_right, residual_right, dresdv_tmp, derivative_right)
        residual_scale = max(1.0_dp, abs(residual_left), abs(residual_right))
        tol_eff = max(abs(tol), 64.0_dp*epsilon(1.0_dp)*residual_scale)
        bracketed = residual_left == 0.0_dp .or. residual_right == 0.0_dp .or. &
            residual_left*residual_right < 0.0_dp
        if (.not. bracketed) then
            write (msg, "(a,g0,a,g0,a,g0)") &
                "driftorbit_root couldn't bracket 0 for v/vth = ", v/vth, &
                LF//TAB//"etamin = ", eta_left, ", etamax = ", eta_right
            call warning(msg)
            driftorbit_root = [ -1.0_dp, 0.0_dp ]
            return
        end if

        if (abs(residual_left) <= tol_eff) then
            driftorbit_root = [eta_left, derivative_left]
            return
        end if
        if (abs(residual_right) <= tol_eff) then
            driftorbit_root = [eta_right, derivative_right]
            return
        end if

        eta = 0.5_dp*(eta_left + eta_right)
        do k = 1, maxit
            residual_old = residual
            call resonance_value(v, eta, residual, dresdv_tmp, derivative)
            driftorbit_root(1) = eta
            if (abs(residual) <= tol_eff) then
                state = 1
                driftorbit_root(2) = derivative
                exit
            end if

            eta_new = eta
            if (ieee_is_finite(derivative)) then
                if (derivative /= 0.0_dp) then
                    eta_new = eta - residual/derivative
                end if
            end if
            if (.not. ieee_is_finite(eta_new) .or. eta_new <= eta_left .or. &
                    eta_new >= eta_right) then
                eta_new = 0.5_dp*(eta_left + eta_right)
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
            eta = eta_new
        end do

        if (state < 0) then
            driftorbit_root(1) = -2.0_dp
            driftorbit_root(2) = 0.0_dp
            write (msg, "(*(g0))") &
                "driftorbit_root: did not converge within ", maxit, " iterations"//LF// &
                TAB//"v/vth = ", v/vth, ", mth = ", mth, ", sign_vpar = ", sign_vpar, LF// &
                TAB//"etamin = ", eta_min, ", etamax = ", eta_max, ", eta = ", eta, LF// &
                TAB//"resmin = ", residual_left, ", resmax = ", residual_right, &
                ", res = ", residual, LF//TAB//"resold = ", residual_old, &
                ", requested tol = ", tol, ", effective tol = ", tol_eff
            call warning(msg)
        end if
    end function driftorbit_root

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
        if (eta_min <= 0.0_dp .or. eta_max <= 0.0_dp) return

        x_min = (eta_min - etatp)/etatp
        x_max = (eta_max - etatp)/etatp
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
            eta = etatp*(1.0_dp + distance)
        else
            eta = etatp*(1.0_dp - distance)
        end if
    end function eta_from_distance

end module neort_resonance
