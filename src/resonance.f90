module neort_resonance
    use iso_fortran_env, only: dp => real64
    use neort_freq, only: Om_th, Om_ph, d_Om_ds
    use driftorbit, only: mth, mph, nlev, vth, sign_vpar
    implicit none

contains
    subroutine driftorbit_coarse(v, eta_min, eta_max, roots, nroots)
        real(dp), intent(in) :: v, eta_min, eta_max
        real(dp), intent(out) :: roots(:, :)
        integer, intent(out) :: nroots
        real(dp) :: deta
        real(dp) :: Omph, dOmphdv, dOmphdeta
        real(dp) :: Omth, dOmthdv, dOmthdeta
        real(dp) :: res, dresdv, dresdeta
        real(dp) :: resold, dresdvold, dresdetaold
        real(dp) :: eta
        integer :: k, ninterv

        ninterv = size(roots, 1)

        deta = (eta_max - eta_min)*1.0_dp/ninterv
        nroots = 0

        do k = 0, ninterv
            eta = eta_min + k*deta
            call Om_th(v, eta, Omth, dOmthdv, dOmthdeta)
            call Om_ph(v, eta, Omph, dOmphdv, dOmphdeta)
            res = mth*Omth + mph*Omph
            dresdv = mth*dOmthdv + mph*dOmphdv
            dresdeta = mth*dOmthdeta + mph*dOmphdeta
            if (k > 0) then
                if (sign(1.0_dp, res) /= sign(1.0_dp, resold)) then
                    nroots = nroots + 1
                    roots(nroots, 1) = eta - deta
                    roots(nroots, 2) = eta
                end if
            end if
            resold = res
            dresdvold = dresdv
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
        use logger, only: warning

        real(dp) :: driftorbit_root(2)
        real(dp), intent(in) :: v, tol, eta_min, eta_max
        real(dp) :: res, res_old, eta0, eta_old
        real(dp) :: Omph, dOmphdv, dOmphdeta
        real(dp) :: Omth, dOmthdv, dOmthdeta
        integer :: maxit, k, state
        real(dp) :: etamin2, etamax2
        logical :: slope_pos
        real(dp) :: resmin, resmax
        real(dp) :: eta

        character(len=1024) :: msg
        character, parameter :: TAB = char(9)
        character, parameter :: LF = char(10)

        maxit = 100
        state = -2
        eta0 = eta
        eta_old = 0.0_dp
        res = 0.0_dp

        etamin2 = eta_min
        etamax2 = eta_max

        eta = etamin2
        call Om_ph(v, eta, Omph, dOmphdv, dOmphdeta)
        call Om_th(v, eta, Omth, dOmthdv, dOmthdeta)
        res = mph*Omph + mth*Omth
        resmin = res

        eta = etamax2
        call Om_ph(v, eta, Omph, dOmphdv, dOmphdeta)
        call Om_th(v, eta, Omth, dOmthdv, dOmthdeta)
        resmax = mph*Omph + mth*Omth
        if (resmax - resmin > 0) then
            slope_pos = .true.
        else
            slope_pos = .false.
        end if

        if (driftorbit_nroot(v, etamin2, etamax2) == 0) then
            write (msg, "(a,g0,a,g0,a,g0)") &
                  "driftorbit_root couldn't bracket 0 for v/vth = ", v / vth, LF // &
                  TAB // "etamin = ", etamin2, ", etamax = ", etamax2
            call warning(msg)
            return
        end if

        do k = 1, maxit
            res_old = res
            call Om_ph(v, eta, Omph, dOmphdv, dOmphdeta)
            call Om_th(v, eta, Omth, dOmthdv, dOmthdeta)
            res = mph*Omph + mth*Omth

            driftorbit_root(1) = eta

            if (abs(res) < tol) then
                state = 1
                driftorbit_root(2) = mph*dOmphdeta + mth*dOmthdeta
                exit
            elseif ((slope_pos .and. res > 0) .or. &
                    ((.not. slope_pos) .and. res < 0)) then
                etamax2 = eta
                eta_old = eta
                eta = (eta + etamin2)/2.0_dp
            else
                etamin2 = eta
                eta_old = eta
                eta = (eta + etamax2)/2.0_dp
            end if
        end do
        if (state < 0) then
            driftorbit_root(2) = mph*dOmphdeta + mth*dOmthdeta
            write (msg, "(a,i0,a,g0,a,i0,a,g0,a,g0,a,g0,a,g0,a,g0,a,g0,a,g0,a,g0,a,g0,a,g0)") &
                  "driftorbit_root: did not converge within ", maxit, " iterations" // LF // &
                  TAB // "v/vth = ", v / vth, ", mth = ", mth, ", sign_vpar = ", sign_vpar, LF // &
                  TAB // "etamin = ", eta_min, ", etamax = ", eta_max, ", eta = ", eta, LF // &
                  TAB // "resmin = ", resmin, ", resmax = ", resmax, ", res = ", res, LF // &
                  TAB // "resold = ", res_old, ", res = ", res, LF // &
                  TAB // "tol = ", tol
            call warning(msg)
        end if
        eta = eta0
    end function driftorbit_root
end module neort_resonance
