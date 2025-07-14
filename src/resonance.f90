module neort_resonance
    use neort_freq, only: Om_th, Om_ph, d_Om_ds
#ifdef USE_THICK_ORBITS
    use neort_freq, only: Om_th_unified, Om_ph_unified, calculate_orbit_width_criterion
#endif
    use driftorbit, only: mth, mph, nlev, vth, sign_vpar
    implicit none

contains
    subroutine driftorbit_coarse(v, eta_min, eta_max, roots, nroots)
        real(8), intent(in) :: v, eta_min, eta_max
        real(8), intent(out) :: roots(:, :)
        integer, intent(out) :: nroots
        real(8) :: deta
        real(8) :: Omph, dOmphdv, dOmphdeta
        real(8) :: Omth, dOmthdv, dOmthdeta
        real(8) :: res, dresdv, dresdeta
        real(8) :: resold, dresdvold, dresdetaold
        real(8) :: eta
        integer :: k, ninterv

        ninterv = size(roots, 1)

        deta = (eta_max - eta_min)*1d0/ninterv
        nroots = 0

        do k = 0, ninterv
            eta = eta_min + k*deta
            call Om_th(v, eta, Omth, dOmthdv, dOmthdeta)
            call Om_ph(v, eta, Omph, dOmphdv, dOmphdeta)
            res = mth*Omth + mph*Omph
            dresdv = mth*dOmthdv + mph*dOmphdv
            dresdeta = mth*dOmthdeta + mph*dOmphdeta
            if (k > 0) then
                if (sign(1d0, res) /= sign(1d0, resold)) then
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
        real(8), intent(in) :: v
        real(8), intent(in) :: eta_min, eta_max
        real(8) :: roots(nlev, 3)

        call driftorbit_coarse(v, eta_min, eta_max, roots, driftorbit_nroot)
    end function driftorbit_nroot

    function driftorbit_root(v, tol, eta_min, eta_max)
        real(8) :: driftorbit_root(2)
        real(8), intent(in) :: v, tol, eta_min, eta_max
        real(8) :: res, res_old, eta0, eta_old
        real(8) :: Omph, dOmphdv, dOmphdeta
        real(8) :: Omth, dOmthdv, dOmthdeta
        integer :: maxit, k, state
        real(8) :: etamin2, etamax2
        logical :: slope_pos
        real(8) :: resmin, resmax
        real(8) :: eta

        maxit = 100
        state = -2
        eta0 = eta
        eta_old = 0d0
        res = 0d0

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
            print *, "ERROR: driftorbit_root couldn't bracket 0 for v/vth = ", v/vth
            print *, "ERROR: etamin = ", etamin2, " etamax = ", etamax2

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
                eta = (eta + etamin2)/2d0
            else
                etamin2 = eta
                eta_old = eta
                eta = (eta + etamax2)/2d0
            end if
        end do
        if (state < 0) then
            driftorbit_root(2) = mph*dOmphdeta + mth*dOmthdeta
            print *, "ERROR: driftorbit_root did not converge in 100 iterations"
            print *, "v/vth  = ", v/vth, "mth    = ", mth, "sign_vpar= ", sign_vpar
            print *, "etamin = ", eta_min, "etamax = ", eta_max, "eta = ", eta
            print *, "resmin = ", resmin, "resmax = ", resmax, "res = ", res
            print *, "resold = ", res_old, "res    = ", res
            print *, "tol = ", tol
        end if
        eta = eta0
    end function driftorbit_root
    
    subroutine driftorbit_coarse_unified(v, eta_min, eta_max, roots, nroots)
        ! Unified resonance finder with thick/thin orbit dispatch
        real(8), intent(in) :: v, eta_min, eta_max
        real(8), intent(out) :: roots(:, :)
        integer, intent(out) :: nroots
        
        logical, parameter :: enable_thick_orbit_resonance = .true.
        real(8), parameter :: orbit_width_threshold = 0.01d0
        
        real(8) :: orbit_width
        
#ifdef USE_THICK_ORBITS
        if (enable_thick_orbit_resonance) then
            orbit_width = calculate_orbit_width_criterion(v, 0.5d0 * (eta_min + eta_max))
            
            if (orbit_width > orbit_width_threshold) then
                ! Use thick orbit resonance calculation
                call driftorbit_coarse_thick(v, eta_min, eta_max, roots, nroots)
                return
            end if
        end if
#endif
        
        ! Use thin orbit calculation (original implementation)
        call driftorbit_coarse(v, eta_min, eta_max, roots, nroots)
        
    end subroutine driftorbit_coarse_unified
    
    subroutine driftorbit_coarse_thick(v, eta_min, eta_max, roots, nroots)
        ! Thick orbit resonance finder with finite orbit width effects
        real(8), intent(in) :: v, eta_min, eta_max
        real(8), intent(out) :: roots(:, :)
        integer, intent(out) :: nroots
        
        real(8) :: deta
        real(8) :: Omph, dOmphdv, dOmphdeta
        real(8) :: Omth, dOmthdv, dOmthdeta
        real(8) :: res, dresdv, dresdeta
        real(8) :: resold, dresdvold, dresdetaold
        real(8) :: eta, orbit_width, resonance_width
        integer :: k, ninterv
        
        ninterv = size(roots, 1)
        
        deta = (eta_max - eta_min)*1d0/ninterv
        nroots = 0
        
        do k = 0, ninterv
            eta = eta_min + k*deta
            
#ifdef USE_THICK_ORBITS
            ! Use unified frequency interface (includes thick orbit dispatch)
            call Om_th_unified(v, eta, Omth, dOmthdv, dOmthdeta)
            call Om_ph_unified(v, eta, Omph, dOmphdv, dOmphdeta)
            
            ! Account for finite orbit width in resonance detection
            orbit_width = calculate_orbit_width_criterion(v, eta)
            resonance_width = abs(Omth) * orbit_width  ! Broadened resonance
#else
            ! Fallback to thin orbit calculation
            call Om_th(v, eta, Omth, dOmthdv, dOmthdeta)
            call Om_ph(v, eta, Omph, dOmphdv, dOmphdeta)
            
            ! No orbit width effects for thin orbits
            resonance_width = 0.0d0
#endif
            
            res = mth*Omth + mph*Omph
            dresdv = mth*dOmthdv + mph*dOmphdv
            dresdeta = mth*dOmthdeta + mph*dOmphdeta
            
            if (k > 0) then
                ! Modified resonance condition: |res| < resonance_width or sign change
                if (abs(res) < resonance_width .or. sign(1d0, res) /= sign(1d0, resold)) then
                    nroots = nroots + 1
                    roots(nroots, 1) = eta - deta
                    roots(nroots, 2) = eta
                end if
            end if
            
            resold = res
            dresdvold = dresdv
            dresdetaold = dresdeta
        end do
        
    end subroutine driftorbit_coarse_thick

end module neort_resonance
