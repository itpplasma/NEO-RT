module neort_resonance
    use neort_freq, only: Om_th, Om_ph, d_Om_ds
    use freq_thick, only: compute_canonical_frequencies_thick
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
        ! Unified resonance finder with runtime thick/thin orbit dispatch
        use runtime_config, only: get_use_thick_orbits
        real(8), intent(in) :: v, eta_min, eta_max
        real(8), intent(out) :: roots(:, :)
        integer, intent(out) :: nroots
        
        logical, parameter :: enable_thick_orbit_resonance = .true.
        real(8), parameter :: orbit_width_threshold = 0.01d0
        
        real(8) :: orbit_width
        
        if (get_use_thick_orbits() .and. enable_thick_orbit_resonance) then
            orbit_width = calculate_orbit_width_parameter(v, 0.5d0 * (eta_min + eta_max))
            
            if (orbit_width > orbit_width_threshold) then
                ! Use thick orbit resonance calculation
                call driftorbit_coarse_thick(v, eta_min, eta_max, roots, nroots)
                return
            end if
        end if
        
        ! Use thin orbit calculation (original implementation)
        call driftorbit_coarse(v, eta_min, eta_max, roots, nroots)
        
    end subroutine driftorbit_coarse_unified
    
    subroutine driftorbit_coarse_thick(v, eta_min, eta_max, roots, nroots)
        ! Thick orbit resonance finder with finite orbit width effects
        real(8), intent(in) :: v, eta_min, eta_max
        real(8), intent(out) :: roots(:, :)
        integer, intent(out) :: nroots
        
        real(8) :: deta
        real(8) :: Omph_thick, Omth_thick
        real(8) :: res, resold
        real(8) :: eta, orbit_width, resonance_width
        integer :: k, ninterv
        logical :: freq_success
        
        ninterv = size(roots, 1)
        
        deta = (eta_max - eta_min)*1d0/ninterv
        nroots = 0
        resold = 0d0
        
        do k = 0, ninterv
            eta = eta_min + k*deta
            
            ! Use thick orbit frequency calculation
            call compute_canonical_frequencies_thick(v, eta, Omth_thick, Omph_thick, freq_success)
            
            if (.not. freq_success) then
                ! Fallback to thin orbit calculation if thick orbit fails
                block
                    real(8) :: dOmthdv_dummy, dOmthdeta_dummy, dOmphdv_dummy, dOmphdeta_dummy
                    call Om_th(v, eta, Omth_thick, dOmthdv_dummy, dOmthdeta_dummy)
                    call Om_ph(v, eta, Omph_thick, dOmphdv_dummy, dOmphdeta_dummy)
                end block
                resonance_width = 0.0d0
            else
                ! Account for finite orbit width in resonance detection
                orbit_width = calculate_orbit_width_parameter(v, eta)
                resonance_width = abs(Omth_thick) * orbit_width  ! Broadened resonance
            end if
            
            res = mth*Omth_thick + mph*Omph_thick
            
            if (k > 0) then
                ! Modified resonance condition: |res| < resonance_width or sign change
                if (abs(res) < resonance_width .or. sign(1d0, res) /= sign(1d0, resold)) then
                    nroots = nroots + 1
                    roots(nroots, 1) = eta - deta
                    roots(nroots, 2) = eta
                end if
            end if
            
            resold = res
        end do
        
    end subroutine driftorbit_coarse_thick

    function calculate_orbit_width_parameter(v, eta) result(orbit_width)
        ! Calculate orbit width parameter for resonance broadening
        real(8), intent(in) :: v, eta
        real(8) :: orbit_width
        
        ! Physical parameters
        real(8), parameter :: B_field = 2.5d0       ! Tesla
        real(8), parameter :: mass_amu = 2.0d0      ! Deuterium
        real(8), parameter :: L_B = 0.5d0           ! Magnetic scale length
        real(8), parameter :: v_thermal_ref = 1.0d6 ! Reference thermal velocity
        real(8) :: rho_gyro
        
        ! Calculate gyroradius
        rho_gyro = (v / v_thermal_ref) * 1.66d-27 * mass_amu * v_thermal_ref / (1.6d-19 * B_field)
        
        ! Orbit width parameter (normalized to magnetic scale length)
        orbit_width = rho_gyro / L_B
        
    end function calculate_orbit_width_parameter

end module neort_resonance
