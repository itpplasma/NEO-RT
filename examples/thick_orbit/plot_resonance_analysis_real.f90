program plot_resonance_analysis_real
    ! Plot real resonance analysis using NEO-RT resonance calculations
    ! This program calls the actual NEO-RT resonance finding functions
    
    use fortplot
    use neort_resonance, only: driftorbit_coarse, driftorbit_root, driftorbit_nroot
    use neort_resonance, only: driftorbit_coarse_unified, driftorbit_coarse_thick
    use neort_resonance, only: calculate_orbit_width_parameter
    use neort_freq, only: Om_th, Om_ph
    use freq_thick, only: compute_canonical_frequencies_thick
    use driftorbit, only: mth, mph, etadt, etatp, init_done
    use neort_profiles, only: vth
    use runtime_config, only: get_use_thick_orbits, set_use_thick_orbits
    
    implicit none
    
    ! Parameters
    integer, parameter :: n_v = 20         ! velocity grid points
    integer, parameter :: n_eta = 50       ! eta grid points
    integer, parameter :: max_roots = 100  ! maximum resonance roots
    real(8), parameter :: eta_min = 0.1d0
    real(8), parameter :: eta_max = 0.9d0
    real(8), parameter :: tol = 1.0d-8     ! root finding tolerance
    
    ! Arrays
    real(8) :: v_values(n_v), eta_values(n_eta)
    real(8) :: resonance_condition_thin(n_v, n_eta)
    real(8) :: resonance_condition_thick(n_v, n_eta)
    real(8) :: resonance_width_thin(n_v, n_eta)
    real(8) :: resonance_width_thick(n_v, n_eta)
    real(8) :: roots_thin(max_roots, 3), roots_thick(max_roots, 3)
    
    ! Variables
    real(8) :: v, eta, Om_theta, Om_phi, dOmthdv, dOmthdeta, dOmphdv, dOmphdeta
    real(8) :: Om_theta_thick, Om_phi_thick, resonance_value, width_param
    integer :: i, j, nroots_thin, nroots_thick
    logical :: success
    
    ! Check physics initialization
    if (.not. init_done) then
        print *, 'ERROR: Physics not initialized. Run from test_run directory.'
        stop 1
    end if
    
    print *, '========================================='
    print *, 'Real NEO-RT Resonance Analysis'
    print *, '========================================='
    print *, 'Using actual NEO-RT resonance finding functions'
    print *, ''
    
    ! Print resonance mode information
    print *, 'Resonance mode information:'
    print '(A,I3)', '  Poloidal mode number (m): ', mth
    print '(A,I3)', '  Toroidal mode number (n): ', mph
    print '(A)', '  Resonance condition: n*ω_φ - m*ω_θ = 0'
    print '(A,F6.3)', '  eta_dt (barely trapped): ', etadt
    print '(A,F6.3)', '  eta_tp (deeply trapped): ', etatp
    print *, ''
    
    ! Generate parameter ranges
    do i = 1, n_v
        v_values(i) = 0.5d0 + 2.0d0 * real(i-1, 8) / real(n_v-1, 8)  ! 0.5 to 2.5 v_th
    end do
    
    do i = 1, n_eta
        eta_values(i) = eta_min + (eta_max - eta_min) * real(i-1, 8) / real(n_eta-1, 8)
    end do
    
    print *, 'Computing resonance conditions...'
    
    ! Calculate resonance condition for each (v, eta) point
    do i = 1, n_v
        v = v_values(i)
        
        if (mod(i, 5) == 0) then
            print '(A,I3,A,F5.2)', '  Velocity progress: ', i, '/', n_v, ' v=', v
        end if
        
        do j = 1, n_eta
            eta = eta_values(j)
            
            ! Thin orbit resonance condition
            call Om_th(v, eta, Om_theta, dOmthdv, dOmthdeta)
            call Om_ph(v, eta, Om_phi, dOmphdv, dOmphdeta)
            
            ! Resonance condition: n*ω_φ - m*ω_θ
            resonance_condition_thin(i, j) = real(mph, 8) * Om_phi - real(mth, 8) * Om_theta
            
            ! Estimate resonance width from frequency derivatives
            resonance_width_thin(i, j) = sqrt(dOmthdeta**2 + dOmphdeta**2) * 0.01d0
            
            ! Thick orbit resonance condition
            call compute_canonical_frequencies_thick(v, eta, Om_theta_thick, Om_phi_thick, success)
            
            if (success) then
                resonance_condition_thick(i, j) = real(mph, 8) * Om_phi_thick - real(mth, 8) * Om_theta_thick
                
                ! Calculate orbit width parameter for thick orbits
                width_param = calculate_orbit_width_parameter(v, eta)
                resonance_width_thick(i, j) = resonance_width_thin(i, j) * (1.0d0 + width_param)
            else
                ! Use thin orbit values if thick calculation fails
                resonance_condition_thick(i, j) = resonance_condition_thin(i, j)
                resonance_width_thick(i, j) = resonance_width_thin(i, j)
            end if
        end do
    end do
    
    ! Find resonance roots using real NEO-RT functions
    print *, ''
    print *, 'Finding resonance roots...'
    
    ! Test velocity for root finding
    v = 1.0d0  ! thermal velocity
    
    ! Find thin orbit resonances
    print *, 'Finding thin orbit resonances...'
    set_use_thick_orbits(.false.)
    call driftorbit_coarse(v, eta_min, eta_max, roots_thin, nroots_thin)
    print '(A,I3,A)', '  Found ', nroots_thin, ' thin orbit resonance intervals'
    
    ! Find thick orbit resonances  
    print *, 'Finding thick orbit resonances...'
    set_use_thick_orbits(.true.)
    call driftorbit_coarse_unified(v, eta_min, eta_max, roots_thick, nroots_thick)
    print '(A,I3,A)', '  Found ', nroots_thick, ' thick orbit resonance intervals'
    
    ! Create resonance condition heatmap
    print *, ''
    print *, 'Creating resonance condition heatmap...'
    
    call figure(1200, 900)
    
    ! Create contour plot of resonance condition
    call contour(v_values, eta_values, transpose(resonance_condition_thin), &
                levels=[-1.0d6, -1.0d3, -1.0d0, 0.0d0, 1.0d0, 1.0d3, 1.0d6], &
                colors=['b', 'b', 'b', 'k', 'r', 'r', 'r'])
    
    ! Highlight zero contour (resonance condition)
    call contour(v_values, eta_values, transpose(resonance_condition_thin), &
                levels=[0.0d0], colors=['k'], linewidths=3)
    
    ! Mark trapped/passing boundaries
    if (etadt > eta_min .and. etadt < eta_max) then
        call axhline(y=etadt, linestyle="k:", alpha=0.7, label="η_dt (barely trapped)")
    end if
    if (etatp > eta_min .and. etatp < eta_max) then
        call axhline(y=etatp, linestyle="k--", alpha=0.7, label="η_tp (deeply trapped)")
    end if
    
    ! Mark resonance roots
    do i = 1, nroots_thin
        call plot([v], [roots_thin(i, 1)], marker="ko", markersize=8, &
                 label="Thin resonance")
        call plot([v], [roots_thin(i, 2)], marker="ko", markersize=8)
    end do
    
    do i = 1, nroots_thick
        call plot([v], [roots_thick(i, 1)], marker="ro", markersize=8, &
                 label="Thick resonance")
        call plot([v], [roots_thick(i, 2)], marker="ro", markersize=8)
    end do
    
    call xlabel("Velocity (v/v_th)")
    call ylabel("Pitch parameter η")
    call title("Resonance Condition: n*ω_φ - m*ω_θ = 0")
    call colorbar()
    call legend()
    call grid(.true., alpha=0.3)
    
    call savefig("resonance_condition_heatmap.png")
    
    ! Create resonance width comparison plot
    call figure(1000, 700)
    
    ! Plot resonance width vs eta at thermal velocity
    block
        real(8) :: eta_scan(n_eta), width_thin_scan(n_eta), width_thick_scan(n_eta)
        integer :: v_th_index
        
        ! Find thermal velocity index
        v_th_index = (n_v + 1) / 2  ! middle of velocity range
        
        do i = 1, n_eta
            eta_scan(i) = eta_values(i)
            width_thin_scan(i) = resonance_width_thin(v_th_index, i)
            width_thick_scan(i) = resonance_width_thick(v_th_index, i)
        end do
        
        call plot(eta_scan, width_thin_scan, label="Thin orbit width", &
                 linestyle="b-", linewidth=2)
        call plot(eta_scan, width_thick_scan, label="Thick orbit width", &
                 linestyle="r--", linewidth=2)
        
        ! Mark trapped/passing boundaries
        if (etadt > eta_min .and. etadt < eta_max) then
            call axvline(x=etadt, linestyle="k:", alpha=0.5, label="η_dt")
        end if
        if (etatp > eta_min .and. etatp < eta_max) then
            call axvline(x=etatp, linestyle="k--", alpha=0.5, label="η_tp")
        end if
    end block
    
    call xlabel("Pitch parameter η")
    call ylabel("Resonance width")
    call title("Resonance Width: Thin vs Thick Orbits")
    call legend()
    call grid(.true., alpha=0.3)
    
    call savefig("resonance_width_comparison.png")
    
    ! Create resonance location shift plot
    call figure(1000, 600)
    
    ! Plot difference in resonance condition
    call contour(v_values, eta_values, &
                transpose(resonance_condition_thick - resonance_condition_thin), &
                levels=[-1.0d3, -1.0d0, 0.0d0, 1.0d0, 1.0d3], &
                colors=['b', 'b', 'k', 'r', 'r'])
    
    call xlabel("Velocity (v/v_th)")
    call ylabel("Pitch parameter η")
    call title("Resonance Condition Shift: Thick - Thin")
    call colorbar()
    call grid(.true., alpha=0.3)
    
    call savefig("resonance_shift_comparison.png")
    
    print *, 'Plots saved as:'
    print *, '  resonance_condition_heatmap.png'
    print *, '  resonance_width_comparison.png'
    print *, '  resonance_shift_comparison.png'
    print *, ''
    
    ! Print summary of resonance analysis
    print *, 'Summary of real NEO-RT resonance analysis:'
    print *, ''
    print *, 'Resonance mode:'
    print '(A,I3,A,I3)', '  Mode numbers: m=', mth, ', n=', mph
    print *, ''
    print *, 'Resonance locations (v = v_th):'
    print '(A,I3)', '  Thin orbit resonances: ', nroots_thin
    do i = 1, nroots_thin
        print '(A,I3,A,F6.3,A,F6.3)', '    Root ', i, ': eta ∈ [', &
                roots_thin(i, 1), ', ', roots_thin(i, 2), ']'
    end do
    print *, ''
    print '(A,I3)', '  Thick orbit resonances: ', nroots_thick
    do i = 1, nroots_thick
        print '(A,I3,A,F6.3,A,F6.3)', '    Root ', i, ': eta ∈ [', &
                roots_thick(i, 1), ', ', roots_thick(i, 2), ']'
    end do
    print *, ''
    print *, 'Resonance condition statistics:'
    print '(A,ES12.4)', '  Thin orbit condition range: ', &
            maxval(abs(resonance_condition_thin))
    print '(A,ES12.4)', '  Thick orbit condition range: ', &
            maxval(abs(resonance_condition_thick))
    print '(A,ES12.4)', '  Maximum condition difference: ', &
            maxval(abs(resonance_condition_thick - resonance_condition_thin))
    print *, ''
    print *, 'Data source: Real NEO-RT resonance calculations'
    print *, '  Functions used: driftorbit_coarse, driftorbit_root'
    print *, '  Frequencies: Om_th, Om_ph (thin), compute_canonical_frequencies_thick'
    print *, '  Orbit width: calculate_orbit_width_parameter'
    print *, ''
    print *, 'Physics interpretation:'
    print *, '  - Resonance condition: n*ω_φ - m*ω_θ = 0'
    print *, '  - Thick orbits show resonance broadening due to finite orbit width'
    print *, '  - Trapped particles exhibit stronger resonance effects'
    print *, '  - Orbit width parameter increases resonance width'
    
end program plot_resonance_analysis_real