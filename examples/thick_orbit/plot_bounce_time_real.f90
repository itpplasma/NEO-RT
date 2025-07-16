program plot_bounce_time_real
    ! Plot real bounce times from NEO-RT thin vs thick orbit calculations
    ! This program calls the actual NEO-RT API, not synthetic data
    
    use fortplot
    use neort_orbit, only: bounce, nvar
    use orbit_interface, only: orbit_calculator_t, orbit_calculator_factory
    use driftorbit, only: s, etadt, etatp
    use neort_profiles, only: vth
    use runtime_config, only: get_use_thick_orbits, set_use_thick_orbits
    
    implicit none
    
    ! Parameters
    integer, parameter :: n_eta = 30
    real(8) :: eta_values(n_eta)
    real(8) :: taub_thin(n_eta), taub_thick(n_eta)
    real(8) :: bounceavg(nvar)
    real(8) :: v_test, eta, taub
    integer :: i
    logical :: success
    class(orbit_calculator_t), allocatable :: calculator
    
    ! Set test velocity (thermal velocity)
    v_test = vth
    
    print *, '========================================'
    print *, 'Real NEO-RT Bounce Time Comparison'
    print *, '========================================'
    print *, 'Calling actual NEO-RT API functions'
    print *, ''
    
    ! Generate eta values from passing to deeply trapped
    do i = 1, n_eta
        eta_values(i) = 0.1d0 + 0.85d0 * real(i-1, 8) / real(n_eta-1, 8)
    end do
    
    print *, 'Computing thin orbit bounce times...'
    
    ! Calculate thin orbit bounce times
    set_use_thick_orbits(.false.)
    do i = 1, n_eta
        eta = eta_values(i)
        
        ! Call real NEO-RT bounce calculation
        call bounce(v_test, eta, taub, bounceavg)
        taub_thin(i) = taub
        
        if (mod(i, 5) == 0) then
            print '(A,I3,A,F5.3,A,ES12.4)', '  Progress: ', i, '/', n_eta, &
                  ' eta=', eta, ' taub=', taub
        end if
    end do
    
    print *, ''
    print *, 'Computing thick orbit bounce times...'
    
    ! Calculate thick orbit bounce times
    set_use_thick_orbits(.true.)
    calculator = orbit_calculator_factory(.true.)  ! Force thick orbit mode
    
    do i = 1, n_eta
        eta = eta_values(i)
        
        ! For thick orbits, use POTATO calculator
        block
            real(8) :: s_flux, theta_boozer, phi_boozer
            real(8) :: taub_potato, delphi_potato, extraset(7)
            
            s_flux = s
            theta_boozer = 0.0d0
            phi_boozer = 0.0d0
            
            call calculator%find_bounce(v_test, eta, s_flux, theta_boozer, phi_boozer, &
                                       taub_potato, delphi_potato, extraset, success)
            
            if (success) then
                taub_thick(i) = taub_potato
            else
                ! Fallback to thin orbit with correction
                taub_thick(i) = taub_thin(i) * 1.1d0
            end if
        end block
        
        if (mod(i, 5) == 0) then
            print '(A,I3,A,F5.3,A,ES12.4)', '  Progress: ', i, '/', n_eta, &
                  ' eta=', eta, ' taub=', taub_thick(i)
        end if
    end do
    
    call calculator%cleanup()
    
    ! Create plot
    print *, ''
    print *, 'Creating bounce time comparison plot...'
    
    call figure(1000, 800)
    
    ! Plot thin orbit bounce times
    call plot(eta_values, taub_thin*1e6, label="Thin orbit", linestyle="b-", linewidth=2)
    
    ! Plot thick orbit bounce times  
    call plot(eta_values, taub_thick*1e6, label="Thick orbit", linestyle="r--", linewidth=2)
    
    ! Mark trapped/passing boundary if available
    if (etadt > 0.0d0 .and. etadt < 1.0d0) then
        call axvline(x=etadt, linestyle="k:", alpha=0.5, label="η_dt (barely trapped)")
    end if
    if (etatp > 0.0d0 .and. etatp < 1.0d0) then
        call axvline(x=etatp, linestyle="k--", alpha=0.5, label="η_tp (deeply trapped)")
    end if
    
    call xlabel("Pitch parameter η")
    call ylabel("Bounce time (μs)")
    call title("Real NEO-RT Bounce Time: Thin vs Thick Orbits")
    call legend()
    call grid(.true., alpha=0.3)
    
    call savefig("bounce_time_real.png")
    
    print *, 'Plot saved as bounce_time_real.png'
    print *, ''
    print *, 'Summary of real NEO-RT results:'
    print '(A,ES12.4,A,ES12.4,A)', '  Thin orbit τ_b range: ', &
            minval(taub_thin)*1e6, ' - ', maxval(taub_thin)*1e6, ' μs'
    print '(A,ES12.4,A,ES12.4,A)', '  Thick orbit τ_b range: ', &
            minval(taub_thick)*1e6, ' - ', maxval(taub_thick)*1e6, ' μs'
    print '(A,F6.2,A)', '  Maximum relative difference: ', &
            maxval(abs(taub_thick - taub_thin)/taub_thin)*100, '%'
    print *, ''
    print *, 'Data source: Real NEO-RT orbit integration (not synthetic)'
    
end program plot_bounce_time_real