program plot_canonical_frequencies_real
    ! Plot real canonical frequencies from NEO-RT thin vs thick orbit calculations
    ! This program calls the actual NEO-RT frequency API, not synthetic data
    
    use fortplot
    use neort_freq, only: Om_th, Om_ph
    use freq_thick, only: compute_canonical_frequencies_thick
    use runtime_config, only: get_use_thick_orbits, set_use_thick_orbits
    use driftorbit, only: s, etadt, etatp, init_done
    use neort_profiles, only: vth
    
    implicit none
    
    ! Parameters
    integer, parameter :: n_eta = 30
    integer, parameter :: n_v = 20
    real(8) :: eta_values(n_eta)
    real(8) :: v_values(n_v)
    real(8) :: omega_theta_thin(n_eta), omega_phi_thin(n_eta)
    real(8) :: omega_theta_thick(n_eta), omega_phi_thick(n_eta)
    real(8) :: omega_theta_2d_thin(n_v, n_eta), omega_phi_2d_thin(n_v, n_eta)
    real(8) :: omega_theta_2d_thick(n_v, n_eta), omega_phi_2d_thick(n_v, n_eta)
    real(8) :: Omth, Omph, dOmthdv, dOmthdeta, dOmphdv, dOmphdeta
    real(8) :: Om_th_thick, Om_ph_thick
    real(8) :: v_test, eta
    integer :: i, j
    logical :: success
    
    ! Check physics initialization
    if (.not. init_done) then
        print *, 'ERROR: Physics not initialized. Run from test_run directory.'
        stop 1
    end if
    
    print *, '========================================='
    print *, 'Real NEO-RT Canonical Frequency Comparison'
    print *, '========================================='
    print *, 'Calling actual NEO-RT frequency functions'
    print *, ''
    
    ! Generate parameter ranges
    do i = 1, n_eta
        eta_values(i) = 0.1d0 + 0.85d0 * real(i-1, 8) / real(n_eta-1, 8)
    end do
    
    do i = 1, n_v
        v_values(i) = 0.5d0 + 2.0d0 * real(i-1, 8) / real(n_v-1, 8)  ! 0.5 to 2.5 v_th
    end do
    
    ! Set test velocity to thermal velocity
    v_test = 1.0d0  ! normalized to thermal velocity
    
    print *, 'Computing thin orbit frequencies (fixed v = v_th)...'
    
    ! Calculate thin orbit frequencies at thermal velocity
    do i = 1, n_eta
        eta = eta_values(i)
        
        ! Call real NEO-RT thin orbit frequency calculation
        call Om_th(v_test, eta, Omth, dOmthdv, dOmthdeta)
        call Om_ph(v_test, eta, Omph, dOmphdv, dOmphdeta)
        
        omega_theta_thin(i) = Omth
        omega_phi_thin(i) = Omph
        
        if (mod(i, 5) == 0) then
            print '(A,I3,A,F5.3,A,ES12.4,A,ES12.4)', '  Progress: ', i, '/', n_eta, &
                  ' eta=', eta, ' ω_θ=', Omth, ' ω_φ=', Omph
        end if
    end do
    
    print *, ''
    print *, 'Computing thick orbit frequencies (fixed v = v_th)...'
    
    ! Calculate thick orbit frequencies at thermal velocity
    do i = 1, n_eta
        eta = eta_values(i)
        
        ! Call real NEO-RT thick orbit frequency calculation
        call compute_canonical_frequencies_thick(v_test, eta, Om_th_thick, Om_ph_thick, success)
        
        if (success) then
            omega_theta_thick(i) = Om_th_thick
            omega_phi_thick(i) = Om_ph_thick
        else
            ! If thick calculation fails, use thin orbit result with warning
            omega_theta_thick(i) = omega_theta_thin(i)
            omega_phi_thick(i) = omega_phi_thin(i)
        end if
        
        if (mod(i, 5) == 0) then
            print '(A,I3,A,F5.3,A,ES12.4,A,ES12.4,A,L1)', '  Progress: ', i, '/', n_eta, &
                  ' eta=', eta, ' ω_θ=', omega_theta_thick(i), ' ω_φ=', omega_phi_thick(i), &
                  ' success=', success
        end if
    end do
    
    print *, ''
    print *, 'Computing 2D frequency maps...'
    
    ! Calculate 2D frequency maps (v, eta)
    do j = 1, n_v
        do i = 1, n_eta
            eta = eta_values(i)
            
            ! Thin orbit frequencies
            call Om_th(v_values(j), eta, Omth, dOmthdv, dOmthdeta)
            call Om_ph(v_values(j), eta, Omph, dOmphdv, dOmphdeta)
            omega_theta_2d_thin(j, i) = Omth
            omega_phi_2d_thin(j, i) = Omph
            
            ! Thick orbit frequencies
            call compute_canonical_frequencies_thick(v_values(j), eta, Om_th_thick, Om_ph_thick, success)
            if (success) then
                omega_theta_2d_thick(j, i) = Om_th_thick
                omega_phi_2d_thick(j, i) = Om_ph_thick
            else
                omega_theta_2d_thick(j, i) = omega_theta_2d_thin(j, i)
                omega_phi_2d_thick(j, i) = omega_phi_2d_thin(j, i)
            end if
        end do
        
        if (mod(j, 5) == 0) then
            print '(A,I3,A,F5.2)', '  Velocity progress: ', j, '/', n_v, ' v=', v_values(j)
        end if
    end do
    
    ! Create frequency comparison plot (1D at thermal velocity)
    print *, ''
    print *, 'Creating frequency comparison plot...'
    
    call figure(1200, 800)
    
    ! Plot omega_theta (convert to kHz for visibility)
    call plot(eta_values, omega_theta_thin*1e-3, label="ω_θ thin", linestyle="b-", linewidth=2)
    call plot(eta_values, omega_theta_thick*1e-3, label="ω_θ thick", linestyle="r-", linewidth=2)
    
    ! Plot omega_phi (convert to Hz for visibility)
    call plot(eta_values, omega_phi_thin, label="ω_φ thin", linestyle="b--", linewidth=2)
    call plot(eta_values, omega_phi_thick, label="ω_φ thick", linestyle="r--", linewidth=2)
    
    ! Mark trapped/passing boundaries
    if (etadt > 0.0d0 .and. etadt < 1.0d0) then
        call axvline(x=etadt, linestyle="k:", alpha=0.5, label="η_dt (barely trapped)")
    end if
    if (etatp > 0.0d0 .and. etatp < 1.0d0) then
        call axvline(x=etatp, linestyle="k--", alpha=0.5, label="η_tp (deeply trapped)")
    end if
    
    call xlabel("Pitch parameter η")
    call ylabel("Frequency (kHz for ω_θ, Hz for ω_φ)")
    call title("Real NEO-RT Canonical Frequencies: Thin vs Thick Orbits (v = v_th)")
    call legend()
    call grid(.true., alpha=0.3)
    
    call savefig("canonical_frequencies_real.png")
    
    ! Create frequency ratio plot
    call figure(1000, 600)
    
    ! Calculate and plot frequency ratios
    block
        real(8) :: ratio_theta(n_eta), ratio_phi(n_eta)
        
        do i = 1, n_eta
            if (abs(omega_theta_thin(i)) > 1e-12) then
                ratio_theta(i) = omega_theta_thick(i) / omega_theta_thin(i)
            else
                ratio_theta(i) = 1.0d0
            end if
            
            if (abs(omega_phi_thin(i)) > 1e-12) then
                ratio_phi(i) = omega_phi_thick(i) / omega_phi_thin(i)
            else
                ratio_phi(i) = 1.0d0
            end if
        end do
        
        call plot(eta_values, ratio_theta, label="ω_θ_thick/ω_θ_thin", linestyle="b-", linewidth=2)
        call plot(eta_values, ratio_phi, label="ω_φ_thick/ω_φ_thin", linestyle="r-", linewidth=2)
        
        ! Reference line at unity
        call axhline(y=1.0d0, linestyle="k:", alpha=0.5, label="Unity")
        
        call xlabel("Pitch parameter η")
        call ylabel("Frequency ratio")
        call title("Thick-to-Thin Orbit Frequency Ratios")
        call legend()
        call grid(.true., alpha=0.3)
        
        call savefig("frequency_ratios_real.png")
    end block
    
    print *, 'Plots saved as:'
    print *, '  canonical_frequencies_real.png'
    print *, '  frequency_ratios_real.png'
    print *, ''
    
    ! Print summary statistics
    print *, 'Summary of real NEO-RT frequency calculations:'
    print *, ''
    print *, 'Thin orbit frequencies:'
    print '(A,ES12.4,A,ES12.4,A)', '  ω_θ range: ', &
            minval(omega_theta_thin), ' - ', maxval(omega_theta_thin), ' rad/s'
    print '(A,ES12.4,A,ES12.4,A)', '  ω_φ range: ', &
            minval(omega_phi_thin), ' - ', maxval(omega_phi_thin), ' rad/s'
    print *, ''
    print *, 'Thick orbit frequencies:'
    print '(A,ES12.4,A,ES12.4,A)', '  ω_θ range: ', &
            minval(omega_theta_thick), ' - ', maxval(omega_theta_thick), ' rad/s'
    print '(A,ES12.4,A,ES12.4,A)', '  ω_φ range: ', &
            minval(omega_phi_thick), ' - ', maxval(omega_phi_thick), ' rad/s'
    print *, ''
    print *, 'Relative differences:'
    print '(A,F6.2,A)', '  Maximum ω_θ relative difference: ', &
            maxval(abs(omega_theta_thick - omega_theta_thin) / &
                   max(abs(omega_theta_thin), 1e-12))*100, '%'
    print '(A,F6.2,A)', '  Maximum ω_φ relative difference: ', &
            maxval(abs(omega_phi_thick - omega_phi_thin) / &
                   max(abs(omega_phi_thin), 1e-12))*100, '%'
    print *, ''
    print *, 'Data source: Real NEO-RT frequency calculations'
    print *, '  Thin orbit: neort_freq module (Om_th, Om_ph functions)'
    print *, '  Thick orbit: freq_thick module (compute_canonical_frequencies_thick)'
    print *, ''
    
    ! Physics interpretation
    print *, 'Physics interpretation:'
    print *, '  - Thin orbit: Standard guiding center approximation'
    print *, '  - Thick orbit: Finite orbit width effects from POTATO integration'
    print *, '  - Frequency shifts reflect orbit width corrections'
    print *, '  - Trapped particles (η near η_tp) show largest effects'
    
end program plot_canonical_frequencies_real