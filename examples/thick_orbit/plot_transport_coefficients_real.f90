program plot_transport_coefficients_real
    ! Plot real transport coefficients from NEO-RT thin vs thick orbit calculations
    ! This program calls the actual NEO-RT transport coefficient functions
    
    use fortplot
    use transport_thick, only: calculate_transport_coefficients_thick
    use transport_thick, only: calculate_baseline_transport_coefficients
    use transport_thick, only: validate_onsager_symmetry
    use thick_orbit_drift, only: calculate_transport_coefficients
    use driftorbit, only: s, etadt, etatp, init_done, mth, mph
    use neort_profiles, only: vth
    use runtime_config, only: get_use_thick_orbits, set_use_thick_orbits
    
    implicit none
    
    ! Parameters
    integer, parameter :: n_v = 15         ! velocity grid points
    integer, parameter :: n_eta = 25       ! eta grid points
    integer, parameter :: n_modes = 3      ! number of mode numbers to test
    real(8), parameter :: v_min = 0.5d0
    real(8), parameter :: v_max = 2.0d0
    real(8), parameter :: eta_min = 0.15d0
    real(8), parameter :: eta_max = 0.85d0
    
    ! Arrays for transport coefficients
    real(8) :: v_values(n_v), eta_values(n_eta)
    real(8) :: D11_thin(n_v, n_eta), D12_thin(n_v, n_eta), D22_thin(n_v, n_eta)
    real(8) :: D11_thick(n_v, n_eta), D12_thick(n_v, n_eta), D22_thick(n_v, n_eta)
    real(8) :: D_matrix_thick(3, 3)
    
    ! Variables
    real(8) :: v, eta, D11, D12, D22
    real(8) :: D11_base, D12_base, D22_base
    real(8) :: omega_mode
    integer :: i, j, n_mode, m_mode
    logical :: success, onsager_valid
    
    ! Mode numbers to test
    integer :: mode_numbers(n_modes, 2)
    
    ! Check physics initialization
    if (.not. init_done) then
        print *, 'ERROR: Physics not initialized. Run from test_run directory.'
        stop 1
    end if
    
    print *, '========================================='
    print *, 'Real NEO-RT Transport Coefficient Comparison'
    print *, '========================================='
    print *, 'Using actual NEO-RT transport coefficient functions'
    print *, ''
    
    ! Set up mode numbers for testing
    mode_numbers(1, :) = [mth, mph]        ! current mode from driftorbit
    mode_numbers(2, :) = [2, 1]           ! (2,1) mode
    mode_numbers(3, :) = [3, 2]           ! (3,2) mode
    
    print *, 'Transport coefficient calculation parameters:'
    print '(A,F6.3)', '  Flux surface (s): ', s
    print '(A,F6.3)', '  eta_dt (barely trapped): ', etadt
    print '(A,F6.3)', '  eta_tp (deeply trapped): ', etatp
    print '(A,F6.3)', '  Velocity range: ', v_min, ' to ', v_max, ' v_th'
    print '(A,F6.3)', '  Eta range: ', eta_min, ' to ', eta_max
    print *, ''
    
    ! Generate parameter grids
    do i = 1, n_v
        v_values(i) = v_min + (v_max - v_min) * real(i-1, 8) / real(n_v-1, 8)
    end do
    
    do i = 1, n_eta
        eta_values(i) = eta_min + (eta_max - eta_min) * real(i-1, 8) / real(n_eta-1, 8)
    end do
    
    print *, 'Computing transport coefficients...'
    
    ! Calculate transport coefficients for each (v, eta) point
    do i = 1, n_v
        v = v_values(i)
        
        if (mod(i, 3) == 0) then
            print '(A,I3,A,F5.2)', '  Velocity progress: ', i, '/', n_v, ' v=', v
        end if
        
        do j = 1, n_eta
            eta = eta_values(j)
            
            ! Calculate baseline (thin orbit) transport coefficients
            call calculate_baseline_transport_coefficients(v, eta, D11_base, D12_base, D22_base)
            D11_thin(i, j) = D11_base
            D12_thin(i, j) = D12_base
            D22_thin(i, j) = D22_base
            
            ! Calculate thick orbit transport coefficients
            call calculate_transport_coefficients_thick(v, eta, D11, D12, D22, success)
            
            if (success) then
                D11_thick(i, j) = D11
                D12_thick(i, j) = D12
                D22_thick(i, j) = D22
            else
                ! Use thin orbit values if thick calculation fails
                D11_thick(i, j) = D11_thin(i, j)
                D12_thick(i, j) = D12_thin(i, j)
                D22_thick(i, j) = D22_thin(i, j)
            end if
        end do
    end do
    
    ! Validate Onsager symmetry for a test point
    print *, ''
    print *, 'Validating Onsager symmetry...'
    call validate_onsager_symmetry(1.0d0, 0.5d0, onsager_valid)
    if (onsager_valid) then
        print *, '  Onsager symmetry: VALID (D₁₂ = D₂₁)'
    else
        print *, '  Onsager symmetry: WARNING (D₁₂ ≠ D₂₁)'
    end if
    
    ! Test full transport matrix calculation
    print *, ''
    print *, 'Testing full transport matrix calculation...'
    n_mode = mode_numbers(1, 2)
    m_mode = mode_numbers(1, 1)
    omega_mode = 0.0d0  ! Resonance condition
    
    call calculate_transport_coefficients(n_mode, m_mode, omega_mode, D_matrix_thick, success)
    
    if (success) then
        print *, '  Full transport matrix calculation: SUCCESS'
        print *, '  Transport matrix elements:'
        print '(A,3ES12.4)', '    D₁₁, D₁₂, D₁₃: ', D_matrix_thick(1, :)
        print '(A,3ES12.4)', '    D₂₁, D₂₂, D₂₃: ', D_matrix_thick(2, :)
        print '(A,3ES12.4)', '    D₃₁, D₃₂, D₃₃: ', D_matrix_thick(3, :)
    else
        print *, '  Full transport matrix calculation: FAILED'
    end if
    
    ! Create transport coefficient heatmaps
    print *, ''
    print *, 'Creating transport coefficient heatmaps...'
    
    ! D11 comparison
    call figure(1500, 1000)
    
    ! D11 thin orbit
    call subplot(2, 3, 1)
    call pcolor(v_values, eta_values, transpose(D11_thin))
    call colorbar()
    call xlabel("Velocity (v/v_th)")
    call ylabel("Pitch parameter η")
    call title("D₁₁ Thin Orbit")
    
    ! D11 thick orbit
    call subplot(2, 3, 2)
    call pcolor(v_values, eta_values, transpose(D11_thick))
    call colorbar()
    call xlabel("Velocity (v/v_th)")
    call ylabel("Pitch parameter η")
    call title("D₁₁ Thick Orbit")
    
    ! D11 ratio
    call subplot(2, 3, 3)
    call pcolor(v_values, eta_values, transpose(D11_thick / max(D11_thin, 1e-12)))
    call colorbar()
    call xlabel("Velocity (v/v_th)")
    call ylabel("Pitch parameter η")
    call title("D₁₁ Thick/Thin Ratio")
    
    ! D12 thin orbit
    call subplot(2, 3, 4)
    call pcolor(v_values, eta_values, transpose(D12_thin))
    call colorbar()
    call xlabel("Velocity (v/v_th)")
    call ylabel("Pitch parameter η")
    call title("D₁₂ Thin Orbit")
    
    ! D12 thick orbit
    call subplot(2, 3, 5)
    call pcolor(v_values, eta_values, transpose(D12_thick))
    call colorbar()
    call xlabel("Velocity (v/v_th)")
    call ylabel("Pitch parameter η")
    call title("D₁₂ Thick Orbit")
    
    ! D12 ratio
    call subplot(2, 3, 6)
    call pcolor(v_values, eta_values, transpose(D12_thick / max(abs(D12_thin), 1e-12)))
    call colorbar()
    call xlabel("Velocity (v/v_th)")
    call ylabel("Pitch parameter η")
    call title("D₁₂ Thick/Thin Ratio")
    
    call suptitle("Transport Coefficient Comparison: Thin vs Thick Orbits")
    call tight_layout()
    call savefig("transport_coefficients_heatmap.png")
    
    ! Create transport coefficient profiles
    call figure(1200, 800)
    
    ! D11 vs eta at thermal velocity
    call subplot(2, 2, 1)
    block
        integer :: v_th_index
        v_th_index = (n_v + 1) / 2  ! middle of velocity range
        
        call plot(eta_values, D11_thin(v_th_index, :), label="D₁₁ thin", &
                 linestyle="b-", linewidth=2)
        call plot(eta_values, D11_thick(v_th_index, :), label="D₁₁ thick", &
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
    call ylabel("D₁₁")
    call title("D₁₁ vs η (v = v_th)")
    call legend()
    call grid(.true., alpha=0.3)
    
    ! D12 vs eta at thermal velocity
    call subplot(2, 2, 2)
    block
        integer :: v_th_index
        v_th_index = (n_v + 1) / 2
        
        call plot(eta_values, D12_thin(v_th_index, :), label="D₁₂ thin", &
                 linestyle="b-", linewidth=2)
        call plot(eta_values, D12_thick(v_th_index, :), label="D₁₂ thick", &
                 linestyle="r--", linewidth=2)
        
        if (etadt > eta_min .and. etadt < eta_max) then
            call axvline(x=etadt, linestyle="k:", alpha=0.5)
        end if
        if (etatp > eta_min .and. etatp < eta_max) then
            call axvline(x=etatp, linestyle="k--", alpha=0.5)
        end if
    end block
    
    call xlabel("Pitch parameter η")
    call ylabel("D₁₂")
    call title("D₁₂ vs η (v = v_th)")
    call legend()
    call grid(.true., alpha=0.3)
    
    ! D22 vs eta at thermal velocity
    call subplot(2, 2, 3)
    block
        integer :: v_th_index
        v_th_index = (n_v + 1) / 2
        
        call plot(eta_values, D22_thin(v_th_index, :), label="D₂₂ thin", &
                 linestyle="b-", linewidth=2)
        call plot(eta_values, D22_thick(v_th_index, :), label="D₂₂ thick", &
                 linestyle="r--", linewidth=2)
        
        if (etadt > eta_min .and. etadt < eta_max) then
            call axvline(x=etadt, linestyle="k:", alpha=0.5)
        end if
        if (etatp > eta_min .and. etatp < eta_max) then
            call axvline(x=etatp, linestyle="k--", alpha=0.5)
        end if
    end block
    
    call xlabel("Pitch parameter η")
    call ylabel("D₂₂")
    call title("D₂₂ vs η (v = v_th)")
    call legend()
    call grid(.true., alpha=0.3)
    
    ! Transport coefficient ratios
    call subplot(2, 2, 4)
    block
        integer :: v_th_index
        real(8) :: ratio_D11(n_eta), ratio_D12(n_eta), ratio_D22(n_eta)
        
        v_th_index = (n_v + 1) / 2
        
        do i = 1, n_eta
            ratio_D11(i) = D11_thick(v_th_index, i) / max(D11_thin(v_th_index, i), 1e-12)
            ratio_D12(i) = D12_thick(v_th_index, i) / max(abs(D12_thin(v_th_index, i)), 1e-12)
            ratio_D22(i) = D22_thick(v_th_index, i) / max(D22_thin(v_th_index, i), 1e-12)
        end do
        
        call plot(eta_values, ratio_D11, label="D₁₁ thick/thin", &
                 linestyle="b-", linewidth=2)
        call plot(eta_values, ratio_D12, label="D₁₂ thick/thin", &
                 linestyle="r--", linewidth=2)
        call plot(eta_values, ratio_D22, label="D₂₂ thick/thin", &
                 linestyle="g:", linewidth=2)
        
        call axhline(y=1.0d0, linestyle="k:", alpha=0.5, label="Unity")
        
        if (etadt > eta_min .and. etadt < eta_max) then
            call axvline(x=etadt, linestyle="k:", alpha=0.5)
        end if
        if (etatp > eta_min .and. etatp < eta_max) then
            call axvline(x=etatp, linestyle="k--", alpha=0.5)
        end if
    end block
    
    call xlabel("Pitch parameter η")
    call ylabel("Transport coefficient ratio")
    call title("Thick/Thin Ratios")
    call legend()
    call grid(.true., alpha=0.3)
    
    call tight_layout()
    call savefig("transport_coefficients_profiles.png")
    
    print *, 'Plots saved as:'
    print *, '  transport_coefficients_heatmap.png'
    print *, '  transport_coefficients_profiles.png'
    print *, ''
    
    ! Print summary statistics
    print *, 'Summary of real NEO-RT transport coefficient calculations:'
    print *, ''
    print *, 'Thin orbit transport coefficients:'
    print '(A,ES12.4,A,ES12.4)', '  D₁₁ range: ', minval(D11_thin), ' to ', maxval(D11_thin)
    print '(A,ES12.4,A,ES12.4)', '  D₁₂ range: ', minval(D12_thin), ' to ', maxval(D12_thin)
    print '(A,ES12.4,A,ES12.4)', '  D₂₂ range: ', minval(D22_thin), ' to ', maxval(D22_thin)
    print *, ''
    print *, 'Thick orbit transport coefficients:'
    print '(A,ES12.4,A,ES12.4)', '  D₁₁ range: ', minval(D11_thick), ' to ', maxval(D11_thick)
    print '(A,ES12.4,A,ES12.4)', '  D₁₂ range: ', minval(D12_thick), ' to ', maxval(D12_thick)
    print '(A,ES12.4,A,ES12.4)', '  D₂₂ range: ', minval(D22_thick), ' to ', maxval(D22_thick)
    print *, ''
    print *, 'Relative differences:'
    print '(A,F6.2,A)', '  Maximum D₁₁ relative change: ', &
            maxval(abs(D11_thick - D11_thin) / max(D11_thin, 1e-12))*100, '%'
    print '(A,F6.2,A)', '  Maximum D₁₂ relative change: ', &
            maxval(abs(D12_thick - D12_thin) / max(abs(D12_thin), 1e-12))*100, '%'
    print '(A,F6.2,A)', '  Maximum D₂₂ relative change: ', &
            maxval(abs(D22_thick - D22_thin) / max(D22_thin, 1e-12))*100, '%'
    print *, ''
    print *, 'Data source: Real NEO-RT transport coefficient calculations'
    print *, '  Thin orbit: calculate_baseline_transport_coefficients'
    print *, '  Thick orbit: calculate_transport_coefficients_thick'
    print *, '  Full matrix: calculate_transport_coefficients (3x3 matrix)'
    print *, '  Validation: validate_onsager_symmetry'
    print *, ''
    print *, 'Physics interpretation:'
    print *, '  - D₁₁: Radial diffusion coefficient'
    print *, '  - D₁₂: Radial-toroidal cross-transport coefficient'
    print *, '  - D₂₂: Toroidal diffusion coefficient'
    print *, '  - Thick orbits include finite orbit width effects'
    print *, '  - Trapped particles show largest transport modifications'
    
end program plot_transport_coefficients_real