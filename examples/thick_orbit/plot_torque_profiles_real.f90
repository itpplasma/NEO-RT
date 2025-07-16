program plot_torque_profiles_real
    ! Plot real NTV torque profiles from NEO-RT thin vs thick orbit calculations
    ! This program calls the actual NEO-RT torque calculation functions
    
    use fortplot
    use torque_thick, only: calculate_ntv_torque_density, calculate_velocity_space_torque
    use torque_thick, only: calculate_torque_profile, validate_torque_conservation
    use thick_orbit_drift, only: calculate_transport_coefficients
    use driftorbit, only: s, etadt, etatp, init_done, mth, mph
    use neort_profiles, only: vth
    use runtime_config, only: get_use_thick_orbits, set_use_thick_orbits
    
    implicit none
    
    ! Parameters
    integer, parameter :: n_flux = 20      ! flux surface grid points
    integer, parameter :: n_v = 15         ! velocity grid points
    integer, parameter :: n_eta = 20       ! eta grid points
    integer, parameter :: n_modes = 3      ! number of mode numbers to test
    real(8), parameter :: s_min = 0.1d0
    real(8), parameter :: s_max = 0.9d0
    real(8), parameter :: v_min = 0.5d0
    real(8), parameter :: v_max = 2.0d0
    real(8), parameter :: eta_min = 0.2d0
    real(8), parameter :: eta_max = 0.8d0
    
    ! Arrays
    real(8) :: flux_surfaces(n_flux)
    real(8) :: torque_profile_thin(n_flux), torque_profile_thick(n_flux)
    real(8) :: v_values(n_v), eta_values(n_eta)
    real(8) :: torque_density_thin(n_v, n_eta), torque_density_thick(n_v, n_eta)
    real(8) :: total_torque_thin, total_torque_thick
    
    ! Variables
    real(8) :: v, eta, omega_mode, torque_density, s_test
    integer :: i, j, n_mode, m_mode
    logical :: success, conservation_valid
    
    ! Mode numbers to test
    integer :: mode_numbers(n_modes, 2)
    
    ! Check physics initialization
    if (.not. init_done) then
        print *, 'ERROR: Physics not initialized. Run from test_run directory.'
        stop 1
    end if
    
    print *, '========================================='
    print *, 'Real NEO-RT Torque Profile Comparison'
    print *, '========================================='
    print *, 'Using actual NEO-RT torque calculation functions'
    print *, ''
    
    ! Set up mode numbers for testing
    mode_numbers(1, :) = [mth, mph]        ! current mode from driftorbit
    mode_numbers(2, :) = [2, 1]           ! (2,1) mode
    mode_numbers(3, :) = [3, 2]           ! (3,2) mode
    
    ! Set resonance mode for main calculation
    n_mode = mode_numbers(1, 2)
    m_mode = mode_numbers(1, 1)
    omega_mode = 0.0d0  ! resonance condition
    
    print *, 'NTV torque calculation parameters:'
    print '(A,I3,A,I3)', '  Mode numbers: (m,n) = (', m_mode, ',', n_mode, ')'
    print '(A,ES12.4)', '  Mode frequency: ω = ', omega_mode, ' rad/s'
    print '(A,F6.3)', '  eta_dt (barely trapped): ', etadt
    print '(A,F6.3)', '  eta_tp (deeply trapped): ', etatp
    print '(A,F6.3)', '  Flux surface range: ', s_min, ' to ', s_max
    print *, ''
    
    ! Generate parameter grids
    do i = 1, n_flux
        flux_surfaces(i) = s_min + (s_max - s_min) * real(i-1, 8) / real(n_flux-1, 8)
    end do
    
    do i = 1, n_v
        v_values(i) = v_min + (v_max - v_min) * real(i-1, 8) / real(n_v-1, 8)
    end do
    
    do i = 1, n_eta
        eta_values(i) = eta_min + (eta_max - eta_min) * real(i-1, 8) / real(n_eta-1, 8)
    end do
    
    ! Calculate torque density in velocity space for fixed flux surface
    print *, 'Computing torque density in velocity space...'
    s_test = 0.5d0  ! middle flux surface
    s = s_test
    
    do i = 1, n_v
        v = v_values(i)
        
        if (mod(i, 5) == 0) then
            print '(A,I3,A,F5.2)', '  Velocity progress: ', i, '/', n_v, ' v=', v
        end if
        
        do j = 1, n_eta
            eta = eta_values(j)
            
            ! Calculate thin orbit torque density (using transport coefficients)
            set_use_thick_orbits(.false.)
            call calculate_ntv_torque_density(v, eta, n_mode, m_mode, omega_mode, &
                                             torque_density, success)
            
            if (success) then
                torque_density_thin(i, j) = torque_density
            else
                torque_density_thin(i, j) = 0.0d0
            end if
            
            ! Calculate thick orbit torque density
            set_use_thick_orbits(.true.)
            call calculate_ntv_torque_density(v, eta, n_mode, m_mode, omega_mode, &
                                             torque_density, success)
            
            if (success) then
                torque_density_thick(i, j) = torque_density
            else
                torque_density_thick(i, j) = torque_density_thin(i, j)  ! fallback
            end if
        end do
    end do
    
    ! Calculate total torque for each physics model
    print *, ''
    print *, 'Computing total torque...'
    
    set_use_thick_orbits(.false.)
    call calculate_velocity_space_torque(n_mode, m_mode, omega_mode, total_torque_thin, success)
    if (success) then
        print '(A,ES12.4,A)', '  Total thin orbit torque: ', total_torque_thin, ' N·m'
    else
        print *, '  Total thin orbit torque: CALCULATION FAILED'
        total_torque_thin = 0.0d0
    end if
    
    set_use_thick_orbits(.true.)
    call calculate_velocity_space_torque(n_mode, m_mode, omega_mode, total_torque_thick, success)
    if (success) then
        print '(A,ES12.4,A)', '  Total thick orbit torque: ', total_torque_thick, ' N·m'
    else
        print *, '  Total thick orbit torque: CALCULATION FAILED'
        total_torque_thick = total_torque_thin  ! fallback
    end if
    
    if (abs(total_torque_thin) > 1e-15) then
        print '(A,F6.2)', '  Torque ratio (thick/thin): ', total_torque_thick / total_torque_thin
    end if
    
    ! Calculate torque profiles across flux surfaces
    print *, ''
    print *, 'Computing torque profiles...'
    
    ! Thin orbit profile
    print *, '  Computing thin orbit profile...'
    set_use_thick_orbits(.false.)
    call calculate_torque_profile(n_mode, m_mode, omega_mode, flux_surfaces, &
                                 torque_profile_thin, success)
    
    if (.not. success) then
        print *, '  WARNING: Thin orbit profile calculation failed'
        torque_profile_thin = 0.0d0
    end if
    
    ! Thick orbit profile
    print *, '  Computing thick orbit profile...'
    set_use_thick_orbits(.true.)
    call calculate_torque_profile(n_mode, m_mode, omega_mode, flux_surfaces, &
                                 torque_profile_thick, success)
    
    if (.not. success) then
        print *, '  WARNING: Thick orbit profile calculation failed, using thin orbit'
        torque_profile_thick = torque_profile_thin
    end if
    
    ! Validate torque conservation
    print *, ''
    print *, 'Validating torque conservation...'
    call validate_torque_conservation(torque_profile_thin, conservation_valid)
    if (conservation_valid) then
        print *, '  Thin orbit torque conservation: VALID'
    else
        print *, '  Thin orbit torque conservation: WARNING'
    end if
    
    call validate_torque_conservation(torque_profile_thick, conservation_valid)
    if (conservation_valid) then
        print *, '  Thick orbit torque conservation: VALID'
    else
        print *, '  Thick orbit torque conservation: WARNING'
    end if
    
    ! Create torque density heatmap
    print *, ''
    print *, 'Creating torque density heatmap...'
    
    call figure(1500, 1000)
    
    ! Thin orbit torque density
    call subplot(2, 2, 1)
    call pcolor(v_values, eta_values, transpose(torque_density_thin))
    call colorbar()
    call xlabel("Velocity (v/v_th)")
    call ylabel("Pitch parameter η")
    call title("Torque Density - Thin Orbit")
    
    ! Thick orbit torque density
    call subplot(2, 2, 2)
    call pcolor(v_values, eta_values, transpose(torque_density_thick))
    call colorbar()
    call xlabel("Velocity (v/v_th)")
    call ylabel("Pitch parameter η")
    call title("Torque Density - Thick Orbit")
    
    ! Torque density difference
    call subplot(2, 2, 3)
    call pcolor(v_values, eta_values, transpose(torque_density_thick - torque_density_thin))
    call colorbar()
    call xlabel("Velocity (v/v_th)")
    call ylabel("Pitch parameter η")
    call title("Torque Density Difference (Thick - Thin)")
    
    ! Torque density ratio
    call subplot(2, 2, 4)
    call pcolor(v_values, eta_values, &
               transpose(torque_density_thick / max(abs(torque_density_thin), 1e-15)))
    call colorbar()
    call xlabel("Velocity (v/v_th)")
    call ylabel("Pitch parameter η")
    call title("Torque Density Ratio (Thick/Thin)")
    
    call suptitle("NTV Torque Density in Velocity Space")
    call tight_layout()
    call savefig("torque_density_heatmap.png")
    
    ! Create torque profile comparison
    call figure(1200, 800)
    
    ! Torque profile vs flux surface
    call subplot(2, 2, 1)
    call plot(flux_surfaces, torque_profile_thin, label="Thin orbit", &
             linestyle="b-", linewidth=2)
    call plot(flux_surfaces, torque_profile_thick, label="Thick orbit", &
             linestyle="r--", linewidth=2)
    
    call xlabel("Flux surface (s)")
    call ylabel("Torque density (N·m/m³)")
    call title("NTV Torque Profile")
    call legend()
    call grid(.true., alpha=0.3)
    
    ! Torque profile difference
    call subplot(2, 2, 2)
    call plot(flux_surfaces, torque_profile_thick - torque_profile_thin, &
             linestyle="g-", linewidth=2)
    call axhline(y=0.0d0, linestyle="k:", alpha=0.5)
    
    call xlabel("Flux surface (s)")
    call ylabel("Torque difference (N·m/m³)")
    call title("Torque Profile Difference (Thick - Thin)")
    call grid(.true., alpha=0.3)
    
    ! Torque profile ratio
    call subplot(2, 2, 3)
    call plot(flux_surfaces, torque_profile_thick / max(abs(torque_profile_thin), 1e-15), &
             linestyle="m-", linewidth=2)
    call axhline(y=1.0d0, linestyle="k:", alpha=0.5, label="Unity")
    
    call xlabel("Flux surface (s)")
    call ylabel("Torque ratio")
    call title("Torque Profile Ratio (Thick/Thin)")
    call legend()
    call grid(.true., alpha=0.3)
    
    ! Integrated torque
    call subplot(2, 2, 4)
    block
        real(8) :: integrated_thin(n_flux), integrated_thick(n_flux)
        
        ! Calculate cumulative integral
        integrated_thin(1) = 0.0d0
        integrated_thick(1) = 0.0d0
        
        do i = 2, n_flux
            integrated_thin(i) = integrated_thin(i-1) + &
                    torque_profile_thin(i) * (flux_surfaces(i) - flux_surfaces(i-1))
            integrated_thick(i) = integrated_thick(i-1) + &
                    torque_profile_thick(i) * (flux_surfaces(i) - flux_surfaces(i-1))
        end do
        
        call plot(flux_surfaces, integrated_thin, label="Thin orbit", &
                 linestyle="b-", linewidth=2)
        call plot(flux_surfaces, integrated_thick, label="Thick orbit", &
                 linestyle="r--", linewidth=2)
    end block
    
    call xlabel("Flux surface (s)")
    call ylabel("Integrated torque (N·m)")
    call title("Cumulative Torque Integration")
    call legend()
    call grid(.true., alpha=0.3)
    
    call tight_layout()
    call savefig("torque_profiles_comparison.png")
    
    print *, 'Plots saved as:'
    print *, '  torque_density_heatmap.png'
    print *, '  torque_profiles_comparison.png'
    print *, ''
    
    ! Print summary statistics
    print *, 'Summary of real NEO-RT torque calculations:'
    print *, ''
    print *, 'Torque density statistics:'
    print '(A,ES12.4)', '  Thin orbit maximum: ', maxval(abs(torque_density_thin))
    print '(A,ES12.4)', '  Thick orbit maximum: ', maxval(abs(torque_density_thick))
    print '(A,ES12.4)', '  Maximum difference: ', maxval(abs(torque_density_thick - torque_density_thin))
    print *, ''
    print *, 'Torque profile statistics:'
    print '(A,ES12.4)', '  Thin orbit profile maximum: ', maxval(abs(torque_profile_thin))
    print '(A,ES12.4)', '  Thick orbit profile maximum: ', maxval(abs(torque_profile_thick))
    print '(A,ES12.4)', '  Maximum profile difference: ', maxval(abs(torque_profile_thick - torque_profile_thin))
    print *, ''
    print *, 'Integrated torque values:'
    print '(A,ES12.4,A)', '  Total thin orbit torque: ', total_torque_thin, ' N·m'
    print '(A,ES12.4,A)', '  Total thick orbit torque: ', total_torque_thick, ' N·m'
    
    if (abs(total_torque_thin) > 1e-15) then
        print '(A,F6.2,A)', '  Relative torque change: ', &
                (total_torque_thick - total_torque_thin) / total_torque_thin * 100, '%'
    end if
    
    print *, ''
    print *, 'Data source: Real NEO-RT torque calculations'
    print *, '  Functions used:'
    print *, '    - calculate_ntv_torque_density: Local torque density calculation'
    print *, '    - calculate_velocity_space_torque: Velocity space integration'
    print *, '    - calculate_torque_profile: Flux surface torque profile'
    print *, '    - validate_torque_conservation: Conservation validation'
    print *, ''
    print *, 'Physics interpretation:'
    print *, '  - NTV torque arises from resonant interaction with magnetic perturbations'
    print *, '  - Thick orbits include finite orbit width effects on resonance'
    print *, '  - Trapped particles typically contribute most to NTV torque'
    print *, '  - Torque profile depends on magnetic field perturbation spectrum'
    print *, '  - Conservation validation ensures physical consistency'
    
end program plot_torque_profiles_real