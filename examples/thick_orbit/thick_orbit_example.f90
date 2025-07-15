program thick_orbit_example
    ! Example demonstrating thick orbit calculations in NEO-RT
    ! This shows how to use the real POTATO integration
    ! to compare thick and thin orbit calculations with plots
    
    use potato_field_bridge, only: real_find_bounce_calculation, &
                                   initialize_potato_field
    use fortplotlib, only: figure_t
    
    implicit none
    
    ! Physical parameters
    real(8) :: v_thermal = 1.0d6  ! Thermal velocity (cm/s)
    real(8) :: v_particle         ! Particle velocity
    real(8) :: eta                ! Pitch parameter
    
    ! Results
    real(8) :: taub_thick, taub_thin, delphi_thick, delphi_thin
    real(8) :: omega_theta_thick, omega_phi_thick
    real(8) :: omega_theta_thin, omega_phi_thin
    logical :: success_thick, success_thin
    
    ! Control parameters
    integer :: i, neta
    real(8) :: eta_min, eta_max, deta
    
    print *, '===================================================='
    print *, 'NEO-RT Thick Orbit Example with Real POTATO'
    print *, '===================================================='
    
    ! Initialize POTATO field
    call initialize_potato_field(success_thick)
    if (.not. success_thick) then
        print *, 'ERROR: Failed to initialize POTATO field'
        stop 1
    end if
    
    ! Set particle velocity
    v_particle = 2.0d0 * v_thermal  ! Fast ion
    
    print *, ''
    print *, 'Particle parameters:'
    print *, '  Velocity =', v_particle / v_thermal, 'v_thermal'
    print *, ''
    
    ! Example 1: Single orbit comparison
    call example_single_orbit_comparison()
    
    ! Example 2: Scan over pitch angles
    call example_pitch_angle_scan()
    
    ! Example 3: Resonance analysis
    call example_resonance_analysis()
    
    ! Example 4: Generate plots
    call example_generate_plots()
    
contains

    subroutine example_single_orbit_comparison()
        print *, '----------------------------------------------------'
        print *, 'Example 1: Single Orbit Comparison'
        print *, '----------------------------------------------------'
        
        eta = 0.5d0
        
        ! Calculate with thick orbit (real POTATO)
        call real_find_bounce_calculation(v_particle, eta, taub_thick, delphi_thick, success_thick)
        
        ! Calculate frequencies from POTATO results
        if (success_thick) then
            omega_phi_thick = delphi_thick / taub_thick  ! Toroidal frequency
            omega_theta_thick = 2.0d0 * 3.14159d0 / taub_thick  ! Poloidal frequency (approximation)
        else
            print *, 'ERROR: Thick orbit calculation failed'
            return
        end if
        
        ! Calculate with thin orbit (stub implementation for comparison)
        call calculate_bounce_time(v_particle, eta, taub_thin, delphi_thin, success_thin)
        
        if (success_thin) then
            omega_phi_thin = delphi_thin / taub_thin
            omega_theta_thin = 2.0d0 * 3.14159d0 / taub_thin
        else
            print *, 'ERROR: Thin orbit calculation failed'
            return
        end if
        
        print *, 'At eta =', eta, ':'
        print *, ''
        print *, '                  Thick Orbit      Thin Orbit      Ratio'
        print *, '  Bounce time:', taub_thick, taub_thin, taub_thick/taub_thin
        print *, '  omega_theta:', omega_theta_thick, omega_theta_thin, &
                 omega_theta_thick/omega_theta_thin
        print *, '  omega_phi:  ', omega_phi_thick, omega_phi_thin, &
                 omega_phi_thick/omega_phi_thin
        print *, ''
        
    end subroutine example_single_orbit_comparison
    
    subroutine example_pitch_angle_scan()
        print *, '----------------------------------------------------'
        print *, 'Example 2: Pitch Angle Scan'
        print *, '----------------------------------------------------'
        
        neta = 11
        eta_min = 0.1d0
        eta_max = 0.9d0
        deta = (eta_max - eta_min) / real(neta - 1, 8)
        
        print *, 'Scanning eta from', eta_min, 'to', eta_max, ' with real POTATO'
        print *, ''
        print *, '   eta      Bounce Time    omega_theta     omega_phi      q_eff'
        print *, '  -----     -----------    -----------     ---------      -----'
        
        do i = 1, neta
            eta = eta_min + (i - 1) * deta
            
            call real_find_bounce_calculation(v_particle, eta, taub_thick, delphi_thick, success_thick)
            
            if (success_thick) then
                omega_phi_thick = delphi_thick / taub_thick
                omega_theta_thick = 2.0d0 * 3.14159d0 / taub_thick
                
                write(*, '(F7.3, 4ES15.6)') eta, taub_thick, omega_theta_thick, &
                        omega_phi_thick, omega_phi_thick/omega_theta_thick
            else
                write(*, '(F7.3, A)') eta, '   CALCULATION FAILED'
            end if
        end do
        print *, ''
        
    end subroutine example_pitch_angle_scan
    
    subroutine example_resonance_analysis()
        real(8) :: omega_mode, delta_omega
        integer :: m, n
        real(8) :: eta_scan(101), resonance_param(101)
        integer :: j, npts
        real(8) :: eta_res
        
        print *, '----------------------------------------------------'
        print *, 'Example 3: Resonance Analysis'
        print *, '----------------------------------------------------'
        
        ! Mode parameters
        m = 2  ! Poloidal mode number
        n = 1  ! Toroidal mode number
        omega_mode = 5.0d4  ! Mode frequency (rad/s)
        
        print *, 'Searching for resonance with (m,n) = (', m, ',', n, ')'
        print *, 'Mode frequency =', omega_mode, 'rad/s'
        print *, ''
        
        ! Scan to find resonance with real POTATO
        npts = 101
        do j = 1, npts
            eta_scan(j) = 0.1d0 + 0.8d0 * real(j-1, 8) / real(npts-1, 8)
            
            call real_find_bounce_calculation(v_particle, eta_scan(j), taub_thick, delphi_thick, success_thick)
            
            if (success_thick) then
                omega_phi_thick = delphi_thick / taub_thick
                omega_theta_thick = 2.0d0 * 3.14159d0 / taub_thick
                
                ! Resonance parameter: how close to n*omega_phi - m*omega_theta = omega_mode
                delta_omega = n * omega_phi_thick - m * omega_theta_thick - omega_mode
                resonance_param(j) = delta_omega
            else
                resonance_param(j) = 1.0d10  ! Large value to indicate failure
            end if
        end do
        
        ! Find approximate resonance location (where parameter changes sign)
        eta_res = 0.0d0
        do j = 2, npts
            if (resonance_param(j-1) * resonance_param(j) < 0.0d0) then
                ! Linear interpolation to find zero crossing
                eta_res = eta_scan(j-1) - resonance_param(j-1) * &
                         (eta_scan(j) - eta_scan(j-1)) / &
                         (resonance_param(j) - resonance_param(j-1))
                exit
            end if
        end do
        
        if (eta_res > 0.0d0) then
            print *, 'Resonance found at eta =', eta_res
            
            ! Analyze resonance width with real POTATO
            call real_find_bounce_calculation(v_particle, eta_res, taub_thick, delphi_thick, success_thick)
            if (success_thick) then
                omega_phi_thick = delphi_thick / taub_thick
                omega_theta_thick = 2.0d0 * 3.14159d0 / taub_thick
                
                print *, 'At resonance:'
                print *, '  omega_theta =', omega_theta_thick
                print *, '  omega_phi =', omega_phi_thick
                print *, '  n*omega_phi - m*omega_theta =', &
                         n * omega_phi_thick - m * omega_theta_thick
                print *, '  Mismatch with mode =', &
                         abs(n * omega_phi_thick - m * omega_theta_thick - omega_mode)
            end if
        else
            print *, 'No resonance found in eta range'
        end if
        
        print *, ''
        
    end subroutine example_resonance_analysis
    
    subroutine example_generate_plots()
        ! Generate plots using fortplotlib
        integer, parameter :: npts = 50
        real(8) :: eta_array(npts), taub_array(npts), omega_phi_array(npts)
        real(8) :: v_array(npts), taub_v_array(npts)
        integer :: j
        real(8) :: eta_plot, v_plot
        type(figure) :: fig
        
        print *, '----------------------------------------------------'
        print *, 'Example 4: Generate Plots with fortplotlib'
        print *, '----------------------------------------------------'
        
        ! Plot 1: Bounce time vs eta
        print *, 'Generating bounce time vs eta plot...'
        do j = 1, npts
            eta_plot = 0.1d0 + 0.8d0 * real(j-1, 8) / real(npts-1, 8)
            eta_array(j) = eta_plot
            
            call real_find_bounce_calculation(v_particle, eta_plot, taub_thick, delphi_thick, success_thick)
            if (success_thick) then
                taub_array(j) = taub_thick
                omega_phi_array(j) = delphi_thick / taub_thick
            else
                taub_array(j) = 0.0d0
                omega_phi_array(j) = 0.0d0
            end if
        end do
        
        call fig%initialize()
        call fig%plot(eta_array, taub_array, linewidth=2.0)
        call fig%xlabel('Pitch Parameter η')
        call fig%ylabel('Bounce Time [s]')
        call fig%title('POTATO Thick Orbit: Bounce Time vs Pitch Parameter')
        call fig%save('potato_bounce_time_vs_eta.png')
        
        ! Plot 2: Bounce time vs velocity
        print *, 'Generating bounce time vs velocity plot...'
        eta_plot = 0.5d0  ! Fixed eta
        do j = 1, npts
            v_plot = 0.5d6 + 2.0d6 * real(j-1, 8) / real(npts-1, 8)  ! 0.5-2.5 * 1e6 cm/s
            v_array(j) = v_plot
            
            call real_find_bounce_calculation(v_plot, eta_plot, taub_thick, delphi_thick, success_thick)
            if (success_thick) then
                taub_v_array(j) = taub_thick
            else
                taub_v_array(j) = 0.0d0
            end if
        end do
        
        call fig%initialize()
        call fig%plot(v_array, taub_v_array, linewidth=2.0)
        call fig%xlabel('Velocity [cm/s]')
        call fig%ylabel('Bounce Time [s]')
        call fig%title('POTATO Thick Orbit: Bounce Time vs Velocity (η=0.5)')
        call fig%save('potato_bounce_time_vs_velocity.png')
        
        print *, 'Plots saved:'
        print *, '  - potato_bounce_time_vs_eta.png'
        print *, '  - potato_bounce_time_vs_velocity.png'
        print *, ''
        
    end subroutine example_generate_plots
    
end program thick_orbit_example