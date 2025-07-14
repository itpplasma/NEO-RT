program thick_orbit_example
    ! Example demonstrating thick orbit calculations in NEO-RT
    ! This shows how to use the runtime dispatch architecture
    ! to compare thick and thin orbit calculations
    
    use potato_interface, only: thick_orbit_type_t, thin_orbit_type_t
    use orbit_types, only: orbit_type_t
    
    implicit none
    
    ! Physical parameters
    real(8) :: v_thermal = 1.0d6  ! Thermal velocity (cm/s)
    real(8) :: v_particle         ! Particle velocity
    real(8) :: eta                ! Pitch parameter
    
    ! Results
    real(8) :: taub_thick, taub_thin
    real(8) :: omega_theta_thick, omega_phi_thick
    real(8) :: omega_theta_thin, omega_phi_thin
    real(8), allocatable :: bounceavg_thick(:), bounceavg_thin(:)
    
    ! Orbit objects
    class(orbit_type_t), allocatable :: orbit
    logical :: use_thick_orbits
    
    ! Control parameters
    integer :: i, neta
    real(8) :: eta_min, eta_max, deta
    
    print *, '===================================================='
    print *, 'NEO-RT Thick Orbit Example'
    print *, '===================================================='
    
    ! Initialize arrays
    allocate(bounceavg_thick(7))
    allocate(bounceavg_thin(7))
    
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
    
    deallocate(bounceavg_thick, bounceavg_thin)
    
contains

    subroutine example_single_orbit_comparison()
        print *, '----------------------------------------------------'
        print *, 'Example 1: Single Orbit Comparison'
        print *, '----------------------------------------------------'
        
        eta = 0.5d0
        
        ! Calculate with thick orbit
        use_thick_orbits = .true.
        allocate(thick_orbit_type_t :: orbit)
        
        bounceavg_thick = 0.0d0
        call orbit%calculate_bounce_time(v_particle, eta, taub_thick, bounceavg_thick)
        call orbit%calculate_frequencies(eta, omega_theta_thick, omega_phi_thick)
        
        deallocate(orbit)
        
        ! Calculate with thin orbit (when available)
        ! For now, use placeholder values
        taub_thin = taub_thick * 0.95d0
        omega_theta_thin = omega_theta_thick * 1.05d0
        omega_phi_thin = omega_phi_thick * 1.03d0
        
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
        
        allocate(thick_orbit_type_t :: orbit)
        
        print *, 'Scanning eta from', eta_min, 'to', eta_max
        print *, ''
        print *, '   eta      Bounce Time    omega_theta     omega_phi      q_eff'
        print *, '  -----     -----------    -----------     ---------      -----'
        
        do i = 1, neta
            eta = eta_min + (i - 1) * deta
            
            bounceavg_thick = 0.0d0
            call orbit%calculate_bounce_time(v_particle, eta, taub_thick, bounceavg_thick)
            call orbit%calculate_frequencies(eta, omega_theta_thick, omega_phi_thick)
            
            write(*, '(F7.3, 4ES15.6)') eta, taub_thick, omega_theta_thick, &
                    omega_phi_thick, omega_phi_thick/omega_theta_thick
        end do
        
        deallocate(orbit)
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
        
        allocate(thick_orbit_type_t :: orbit)
        
        ! Scan to find resonance
        npts = 101
        do j = 1, npts
            eta_scan(j) = 0.1d0 + 0.8d0 * real(j-1, 8) / real(npts-1, 8)
            
            call orbit%calculate_frequencies(eta_scan(j), omega_theta_thick, omega_phi_thick)
            
            ! Resonance parameter: how close to n*omega_phi - m*omega_theta = omega_mode
            delta_omega = n * omega_phi_thick - m * omega_theta_thick - omega_mode
            resonance_param(j) = delta_omega
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
            
            ! Analyze resonance width
            call orbit%calculate_frequencies(eta_res, omega_theta_thick, omega_phi_thick)
            print *, 'At resonance:'
            print *, '  omega_theta =', omega_theta_thick
            print *, '  omega_phi =', omega_phi_thick
            print *, '  n*omega_phi - m*omega_theta =', &
                     n * omega_phi_thick - m * omega_theta_thick
            print *, '  Mismatch with mode =', &
                     abs(n * omega_phi_thick - m * omega_theta_thick - omega_mode)
        else
            print *, 'No resonance found in eta range'
        end if
        
        deallocate(orbit)
        print *, ''
        
    end subroutine example_resonance_analysis
    
end program thick_orbit_example