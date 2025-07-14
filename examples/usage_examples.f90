program usage_examples
    ! Practical examples showing how to use the thick orbit framework
    ! for immediate physics research applications
    
    use potato_interface, only: thick_orbit_type_t, thin_orbit_type_t
    use orbit_types, only: orbit_type_t
    
    implicit none
    
    ! Example 1: Compare orbit types for single particle
    call example_orbit_comparison()
    
    ! Example 2: Parameter scan with automatic orbit selection
    call example_parameter_scan()
    
    ! Example 3: Resonance width study
    call example_resonance_study()
    
contains

    subroutine example_orbit_comparison()
        ! Direct comparison of thin vs thick orbit physics
        
        type(thick_orbit_type_t) :: thick_orbit
        type(thin_orbit_type_t) :: thin_orbit
        
        real(8) :: v, eta
        real(8) :: taub_thick, taub_thin
        real(8) :: omega_theta_thick, omega_theta_thin
        real(8) :: omega_phi_thick, omega_phi_thin
        real(8), allocatable :: bounceavg_thick(:), bounceavg_thin(:)
        
        print *, '=== Orbit Type Comparison ==='
        
        v = 2.0d0  ! 2 * v_thermal
        eta = 0.5d0
        
        allocate(bounceavg_thick(7), bounceavg_thin(7))
        
        ! Calculate with thick orbit
        call thick_orbit%calculate_bounce_time(v, eta, taub_thick, bounceavg_thick)
        call thick_orbit%calculate_frequencies(eta, omega_theta_thick, omega_phi_thick)
        
        ! Calculate with thin orbit (when available)
        ! Currently using stub, but interface is ready
        ! call thin_orbit%calculate_bounce_time(v, eta, taub_thin, bounceavg_thin)
        ! call thin_orbit%calculate_frequencies(eta, omega_theta_thin, omega_phi_thin)
        
        print *, 'Thick orbit results:'
        print *, '  Bounce time:', taub_thick
        print *, '  omega_theta:', omega_theta_thick
        print *, '  omega_phi:', omega_phi_thick
        print *, '  q_eff:', omega_phi_thick/omega_theta_thick
        
        deallocate(bounceavg_thick, bounceavg_thin)
        
    end subroutine example_orbit_comparison
    
    subroutine example_parameter_scan()
        ! Scan over velocity and pitch angle
        
        class(orbit_type_t), allocatable :: orbit
        real(8) :: v, eta, taub, omega_theta, omega_phi
        real(8), allocatable :: bounceavg(:)
        integer :: i, j, nv, neta
        logical :: use_thick_orbits
        
        print *, '=== Parameter Scan Example ==='
        
        nv = 5
        neta = 10
        allocate(bounceavg(7))
        
        ! For this example, always use thick orbits
        use_thick_orbits = .true.
        allocate(thick_orbit_type_t :: orbit)
        
        print *, 'v/v_th    eta      taub        omega_theta   omega_phi    q_eff'
        print *, '------    ---      ----        -----------   ---------    -----'
        
        do i = 1, nv
            v = 1.0d0 + real(i-1, 8) * 2.0d0 / real(nv-1, 8)  ! 1 to 3 v_thermal
            
            do j = 1, neta
                eta = 0.1d0 + real(j-1, 8) * 0.8d0 / real(neta-1, 8)  ! 0.1 to 0.9
                
                call orbit%calculate_bounce_time(v, eta, taub, bounceavg)
                call orbit%calculate_frequencies(eta, omega_theta, omega_phi)
                
                write(*, '(F6.2, 4X, F6.3, 2X, ES10.3, 2X, ES11.4, 2X, ES9.2, 2X, F8.5)') &
                    v, eta, taub, omega_theta, omega_phi, omega_phi/omega_theta
            end do
        end do
        
        deallocate(bounceavg)
        
    end subroutine example_parameter_scan
    
    subroutine example_resonance_study()
        ! Study resonance conditions across parameter space
        
        type(thick_orbit_type_t) :: thick_orbit
        real(8) :: v, eta, omega_theta, omega_phi
        real(8) :: omega_mode, resonance_parameter
        integer :: m, n, i
        integer, parameter :: neta = 50
        real(8), parameter :: eta_min = 0.1d0, eta_max = 0.9d0
        integer :: resonance_count
        
        print *, '=== Resonance Study Example ==='
        
        ! Mode parameters (example: 2/1 mode)
        m = 2
        n = 1
        omega_mode = 5.0d4  ! rad/s
        v = 2.5d0  ! Fast ion
        
        print *, 'Scanning for', m, '/', n, 'resonances at v =', v, 'v_thermal'
        print *, 'Mode frequency =', omega_mode, 'rad/s'
        print *, ''
        
        resonance_count = 0
        
        do i = 1, neta
            eta = eta_min + real(i-1, 8) * (eta_max - eta_min) / real(neta-1, 8)
            
            call thick_orbit%calculate_frequencies(eta, omega_theta, omega_phi)
            
            ! Resonance condition: n*omega_phi - m*omega_theta = omega_mode
            resonance_parameter = n * omega_phi - m * omega_theta - omega_mode
            
            ! Check if close to resonance (within 1% of mode frequency)
            if (abs(resonance_parameter) < 0.01d0 * omega_mode) then
                resonance_count = resonance_count + 1
                print *, 'Resonance found:'
                print *, '  eta =', eta
                print *, '  omega_theta =', omega_theta
                print *, '  omega_phi =', omega_phi
                print *, '  resonance parameter =', resonance_parameter
                print *, '  relative error =', abs(resonance_parameter)/omega_mode
                print *, ''
            endif
        end do
        
        print *, 'Total resonances found:', resonance_count
        
    end subroutine example_resonance_study
    
end program usage_examples