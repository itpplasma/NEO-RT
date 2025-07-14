program test_orbit_comparison
    implicit none
    
    call test_orbit_type_polymorphism
    call test_frequency_comparison_framework
    call test_bounce_time_comparison_framework
    
    contains
    
    subroutine test_orbit_type_polymorphism
        use potato_interface, only: thick_orbit_type_t, thin_orbit_type_t
        use orbit_types, only: orbit_type_t
        
        class(orbit_type_t), allocatable :: orbit
        real(8) :: v, eta, taub, bounceavg(7)
        real(8) :: omega_theta, omega_phi
        integer :: orbit_type_choice
        
        print *, 'Testing orbit type polymorphism...'
        
        v = 1.0d6
        eta = 0.5d0
        bounceavg = 0.0d0
        
        ! Test polymorphic orbit selection - thick orbit only for now
        print *, '  Testing thick orbit type...'
        allocate(thick_orbit_type_t :: orbit)
        
        ! Same interface for both
        call orbit%calculate_bounce_time(v, eta, taub, bounceavg)
        call orbit%calculate_frequencies(eta, omega_theta, omega_phi)
        
        if (taub <= 0.0d0 .or. omega_theta <= 0.0d0 .or. omega_phi <= 0.0d0) then
            print *, 'test_orbit_type_polymorphism FAILED - invalid results'
            error stop
        end if
        
        deallocate(orbit)
        
        ! Thin orbit testing would require full initialization
        print *, '  (Thin orbit type testing skipped - requires full initialization)'
        
        print *, 'test_orbit_type_polymorphism OK - polymorphism works correctly'
        
    end subroutine test_orbit_type_polymorphism
    
    subroutine test_frequency_comparison_framework
        use potato_interface, only: thick_orbit_type_t
        
        type(thick_orbit_type_t) :: thick_orbit
        real(8) :: eta
        real(8) :: omega_theta_thick, omega_phi_thick
        real(8) :: omega_theta_expected, omega_phi_expected
        real(8) :: relative_error_theta, relative_error_phi
        
        print *, 'Testing frequency comparison framework...'
        
        eta = 0.3d0
        
        call thick_orbit%calculate_frequencies(eta, omega_theta_thick, omega_phi_thick)
        
        ! Expected values from stub implementation
        ! omega_bounce = 2*pi / taub, omega_toroidal = delphi / taub
        ! taub = 1e-6 * (1 + 0.1*eta) = 1e-6 * 1.03
        ! delphi = 0.05 * (1 + 0.2*eta) = 0.05 * 1.06
        omega_theta_expected = 2.0d0 * 3.14159265358979323846d0 / (1.0d-6 * 1.03d0)
        omega_phi_expected = 0.05d0 * 1.06d0 / (1.0d-6 * 1.03d0)
        
        relative_error_theta = abs(omega_theta_thick - omega_theta_expected) / omega_theta_expected
        relative_error_phi = abs(omega_phi_thick - omega_phi_expected) / omega_phi_expected
        
        print *, '  Thick orbit frequency calculation:'
        print *, '    omega_theta =', omega_theta_thick, ', expected =', omega_theta_expected
        print *, '    omega_phi =', omega_phi_thick, ', expected =', omega_phi_expected
        print *, '    Relative errors: theta =', relative_error_theta, ', phi =', relative_error_phi
        
        ! Framework should calculate expected values correctly
        if (relative_error_theta > 1.0d-12 .or. relative_error_phi > 1.0d-12) then
            print *, 'test_frequency_comparison_framework FAILED - unexpected results'
            error stop
        end if
        
        print *, 'test_frequency_comparison_framework OK - comparison framework ready'
        print *, '  (Full comparison with thin orbits requires initialization)'
        
    end subroutine test_frequency_comparison_framework
    
    subroutine test_bounce_time_comparison_framework
        use potato_interface, only: thick_orbit_type_t
        
        type(thick_orbit_type_t) :: thick_orbit
        real(8) :: v, eta
        real(8) :: taub_thick, bounceavg_thick(7)
        real(8) :: taub_expected
        real(8) :: relative_error
        
        print *, 'Testing bounce time comparison framework...'
        
        v = 1.0d6
        eta = 0.4d0
        bounceavg_thick = 0.0d0
        
        call thick_orbit%calculate_bounce_time(v, eta, taub_thick, bounceavg_thick)
        
        ! Expected value from stub implementation
        ! taub = 1.0d-6 * (1.0d0 + 0.1d0 * eta)
        taub_expected = 1.0d-6 * (1.0d0 + 0.1d0 * 0.4d0)
        
        relative_error = abs(taub_thick - taub_expected) / taub_expected
        
        print *, '  Bounce time calculation:'
        print *, '    taub_thick =', taub_thick
        print *, '    taub_expected =', taub_expected
        print *, '    Relative error =', relative_error
        
        ! Framework should calculate expected value correctly
        if (relative_error > 1.0d-12) then
            print *, 'test_bounce_time_comparison_framework FAILED - unexpected result'
            error stop
        end if
        
        ! Verify bounce averages were modified (stub adds 1.0 to each element)
        if (all(bounceavg_thick == 0.0d0)) then
            print *, 'test_bounce_time_comparison_framework FAILED - bounceavg not modified'
            error stop
        end if
        
        print *, 'test_bounce_time_comparison_framework OK - comparison framework ready'
        print *, '  (Full comparison with thin orbits requires initialization)'
        
    end subroutine test_bounce_time_comparison_framework
    
end program test_orbit_comparison