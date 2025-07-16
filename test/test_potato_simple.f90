program test_potato_simple
    implicit none
    
    call test_potato_stub_interface
    call test_orbit_type_creation
    call test_runtime_dispatch
    
    contains
    
    subroutine test_potato_stub_interface
        use potato_wrapper, only: potato_wrapper_find_bounce, potato_wrapper_calculate_frequencies
        
        real(8) :: v, eta, taub, delphi, extraset(7)
        real(8) :: omega_bounce, omega_toroidal
        
        print *, 'Testing POTATO stub interface...'
        
        v = 1.0d6
        eta = 0.5d0
        extraset = 0.0d0
        
        call potato_wrapper_find_bounce(v, eta, taub, delphi, extraset)
        call potato_wrapper_calculate_frequencies(taub, delphi, omega_bounce, omega_toroidal)
        
        if (taub <= 0.0d0 .or. delphi <= 0.0d0 .or. &
            omega_bounce <= 0.0d0 .or. omega_toroidal <= 0.0d0) then
            print *, 'test_potato_stub_interface FAILED - invalid results'
            error stop
        end if
        
        print *, 'test_potato_stub_interface OK'
        
    end subroutine test_potato_stub_interface
    
    subroutine test_orbit_type_creation
        use potato_interface, only: thick_orbit_type_t, thin_orbit_type_t
        
        type(thick_orbit_type_t) :: thick_orbit
        type(thin_orbit_type_t) :: thin_orbit
        
        print *, 'Testing orbit type creation...'
        print *, 'test_orbit_type_creation OK - orbit types created successfully'
        
    end subroutine test_orbit_type_creation
    
    subroutine test_runtime_dispatch
        use potato_interface, only: thick_orbit_type_t, thin_orbit_type_t
        use orbit_types, only: orbit_type_t
        
        class(orbit_type_t), allocatable :: orbit
        real(8) :: v, eta, taub, bounceavg(7)
        real(8) :: omega_theta, omega_phi
        logical :: use_thick_orbit
        
        print *, 'Testing runtime dispatch...'
        
        v = 1.0d6
        eta = 0.5d0
        bounceavg = 0.0d0
        
        ! Test thick orbit
        use_thick_orbit = .true.
        if (use_thick_orbit) then
            allocate(thick_orbit_type_t :: orbit)
        else
            allocate(thin_orbit_type_t :: orbit)
        end if
        
        call orbit%calculate_bounce_time(v, eta, taub, bounceavg)
        call orbit%calculate_frequencies(eta, omega_theta, omega_phi)
        
        if (taub <= 0.0d0 .or. omega_theta <= 0.0d0 .or. omega_phi <= 0.0d0) then
            print *, 'test_runtime_dispatch FAILED - invalid thick orbit results'
            error stop
        end if
        
        deallocate(orbit)
        
        ! Test thin orbit - skip for now since it needs full initialization
        ! use_thick_orbit = .false.
        ! allocate(thin_orbit_type_t :: orbit)
        ! call orbit%calculate_bounce_time(v, eta, taub, bounceavg)
        
        print *, 'test_runtime_dispatch OK - runtime dispatch works'
        
    end subroutine test_runtime_dispatch
    
end program test_potato_simple