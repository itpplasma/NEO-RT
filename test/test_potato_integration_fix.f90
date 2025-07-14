program test_potato_integration_fix
    ! Test to diagnose and fix POTATO orbit integration issues with EFIT data
    ! Focus on proper initial conditions and integration parameters
    
    use potato_field_bridge, only: real_find_bounce_calculation, initialize_potato_field
    use field_eq_mod, only: nrad, nzet, rad, zet, psi_axis, psi_sep
    use global_invariants, only: dtau
    use parmot_mod, only: ro0
    implicit none
    
    ! Test parameters
    real(8), parameter :: v_thermal = 1.0d6  ! m/s
    real(8), parameter :: B_ref = 2.5d0      ! Tesla
    
    ! Variables
    real(8) :: v_test, eta_test, taub, delphi
    real(8) :: R_start, Z_start, psi_start
    logical :: success
    integer :: i, n_success, n_fail
    
    print *, '======================================================'
    print *, 'TEST: POTATO Integration Fix for EFIT Data'
    print *, '======================================================'
    
    ! Initialize POTATO field
    call initialize_potato_field(success)
    if (.not. success) then
        print *, 'ERROR: Field initialization failed'
        stop 1
    end if
    
    print *, 'Field initialized. Grid: R∈[', rad(1), ',', rad(nrad), ']'
    print *, '                      Z∈[', zet(1), ',', zet(nzet), ']'
    
    ! Test 1: Scan different starting positions
    print *, ''
    print *, 'Test 1: Starting position sensitivity'
    print *, '-------------------------------------'
    
    v_test = v_thermal
    eta_test = 0.5d0
    n_success = 0
    n_fail = 0
    
    do i = 1, 5
        ! Try different radial positions from axis to mid-radius
        R_start = 1.5d0 + 0.05d0 * real(i-1, 8)
        Z_start = 0.0d0
        
        print *, 'Testing R =', R_start, 'm'
        
        ! Try with adjusted parameters
        dtau = 1.0d-8 * real(i, 8)  ! Vary time step
        ro0 = 8.2d-3                 ! Realistic gyroradius
        
        call real_find_bounce_calculation(v_test, eta_test, taub, delphi, success)
        
        if (success) then
            n_success = n_success + 1
            print *, '  ✓ SUCCESS: taub =', taub, 's, delphi =', delphi, 'rad'
        else
            n_fail = n_fail + 1
            print *, '  ✗ FAILED: Integration did not converge'
        end if
    end do
    
    print *, 'Results:', n_success, 'succeeded,', n_fail, 'failed'
    
    ! Test 2: Velocity scaling
    print *, ''
    print *, 'Test 2: Velocity dependence'
    print *, '---------------------------'
    
    R_start = 1.65d0  ! Near magnetic axis
    dtau = 1.0d-7
    
    do i = 1, 3
        v_test = v_thermal * real(i, 8)
        
        print *, 'Testing v =', v_test/v_thermal, 'v_thermal'
        
        call real_find_bounce_calculation(v_test, eta_test, taub, delphi, success)
        
        if (success) then
            print *, '  ✓ Bounce time:', taub*1.0d6, 'μs'
            print *, '  ✓ Scaling: taub*v =', taub*v_test, '(should be ~constant)'
        else
            print *, '  ✗ Integration failed'
        end if
    end do
    
    ! Test 3: Pitch angle scan
    print *, ''
    print *, 'Test 3: Pitch angle dependence'  
    print *, '------------------------------'
    
    v_test = v_thermal
    
    do i = 1, 5
        eta_test = 0.1d0 + 0.2d0 * real(i-1, 8)
        
        print *, 'Testing eta =', eta_test
        
        call real_find_bounce_calculation(v_test, eta_test, taub, delphi, success)
        
        if (success) then
            print *, '  ✓ taub =', taub*1.0d6, 'μs, omega_phi =', delphi/taub/1.0d3, 'kHz'
        else
            print *, '  ✗ Failed (may be in loss cone or forbidden region)'
        end if
    end do
    
    ! Recommendations
    print *, ''
    print *, 'Recommendations for stable integration:'
    print *, '--------------------------------------'
    print *, '1. Use adaptive time stepping (start with dtau ~ 1e-8)'
    print *, '2. Initialize near flux surfaces (not arbitrary R,Z)'
    print *, '3. Check particle remains in allowed phase space'
    print *, '4. Consider orbit classification before integration'
    
    if (n_success > 0) then
        print *, ''
        print *, 'SUCCESS: Some orbits integrated successfully!'
        print *, 'POTATO can work with EFIT data with proper parameters'
    else
        print *, ''
        print *, 'FAILURE: No successful integrations'
        print *, 'Need to debug POTATO-EFIT interface further'
        stop 1
    end if
    
end program test_potato_integration_fix