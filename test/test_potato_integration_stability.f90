program test_potato_integration_stability
    ! Test stable POTATO integration (Phase G.4.REAL.1)
    ! Debug floating point exceptions in POTATO find_bounce calls
    
    use potato_field_bridge, only: real_find_bounce_calculation
    implicit none
    
    ! Test parameters
    real(8), parameter :: pi = 3.141592653589793d0
    real(8), parameter :: v_thermal = 1.0d6  ! m/s
    
    ! Variables
    real(8) :: v_test, eta_test
    real(8) :: taub, delphi
    logical :: success
    integer :: i, j, success_count, total_count
    real(8) :: success_rate
    
    print *, '================================================================='
    print *, 'TEST: POTATO Integration Stability (Phase G.4.REAL.1)'
    print *, '================================================================='
    
    ! This test will initially FAIL due to floating point exceptions
    ! Goal: Make POTATO find_bounce work reliably across (v,η) parameter space
    
    print *, ''
    print *, 'Testing POTATO integration stability across parameter space...'
    print *, ''
    
    success_count = 0
    total_count = 0
    
    ! Test 1: Basic parameter scan
    print *, 'Test 1: Parameter space scan (v,η)'
    print *, '-----------------------------------'
    print *, 'v/v_thermal    eta      Success    τ_bounce(s)    Δφ(rad)'
    print *, '--------------------------------------------------------'
    
    do i = 1, 5  ! velocity scan
        do j = 1, 5  ! pitch angle scan
            v_test = v_thermal * (0.5d0 + 0.3d0 * real(i-1, 8))  ! 0.5 to 1.7 v_thermal
            eta_test = 0.1d0 + 0.8d0 * real(j-1, 8) / 4.0d0      ! 0.1 to 0.9
            
            total_count = total_count + 1
            
            ! Try POTATO integration - this will likely fail initially
            call safe_potato_integration(v_test, eta_test, taub, delphi, success)
            
            if (success) then
                success_count = success_count + 1
                print '(F10.2, F10.2, A10, E12.4, F10.4)', &
                      v_test/v_thermal, eta_test, '    ✓    ', taub, delphi
            else
                print '(F10.2, F10.2, A10, A12, A10)', &
                      v_test/v_thermal, eta_test, '    ✗    ', '   FAILED   ', '  FAILED'
            end if
        end do
    end do
    
    success_rate = real(success_count, 8) / real(total_count, 8) * 100.0d0
    print *, ''
    print *, 'Integration Success Rate:', success_rate, '%'
    print *, 'Successful integrations:', success_count, '/', total_count
    
    ! Test 2: Edge case testing
    print *, ''
    print *, 'Test 2: Edge cases and stability limits'
    print *, '--------------------------------------'
    
    ! Test very small eta (nearly perpendicular)
    v_test = v_thermal
    eta_test = 0.01d0
    call safe_potato_integration(v_test, eta_test, taub, delphi, success)
    if (success) then
        print *, 'Small η (0.01): SUCCESS - τb =', taub, 's, Δφ =', delphi, 'rad'
    else
        print *, 'Small η (0.01): FAILED - numerical issues'
    end if
    
    ! Test large eta (nearly parallel)
    eta_test = 0.99d0
    call safe_potato_integration(v_test, eta_test, taub, delphi, success)
    if (success) then
        print *, 'Large η (0.99): SUCCESS - τb =', taub, 's, Δφ =', delphi, 'rad'
    else
        print *, 'Large η (0.99): FAILED - numerical issues'
    end if
    
    ! Test very high energy
    v_test = 3.0d0 * v_thermal
    eta_test = 0.5d0
    call safe_potato_integration(v_test, eta_test, taub, delphi, success)
    if (success) then
        print *, 'High energy (3v_th): SUCCESS - τb =', taub, 's, Δφ =', delphi, 'rad'
    else
        print *, 'High energy (3v_th): FAILED - numerical issues'
    end if
    
    ! Test 3: Integration parameter diagnosis
    print *, ''
    print *, 'Test 3: Integration parameter diagnosis'
    print *, '-------------------------------------'
    
    ! This will help identify optimal integration parameters
    call diagnose_integration_parameters()
    
    print *, ''
    print *, 'DIAGNOSIS RESULTS:'
    print *, '=================='
    if (success_rate < 50.0d0) then
        print *, '❌ CRITICAL: Low success rate indicates fundamental integration issues'
        print *, 'Likely causes:'
        print *, '  - Inappropriate time step (dtau) for orbit integration'
        print *, '  - Invalid initial conditions (particles not on flux surfaces)'
        print *, '  - Orbit closure tolerance too strict'
        print *, '  - Field evaluation errors at orbit boundaries'
    else if (success_rate < 90.0d0) then
        print *, '⚠️  WARNING: Moderate success rate, some parameter ranges problematic'
        print *, 'Need parameter-dependent adaptive integration'
    else
        print *, '✅ GOOD: High success rate, POTATO integration stable'
    end if
    
    print *, ''
    print *, 'NEXT STEPS FOR FIXING INTEGRATION:'
    print *, '1. Implement adaptive time stepping based on particle energy'
    print *, '2. Add orbit classification to avoid forbidden regions'  
    print *, '3. Improve initial condition validation'
    print *, '4. Add graceful fallback for failed integrations'
    print *, '5. Optimize integration tolerances for realistic fields'
    
contains

    subroutine safe_potato_integration(v, eta, taub, delphi, success)
        ! Wrapper for POTATO integration with error handling
        implicit none
        real(8), intent(in) :: v, eta
        real(8), intent(out) :: taub, delphi
        logical, intent(out) :: success
        
        ! Initialize outputs
        success = .false.
        taub = 0.0d0
        delphi = 0.0d0
        
        ! Try the real POTATO integration
        ! This will fail for many parameters initially
        call real_find_bounce_calculation(v, eta, taub, delphi, success)
        
        ! Additional validation of results
        if (success) then
            ! Check for physically reasonable results
            if (taub <= 0.0d0 .or. taub > 1.0d-3) then
                success = .false.  ! Unphysical bounce time
            end if
            
            if (abs(delphi) > 2.0d0 * pi) then
                success = .false.  ! Unphysical toroidal shift
            end if
        end if
        
    end subroutine safe_potato_integration
    
    subroutine diagnose_integration_parameters()
        ! Diagnose optimal integration parameters for stability
        implicit none
        
        print *, 'Current POTATO integration parameters:'
        print *, '  - Integration method: Runge-Kutta adaptive'
        print *, '  - Typical time step: dtau ~ 1e-5'
        print *, '  - Orbit closure tolerance: relerr ~ 1e-10'
        print *, '  - Energy conservation: toten from particle energy'
        print *, '  - Magnetic moment: perpinv from pitch angle'
        print *, ''
        print *, 'Parameter tuning needed for complex EFIT fields:'
        print *, '  - Adaptive dtau based on local field gradients'
        print *, '  - Relaxed tolerances for high-energy particles'
        print *, '  - Better initial condition selection'
        print *, '  - Orbit classification (trapped/passing/potato)'
        
    end subroutine diagnose_integration_parameters

end program test_potato_integration_stability