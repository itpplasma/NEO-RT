program test_drift_standalone
    ! Test the thick orbit drift calculations without POTATO dependency
    ! This tests the drift physics implementation directly
    
    use thick_orbit_drift, only: gradB_drift_velocity, curvature_drift_velocity, magnetic_moment_conservation
    implicit none
    
    ! Test parameters
    real(8), parameter :: pi = 3.141592653589793d0
    real(8), parameter :: v_thermal = 1.0d6  ! m/s
    
    ! Variables
    real(8) :: v_test, eta_test
    real(8) :: R_test, Z_test, phi_test, B_test
    real(8) :: v_gradB(3), v_curv(3), mu
    logical :: success
    integer :: test_count, success_count
    
    print *, '================================================'
    print *, 'TEST: Drift Calculations Standalone (G.4.REAL.2)'
    print *, '================================================'
    
    test_count = 0
    success_count = 0
    
    ! Set up test parameters
    v_test = v_thermal
    eta_test = 0.5d0
    R_test = 1.5d0
    Z_test = 0.0d0
    phi_test = 0.0d0
    B_test = 2.5d0
    
    print *, ''
    print *, 'Test 1: Grad-B drift velocity'
    print *, '-----------------------------'
    print *, 'Input parameters:'
    print *, '  v =', v_test, 'm/s'
    print *, '  eta =', eta_test
    print *, '  R =', R_test, 'm'
    print *, '  Z =', Z_test, 'm'
    print *, '  B =', B_test, 'T'
    
    test_count = test_count + 1
    call gradB_drift_velocity(v_test, eta_test, R_test, Z_test, phi_test, B_test, v_gradB, success)
    
    if (success) then
        success_count = success_count + 1
        print *, '✓ Grad-B drift calculated successfully:'
        print *, '  v_gradB_R   =', v_gradB(1), 'm/s'
        print *, '  v_gradB_phi =', v_gradB(2), 'm/s'
        print *, '  v_gradB_Z   =', v_gradB(3), 'm/s'
        print *, '  |v_gradB|   =', sqrt(sum(v_gradB**2)), 'm/s'
    else
        print *, '✗ Grad-B drift calculation failed'
    end if
    
    print *, ''
    print *, 'Test 2: Curvature drift velocity'
    print *, '-------------------------------'
    
    test_count = test_count + 1
    call curvature_drift_velocity(v_test, eta_test, R_test, Z_test, phi_test, B_test, v_curv, success)
    
    if (success) then
        success_count = success_count + 1
        print *, '✓ Curvature drift calculated successfully:'
        print *, '  v_curv_R   =', v_curv(1), 'm/s'
        print *, '  v_curv_phi =', v_curv(2), 'm/s'
        print *, '  v_curv_Z   =', v_curv(3), 'm/s'
        print *, '  |v_curv|   =', sqrt(sum(v_curv**2)), 'm/s'
    else
        print *, '✗ Curvature drift calculation failed'
    end if
    
    print *, ''
    print *, 'Test 3: Magnetic moment conservation'
    print *, '----------------------------------'
    
    test_count = test_count + 1
    call magnetic_moment_conservation(v_test, eta_test, B_test, mu, success)
    
    if (success) then
        success_count = success_count + 1
        print *, '✓ Magnetic moment calculated successfully:'
        print *, '  μ =', mu, 'J/T'
        print *, '  Expected thermal: ~1.7e-20 J/T'
    else
        print *, '✗ Magnetic moment calculation failed'
    end if
    
    print *, ''
    print *, 'Test 4: Combined drift physics'
    print *, '-----------------------------'
    
    if (success) then
        print *, 'Total drift velocity:'
        print *, '  v_total_R   =', v_gradB(1) + v_curv(1), 'm/s'
        print *, '  v_total_phi =', v_gradB(2) + v_curv(2), 'm/s'  
        print *, '  v_total_Z   =', v_gradB(3) + v_curv(3), 'm/s'
        print *, '  |v_total|   =', sqrt(sum((v_gradB + v_curv)**2)), 'm/s'
    end if
    
    print *, ''
    print *, 'SUMMARY:'
    print *, '========'
    print *, 'Tests passed:', success_count, '/', test_count
    print *, 'Success rate:', real(success_count)/real(test_count)*100.0, '%'
    
    if (success_count == test_count) then
        print *, ''
        print *, '🎯 PHASE G.4.REAL.2 COMPONENT STATUS: ✅ WORKING'
        print *, 'Drift velocity calculations functional'
        print *, 'Ready for POTATO orbit integration'
    else
        print *, ''
        print *, '❌ PHASE G.4.REAL.2 COMPONENT STATUS: INCOMPLETE'
        print *, 'Some drift calculations failing'
    end if
    
end program test_drift_standalone