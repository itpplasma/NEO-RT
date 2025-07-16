program test_thick_orbit_torque
    ! Test for NTV torque calculation with thick orbit corrections
    ! Following TDD methodology - write failing test first
    
    use neort_resonance
    use transport_thick
    use torque_thick
    use runtime_config, only: set_use_thick_orbits
    implicit none
    
    ! Test parameters
    real(8), parameter :: pi = 3.141592653589793d0
    real(8), parameter :: v_thermal = 1.0d6  ! m/s
    real(8), parameter :: n_density = 1.0d19  ! m^-3
    real(8), parameter :: B_field = 2.5d0     ! Tesla
    
    ! Variables
    real(8) :: torque_density_thin, torque_density_thick
    real(8) :: torque_difference, relative_difference
    real(8) :: v_test, eta_test
    integer :: n_mode, m_mode
    real(8) :: omega_mode
    logical :: test_passed
    
    print *, "======================================================="
    print *, "TEST: NTV Torque Calculation with Thick Orbit Effects"
    print *, "======================================================="
    print *, ""
    print *, "Following TDD methodology - this test will initially fail"
    print *, "until the torque calculation module is implemented."
    print *, ""
    
    ! Initialize test parameters
    v_test = v_thermal
    eta_test = 0.5d0
    n_mode = 2     ! Toroidal mode number
    m_mode = 4     ! Poloidal mode number  
    omega_mode = 0.0d0  ! Static RMP
    
    test_passed = .false.
    
    ! Test 1: Thin orbit torque calculation (baseline)
    print *, "Test 1: Thin orbit torque calculation"
    print *, "------------------------------------"
    call set_use_thick_orbits(.false.)
    call calculate_ntv_torque_density(v_test, eta_test, n_mode, m_mode, omega_mode, &
                                     torque_density_thin, success=test_passed)
    
    if (test_passed) then
        print *, "✓ Thin orbit torque calculated:", torque_density_thin, "N/m³"
    else
        print *, "✗ Thin orbit torque calculation failed"
        print *, "  Expected: Function not implemented yet (TDD)"
        print *, "  This is the correct behavior for TDD - implement next!"
    end if
    
    ! Test 2: Thick orbit torque calculation
    print *, ""
    print *, "Test 2: Thick orbit torque calculation"
    print *, "-------------------------------------"
    call set_use_thick_orbits(.true.)
    call calculate_ntv_torque_density(v_test, eta_test, n_mode, m_mode, omega_mode, &
                                     torque_density_thick, success=test_passed)
    
    if (test_passed) then
        print *, "✓ Thick orbit torque calculated:", torque_density_thick, "N/m³"
    else
        print *, "✗ Thick orbit torque calculation failed"
        print *, "  Expected: Function not implemented yet (TDD)"
        print *, "  This is the correct behavior for TDD - implement next!"
    end if
    
    ! Test 3: Compare thin vs thick orbit results
    print *, ""
    print *, "Test 3: Torque comparison (thin vs thick orbits)"
    print *, "-----------------------------------------------"
    
    if (test_passed) then
        torque_difference = torque_density_thick - torque_density_thin
        relative_difference = abs(torque_difference) / abs(torque_density_thin) * 100.0d0
        
        print *, "Thin orbit torque:   ", torque_density_thin, "N/m³"
        print *, "Thick orbit torque:  ", torque_density_thick, "N/m³"
        print *, "Absolute difference: ", torque_difference, "N/m³"
        print *, "Relative difference: ", relative_difference, "%"
        
        if (relative_difference > 5.0d0) then
            print *, "✓ Significant thick orbit correction detected (>5%)"
        else
            print *, "⚠ Small thick orbit correction (<5%)"
        end if
    else
        print *, "Cannot compare - torque calculation not implemented"
        print *, "Expected behavior for TDD approach"
    end if
    
    ! Test 4: Velocity space integration
    print *, ""
    print *, "Test 4: Velocity space integration"
    print *, "----------------------------------"
    
    call test_velocity_space_integration()
    
    ! Test 5: Resonance-torque coupling
    print *, ""
    print *, "Test 5: Resonance-torque coupling"
    print *, "--------------------------------"
    
    call test_resonance_torque_coupling(n_mode, m_mode, omega_mode)
    
    ! Test Summary
    print *, ""
    print *, "======================================================="
    print *, "TEST SUMMARY"
    print *, "======================================================="
    
    if (.not. test_passed) then
        print *, "✓ TDD Test Status: CORRECTLY FAILING"
        print *, "  This is expected behavior - implement torque module next"
        print *, ""
        print *, "NEXT IMPLEMENTATION STEPS:"
        print *, "1. Create src/torque_thick.f90 module"
        print *, "2. Implement calculate_ntv_torque_density() function"
        print *, "3. Integrate: field → orbit → frequency → resonance → transport → torque"
        print *, "4. Add velocity space integration with orbit averaging"
        print *, "5. Run test again to verify GREEN state"
        print *, ""
        print *, "Physics Pipeline Required:"
        print *, "- Resonance identification (✓ completed)"
        print *, "- Transport coefficients (✓ completed)"
        print *, "- Velocity space integration (pending)"
        print *, "- Torque density calculation (pending)"
    else
        print *, "✗ Unexpected: Test should fail initially (TDD violation)"
    end if
    
contains

    ! Note: calculate_ntv_torque_density is now provided by torque_thick module
    
    subroutine test_velocity_space_integration()
        ! Test velocity space integration with orbit averaging
        use torque_thick, only: calculate_velocity_space_torque
        implicit none
        
        real(8) :: total_torque
        logical :: calc_success
        
        print *, "Testing velocity space integration framework..."
        print *, "  Required: ∫∫ f₀(v,η) × resonance_condition × transport_coeffs dv dη"
        
        ! Test velocity space integration
        call calculate_velocity_space_torque(2, 4, 0.0d0, total_torque, calc_success)
        
        if (calc_success) then
            print *, "  ✓ Velocity space integration successful"
            print *, "  Total torque:", total_torque, "N·m"
        else
            print *, "  ✗ Velocity space integration failed"
        end if
        
    end subroutine test_velocity_space_integration
    
    subroutine test_resonance_torque_coupling(n_mode, m_mode, omega_mode)
        ! Test coupling between resonance analysis and torque calculation
        implicit none
        integer, intent(in) :: n_mode, m_mode
        real(8), intent(in) :: omega_mode
        
        print *, "Testing resonance-torque coupling..."
        print *, "  Mode: n =", n_mode, ", m =", m_mode
        print *, "  Coupling: resonance finder → transport matrix → torque"
        print *, "  Status: NOT IMPLEMENTED (TDD)"
        print *, "  Next: Integrate resonance.f90 with torque calculation"
        
    end subroutine test_resonance_torque_coupling

end program test_thick_orbit_torque