program test_frequency_dispatch
    ! Test unified frequency interface (Phase G.2.2)
    ! Test runtime selection between thin/thick calculations
    
    implicit none
    
    ! Test parameters
    real(8), parameter :: v_thermal = 1.0d6  ! m/s
    real(8), parameter :: pi = 3.141592653589793d0
    
    ! Variables
    real(8) :: v_test, eta_test
    real(8) :: Om_th_unified, Om_ph_unified
    real(8) :: Om_th_thin, Om_ph_thin
    real(8) :: Om_th_thick, Om_ph_thick
    real(8) :: orbit_width_criterion, frequency_difference
    logical :: success, use_thick_orbits
    integer :: i
    
    print *, '======================================================='
    print *, 'TEST: Unified Frequency Interface (Phase G.2.2)'
    print *, '======================================================='
    
    ! This test MUST FAIL initially - we haven't implemented unified interface yet
    
    ! Test 1: Orbit width criterion for thin/thick selection
    print *, ''
    print *, 'Test 1: Orbit width criterion calculation'
    print *, '----------------------------------------'
    
    do i = 1, 5
        v_test = v_thermal * real(i, 8) * 0.5d0  ! 0.5, 1.0, 1.5, 2.0, 2.5 v_thermal
        eta_test = 0.5d0
        
        call calculate_orbit_width_criterion(v_test, eta_test, orbit_width_criterion)
        use_thick_orbits = orbit_width_criterion > 0.01d0  ! 1% threshold
        
        print *, 'v =', v_test/v_thermal, 'v_thermal:'
        print *, '  Orbit width parameter δr/L_B =', orbit_width_criterion
        print *, '  → Use thick orbits:', use_thick_orbits
    end do
    
    ! Test 2: Unified frequency interface (THIS SHOULD FAIL)
    print *, ''
    print *, 'Test 2: Unified frequency interface calls'
    print *, '-----------------------------------------'
    
    v_test = v_thermal
    eta_test = 0.5d0
    
    ! Call unified interface (THIS SHOULD FAIL)
    call unified_frequency_calculation(v_test, eta_test, Om_th_unified, Om_ph_unified, success)
    
    if (.not. success) then
        print *, '✗ EXPECTED FAILURE: Unified frequency interface not implemented yet'
        print *, '  This test will pass once freq.f90 is modified for dispatch'
    else
        print *, '✓ Unified frequencies calculated:'
        print *, '  Om_th =', Om_th_unified, 'rad/s'
        print *, '  Om_ph =', Om_ph_unified, 'rad/s'
    end if
    
    ! Test 3: Smooth transition validation
    print *, ''
    print *, 'Test 3: Transition smoothness near threshold'
    print *, '-------------------------------------------'
    
    eta_test = 0.5d0
    
    ! Test near the transition threshold
    do i = 1, 3
        v_test = v_thermal * (0.8d0 + 0.2d0 * real(i-1, 8))  ! 0.8, 1.0, 1.2 v_thermal
        
        call calculate_orbit_width_criterion(v_test, eta_test, orbit_width_criterion)
        
        print *, 'v =', v_test/v_thermal, 'v_thermal, δr/L_B =', orbit_width_criterion
        
        ! Calculate frequencies with both methods for comparison
        call calculate_reference_frequencies(v_test, eta_test, Om_th_thin, Om_ph_thin, &
                                           Om_th_thick, Om_ph_thick)
        
        frequency_difference = abs(Om_th_thick - Om_th_thin) / Om_th_thin
        
        print *, '  Thin  Om_th =', Om_th_thin
        print *, '  Thick Om_th =', Om_th_thick  
        print *, '  Relative difference =', frequency_difference * 100.0d0, '%'
        
        if (frequency_difference > 0.05d0) then
            print *, '  → Significant difference detected'
        else
            print *, '  → Frequencies similar, transition smooth'
        end if
    end do
    
    print *, ''
    print *, 'NEXT STEPS:'
    print *, '1. Modify src/freq.f90 to dispatch between thin/thick'
    print *, '2. Implement orbit width criterion in runtime'
    print *, '3. Ensure smooth transition between regimes'
    print *, '4. Test backwards compatibility with existing code'
    
contains

    subroutine unified_frequency_calculation(v, eta, Om_th, Om_ph, success)
        ! Unified frequency interface with runtime dispatch
        use neort_freq, only: Om_th_unified, Om_ph_unified
        implicit none
        real(8), intent(in) :: v, eta
        real(8), intent(out) :: Om_th, Om_ph
        logical, intent(out) :: success
        
        real(8) :: dOmthdv, dOmthdeta, dOmphdv, dOmphdeta
        
        ! Call the new unified interface
        call Om_th_unified(v, eta, Om_th, dOmthdv, dOmthdeta)
        call Om_ph_unified(v, eta, Om_ph, dOmphdv, dOmphdeta)
        
        success = .true.
        
    end subroutine unified_frequency_calculation
    
    subroutine calculate_orbit_width_criterion(v, eta, criterion)
        ! Calculate orbit width parameter δr/L_B for thin/thick selection
        implicit none
        real(8), intent(in) :: v, eta
        real(8), intent(out) :: criterion
        
        ! Physical parameters for ASDEX Upgrade-like conditions
        real(8), parameter :: B_field = 2.5d0       ! Tesla
        real(8), parameter :: mass_amu = 2.0d0      ! Deuterium mass (amu)
        real(8), parameter :: charge = 1.0d0        ! Elementary charge
        real(8), parameter :: L_B = 0.5d0           ! Magnetic scale length (m)
        real(8), parameter :: v_thermal = 1.0d6     ! Thermal velocity (m/s)
        
        real(8) :: rho_gyro
        
        ! Simplified gyroradius calculation: ρ = mv/(qB)
        ! For thermal particles: ρ ~ v_thermal * m_deuterium / (e * B)
        rho_gyro = (v / v_thermal) * 1.66d-27 * 2.0d0 * v_thermal / (1.6d-19 * B_field)
        
        ! Orbit width parameter (should be ~ 0.001-0.01 for typical conditions)
        criterion = rho_gyro / L_B
        
    end subroutine calculate_orbit_width_criterion
    
    subroutine calculate_reference_frequencies(v, eta, Om_th_thin, Om_ph_thin, &
                                             Om_th_thick, Om_ph_thick)
        ! Calculate frequencies with both methods for comparison
        implicit none
        real(8), intent(in) :: v, eta
        real(8), intent(out) :: Om_th_thin, Om_ph_thin, Om_th_thick, Om_ph_thick
        
        real(8) :: taub, delphi
        logical :: success
        
        ! Thin orbit calculation (simplified)
        taub = 1.0d-4 / v * 1.0d6    ! Bounce time ~ 1/v
        delphi = 0.1d0 * eta         ! Toroidal shift proportional to pitch
        
        Om_th_thin = 2.0d0 * pi / taub
        Om_ph_thin = delphi / taub
        
        ! Thick orbit calculation (approximation for testing)
        call approximate_thick_frequencies(v, eta, Om_th_thick, Om_ph_thick, success)
        
        if (.not. success) then
            Om_th_thick = Om_th_thin
            Om_ph_thick = Om_ph_thin
        end if
        
    end subroutine calculate_reference_frequencies
    
    subroutine approximate_thick_frequencies(v, eta, Om_th, Om_ph, success)
        ! Approximate thick orbit frequencies for testing
        implicit none
        real(8), intent(in) :: v, eta
        real(8), intent(out) :: Om_th, Om_ph
        logical, intent(out) :: success
        
        real(8) :: taub_thin, delphi_thin, correction_factor
        real(8) :: orbit_width_param
        
        ! Calculate orbit width parameter
        call calculate_orbit_width_criterion(v, eta, orbit_width_param)
        
        ! Thin orbit baseline
        taub_thin = 1.0d-4 / v * 1.0d6
        delphi_thin = 0.1d0 * eta
        
        ! Finite orbit width correction
        correction_factor = 1.0d0 + orbit_width_param**2
        
        Om_th = 2.0d0 * pi / taub_thin * correction_factor
        Om_ph = delphi_thin / taub_thin * (1.0d0 + 0.5d0 * orbit_width_param**2)
        
        success = .true.
        
    end subroutine approximate_thick_frequencies

end program test_frequency_dispatch