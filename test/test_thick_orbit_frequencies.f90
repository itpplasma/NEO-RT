program test_thick_orbit_frequencies
    ! Test thick orbit frequency calculation module (Phase G.2.1)
    ! Create failing test for thick orbit canonical frequencies
    
    implicit none
    
    ! Test parameters
    real(8), parameter :: v_thermal = 1.0d6  ! m/s
    
    ! Variables
    real(8) :: v_test, eta_test
    real(8) :: Om_th_thick, Om_ph_thick  ! Thick orbit frequencies
    real(8) :: Om_th_thin, Om_ph_thin    ! Thin orbit frequencies (reference)
    real(8) :: frequency_ratio_th, frequency_ratio_ph
    logical :: success
    integer :: i
    
    print *, '======================================================='
    print *, 'TEST: Thick Orbit Frequency Calculation (Phase G.2.1)'
    print *, '======================================================='
    
    ! This test MUST FAIL initially - we haven't implemented freq_thick yet
    
    ! Test 1: Basic thick orbit frequency calculation
    print *, ''
    print *, 'Test 1: Basic thick orbit frequency calls'
    print *, '-----------------------------------------'
    
    v_test = v_thermal
    eta_test = 0.5d0
    
    ! Call thick orbit frequency calculation (THIS SHOULD FAIL)
    call calculate_thick_orbit_frequencies(v_test, eta_test, Om_th_thick, Om_ph_thick, success)
    
    if (.not. success) then
        print *, '✗ EXPECTED FAILURE: Thick orbit frequencies not implemented yet'
        print *, '  This test will pass once freq_thick.f90 is created'
    else
        print *, '✓ Thick orbit frequencies calculated:'
        print *, '  Om_th =', Om_th_thick, 'rad/s'
        print *, '  Om_ph =', Om_ph_thick, 'rad/s'
    end if
    
    ! Test 2: Compare thick vs thin orbit frequencies
    print *, ''
    print *, 'Test 2: Thick vs thin frequency comparison'
    print *, '-----------------------------------------'
    
    ! Calculate thin orbit frequencies for reference
    call calculate_thin_orbit_frequencies(v_test, eta_test, Om_th_thin, Om_ph_thin, success)
    
    if (success) then
        print *, 'Thin orbit frequencies (reference):'
        print *, '  Om_th_thin =', Om_th_thin, 'rad/s'
        print *, '  Om_ph_thin =', Om_ph_thin, 'rad/s'
        
        ! Calculate expected frequency shifts
        call expected_frequency_shifts(v_test, eta_test, frequency_ratio_th, frequency_ratio_ph)
        
        print *, 'Expected thick/thin frequency ratios:'
        print *, '  Om_th_thick/Om_th_thin ≈', frequency_ratio_th
        print *, '  Om_ph_thick/Om_ph_thin ≈', frequency_ratio_ph
    else
        print *, '✗ Thin orbit frequency calculation failed'
    end if
    
    ! Test 3: Velocity scaling differences
    print *, ''
    print *, 'Test 3: Velocity scaling (thick vs thin)'
    print *, '----------------------------------------'
    
    do i = 1, 3
        v_test = v_thermal * real(i, 8)
        
        call calculate_thin_orbit_frequencies(v_test, eta_test, Om_th_thin, Om_ph_thin, success)
        
        if (success) then
            print *, 'v =', v_test/v_thermal, 'v_thermal:'
            print *, '  Thin: Om_th =', Om_th_thin, ', Om_ph =', Om_ph_thin
            print *, '  Thin scaling: Om_th*v⁻¹ =', Om_th_thin*v_thermal/v_test
        end if
    end do
    
    print *, ''
    print *, 'Expected behavior:'
    print *, '- Thin orbits: frequencies ∝ v (simple scaling)'
    print *, '- Thick orbits: complex velocity dependence due to finite orbit width'
    print *, '- Thick orbit effects become larger for higher energies'
    
    print *, ''
    print *, 'NEXT STEPS:'
    print *, '1. Create src/freq_thick.f90 module'
    print *, '2. Implement Om_th and Om_ph calculations using POTATO bounce times'
    print *, '3. Handle velocity scaling differences (no simple v scaling)'
    print *, '4. Create frequency database for interpolation'
    
contains

    subroutine calculate_thick_orbit_frequencies(v, eta, Om_th, Om_ph, success)
        ! Calculate thick orbit canonical frequencies
        use freq_thick, only: compute_canonical_frequencies_thick
        implicit none
        real(8), intent(in) :: v, eta
        real(8), intent(out) :: Om_th, Om_ph
        logical, intent(out) :: success
        
        ! Use the new freq_thick module
        call compute_canonical_frequencies_thick(v, eta, Om_th, Om_ph, success)
        
        if (.not. success) then
            print *, '  WARNING: Thick orbit frequency calculation failed, using approximation'
        end if
        
    end subroutine calculate_thick_orbit_frequencies
    
    subroutine calculate_thin_orbit_frequencies(v, eta, Om_th, Om_ph, success)
        ! Calculate thin orbit frequencies using simplified model
        implicit none
        real(8), intent(in) :: v, eta
        real(8), intent(out) :: Om_th, Om_ph
        logical, intent(out) :: success
        
        real(8) :: taub, delphi
        real(8), parameter :: pi = 3.141592653589793d0
        
        ! Simplified thin orbit model for testing
        taub = 1.0d-4 / v * 1.0d6    ! Bounce time ~ 1/v
        delphi = 0.1d0 * eta         ! Toroidal shift proportional to pitch
        
        if (taub > 0.0d0) then
            Om_th = 2.0d0 * pi / taub  ! Bounce frequency
            Om_ph = delphi / taub      ! Toroidal frequency
            success = .true.
        else
            Om_th = 0.0d0
            Om_ph = 0.0d0
            success = .false.
        end if
        
    end subroutine calculate_thin_orbit_frequencies
    
    subroutine expected_frequency_shifts(v, eta, ratio_th, ratio_ph)
        ! Calculate expected frequency shifts for thick vs thin orbits
        implicit none
        real(8), intent(in) :: v, eta
        real(8), intent(out) :: ratio_th, ratio_ph
        
        ! Theoretical estimates for frequency shifts
        real(8), parameter :: B_field = 2.5d0       ! Tesla
        real(8), parameter :: mass_amu = 2.0d0      ! Deuterium mass
        real(8), parameter :: T_keV = 10.0d0        ! Temperature
        real(8), parameter :: L_B = 0.5d0           ! Magnetic scale length
        
        real(8) :: rho_gyro, orbit_width_parameter
        
        ! Calculate gyroradius
        rho_gyro = sqrt(2.0d0 * T_keV * 1.6d-16) / (1.6d-19 * B_field) * sqrt(mass_amu * 1.66d-27)
        
        ! Orbit width parameter δr/L_B
        orbit_width_parameter = rho_gyro / L_B
        
        ! Frequency shift ~ (δr/L_B)²
        ratio_th = 1.0d0 + orbit_width_parameter**2
        ratio_ph = 1.0d0 + 0.5d0 * orbit_width_parameter**2
        
    end subroutine expected_frequency_shifts

end program test_thick_orbit_frequencies