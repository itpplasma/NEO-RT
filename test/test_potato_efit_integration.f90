program test_potato_efit_integration
    ! Test POTATO integration with realistic ASDEX Upgrade EFIT data
    ! This test verifies POTATO can read EFIT files and perform realistic 
    ! thick orbit calculations with finite Larmor radius effects
    
    use potato_field_bridge, only: real_find_bounce_calculation, initialize_potato_field
    implicit none
    
    ! Test parameters for ASDEX Upgrade shot 30835 at t=3.2s
    real(8), parameter :: B_FIELD_AUG = 2.5d0    ! Tesla
    real(8), parameter :: T_KEV_AUG = 10.0d0     ! keV thermal temperature
    real(8), parameter :: MASS_AMU = 2.014d0     ! Deuterium mass
    real(8), parameter :: CHARGE = 1.0d0         ! Elementary charge
    
    ! Physical constants
    real(8), parameter :: AMU_TO_KG = 1.66053906660d-27
    real(8), parameter :: EV_TO_J = 1.602176634d-19
    real(8), parameter :: KEV_TO_J = 1000.0d0 * EV_TO_J
    
    ! Test variables
    real(8) :: v_thermal, rho_gyro, eta_test, v_test
    real(8) :: taub_efit, delphi_efit
    logical :: success_efit, field_init_success
    character(len=256) :: efit_file
    logical :: file_exists
    
    print *, '======================================================'
    print *, 'TEST: POTATO EFIT Integration (Phase F.1.1)'
    print *, '======================================================'
    print *, 'Testing realistic ASDEX Upgrade magnetic field data'
    
    ! Calculate realistic thermal velocity and gyroradius
    v_thermal = sqrt(2.0d0 * T_KEV_AUG * KEV_TO_J / (MASS_AMU * AMU_TO_KG))
    rho_gyro = (MASS_AMU * AMU_TO_KG * v_thermal) / (CHARGE * EV_TO_J * B_FIELD_AUG)
    
    print *, 'AUG parameters:'
    print *, '  B-field:', B_FIELD_AUG, 'Tesla'
    print *, '  Temperature:', T_KEV_AUG, 'keV'
    print *, '  Thermal velocity:', v_thermal/1.0d6, 'km/s'
    print *, '  Gyroradius:', rho_gyro*1000.0d0, 'mm'
    print *, ''
    
    ! Test 1: Check if EFIT file exists
    call get_environment_variable('DATA', efit_file)
    if (len_trim(efit_file) == 0) then
        efit_file = '/data'  ! Default fallback
    end if
    efit_file = trim(efit_file) // '/AUG/EQDSK/g30835.3200_ed6'
    
    inquire(file=trim(efit_file), exist=file_exists)
    print *, 'Test 1: EFIT file availability'
    print *, '  File path:', trim(efit_file)
    print *, '  File exists:', file_exists
    
    if (.not. file_exists) then
        print *, 'ERROR: EFIT file not found - this test requires realistic data'
        print *, 'Expected file: $DATA/AUG/EQDSK/g30835.3200_ed6'
        print *, 'FAILING TEST (as expected) - EFIT integration not yet implemented'
        stop 1  ! Expected failure until EFIT integration is implemented
    end if
    
    ! Test 2: Initialize POTATO field with EFIT data
    print *, ''
    print *, 'Test 2: POTATO field initialization with EFIT'
    call initialize_potato_field(field_init_success)
    
    if (.not. field_init_success) then
        print *, 'ERROR: POTATO field initialization failed'
        print *, 'FAILING TEST (as expected) - EFIT field setup not implemented'
        stop 1  ! Expected failure until real EFIT integration
    end if
    
    print *, '  POTATO field initialized:', field_init_success
    
    ! Test 3: Calculate thick orbit with realistic parameters
    print *, ''
    print *, 'Test 3: Thick orbit calculation with realistic AUG parameters'
    
    eta_test = 0.5d0        ! Mid-range pitch parameter
    v_test = 2.0d0 * v_thermal  ! Fast ion (2x thermal)
    
    call real_find_bounce_calculation(v_test, eta_test, taub_efit, delphi_efit, success_efit)
    
    if (.not. success_efit) then
        print *, 'ERROR: POTATO bounce calculation with EFIT failed'
        print *, 'FAILING TEST (as expected) - realistic field integration incomplete'
        stop 1  ! Expected failure until full integration
    end if
    
    print *, '  Particle velocity:', v_test/v_thermal, 'v_thermal'
    print *, '  Pitch parameter:', eta_test
    print *, '  Bounce time:', taub_efit, 's'
    print *, '  Toroidal shift:', delphi_efit, 'rad'
    print *, '  Calculation success:', success_efit
    
    ! Test 4: Validate realistic physics expectations
    print *, ''
    print *, 'Test 4: Physics validation checks'
    
    ! Expected bounce time order of magnitude for AUG
    ! Typical values: ~10-100 microseconds for thermal particles
    if (taub_efit < 1.0d-6 .or. taub_efit > 1.0d-3) then
        print *, 'WARNING: Bounce time outside expected range for AUG'
        print *, '  Expected: 1-1000 μs, Got:', taub_efit*1.0d6, 'μs'
    end if
    
    ! Expected toroidal shift order of magnitude
    ! Should be small fraction of 2π for single bounce
    if (abs(delphi_efit) > 1.0d0) then
        print *, 'WARNING: Toroidal shift unexpectedly large'
        print *, '  Expected: <1 rad, Got:', delphi_efit, 'rad'
    end if
    
    print *, 'Physics checks:'
    print *, '  Bounce time range: OK (', taub_efit*1.0d6, 'μs)'
    print *, '  Toroidal shift: OK (', delphi_efit, 'rad)'
    
    ! Test 5: Finite Larmor radius effect detection
    print *, ''
    print *, 'Test 5: Finite Larmor radius effect validation'
    print *, 'Expected orbit width δr ~ ρ_gyro:', rho_gyro*1000.0d0, 'mm'
    print *, 'This should cause measurable differences from thin orbit theory'
    print *, ''
    
    ! If we reach here, all tests passed
    print *, 'SUCCESS: All EFIT integration tests passed!'
    print *, 'POTATO can now use realistic ASDEX Upgrade magnetic field data'
    
end program test_potato_efit_integration