program test_efit_boozer_consistency
    ! Test coordinate system consistency between EFIT (POTATO) and Boozer (NEO-RT)
    ! This ensures both backends use compatible magnetic field representations
    
    use potato_field_bridge, only: initialize_potato_field
    implicit none
    
    ! Test parameters for coordinate consistency
    real(8), parameter :: R_TEST = 1.65d0     ! Major radius (m) - typical AUG
    real(8), parameter :: Z_TEST = 0.1d0      ! Vertical position (m)
    real(8), parameter :: PSI_TOLERANCE = 1.0d-3  ! Relative tolerance for psi
    
    ! Test variables
    logical :: potato_init_success, boozer_available
    real(8) :: psi_efit, psi_boozer, psi_difference
    character(len=256) :: efit_file, boozer_file
    logical :: efit_exists, boozer_exists
    
    print *, '======================================================'
    print *, 'TEST: EFIT-Boozer Coordinate Consistency (Phase F.1.2)'
    print *, '======================================================'
    print *, 'Verifying magnetic field coordinate compatibility'
    
    ! Test 1: Check file availability
    call get_environment_variable('DATA', efit_file)
    if (len_trim(efit_file) == 0) then
        efit_file = '/home/ert/data'
    end if
    
    efit_file = trim(efit_file) // '/AUG/EQDSK/g30835.3200_ed6'
    boozer_file = trim(efit_file(1:len_trim(efit_file)-19)) // &
                  'BOOZER/30835_micdu_eqb_6_t3.2/out_neo-2_rmp_90-n0'
    
    inquire(file=trim(efit_file), exist=efit_exists)
    inquire(file=trim(boozer_file), exist=boozer_exists)
    
    print *, 'Test 1: File availability check'
    print *, '  EFIT file:', trim(efit_file)
    print *, '  EFIT exists:', efit_exists
    print *, '  Boozer file:', trim(boozer_file)
    print *, '  Boozer exists:', boozer_exists
    print *, ''
    
    if (.not. efit_exists) then
        print *, 'ERROR: EFIT file not found'
        print *, 'FAILING TEST (as expected) - need realistic data files'
        stop 1
    end if
    
    ! Test 2: POTATO field initialization with EFIT
    print *, 'Test 2: POTATO EFIT field initialization'
    call initialize_potato_field(potato_init_success)
    
    if (.not. potato_init_success) then
        print *, 'ERROR: POTATO EFIT initialization failed'
        print *, 'FAILING TEST (as expected) - EFIT integration incomplete'
        stop 1
    end if
    
    print *, '  POTATO EFIT field initialized:', potato_init_success
    
    ! Test 3: Field evaluation at test point
    print *, ''
    print *, 'Test 3: Field evaluation consistency'
    print *, '  Test point: R =', R_TEST, 'm, Z =', Z_TEST, 'm'
    
    ! Get psi from POTATO (EFIT-based)
    call evaluate_potato_psi(R_TEST, Z_TEST, psi_efit)
    print *, '  POTATO (EFIT) psi:', psi_efit
    
    ! Note: Boozer coordinate comparison would require NEO-RT magfie interface
    ! For now, document the expected approach
    psi_boozer = 0.0d0  ! Placeholder - would call NEO-RT magfie
    boozer_available = .false.
    
    if (boozer_available) then
        psi_difference = abs(psi_efit - psi_boozer) / abs(psi_efit)
        print *, '  Boozer psi:', psi_boozer
        print *, '  Relative difference:', psi_difference
        
        if (psi_difference > PSI_TOLERANCE) then
            print *, 'ERROR: Coordinate systems inconsistent'
            print *, 'FAILING TEST - psi difference exceeds tolerance'
            stop 1
        end if
        
        print *, '  SUCCESS: Coordinate systems consistent'
    else
        print *, '  Boozer comparison: NOT IMPLEMENTED YET'
        print *, '  FAILING TEST (as expected) - Boozer interface needed'
        stop 1
    end if
    
    ! Test 4: Flux surface consistency
    print *, ''
    print *, 'Test 4: Flux surface topology consistency'
    print *, 'Expected: Both EFIT and Boozer should give same flux surfaces'
    print *, 'This test validates same magnetic axis and separatrix definitions'
    print *, ''
    print *, 'FAILING TEST (as expected) - full comparison not implemented'
    stop 1
    
contains

    subroutine evaluate_potato_psi(R, Z, psi)
        ! Evaluate psi using POTATO's EFIT-based field
        use field_eq_mod, only: psif
        implicit none
        real(8), intent(in) :: R, Z
        real(8), intent(out) :: psi
        
        ! Call POTATO's field evaluation (sets module variable psif)
        call field_eq(R, Z)
        psi = psif
    end subroutine evaluate_potato_psi

end program test_efit_boozer_consistency