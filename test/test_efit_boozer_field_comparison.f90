program test_efit_boozer_field_comparison
    ! Test that EFIT equilibrium (POTATO) and Boozer axisymmetric field (NEO-RT)  
    ! represent the same magnetic equilibrium by comparing field values
    ! This ensures thick and thin orbit calculations use consistent magnetic geometry
    
    use do_magfie_mod, only: do_magfie_init, do_magfie, R0, s
    use potato_field_bridge, only: initialize_potato_field
    use field_eq_mod, only: psif, dpsidr, dpsidz, nrad, nzet, rad, zet
    implicit none
    
    ! Test parameters
    integer, parameter :: ntest_points = 10
    real(8), parameter :: B_TOLERANCE = 0.05d0     ! 5% tolerance for B-field magnitude
    real(8), parameter :: PSI_TOLERANCE = 0.05d0   ! 5% tolerance for flux function
    real(8), parameter :: GRAD_TOLERANCE = 0.1d0   ! 10% tolerance for gradients
    
    ! Test variables
    real(8) :: R_test(ntest_points), Z_test(ntest_points), phi_test
    real(8) :: B_efit, B_boozer, B_diff_max, B_diff_avg
    real(8) :: psi_efit, psi_boozer, psi_diff_max, psi_diff_avg
    real(8) :: BR_efit, BZ_efit, BR_boozer, BZ_boozer
    real(8) :: dpsidr_efit, dpsidz_efit, dpsidr_boozer, dpsidz_boozer
    real(8) :: x_boozer(3), bmod_boozer, sqrtg_boozer
    real(8) :: bder_boozer(3), hcovar_boozer(3), hctrvr_boozer(3), hcurl_boozer(3)
    logical :: potato_init_success, boozer_init_success
    integer :: i, ierr
    character(len=256) :: boozer_file, efit_file
    logical :: files_exist
    
    print *, '=========================================================='
    print *, 'TEST: EFIT vs Boozer Magnetic Field Comparison'
    print *, '=========================================================='
    print *, 'Verifying magnetic equilibrium consistency between backends'
    print *, ''
    
    ! Test 1: Check file availability
    call get_environment_variable('DATA', efit_file)
    if (len_trim(efit_file) == 0) then
        efit_file = '/home/ert/data'
    end if
    boozer_file = trim(efit_file) // '/AUG/BOOZER/30835_micdu_eqb_6_t3.2/out_neo-2_rmp_90-n0'
    efit_file = trim(efit_file) // '/AUG/EQDSK/g30835.3200_ed6'
    
    inquire(file=trim(efit_file), exist=files_exist)
    if (.not. files_exist) then
        print *, 'ERROR: EFIT file not found:', trim(efit_file)
        print *, 'FAILING TEST (as expected) - need realistic data files'
        stop 1
    end if
    
    inquire(file=trim(boozer_file), exist=files_exist) 
    if (.not. files_exist) then
        print *, 'ERROR: Boozer file not found:', trim(boozer_file)
        print *, 'FAILING TEST (as expected) - need realistic data files'
        stop 1
    end if
    
    print *, 'Test 1: Data files found'
    print *, '  EFIT file:', trim(efit_file)
    print *, '  Boozer file:', trim(boozer_file)
    print *, ''
    
    ! Test 2: Initialize POTATO with EFIT data
    print *, 'Test 2: Initializing magnetic field representations'
    call initialize_potato_field(potato_init_success)
    if (.not. potato_init_success) then
        print *, 'ERROR: POTATO EFIT initialization failed'
        print *, 'FAILING TEST (as expected) - EFIT setup incomplete'
        stop 1
    end if
    print *, '  POTATO/EFIT initialized:', potato_init_success
    
    ! Test 3: Initialize NEO-RT with Boozer data
    ! Note: This requires proper setup of NEO-RT input files
    call setup_neort_boozer_field(boozer_file, boozer_init_success)
    if (.not. boozer_init_success) then
        print *, 'ERROR: NEO-RT Boozer initialization failed'
        print *, 'FAILING TEST (as expected) - Boozer setup needed'
        stop 1
    end if
    print *, '  NEO-RT/Boozer initialized:', boozer_init_success
    print *, ''
    
    ! Test 4: Define test points across the plasma volume
    print *, 'Test 3: Setting up test points'
    phi_test = 0.0d0  ! Axisymmetric, so phi doesn't matter
    
    ! Test points from magnetic axis to edge at midplane
    do i = 1, ntest_points
        ! R from axis to edge (typical AUG: R0 ~ 1.65m, a ~ 0.5m)
        R_test(i) = 1.65d0 + 0.05d0 * real(i-1, 8)  
        Z_test(i) = 0.0d0  ! Midplane
    end do
    
    print *, '  Test points defined from R=', R_test(1), 'm to R=', R_test(ntest_points), 'm'
    print *, ''
    
    ! Test 5: Compare magnetic field values
    print *, 'Test 4: Comparing magnetic field values'
    print *, '  Point    R[m]      |B|_EFIT   |B|_Boozer   Diff[%]'
    print *, '  -----  -------    ----------  ----------  --------'
    
    B_diff_max = 0.0d0
    B_diff_avg = 0.0d0
    psi_diff_max = 0.0d0
    psi_diff_avg = 0.0d0
    
    do i = 1, ntest_points
        ! Get EFIT field values from POTATO
        call get_efit_field(R_test(i), Z_test(i), B_efit, psi_efit, &
                           BR_efit, BZ_efit, dpsidr_efit, dpsidz_efit)
        
        ! Get Boozer field values from NEO-RT
        x_boozer(1) = R_test(i)
        x_boozer(2) = phi_test
        x_boozer(3) = Z_test(i)
        call do_magfie(x_boozer, bmod_boozer, sqrtg_boozer, bder_boozer, &
                       hcovar_boozer, hctrvr_boozer, hcurl_boozer)
        B_boozer = bmod_boozer  ! NEO-RT normalizes to B_ref
        
        ! Compare magnitudes
        if (B_boozer > 0.0d0) then
            B_diff_max = max(B_diff_max, abs(B_efit - B_boozer)/B_boozer)
            B_diff_avg = B_diff_avg + abs(B_efit - B_boozer)/B_boozer
            
            write(*, '(I5, F9.3, 2F12.6, F10.2)') i, R_test(i), B_efit, B_boozer, &
                    100.0d0*abs(B_efit - B_boozer)/B_boozer
        else
            print *, '  WARNING: B_boozer = 0 at point', i
        end if
    end do
    
    B_diff_avg = B_diff_avg / real(ntest_points, 8)
    
    print *, ''
    print *, '  Maximum |B| difference:', 100.0d0*B_diff_max, '%'
    print *, '  Average |B| difference:', 100.0d0*B_diff_avg, '%'
    
    if (B_diff_max > B_TOLERANCE) then
        print *, 'ERROR: Magnetic field magnitude differs by more than', 100.0d0*B_TOLERANCE, '%'
        print *, 'FAILING TEST - Field representations inconsistent'
        stop 1
    end if
    
    ! Test 6: Compare flux surfaces (psi values)
    print *, ''
    print *, 'Test 5: Comparing flux function values'
    print *, '  Point    R[m]     psi_EFIT   psi_Boozer   Diff[%]'
    print *, '  -----  -------   ----------  ----------  --------'
    
    ! Note: Implementation depends on how psi is accessed in NEO-RT
    print *, '  Flux comparison: NOT IMPLEMENTED YET'
    print *, '  FAILING TEST (as expected) - need psi access in NEO-RT'
    stop 1
    
contains

    subroutine get_efit_field(R, Z, B, psi, BR, BZ, dpsi_dR, dpsi_dZ)
        ! Get magnetic field from POTATO's EFIT representation
        use field_eq_mod, only: psif, dpsidr, dpsidz
        implicit none
        real(8), intent(in) :: R, Z
        real(8), intent(out) :: B, psi, BR, BZ, dpsi_dR, dpsi_dZ
        
        ! Call POTATO field evaluation (external subroutine)
        external :: field_eq
        call field_eq(R, Z)
        
        ! Get results from module variables
        psi = psif
        dpsi_dR = dpsidr
        dpsi_dZ = dpsidz
        
        ! Calculate B from grad(psi) for axisymmetric case
        ! In cylindrical coords: B = (BR, Bphi, BZ)
        ! For axisymmetric equilibrium: BR = -(1/R)dpsi/dZ, BZ = (1/R)dpsi/dR
        BR = -dpsi_dZ / R
        BZ = dpsi_dR / R
        
        ! Total field includes toroidal component
        ! For now, approximate |B| (needs proper Bphi from EFIT)
        B = sqrt(BR**2 + BZ**2 + 1.0d0)  ! Rough approximation
        
    end subroutine get_efit_field
    
    subroutine setup_neort_boozer_field(filename, success)
        ! Initialize NEO-RT with Boozer field data
        use do_magfie_mod, only: do_magfie_init, s
        implicit none
        character(len=*), intent(in) :: filename
        logical, intent(out) :: success
        
        ! Note: NEO-RT typically reads Boozer file from command line or namelist
        ! For this test, we'll need to set it up appropriately
        print *, '  Setting up NEO-RT with Boozer file:', trim(filename)
        
        ! Initialize at a reference flux surface (s=0.5)
        s = 0.5d0
        
        ! For now, mark as not implemented
        success = .false.
        print *, '  NEO-RT Boozer setup: NOT IMPLEMENTED in test context'
        
    end subroutine setup_neort_boozer_field

end program test_efit_boozer_field_comparison