program test_efit_boozer_simple_comparison
    ! Simple test to verify EFIT and Boozer files represent same equilibrium
    ! by checking basic properties like magnetic axis location and field strength
    
    use potato_field_bridge, only: initialize_potato_field
    use field_eq_mod, only: nrad, nzet, rad, zet, psi_axis, psi_sep
    implicit none
    
    ! Test parameters  
    real(8), parameter :: AXIS_TOLERANCE = 0.05d0  ! 5cm tolerance for axis location
    
    ! Test variables
    logical :: potato_init_success
    real(8) :: R_axis_efit, Z_axis_efit, R_axis_expected, Z_axis_expected
    real(8) :: psi_range_efit, psi_range_expected
    character(len=256) :: efit_file, data_dir
    logical :: file_exists
    integer :: i
    
    print *, '=========================================================='
    print *, 'TEST: Simple EFIT vs Boozer Equilibrium Verification' 
    print *, '=========================================================='
    print *, 'Basic consistency check between magnetic field representations'
    print *, ''
    
    ! Test 1: Verify EFIT file exists
    call get_environment_variable('DATA', data_dir)
    if (len_trim(data_dir) == 0) then
        data_dir = '/home/ert/data'
    end if
    efit_file = trim(data_dir) // '/AUG/EQDSK/g30835.3200_ed6'
    
    inquire(file=trim(efit_file), exist=file_exists)
    if (.not. file_exists) then
        print *, 'ERROR: EFIT file not found:', trim(efit_file)
        print *, 'FAILING TEST (as expected) - need EFIT data file'
        stop 1
    end if
    
    print *, 'Test 1: EFIT file found'
    print *, '  File:', trim(efit_file)
    print *, ''
    
    ! Test 2: Initialize POTATO with EFIT data
    print *, 'Test 2: Loading EFIT equilibrium into POTATO'
    call initialize_potato_field(potato_init_success)
    if (.not. potato_init_success) then
        print *, 'ERROR: POTATO EFIT initialization failed'
        print *, 'FAILING TEST (as expected) - EFIT loading incomplete'
        stop 1
    end if
    
    ! Display EFIT grid information
    print *, '  EFIT grid loaded successfully'
    print *, '  Grid dimensions: nR =', nrad, ', nZ =', nzet
    print *, '  R range:', rad(1), 'to', rad(nrad), 'm'
    print *, '  Z range:', zet(1), 'to', zet(nzet), 'm'
    print *, '  Psi range:', psi_axis, 'to', psi_sep
    print *, ''
    
    ! Test 3: Find magnetic axis from EFIT data
    print *, 'Test 3: Locating magnetic axis from EFIT'
    call find_magnetic_axis_efit(R_axis_efit, Z_axis_efit)
    print *, '  EFIT magnetic axis: R =', R_axis_efit, 'm, Z =', Z_axis_efit, 'm'
    
    ! Expected values for AUG shot 30835 at t=3.2s
    R_axis_expected = 1.65d0  ! Typical AUG magnetic axis
    Z_axis_expected = 0.0d0   ! Usually near midplane
    
    print *, '  Expected axis (approx): R =', R_axis_expected, 'm, Z =', Z_axis_expected, 'm'
    
    if (abs(R_axis_efit - R_axis_expected) > AXIS_TOLERANCE) then
        print *, 'WARNING: Magnetic axis R position differs by', &
                 abs(R_axis_efit - R_axis_expected)*100.0d0, 'cm'
    end if
    
    if (abs(Z_axis_efit - Z_axis_expected) > AXIS_TOLERANCE) then
        print *, 'WARNING: Magnetic axis Z position differs by', &
                 abs(Z_axis_efit - Z_axis_expected)*100.0d0, 'cm'
    end if
    print *, ''
    
    ! Test 4: Check flux surface topology
    print *, 'Test 4: Flux surface topology check'
    psi_range_efit = psi_sep - psi_axis
    print *, '  EFIT psi range:', psi_range_efit
    print *, '  Normalized psi_axis:', psi_axis
    print *, '  Normalized psi_sep:', psi_sep
    
    ! For properly loaded EFIT data, we expect:
    ! - psi_axis should be minimum (usually normalized to 0)
    ! - psi_sep should be maximum (usually normalized to 1)
    if (psi_axis > psi_sep) then
        print *, 'ERROR: psi_axis > psi_sep - flux surface ordering incorrect'
        print *, 'FAILING TEST - EFIT data not properly loaded'
        stop 1
    end if
    
    if (abs(psi_range_efit - 1.0d0) > 0.1d0) then
        print *, 'WARNING: Flux normalization unexpected (should be ~1.0)'
    end if
    print *, ''
    
    ! Test 5: Basic field strength check
    print *, 'Test 5: Magnetic field strength verification'
    print *, 'For ASDEX Upgrade, expect |B| ~ 2-3 Tesla'
    print *, 'Detailed field comparison requires Boozer data loading'
    print *, ''
    
    ! Since we can't easily load Boozer data in test context without
    ! full NEO-RT initialization, we mark this as expected failure
    print *, 'Full EFIT vs Boozer comparison: NOT IMPLEMENTED'
    print *, 'FAILING TEST (as expected) - need Boozer field access'
    print *, ''
    print *, 'Recommendation: Run NEO-RT with both EFIT and Boozer data'
    print *, 'to verify consistency in production calculations'
    stop 1
    
contains

    subroutine find_magnetic_axis_efit(R_axis, Z_axis)
        ! Find approximate magnetic axis location from EFIT grid
        ! Look for minimum |grad(psi)| location
        use field_eq_mod, only: psi
        implicit none
        real(8), intent(out) :: R_axis, Z_axis
        integer :: i, j, i_axis, j_axis
        real(8) :: psi_min
        
        ! Simple search for psi minimum (magnetic axis)
        psi_min = 1.0d10
        i_axis = nrad/2
        j_axis = nzet/2
        
        do i = 2, nrad-1
            do j = 2, nzet-1
                if (psi(i,j) < psi_min) then
                    psi_min = psi(i,j)
                    i_axis = i
                    j_axis = j
                end if
            end do
        end do
        
        R_axis = rad(i_axis)
        Z_axis = zet(j_axis)
        
    end subroutine find_magnetic_axis_efit

end program test_efit_boozer_simple_comparison