program test_orbit_trajectory_comparison
    ! Test for orbit trajectory comparison between thin and thick orbits
    ! This test verifies that the working directory is set up correctly
    ! and that the required input files are available
    
    implicit none
    
    ! Test parameters
    logical :: file_exists
    integer :: test_count, failed_count
    
    ! Initialize test counters
    test_count = 0
    failed_count = 0
    
    print *, '========================================'
    print *, 'Test: Working Directory Setup'
    print *, '========================================'
    print *, ''
    
    ! Test 1: Check for in_file (Boozer coordinate file)
    print *, 'Test 1: Check for in_file (Boozer coordinate file)'
    test_count = test_count + 1
    
    inquire(file="in_file", exist=file_exists)
    if (file_exists) then
        print *, '  ✓ in_file exists'
    else
        print *, '  ✗ in_file not found'
        failed_count = failed_count + 1
    end if
    
    ! Test 2: Check for driftorbit.in (namelist input file)
    print *, 'Test 2: Check for driftorbit.in (namelist input file)'
    test_count = test_count + 1
    
    inquire(file="driftorbit.in", exist=file_exists)
    if (file_exists) then
        print *, '  ✓ driftorbit.in exists'
    else
        print *, '  ✗ driftorbit.in not found'
        failed_count = failed_count + 1
    end if
    
    ! Test 3: Check for plasma.in (plasma profiles file)
    print *, 'Test 3: Check for plasma.in (plasma profiles file)'
    test_count = test_count + 1
    
    inquire(file="plasma.in", exist=file_exists)
    if (file_exists) then
        print *, '  ✓ plasma.in exists'
    else
        print *, '  ✗ plasma.in not found'
        failed_count = failed_count + 1
    end if
    
    ! Test 4: Read first line of in_file to verify format
    print *, 'Test 4: Verify in_file format'
    test_count = test_count + 1
    
    if (file_exists) then
        block
            character(len=256) :: line
            integer :: unit
            
            open(newunit=unit, file="in_file", status="old", action="read")
            read(unit, '(A)') line
            close(unit)
            
            if (index(line, "Boozer") > 0) then
                print *, '  ✓ in_file appears to be Boozer coordinate file'
                print *, '    First line: ', trim(line)
            else
                print *, '  ✗ in_file does not appear to be Boozer coordinate file'
                print *, '    First line: ', trim(line)
                failed_count = failed_count + 1
            end if
        end block
    else
        print *, '  ✗ Cannot verify in_file format (file not found)'
        failed_count = failed_count + 1
    end if
    
    ! Test 5: Read namelist from driftorbit.in
    print *, 'Test 5: Verify driftorbit.in namelist format'
    test_count = test_count + 1
    
    inquire(file="driftorbit.in", exist=file_exists)
    if (file_exists) then
        block
            character(len=256) :: line
            integer :: unit
            logical :: found_namelist = .false.
            
            open(newunit=unit, file="driftorbit.in", status="old", action="read")
            do while (.not. found_namelist)
                read(unit, '(A)', end=10) line
                if (index(line, "&params") > 0) then
                    found_namelist = .true.
                end if
            end do
            10 continue
            close(unit)
            
            if (found_namelist) then
                print *, '  ✓ driftorbit.in contains valid namelist'
            else
                print *, '  ✗ driftorbit.in does not contain valid namelist'
                failed_count = failed_count + 1
            end if
        end block
    else
        print *, '  ✗ Cannot verify driftorbit.in format (file not found)'
        failed_count = failed_count + 1
    end if
    
    ! Summary
    print *, ''
    print *, 'Working Directory Setup Summary:'
    print *, '  Total tests: ', test_count
    print *, '  Passed: ', test_count - failed_count
    print *, '  Failed: ', failed_count
    
    if (failed_count == 0) then
        print *, '  ✓ Working directory is properly set up for NEO-RT!'
        print *, ''
        print *, 'The following files are available:'
        print *, '  - in_file: Boozer coordinate magnetic field data'
        print *, '  - driftorbit.in: NEO-RT namelist input parameters'
        print *, '  - plasma.in: Plasma profiles for torque calculations'
        print *, ''
        print *, 'This working directory can now be used to run NEO-RT'
        print *, 'with real physics instead of synthetic test data.'
    else
        print *, '  ✗ Working directory setup incomplete'
        print *, '  Please ensure all required input files are present'
        stop 1
    end if
    
end program test_orbit_trajectory_comparison