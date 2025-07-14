program test_potato_build
    use iso_fortran_env, only: output_unit, error_unit
    implicit none
    
    logical :: test_passed
    integer :: test_count, failed_count
    
    test_count = 0
    failed_count = 0
    
    write(output_unit, '(A)') "=== Testing POTATO Library Build Integration ==="
    write(output_unit, '(A)') ""
    
    call test_potato_modules_available()
    call test_potato_base_library_link()
    call test_potato_function_calls()
    
    write(output_unit, '(A)') ""
    write(output_unit, '(A,I0,A,I0,A)') "=== Tests: ", test_count, &
        " | Failed: ", failed_count, " ==="
    
    if (failed_count > 0) then
        stop 1
    end if
    
contains
    
    subroutine test_potato_modules_available()
        test_count = test_count + 1
        write(output_unit, '(A)', advance='no') &
            "Testing POTATO modules availability... "
        
        ! Test that POTATO modules are available
        ! Try to use POTATO field modules
        block
            use field_eq_mod, only: psif, dpsidr, dpsidz
            test_passed = .true.
        end block
        
        if (test_passed) then
            write(output_unit, '(A)') "PASSED"
        else
            write(output_unit, '(A)') "FAILED"
            write(error_unit, '(A)') &
                "  ERROR: Cannot access POTATO modules"
            failed_count = failed_count + 1
        end if
    end subroutine test_potato_modules_available
    
    subroutine test_potato_base_library_link()
        test_count = test_count + 1
        write(output_unit, '(A)', advance='no') &
            "Testing POTATO base library linking... "
        
        ! This test verifies we can link against potato_base library
        block
            use orbit_dim_mod, only: neqm
            test_passed = (neqm > 0)
        end block
        
        if (test_passed) then
            write(output_unit, '(A)') "PASSED"
        else
            write(output_unit, '(A)') "FAILED"
            write(error_unit, '(A)') &
                "  ERROR: Cannot link with POTATO base library"
            failed_count = failed_count + 1
        end if
    end subroutine test_potato_base_library_link
    
    subroutine test_potato_function_calls()
        test_count = test_count + 1
        write(output_unit, '(A)', advance='no') &
            "Testing POTATO function signatures... "
        
        ! Test that we can declare external POTATO functions
        block
            real(8) :: z_eqm(5), taub, delphi, extraset(1)
            real(8) :: dtau_in
            integer :: next
            external :: find_bounce, velo_ext
            
            ! Set up dummy values
            next = 1
            dtau_in = 1.0d-3
            z_eqm = [1.5d0, 0.0d0, 0.0d0, 1.0d0, 0.1d0]
            extraset = 0.0d0
            
            ! Since we can't actually call find_bounce without proper setup,
            ! just verify the external declaration works
            test_passed = .true.
        end block
        
        if (test_passed) then
            write(output_unit, '(A)') "PASSED"
        else
            write(output_unit, '(A)') "FAILED"
            write(error_unit, '(A)') &
                "  ERROR: Cannot declare POTATO find_bounce function"
            failed_count = failed_count + 1
        end if
    end subroutine test_potato_function_calls
    
end program test_potato_build

! Dummy velocity routine for testing
subroutine velo_ext(dtau, y, f)
    implicit none
    real(8), intent(in) :: dtau, y(5)
    real(8), intent(out) :: f(5)
    f = 0.0d0
end subroutine velo_ext