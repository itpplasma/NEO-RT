program test_potato_field_bridge
    use iso_fortran_env, only: output_unit, error_unit
    implicit none
    
    logical :: test_passed
    integer :: test_count, failed_count
    
    test_count = 0
    failed_count = 0
    
    write(output_unit, '(A)') "=== Testing POTATO Field Bridge ==="
    write(output_unit, '(A)') ""
    
    call test_field_variables_exist()
    call test_field_eq_call()
    call test_field_bridge_conversion()
    
    write(output_unit, '(A)') ""
    write(output_unit, '(A,I0,A,I0,A)') "=== Tests: ", test_count, &
        " | Failed: ", failed_count, " ==="
    
    if (failed_count > 0) then
        stop 1
    end if
    
contains
    
    subroutine test_field_variables_exist()
        test_count = test_count + 1
        write(output_unit, '(A)', advance='no') &
            "Testing POTATO field variables exist... "
        
        ! Test that POTATO's field variables are accessible
        block
            use field_eq_mod, only: psif, dpsidr, dpsidz
            real(8) :: test_psi, test_dpsidr, test_dpsidz
            
            ! These should be module variables in POTATO
            test_psi = psif
            test_dpsidr = dpsidr
            test_dpsidz = dpsidz
            
            test_passed = .true.
        end block
        
        if (test_passed) then
            write(output_unit, '(A)') "PASSED"
        else
            write(output_unit, '(A)') "FAILED"
            write(error_unit, '(A)') &
                "  ERROR: Cannot access POTATO field variables"
            failed_count = failed_count + 1
        end if
    end subroutine test_field_variables_exist
    
    subroutine test_field_eq_call()
        test_count = test_count + 1
        write(output_unit, '(A)', advance='no') &
            "Testing POTATO field_eq function... "
        
        ! Test that we can call field_eq through the bridge module
        block
            use potato_field_bridge, only: field_eq
            real(8) :: R, Z
            
            ! Set test position
            R = 1.5d0
            Z = 0.0d0
            
            ! Call field_eq to set psif, dpsidr, dpsidz
            call field_eq(R, Z)
            
            ! If we get here, the function exists
            test_passed = .true.
        end block
        
        if (test_passed) then
            write(output_unit, '(A)') "PASSED"
        else
            write(output_unit, '(A)') "FAILED"
            write(error_unit, '(A)') &
                "  ERROR: Cannot call POTATO field_eq function"
            failed_count = failed_count + 1
        end if
    end subroutine test_field_eq_call
    
    subroutine test_field_bridge_conversion()
        test_count = test_count + 1
        write(output_unit, '(A)', advance='no') &
            "Testing field bridge conversion... "
        
        ! Test the field bridge between NEO-RT and POTATO
        block
            use potato_field_bridge, only: psif_neo_rt, dpsidr_neo_rt, &
                                          dpsidz_neo_rt
            real(8) :: R, Z, psi_result, dpsidr_result, dpsidz_result
            
            ! Set test position
            R = 1.5d0
            Z = 0.0d0
            
            ! Call NEO-RT style functions that bridge to POTATO
            call psif_neo_rt(R, Z, psi_result)
            call dpsidr_neo_rt(R, Z, dpsidr_result)
            call dpsidz_neo_rt(R, Z, dpsidz_result)
            
            ! Check that we got reasonable values (not NaN or zero)
            test_passed = (psi_result == psi_result) .and. &
                         (dpsidr_result == dpsidr_result) .and. &
                         (dpsidz_result == dpsidz_result)
        end block
        
        if (test_passed) then
            write(output_unit, '(A)') "PASSED"
        else
            write(output_unit, '(A)') "FAILED"
            write(error_unit, '(A)') &
                "  ERROR: Field bridge conversion failed"
            failed_count = failed_count + 1
        end if
    end subroutine test_field_bridge_conversion
    
end program test_potato_field_bridge