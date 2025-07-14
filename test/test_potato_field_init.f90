program test_potato_field_init
    use iso_fortran_env, only: output_unit, error_unit
    implicit none
    
    logical :: test_passed
    integer :: test_count, failed_count
    
    test_count = 0
    failed_count = 0
    
    write(output_unit, '(A)') "=== Testing POTATO Field Initialization ==="
    write(output_unit, '(A)') ""
    
    ! Initialize field first
    block
        use potato_field_bridge, only: initialize_potato_field
        logical :: init_success
        call initialize_potato_field(init_success)
        if (.not. init_success) then
            write(error_unit, '(A)') "FATAL: Could not initialize POTATO field for testing"
            stop 1
        end if
    end block
    
    call test_field_data_structures()
    call test_field_initialization()
    call test_spline_setup()
    call test_field_evaluation_accuracy()
    
    write(output_unit, '(A)') ""
    write(output_unit, '(A,I0,A,I0,A)') "=== Tests: ", test_count, &
        " | Failed: ", failed_count, " ==="
    
    if (failed_count > 0) then
        stop 1
    end if
    
contains
    
    subroutine test_field_data_structures()
        test_count = test_count + 1
        write(output_unit, '(A)', advance='no') &
            "Testing POTATO field data structures... "
        
        ! Test that POTATO field data arrays are properly allocated
        block
            use field_eq_mod, only: nrad, nzet, rad, zet, psi, splpsi
            
            ! Check if field data structures are allocated
            ! This should fail initially because field data isn't initialized
            test_passed = allocated(rad) .and. allocated(zet) .and. &
                         allocated(psi) .and. allocated(splpsi)
            
            if (test_passed) then
                test_passed = (nrad > 0) .and. (nzet > 0)
            end if
        end block
        
        if (test_passed) then
            write(output_unit, '(A)') "PASSED"
        else
            write(output_unit, '(A)') "FAILED"
            write(error_unit, '(A)') &
                "  ERROR: POTATO field data structures not initialized"
            failed_count = failed_count + 1
        end if
    end subroutine test_field_data_structures
    
    subroutine test_field_initialization()
        test_count = test_count + 1
        write(output_unit, '(A)', advance='no') &
            "Testing POTATO field initialization routine... "
        
        ! Test that we can call POTATO's field initialization
        block
            use potato_field_bridge, only: initialize_potato_field
            logical :: init_success
            
            ! Try to initialize POTATO field from NEO-RT data
            call initialize_potato_field(init_success)
            test_passed = init_success
        end block
        
        if (test_passed) then
            write(output_unit, '(A)') "PASSED"
        else
            write(output_unit, '(A)') "FAILED"
            write(error_unit, '(A)') &
                "  ERROR: Cannot initialize POTATO field"
            failed_count = failed_count + 1
        end if
    end subroutine test_field_initialization
    
    subroutine test_spline_setup()
        test_count = test_count + 1
        write(output_unit, '(A)', advance='no') &
            "Testing POTATO spline interpolation setup... "
        
        ! Test that spline coefficients are properly set up
        block
            use field_eq_mod, only: splpsi, ipoint
            
            ! Check if spline data is allocated and contains data
            test_passed = allocated(splpsi) .and. allocated(ipoint)
            
            if (test_passed) then
                ! Check if spline contains non-zero data
                test_passed = any(splpsi /= 0.0d0)
            end if
        end block
        
        if (test_passed) then
            write(output_unit, '(A)') "PASSED"
        else
            write(output_unit, '(A)') "FAILED"
            write(error_unit, '(A)') &
                "  ERROR: POTATO spline setup failed"
            failed_count = failed_count + 1
        end if
    end subroutine test_spline_setup
    
    subroutine test_field_evaluation_accuracy()
        test_count = test_count + 1
        write(output_unit, '(A)', advance='no') &
            "Testing field evaluation accuracy... "
        
        ! Test that field evaluation gives physically reasonable results
        block
            use potato_field_bridge, only: psif_neo_rt, dpsidr_neo_rt, &
                                          dpsidz_neo_rt
            real(8) :: R, Z, psi1, psi2, dpsidr_result, dpsidz_result
            real(8) :: dR, expected_dpsidr
            
            ! Test at a typical tokamak position
            R = 1.5d0
            Z = 0.0d0
            dR = 0.01d0
            
            ! Get field values
            call psif_neo_rt(R, Z, psi1)
            call psif_neo_rt(R + dR, Z, psi2)
            call dpsidr_neo_rt(R, Z, dpsidr_result)
            call dpsidz_neo_rt(R, Z, dpsidz_result)
            
            ! Check finite difference vs analytical derivative
            expected_dpsidr = (psi2 - psi1) / dR
            test_passed = abs(dpsidr_result - expected_dpsidr) < 0.1d0 * abs(expected_dpsidr)
            
            ! Also check that field values are finite and reasonable
            test_passed = test_passed .and. &
                         (psi1 == psi1) .and. &  ! Not NaN
                         (dpsidr_result == dpsidr_result) .and. &
                         (dpsidz_result == dpsidz_result) .and. &
                         (abs(psi1) < 100.0d0) .and. &  ! Reasonable magnitude
                         (abs(dpsidr_result) < 100.0d0) .and. &
                         (abs(dpsidz_result) < 100.0d0)
        end block
        
        if (test_passed) then
            write(output_unit, '(A)') "PASSED"
        else
            write(output_unit, '(A)') "FAILED"
            write(error_unit, '(A)') &
                "  ERROR: Field evaluation accuracy test failed"
            failed_count = failed_count + 1
        end if
    end subroutine test_field_evaluation_accuracy
    
end program test_potato_field_init