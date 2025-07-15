program test_potato_real_integration
    use iso_fortran_env, only: output_unit, error_unit
    implicit none
    
    logical :: test_passed
    integer :: test_count, failed_count
    
    test_count = 0
    failed_count = 0
    
    write(output_unit, '(A)') "=== Testing Real POTATO Integration ===" 
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
    
    call test_real_find_bounce_call()
    call test_real_orbit_integration()
    call test_real_physics_validation()
    
    write(output_unit, '(A)') ""
    write(output_unit, '(A,I0,A,I0,A)') "=== Tests: ", test_count, &
        " | Failed: ", failed_count, " ==="
    
    if (failed_count > 0) then
        stop 1
    end if
    
contains
    
    subroutine test_real_find_bounce_call()
        test_count = test_count + 1
        write(output_unit, '(A)', advance='no') &
            "Testing real POTATO find_bounce function call... "
        
        ! Test calling actual POTATO find_bounce with proper velo_ext
        block
            use potato_field_bridge, only: real_find_bounce_calculation
            real(8) :: v, eta, taub, delphi
            logical :: success
            
            v = 1.0d5
            eta = 0.5d0
            
            call real_find_bounce_calculation(v, eta, taub, delphi, success)
            
            test_passed = success .and. (taub > 0.0d0) .and. (taub < 1.0d-3) .and. &
                         (abs(delphi) < 1.0d0)
        end block
        
        if (test_passed) then
            write(output_unit, '(A)') "PASSED"
        else
            write(output_unit, '(A)') "FAILED"
            write(error_unit, '(A)') &
                "  ERROR: Real POTATO find_bounce call failed"
            failed_count = failed_count + 1
        end if
    end subroutine test_real_find_bounce_call
    
    subroutine test_real_orbit_integration()
        test_count = test_count + 1
        write(output_unit, '(A)', advance='no') &
            "Testing real POTATO orbit integration... "
        
        ! Test that orbit integration produces physically reasonable results
        block
            use potato_field_bridge, only: real_find_bounce_calculation
            real(8) :: v1, v2, eta_same, taub1, taub2, delphi1, delphi2
            logical :: success1, success2
            
            ! Test velocity scaling - higher velocity should give shorter bounce time
            v1 = 5.0d4
            v2 = 1.0d5
            eta_same = 0.5d0
            
            call real_find_bounce_calculation(v1, eta_same, taub1, delphi1, success1)
            call real_find_bounce_calculation(v2, eta_same, taub2, delphi2, success2)
            
            test_passed = success1 .and. success2 .and. &
                         (taub1 > taub2) .and. &  ! Slower particles take longer
                         (taub1 > 0.0d0) .and. (taub2 > 0.0d0)
        end block
        
        if (test_passed) then
            write(output_unit, '(A)') "PASSED"
        else
            write(output_unit, '(A)') "FAILED"
            write(error_unit, '(A)') &
                "  ERROR: Real orbit integration physics incorrect"
            failed_count = failed_count + 1
        end if
    end subroutine test_real_orbit_integration
    
    subroutine test_real_physics_validation()
        test_count = test_count + 1
        write(output_unit, '(A)', advance='no') &
            "Testing real POTATO physics validation... "
        
        ! Test that real POTATO calculations give different results from stub
        block
            use potato_field_bridge, only: real_find_bounce_calculation
            use potato_stub, only: potato_find_bounce
            real(8) :: v, eta, taub_real, delphi_real, taub_stub, delphi_stub
            real(8) :: extraset_stub(7)
            logical :: success_real, success_stub
            
            v = 1.0d5
            eta = 0.3d0
            
            call real_find_bounce_calculation(v, eta, taub_real, delphi_real, success_real)
            call potato_find_bounce(v, eta, taub_stub, delphi_stub, extraset_stub)
            success_stub = .true.
            
            test_passed = success_real .and. success_stub .and. &
                         (abs(taub_real - taub_stub) / taub_stub > 0.01d0)  ! Should be different
        end block
        
        if (test_passed) then
            write(output_unit, '(A)') "PASSED"
        else
            write(output_unit, '(A)') "FAILED"
            write(error_unit, '(A)') &
                "  ERROR: Real POTATO physics not distinguishable from stub"
            failed_count = failed_count + 1
        end if
    end subroutine test_real_physics_validation
    
end program test_potato_real_integration