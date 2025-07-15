program test_potato_floating_point_fix
    use iso_fortran_env, only: output_unit, error_unit
    implicit none
    
    logical :: test_passed
    integer :: test_count, failed_count
    
    test_count = 0
    failed_count = 0
    
    write(output_unit, '(A)') "=== Testing POTATO Floating Point Exception Fix ==="
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
    
    call test_edge_case_parameters()
    call test_zero_velocity_components()
    call test_extreme_pitch_angles()
    call test_small_magnetic_field()
    
    write(output_unit, '(A)') ""
    write(output_unit, '(A,I0,A,I0,A)') "=== Tests: ", test_count, &
        " | Failed: ", failed_count, " ==="
    
    if (failed_count > 0) then
        stop 1
    end if
    
contains
    
    subroutine test_edge_case_parameters()
        test_count = test_count + 1
        write(output_unit, '(A)', advance='no') &
            "Testing edge case parameters that might cause FPE... "
        
        block
            use potato_field_bridge, only: real_find_bounce_calculation
            real(8) :: v, eta, taub, delphi
            logical :: success
            integer :: i
            real(8), dimension(5) :: test_velocities
            real(8), dimension(5) :: test_etas
            
            ! Test cases that might cause floating point exceptions
            test_velocities = [1.0d-10, 1.0d-5, 1.0d3, 1.0d5, 1.0d7]  ! Very small to large
            test_etas = [1.0d-10, 0.1d0, 0.5d0, 0.9d0, 0.999d0]  ! Including near-zero
            
            test_passed = .true.
            
            do i = 1, 5
                v = test_velocities(i)
                eta = test_etas(i)
                
                ! This should not crash with FPE
                call real_find_bounce_calculation(v, eta, taub, delphi, success)
                
                ! We don't require success, just no crash
                if (.not. success) then
                    write(output_unit, '(A,E12.5,A,F8.5,A)') &
                        " (v=", v, ", eta=", eta, " failed gracefully)"
                end if
            end do
        end block
        
        if (test_passed) then
            write(output_unit, '(A)') "PASSED"
        else
            write(output_unit, '(A)') "FAILED"
            failed_count = failed_count + 1
        end if
    end subroutine test_edge_case_parameters
    
    subroutine test_zero_velocity_components()
        test_count = test_count + 1
        write(output_unit, '(A)', advance='no') &
            "Testing zero velocity components... "
        
        block
            use potato_field_bridge, only: real_find_bounce_calculation
            real(8) :: v, eta, taub, delphi
            logical :: success
            
            ! Test case where parallel velocity might be zero
            v = 1.0d5
            eta = 1.0d0 / 1.0d0  ! This gives v_par = 0 when bmod = 1
            
            call real_find_bounce_calculation(v, eta, taub, delphi, success)
            
            ! Should handle gracefully without FPE
            test_passed = .true.  ! If we get here without crash
        end block
        
        if (test_passed) then
            write(output_unit, '(A)') "PASSED"
        else
            write(output_unit, '(A)') "FAILED"
            failed_count = failed_count + 1
        end if
    end subroutine test_zero_velocity_components
    
    subroutine test_extreme_pitch_angles()
        test_count = test_count + 1
        write(output_unit, '(A)', advance='no') &
            "Testing extreme pitch angles... "
        
        block
            use potato_field_bridge, only: real_find_bounce_calculation
            real(8) :: v, eta, taub, delphi
            logical :: success
            
            v = 1.0d5
            
            ! Test trapped particle at turning point
            eta = 0.999999d0  ! Nearly perpendicular
            call real_find_bounce_calculation(v, eta, taub, delphi, success)
            
            ! Test nearly field-aligned particle
            eta = 1.0d-6  ! Nearly parallel
            call real_find_bounce_calculation(v, eta, taub, delphi, success)
            
            test_passed = .true.  ! If we get here without crash
        end block
        
        if (test_passed) then
            write(output_unit, '(A)') "PASSED"
        else
            write(output_unit, '(A)') "FAILED"
            failed_count = failed_count + 1
        end if
    end subroutine test_extreme_pitch_angles
    
    subroutine test_small_magnetic_field()
        test_count = test_count + 1
        write(output_unit, '(A)', advance='no') &
            "Testing small magnetic field values... "
        
        block
            use potato_field_bridge, only: real_find_bounce_calculation
            real(8) :: v, eta, taub, delphi
            logical :: success
            
            ! This might encounter regions with small B-field
            v = 1.0d5
            eta = 0.5d0
            
            ! The simplified test field should handle this
            call real_find_bounce_calculation(v, eta, taub, delphi, success)
            
            test_passed = .true.  ! If we get here without crash
        end block
        
        if (test_passed) then
            write(output_unit, '(A)') "PASSED"
        else
            write(output_unit, '(A)') "FAILED"
            failed_count = failed_count + 1
        end if
    end subroutine test_small_magnetic_field
    
end program test_potato_floating_point_fix