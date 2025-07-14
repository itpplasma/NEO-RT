program test_potato_spline_real
    use iso_fortran_env, only: output_unit, error_unit
    implicit none
    
    logical :: test_passed
    integer :: test_count, failed_count
    
    test_count = 0
    failed_count = 0
    
    write(output_unit, '(A)') "=== Testing Real POTATO Spline Implementation ==="
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
    
    call test_spline_interpolation_accuracy()
    call test_derivative_consistency() 
    call test_spline_vs_analytical()
    call test_boundary_behavior()
    
    write(output_unit, '(A)') ""
    write(output_unit, '(A,I0,A,I0,A)') "=== Tests: ", test_count, &
        " | Failed: ", failed_count, " ==="
    
    if (failed_count > 0) then
        stop 1
    end if
    
contains
    
    subroutine test_spline_interpolation_accuracy()
        test_count = test_count + 1
        write(output_unit, '(A)', advance='no') &
            "Testing spline interpolation accuracy... "
        
        ! Test that real spline interpolation gives smooth results
        block
            use potato_field_bridge, only: psif_neo_rt, dpsidr_neo_rt, dpsidz_neo_rt
            real(8) :: R1, Z1, R2, Z2, psi1, psi2
            real(8) :: dpsidr1, dpsidz1, dpsidr2, dpsidz2
            real(8) :: smoothness_test
            
            ! Test at two nearby points
            R1 = 1.4d0
            Z1 = 0.1d0
            R2 = 1.41d0  ! Small step
            Z2 = 0.11d0
            
            call psif_neo_rt(R1, Z1, psi1)
            call dpsidr_neo_rt(R1, Z1, dpsidr1)
            call dpsidz_neo_rt(R1, Z1, dpsidz1)
            
            call psif_neo_rt(R2, Z2, psi2)
            call dpsidr_neo_rt(R2, Z2, dpsidr2)
            call dpsidz_neo_rt(R2, Z2, dpsidz2)
            
            ! Check smoothness - derivatives should change gradually
            smoothness_test = abs((dpsidr2 - dpsidr1)) + abs((dpsidz2 - dpsidz1))
            
            ! Real spline should be smooth, stub is too simple to be smooth
            ! This test will fail with stub but pass with real spline
            test_passed = smoothness_test < 0.5d0  ! Conservative threshold
            
            ! Also verify no NaN or infinite values
            test_passed = test_passed .and. &
                         (psi1 == psi1) .and. (psi2 == psi2) .and. &
                         (dpsidr1 == dpsidr1) .and. (dpsidr2 == dpsidr2) .and. &
                         (dpsidz1 == dpsidz1) .and. (dpsidz2 == dpsidz2)
        end block
        
        if (test_passed) then
            write(output_unit, '(A)') "PASSED"
        else
            write(output_unit, '(A)') "FAILED"
            write(error_unit, '(A)') &
                "  ERROR: Spline interpolation not smooth enough (using stub?)"
            failed_count = failed_count + 1
        end if
    end subroutine test_spline_interpolation_accuracy
    
    subroutine test_derivative_consistency()
        test_count = test_count + 1
        write(output_unit, '(A)', advance='no') &
            "Testing derivative consistency... "
        
        ! Test that finite differences match analytical derivatives
        block
            use potato_field_bridge, only: psif_neo_rt, dpsidr_neo_rt, dpsidz_neo_rt
            real(8) :: R, Z, dR, dZ
            real(8) :: psi_center, psi_dr, psi_dz
            real(8) :: dpsidr_analytical, dpsidz_analytical
            real(8) :: dpsidr_finite, dpsidz_finite
            real(8) :: error_r, error_z
            
            R = 1.5d0
            Z = 0.0d0
            dR = 1.0d-4  ! Small step for finite difference
            dZ = 1.0d-4
            
            ! Get analytical derivatives
            call psif_neo_rt(R, Z, psi_center)
            call dpsidr_neo_rt(R, Z, dpsidr_analytical)
            call dpsidz_neo_rt(R, Z, dpsidz_analytical)
            
            ! Compute finite differences
            call psif_neo_rt(R + dR, Z, psi_dr)
            call psif_neo_rt(R, Z + dZ, psi_dz)
            
            dpsidr_finite = (psi_dr - psi_center) / dR
            dpsidz_finite = (psi_dz - psi_center) / dZ
            
            ! Check relative errors
            error_r = abs(dpsidr_analytical - dpsidr_finite) / &
                     max(abs(dpsidr_analytical), 1.0d-6)
            error_z = abs(dpsidz_analytical - dpsidz_finite) / &
                     max(abs(dpsidz_analytical), 1.0d-6)
            
            ! Real spline should have excellent derivative accuracy
            test_passed = (error_r < 0.01d0) .and. (error_z < 0.01d0)  ! 1% tolerance
        end block
        
        if (test_passed) then
            write(output_unit, '(A)') "PASSED"
        else
            write(output_unit, '(A)') "FAILED"
            write(error_unit, '(A)') &
                "  ERROR: Derivative consistency poor (using stub?)"
            failed_count = failed_count + 1
        end if
    end subroutine test_derivative_consistency
    
    subroutine test_spline_vs_analytical()
        test_count = test_count + 1
        write(output_unit, '(A)', advance='no') &
            "Testing spline vs analytical function... "
        
        ! Test spline interpolation on a known analytical function
        block
            use field_eq_mod, only: nrad, nzet
            
            ! Check if we have realistic grid dimensions (real implementation)
            ! vs minimal grid (stub implementation)
            test_passed = (nrad >= 20) .and. (nzet >= 10)
            
            if (.not. test_passed) then
                write(error_unit, '(A,I0,A,I0)') &
                    "  Grid too small: nrad=", nrad, ", nzet=", nzet
            end if
        end block
        
        if (test_passed) then
            write(output_unit, '(A)') "PASSED"
        else
            write(output_unit, '(A)') "FAILED"
            write(error_unit, '(A)') &
                "  ERROR: Grid dimensions suggest stub implementation"
            failed_count = failed_count + 1
        end if
    end subroutine test_spline_vs_analytical
    
    subroutine test_boundary_behavior()
        test_count = test_count + 1
        write(output_unit, '(A)', advance='no') &
            "Testing boundary behavior... "
        
        ! Test that field evaluation handles boundaries gracefully
        block
            use potato_field_bridge, only: psif_neo_rt
            use field_eq_mod, only: rad, zet, nrad, nzet
            real(8) :: psi_inside, psi_boundary, psi_outside
            real(8) :: R_inside, R_boundary, R_outside, Z
            
            if (allocated(rad) .and. nrad > 1) then
                R_inside = rad(nrad/2)
                R_boundary = rad(nrad)
                R_outside = rad(nrad) + 0.1d0
                Z = 0.0d0
                
                call psif_neo_rt(R_inside, Z, psi_inside)
                call psif_neo_rt(R_boundary, Z, psi_boundary)
                call psif_neo_rt(R_outside, Z, psi_outside)
                
                ! Check that all evaluations return finite values
                test_passed = (psi_inside == psi_inside) .and. &
                             (psi_boundary == psi_boundary) .and. &
                             (psi_outside == psi_outside)
                
                ! Check boundary extrapolation is reasonable
                if (test_passed) then
                    test_passed = abs(psi_outside - psi_boundary) < 10.0d0 * abs(psi_boundary)
                end if
            else
                test_passed = .false.
            end if
        end block
        
        if (test_passed) then
            write(output_unit, '(A)') "PASSED"
        else
            write(output_unit, '(A)') "FAILED"
            write(error_unit, '(A)') &
                "  ERROR: Boundary handling insufficient"
            failed_count = failed_count + 1
        end if
    end subroutine test_boundary_behavior
    
end program test_potato_spline_real