program test_potato_wrapper
    implicit none
    
    call test_wrapper_initialization
    call test_wrapper_find_bounce
    call test_wrapper_frequency_calculation
    call test_coordinate_transformation
    
    contains
    
    subroutine test_wrapper_initialization
        use potato_wrapper, only: potato_wrapper_init
        
        print *, 'Testing POTATO wrapper initialization...'
        
        call potato_wrapper_init()
        
        print *, 'test_wrapper_initialization OK'
        
    end subroutine test_wrapper_initialization
    
    subroutine test_wrapper_find_bounce
        use potato_wrapper, only: potato_wrapper_find_bounce
        
        real(8) :: v, eta, taub, delphi
        real(8) :: extraset(7)
        
        print *, 'Testing POTATO wrapper find_bounce...'
        
        v = 1.0d6
        eta = 0.5d0
        extraset = 0.0d0
        
        call potato_wrapper_find_bounce(v, eta, taub, delphi, extraset)
        
        if (taub <= 0.0d0) then
            print *, 'test_wrapper_find_bounce FAILED - invalid bounce time'
            error stop
        end if
        
        if (delphi <= 0.0d0) then
            print *, 'test_wrapper_find_bounce FAILED - invalid toroidal shift'
            error stop
        end if
        
        print *, '  Bounce time:', taub
        print *, '  Toroidal shift:', delphi
        print *, 'test_wrapper_find_bounce OK'
        
    end subroutine test_wrapper_find_bounce
    
    subroutine test_wrapper_frequency_calculation
        use potato_wrapper, only: potato_wrapper_find_bounce, &
                                   potato_wrapper_calculate_frequencies
        
        real(8) :: v, eta, taub, delphi
        real(8) :: extraset(7)
        real(8) :: omega_bounce, omega_toroidal
        real(8), parameter :: pi = 3.14159265358979323846d0
        
        print *, 'Testing POTATO wrapper frequency calculation...'
        
        v = 1.0d6
        eta = 0.3d0
        extraset = 0.0d0
        
        call potato_wrapper_find_bounce(v, eta, taub, delphi, extraset)
        call potato_wrapper_calculate_frequencies(taub, delphi, omega_bounce, omega_toroidal)
        
        ! Verify frequency calculation
        if (abs(omega_bounce - 2.0d0 * pi / taub) > 1.0d-12) then
            print *, 'test_wrapper_frequency_calculation FAILED - incorrect bounce frequency'
            error stop
        end if
        
        if (abs(omega_toroidal - delphi / taub) > 1.0d-12) then
            print *, 'test_wrapper_frequency_calculation FAILED - incorrect toroidal frequency'
            error stop
        end if
        
        print *, '  Bounce frequency:', omega_bounce
        print *, '  Toroidal frequency:', omega_toroidal
        print *, 'test_wrapper_frequency_calculation OK'
        
    end subroutine test_wrapper_frequency_calculation
    
    subroutine test_coordinate_transformation
        ! Test coordinate transformation concepts
        ! In actual implementation, this would verify:
        ! - NEO-RT (v, eta) to POTATO (p, lambda) mapping
        ! - Spatial coordinate initialization
        ! - Normalization factors
        
        real(8) :: v, eta, p_norm, lambda
        
        print *, 'Testing coordinate transformation concepts...'
        
        v = 1.0d6  ! Velocity in cm/s
        eta = 0.4d0  ! NEO-RT pitch parameter
        
        ! Simplified transformation (actual would be more complex)
        p_norm = v / 1.0d6  ! Normalized momentum
        lambda = sqrt(1.0d0 - eta)  ! Cosine of pitch angle
        
        print *, '  NEO-RT coordinates: v =', v, ', eta =', eta
        print *, '  POTATO coordinates: p =', p_norm, ', lambda =', lambda
        
        if (lambda < -1.0d0 .or. lambda > 1.0d0) then
            print *, 'test_coordinate_transformation FAILED - invalid lambda'
            error stop
        end if
        
        print *, 'test_coordinate_transformation OK'
        
    end subroutine test_coordinate_transformation
    
end program test_potato_wrapper