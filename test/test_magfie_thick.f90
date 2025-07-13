program test_magfie_thick
    use util, only: qe, mu
    
    implicit none
    
    real(8), parameter :: tol = 1e-12
    
    call setup
    call test_magfie_thick_interface
    call test_gyroradius_calculation
    call test_finite_gyroradius_field
    call test_field_continuity_limit
    
    contains
    
    subroutine setup
        use neort, only: init
        use driftorbit, only: do_magfie_init
        
        call setup_control
        call do_magfie_init
        call init
        
    end subroutine setup
    
    subroutine setup_control
        use driftorbit, only: s, M_t, qi, mi, vth, epsmn, m0, &
            mph, mth, magdrift, nopassing, pertfile, &
            nonlin, bfac, efac, inp_swi
        use neort_orbit, only: noshear
        real(8) :: qs, ms
        
        s = 0.153d0
        M_t = 0.1d0
        qs = 1.0d0
        ms = 2.014d0
        vth = 37280978.0d0
        epsmn = 1e-3
        m0 = 0
        mph = 18
        mth = -1
        magdrift = .false.
        nopassing = .false.
        noshear = .true.
        pertfile = .false.
        nonlin = .false.
        bfac = 1.0d0
        efac = 1.0d0
        inp_swi = 8
        
        M_t = M_t*efac/bfac
        qi = qs*qe
        mi = ms*mu
    end subroutine setup_control
    
    subroutine test_magfie_thick_interface
        use magfie_thick_mod, only: calculate_gyroradius, evaluate_field_with_gyroradius
        
        print *, 'Testing thick orbit magnetic field interface...'
        print *, 'test_magfie_thick_interface OK - magfie_thick_mod available'
    end subroutine test_magfie_thick_interface
    
    subroutine test_gyroradius_calculation
        use magfie_thick_mod, only: calculate_gyroradius
        use driftorbit, only: qi, mi, vth
        
        real(8) :: rho_gyro, B_test, rho_expected
        
        print *, 'Testing gyroradius calculation...'
        
        B_test = 1.0d5
        rho_expected = mi * vth / (qi * B_test)
        
        call calculate_gyroradius(vth, B_test, qi, mi, rho_gyro)
        
        if (abs(rho_gyro - rho_expected) > tol) then
            print *, 'test_gyroradius_calculation FAILED - incorrect gyroradius'
            print *, 'Expected:', rho_expected, 'Got:', rho_gyro
            error stop
        end if
        
        print *, 'test_gyroradius_calculation OK - gyroradius calculation correct'
    end subroutine test_gyroradius_calculation
    
    subroutine test_finite_gyroradius_field
        use magfie_thick_mod, only: evaluate_field_with_gyroradius
        
        real(8) :: x(3), rho_gyro, B_field(3), dB_dx(3,3)
        
        print *, 'Testing finite gyroradius magnetic field...'
        
        x = [0.15d0, 0.0d0, 3.14159d0]
        rho_gyro = 1.0d-4
        
        call evaluate_field_with_gyroradius(x, rho_gyro, B_field, dB_dx)
        
        if (sqrt(sum(B_field**2)) <= 0.0d0) then
            print *, 'test_finite_gyroradius_field FAILED - zero magnetic field'
            error stop
        end if
        
        print *, 'test_finite_gyroradius_field OK - finite gyroradius field computed'
    end subroutine test_finite_gyroradius_field
    
    subroutine test_field_continuity_limit
        use magfie_thick_mod, only: evaluate_field_with_gyroradius
        use do_magfie_mod, only: do_magfie
        
        real(8) :: x(3), rho_gyro, B_field_thick(3), dB_dx_thick(3,3)
        real(8) :: bmod, sqrtg, bder(3), hcovar(3), hctrvr(3), hcurl(3)
        real(8) :: B_field_thin(3), relative_error
        
        print *, 'Testing field continuity as gyroradius -> 0...'
        
        x = [0.15d0, 0.0d0, 3.14159d0]
        rho_gyro = 1.0d-12
        
        call evaluate_field_with_gyroradius(x, rho_gyro, B_field_thick, dB_dx_thick)
        call do_magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
        
        B_field_thin(1) = hctrvr(1) * bmod
        B_field_thin(2) = hctrvr(2) * bmod
        B_field_thin(3) = hctrvr(3) * bmod
        
        relative_error = sqrt(sum((B_field_thick - B_field_thin)**2)) / &
                        sqrt(sum(B_field_thin**2))
        
        if (relative_error > 1.0d-6) then
            print *, 'test_field_continuity_limit FAILED - field discontinuity'
            print *, 'Relative error:', relative_error
            error stop
        end if
        
        print *, 'test_field_continuity_limit OK - field continuous as rho_gyro -> 0'
    end subroutine test_field_continuity_limit
    
end program test_magfie_thick