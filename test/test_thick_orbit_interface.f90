program test_thick_orbit_interface
    use util, only: qe, mu
    
    implicit none
    
    real(8), parameter :: tol = 1e-12
    real(8) :: v, eta
    
    call setup
    call test_orbit_type_interface
    call test_thin_orbit_wrapper
    call test_integration_thin_wrapper
    
    contains
    
    subroutine setup
        use neort, only: init
        use driftorbit, only: do_magfie_init, etamin, etamax, &
            Om_tE, dOm_tEds, etatp, etadt, epst, sign_vpar, vth, M_t, dM_tds
        use do_magfie_mod, only: R0
        
        call setup_control
        call do_magfie_init
        call init
        
        Om_tE = vth*M_t/R0
        dOm_tEds = vth*dM_tds/R0
        
        etamin = (1+epst)*etatp
        etamax = (1-epst)*etadt
        sign_vpar = 1
        
        v = vth
        eta = 0.5d0*(etamin + etamax)
        
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
    
    subroutine test_orbit_type_interface
        use orbit_types, only: orbit_type_t
        
        print *, 'Testing orbit type interface...'
        print *, 'test_orbit_type_interface OK - orbit_types module available'
    end subroutine test_orbit_type_interface
    
    subroutine test_thin_orbit_wrapper
        use orbit_types, only: thin_orbit_type_t
        use neort_orbit, only: nvar
        
        type(thin_orbit_type_t) :: thin_orbit
        real(8) :: taub, bounceavg(nvar)
        real(8) :: omega_theta, omega_phi
        
        print *, 'Testing thin orbit wrapper...'
        
        call thin_orbit%calculate_bounce_time(v, eta, taub, bounceavg)
        call thin_orbit%calculate_frequencies(eta, omega_theta, omega_phi)
        
        if (taub <= 0.0d0) then
            print *, 'test_thin_orbit_wrapper FAILED - invalid bounce time'
            error stop
        end if
        
        print *, 'test_thin_orbit_wrapper OK - thin orbit wrapper works'
    end subroutine test_thin_orbit_wrapper
    
    subroutine test_integration_thin_wrapper
        use orbit_types, only: thin_orbit_type_t
        use neort_orbit, only: nvar, bounce
        
        type(thin_orbit_type_t) :: thin_orbit
        real(8) :: taub_wrapper, bounceavg_wrapper(nvar)
        real(8) :: taub_direct, bounceavg_direct(nvar)
        real(8) :: relative_error
        
        print *, 'Testing integration: thin wrapper vs direct calls...'
        
        call thin_orbit%calculate_bounce_time(v, eta, taub_wrapper, bounceavg_wrapper)
        call bounce(v, eta, taub_direct, bounceavg_direct)
        
        relative_error = abs(taub_wrapper - taub_direct) / taub_direct
        if (relative_error > tol) then
            print *, 'test_integration_thin_wrapper FAILED - bounce time mismatch'
            print *, 'Wrapper:', taub_wrapper, 'Direct:', taub_direct
            print *, 'Relative error:', relative_error
            error stop
        end if
        
        relative_error = maxval(abs(bounceavg_wrapper - bounceavg_direct) / &
                               (abs(bounceavg_direct) + 1e-16))
        if (relative_error > tol) then
            print *, 'test_integration_thin_wrapper FAILED - bounce averages mismatch'
            print *, 'Max relative error:', relative_error
            error stop
        end if
        
        print *, 'test_integration_thin_wrapper OK - wrapper matches direct calls'
    end subroutine test_integration_thin_wrapper
    
end program test_thick_orbit_interface