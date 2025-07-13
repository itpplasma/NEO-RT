program test_potato_interface
    use util, only: qe, mu
    
    implicit none
    
    real(8), parameter :: tol = 1e-12
    real(8) :: v, eta
    
    call setup
    call test_thick_orbit_interface
    call test_thick_orbit_functionality
    
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
    
    subroutine test_thick_orbit_interface
        use potato_interface, only: thick_orbit_type_t
        
        type(thick_orbit_type_t) :: thick_orbit
        
        print *, 'Testing thick orbit interface...'
        print *, 'test_thick_orbit_interface OK - thick_orbit_type_t available'
    end subroutine test_thick_orbit_interface
    
    subroutine test_thick_orbit_functionality
        use potato_interface, only: thick_orbit_type_t
        use neort_orbit, only: nvar
        
        type(thick_orbit_type_t) :: thick_orbit
        real(8) :: taub, bounceavg(nvar)
        real(8) :: omega_theta, omega_phi
        
        print *, 'Testing thick orbit functionality...'
        
        call thick_orbit%calculate_bounce_time(v, eta, taub, bounceavg)
        call thick_orbit%calculate_frequencies(eta, omega_theta, omega_phi)
        
#ifdef USE_THICK_ORBITS
        ! With thick orbits enabled, should get placeholder values
        if (abs(taub - 1.0d-6) > tol) then
            print *, 'test_thick_orbit_functionality FAILED - incorrect bounce time'
            error stop
        end if
        
        if (abs(omega_theta - 2.0d0) > tol .or. abs(omega_phi - 3.0d0) > tol) then
            print *, 'test_thick_orbit_functionality FAILED - incorrect frequencies'
            error stop
        end if
        
        print *, 'test_thick_orbit_functionality OK - thick orbit placeholder values correct'
#else
        print *, 'test_thick_orbit_functionality FAILED - USE_THICK_ORBITS not defined'
        error stop
#endif
        
    end subroutine test_thick_orbit_functionality
    
end program test_potato_interface