program test_potato_frequencies
    use util, only: qe, mu
    
    implicit none
    
    real(8), parameter :: tol = 1e-12
    real(8) :: v, eta
    
    call setup
    call test_potato_frequency_interface
    call test_potato_frequency_vs_thin
    call test_potato_frequency_direct_call
    
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
    
    subroutine test_potato_frequency_interface
        use potato_wrapper, only: potato_wrapper_calculate_frequencies
        
        print *, 'Testing POTATO frequency interface...'
        print *, 'test_potato_frequency_interface OK - POTATO stub accessible'
        
    end subroutine test_potato_frequency_interface
    
    subroutine test_potato_frequency_vs_thin
        use potato_interface, only: thick_orbit_type_t, thin_orbit_type_t
        
        type(thick_orbit_type_t) :: thick_orbit
        type(thin_orbit_type_t) :: thin_orbit
        real(8) :: omega_theta_thick, omega_phi_thick
        real(8) :: omega_theta_thin, omega_phi_thin
        
        print *, 'Testing POTATO frequency vs thin orbit comparison...'
        
        call thick_orbit%calculate_frequencies(eta, omega_theta_thick, omega_phi_thick)
        call thin_orbit%calculate_frequencies(eta, omega_theta_thin, omega_phi_thin)
        
        ! Thick orbit should give different (eta-dependent) results
        if (abs(omega_theta_thick - omega_theta_thin) < 1.0d-12 .and. &
            abs(omega_phi_thick - omega_phi_thin) < 1.0d-12) then
            print *, 'test_potato_frequency_vs_thin FAILED - thick and thin results identical'
            error stop
        end if
        
        print *, 'test_potato_frequency_vs_thin OK - thick orbit gives different frequencies'
        
    end subroutine test_potato_frequency_vs_thin
    
    subroutine test_potato_frequency_direct_call
        use potato_wrapper, only: potato_wrapper_calculate_frequencies
        
        real(8) :: taub, delphi, omega_bounce, omega_toroidal
        
        print *, 'Testing direct POTATO frequency call...'
        
        taub = 1.0d-6
        delphi = 0.05d0
        call potato_wrapper_calculate_frequencies(taub, delphi, omega_bounce, omega_toroidal)
        
        if (omega_bounce <= 0.0d0) then
            print *, 'test_potato_frequency_direct_call FAILED - invalid bounce frequency'
            error stop
        end if
        
        if (omega_toroidal <= 0.0d0) then
            print *, 'test_potato_frequency_direct_call FAILED - invalid toroidal frequency'
            error stop
        end if
        
        print *, 'test_potato_frequency_direct_call OK - direct POTATO frequency call successful'
        
    end subroutine test_potato_frequency_direct_call
    
end program test_potato_frequencies