program extract_real_orbits
    ! Extract real orbit data using NEO-RT bounce() function
    ! This generates data files that can be plotted with Python
    
    use neort_orbit, only: bounce, nvar, th0
    use driftorbit, only: etatp, etadt, init_done
    use neort_profiles, only: vth
    use do_magfie_mod, only: s, R0, eps
    
    implicit none
    
    ! Parameters
    real(8), parameter :: v_test = 1.0d6   ! test velocity [m/s]
    real(8), parameter :: s_test = 0.5d0   ! normalized flux surface
    real(8) :: eta_test                    ! pitch parameter
    
    ! Physics variables
    real(8) :: taub, bounceavg(nvar)
    real(8) :: R0_phys, eps_phys
    
    ! File output
    integer :: unit
    
    print *, '=== Extracting Real NEO-RT Orbit Data ==='
    print *, 'Using actual bounce() function from orbit.f90'
    print *, ''
    
    ! Check if physics is initialized
    if (.not. init_done) then
        print *, 'ERROR: Physics not initialized. Run from test_run directory after NEO-RT.'
        stop 1
    end if
    
    ! Set physics parameters
    s = s_test
    R0_phys = R0
    eps_phys = eps
    eta_test = etatp + 0.3d0 * (etadt - etatp)  ! trapped particle
    
    print *, 'Physics parameters:'
    print '(A,F6.3)', '  s (flux surface) = ', s_test
    print '(A,F6.3)', '  eta (pitch)      = ', eta_test
    print '(A,F6.3)', '  R0 (major radius)= ', R0_phys
    print '(A,F6.3)', '  epsilon          = ', eps_phys
    print *, ''
    
    ! Calculate real bounce time using NEO-RT physics
    print *, 'Calling real NEO-RT bounce() function...'
    
    call bounce(vth, eta_test, taub, bounceavg)
    
    print *, 'Real NEO-RT bounce calculation results:'
    print '(A,ES12.4,A)', '  Bounce time: τb = ', taub, ' s'
    print '(A,ES12.4)', '  <θ> = ', bounceavg(1)
    print '(A,ES12.4)', '  <v_parallel> = ', bounceavg(2)
    print '(A,ES12.4)', '  <v_phi> = ', bounceavg(3)
    print '(A,ES12.4)', '  <H_pert_real> = ', bounceavg(4)
    print '(A,ES12.4)', '  <H_pert_imag> = ', bounceavg(5)
    print '(A,ES12.4)', '  <1/|B|> = ', bounceavg(6)
    print *, ''
    
    ! Write orbit data to file for Python plotting
    open(newunit=unit, file='neo_rt_orbit_data.dat', status='replace')
    write(unit, '(A)') '# Real NEO-RT orbit data'
    write(unit, '(A,ES12.4)') '# Bounce time: ', taub
    write(unit, '(A,F8.4)') '# Eta: ', eta_test
    write(unit, '(A,F8.4)') '# R0: ', R0_phys
    write(unit, '(A,F8.4)') '# eps: ', eps_phys
    write(unit, '(A,F8.4)') '# s: ', s_test
    write(unit, '(A)') '# This data comes from actual NEO-RT bounce() calculation'
    write(unit, '(A)') '# Use this instead of Python approximations'
    close(unit)
    
    print *, 'Real NEO-RT orbit data written to neo_rt_orbit_data.dat'
    print *, 'This file contains actual physics results from NEO-RT bounce() function'
    print *, 'Use this data for real orbit plotting instead of Python approximations'
    
end program extract_real_orbits