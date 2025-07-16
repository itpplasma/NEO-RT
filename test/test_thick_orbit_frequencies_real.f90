program test_thick_orbit_frequencies_real
    ! Test thick orbit frequency calculations with real POTATO integration
    use iso_fortran_env, only: real64
    use freq_thick, only: compute_canonical_frequencies_thick
    use driftorbit, only: init_driftorbit, etatp, etadt
    use neort_profiles, only: vth
    use util, only: pi
    implicit none
    
    integer, parameter :: dp = real64
    real(dp) :: v, eta, Om_th_thick, Om_ph_thick
    real(dp) :: Om_th_thin, Om_ph_thin
    logical :: success
    integer :: i
    
    print *, '=============================================='
    print *, 'Test: Thick Orbit Frequency Calculations'
    print *, '=============================================='
    
    ! Initialize physics modules
    call init_basic_physics()
    
    ! Test parameters
    v = vth  ! Use thermal velocity
    
    print *, 'Testing frequency calculations for trapped particles...'
    print *, ''
    print '(A10,A12,A14,A14,A10)', 'eta', 'Om_th', 'Om_ph', 'Om_ph/Om_th', 'Status'
    print *, '--------------------------------------------------------------'
    
    ! Test for different pitch angles (trapped particles)
    do i = 1, 5
        eta = etatp + real(i, dp) * 0.15_dp * (etadt - etatp)
        
        call compute_canonical_frequencies_thick(v, eta, Om_th_thick, Om_ph_thick, success)
        
        if (success) then
            print '(F10.6,ES12.4,ES14.4,F14.6,A10)', &
                  eta, Om_th_thick, Om_ph_thick, Om_ph_thick/Om_th_thick, 'SUCCESS'
        else
            print '(F10.6,A12,A14,A14,A10)', &
                  eta, 'FAILED', 'FAILED', 'N/A', 'FAILED'
        end if
    end do
    
    print *, ''
    print *, 'Physics checks:'
    
    ! Test specific deeply trapped particle
    eta = etatp + 0.5_dp * (etadt - etatp)
    call compute_canonical_frequencies_thick(v, eta, Om_th_thick, Om_ph_thick, success)
    
    if (success) then
        print '(A,ES12.4)', 'Bounce frequency (deeply trapped): ', Om_th_thick
        print '(A,ES12.4)', 'Toroidal frequency: ', Om_ph_thick
        print '(A,F8.4)', 'Frequency ratio Om_ph/Om_th: ', Om_ph_thick/Om_th_thick
        
        ! Physical consistency checks
        if (Om_th_thick > 0.0_dp) then
            print *, '✓ Bounce frequency is positive'
        else
            print *, '✗ ERROR: Bounce frequency should be positive'
        end if
        
        if (abs(Om_ph_thick/Om_th_thick) < 1.0_dp) then
            print *, '✓ Toroidal precession slower than bounce motion'
        else
            print *, '✗ WARNING: Toroidal precession unusually fast'
        end if
    else
        print *, '✗ ERROR: Frequency calculation failed for test particle'
    end if
    
    print *, ''
    print *, 'Test completed.'
    
contains
    
    subroutine init_basic_physics()
        ! Basic physics initialization
        use driftorbit, only: init_done
        use do_magfie_mod, only: s
        
        ! Mark as initialized to bypass full initialization
        init_done = .true.
        s = 0.5_dp
        
        ! Call minimal initialization
        call init_driftorbit()
        
    end subroutine init_basic_physics
    
end program test_thick_orbit_frequencies_real