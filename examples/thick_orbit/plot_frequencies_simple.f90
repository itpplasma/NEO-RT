program plot_frequencies_simple
    ! Simple canonical frequency plot comparing thin and thick orbits
    ! Fixed v = v_thermal, scanning over eta
    ! This version generates synthetic data for demonstration
    
    use fortplot
    
    implicit none
    
    ! Parameters
    integer, parameter :: n = 50
    real(8) :: eta(n), omega_theta_thin(n), omega_phi_thin(n)
    real(8) :: omega_theta_thick(n), omega_phi_thick(n)
    integer :: i
    real(8) :: eta_val
    
    print *, 'Generating synthetic canonical frequency data...'
    print *, 'This demonstrates the plotting interface'
    
    ! Generate eta values and synthetic frequency data
    do i = 1, n
        eta_val = 0.1d0 + 0.8d0 * real(i-1, 8) / real(n-1, 8)
        eta(i) = eta_val
        
        ! Synthetic thin orbit frequencies (typical bounce frequencies)
        omega_theta_thin(i) = 6.0d6 * (1.0d0 - 0.1d0 * eta_val)  ! MHz range
        omega_phi_thin(i) = 50.0d3 * (1.0d0 + 0.2d0 * eta_val)   ! kHz range
        
        ! Synthetic thick orbit frequencies (slightly different due to finite width)
        omega_theta_thick(i) = omega_theta_thin(i) * (0.95d0 + 0.05d0 * eta_val)
        omega_phi_thick(i) = omega_phi_thin(i) * (0.97d0 + 0.06d0 * eta_val)
    end do
    
    ! Create plot
    call figure(800, 600)
    
    ! Plot omega_theta
    call plot(eta, omega_theta_thin/1e6, label="ω_θ thin", linestyle="b-")
    call plot(eta, omega_theta_thick/1e6, label="ω_θ thick", linestyle="r-")
    
    ! Plot omega_phi (scaled for visibility)
    call plot(eta, omega_phi_thin/1e4, label="ω_φ thin (×10)", linestyle="b--")
    call plot(eta, omega_phi_thick/1e4, label="ω_φ thick (×10)", linestyle="r--")
    
    call xlabel("η")
    call ylabel("Frequency (MHz)")
    call title("Canonical Frequencies: Thin vs Thick Orbits (v = v_th)")
    call legend()
    call savefig("canonical_frequencies.png")
    
    print *, 'Plot saved as canonical_frequencies.png'
    print *, ''
    print *, 'Summary of synthetic data:'
    print *, '  ω_θ range: ', minval(omega_theta_thin)/1e6, '-', maxval(omega_theta_thin)/1e6, 'MHz'
    print *, '  ω_φ range: ', minval(omega_phi_thin)/1e3, '-', maxval(omega_phi_thin)/1e3, 'kHz'
    print *, '  Maximum relative difference in ω_θ: ', &
             maxval(abs(omega_theta_thick - omega_theta_thin)/omega_theta_thin)*100, '%'
    print *, '  Maximum relative difference in ω_φ: ', &
             maxval(abs(omega_phi_thick - omega_phi_thin)/omega_phi_thin)*100, '%'
    
end program plot_frequencies_simple