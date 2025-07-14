program plot_canonical_frequencies_demo
    ! Comprehensive plotting demo of canonical frequencies omega_theta and omega_phi vs eta
    ! for both thin and thick orbit calculations at fixed v = v_thermal
    ! This version uses synthetic data for demonstration
    
    use fortplot
    
    implicit none
    
    ! Parameters
    integer, parameter :: neta = 100
    real(8), parameter :: eta_min = 0.05d0, eta_max = 0.95d0
    real(8), parameter :: v_thermal = 1.0d0  ! thermal velocity
    
    ! Arrays for plotting
    real(8) :: eta_array(neta)
    real(8) :: omega_theta_thin(neta), omega_phi_thin(neta)
    real(8) :: omega_theta_thick(neta), omega_phi_thick(neta)
    real(8) :: q_eff_thin(neta), q_eff_thick(neta)
    
    ! Local variables
    real(8) :: eta, deta
    integer :: i
    
    print *, '======================================================'
    print *, 'Canonical Frequency Comparison: Thin vs Thick Orbits'
    print *, '======================================================'
    print *, 'v = v_thermal (synthetic data for demonstration)'
    print *, 'Generating frequency data across eta range...'
    
    deta = (eta_max - eta_min) / real(neta - 1, 8)
    
    ! Generate synthetic frequency data
    do i = 1, neta
        eta = eta_min + real(i - 1, 8) * deta
        eta_array(i) = eta
        
        ! Synthetic thin orbit frequencies (realistic values from literature)
        omega_theta_thin(i) = 6.2d6 * (1.0d0 - 0.15d0 * eta**2)  ! MHz range, decreases with eta
        omega_phi_thin(i) = 50.0d3 * (1.0d0 + 0.3d0 * eta + 0.1d0 * eta**2)  ! kHz range
        
        ! Synthetic thick orbit frequencies (finite Larmor radius effects)
        omega_theta_thick(i) = omega_theta_thin(i) * (0.94d0 + 0.08d0 * eta)
        omega_phi_thick(i) = omega_phi_thin(i) * (0.96d0 + 0.08d0 * eta)
        
        ! Calculate effective safety factors
        q_eff_thin(i) = omega_phi_thin(i) / omega_theta_thin(i)
        q_eff_thick(i) = omega_phi_thick(i) / omega_theta_thick(i)
        
        ! Progress indicator
        if (mod(i, 10) == 0) then
            print *, '  Progress:', i, '/', neta
        endif
    end do
    
    print *, 'Data generation complete. Creating plots...'
    
    ! Plot 1: Poloidal frequency omega_theta vs eta
    call figure(800, 600)
    call plot(eta_array, omega_theta_thin/1e6, label="Thin orbit", linestyle="b-")
    call plot(eta_array, omega_theta_thick/1e6, label="Thick orbit", linestyle="r--")
    call xlabel("Pitch parameter \\eta")
    call ylabel("\\omega_\\theta (MHz)")
    call title("Poloidal Frequency vs Pitch Parameter (v = v_{thermal})")
    call legend()
    call savefig("omega_theta_comparison.png")
    print *, 'Saved: omega_theta_comparison.png'
    
    ! Plot 2: Toroidal frequency omega_phi vs eta
    call figure(800, 600)
    call plot(eta_array, omega_phi_thin/1e3, label="Thin orbit", linestyle="b-")
    call plot(eta_array, omega_phi_thick/1e3, label="Thick orbit", linestyle="r--")
    call xlabel("Pitch parameter \\eta")
    call ylabel("\\omega_\\phi (kHz)")
    call title("Toroidal Frequency vs Pitch Parameter (v = v_{thermal})")
    call legend()
    call savefig("omega_phi_comparison.png")
    print *, 'Saved: omega_phi_comparison.png'
    
    ! Plot 3: Effective safety factor q_eff = omega_phi/omega_theta
    call figure(800, 600)
    call plot(eta_array, q_eff_thin*1000, label="Thin orbit", linestyle="b-")
    call plot(eta_array, q_eff_thick*1000, label="Thick orbit", linestyle="r--")
    call xlabel("Pitch parameter \\eta")
    call ylabel("q_{eff} = \\omega_\\phi/\\omega_\\theta (x10^{-3})")
    call title("Effective Safety Factor vs Pitch Parameter (v = v_{thermal})")
    call legend()
    call savefig("q_eff_comparison.png")
    print *, 'Saved: q_eff_comparison.png'
    
    ! Plot 4: Relative difference in frequencies
    call figure(800, 600)
    call plot(eta_array, 100*(omega_theta_thick - omega_theta_thin)/omega_theta_thin, &
              label="\\Delta\\omega_\\theta", linestyle="g-")
    call plot(eta_array, 100*(omega_phi_thick - omega_phi_thin)/omega_phi_thin, &
              label="\\Delta\\omega_\\phi", linestyle="m-")
    call xlabel("Pitch parameter \\eta")
    call ylabel("Relative difference (%)")
    call title("Orbit Width Effects on Canonical Frequencies (v = v_{thermal})")
    call legend()
    call savefig("frequency_differences.png")
    print *, 'Saved: frequency_differences.png'
    
    ! Plot 5: All frequencies on one plot with dual y-axis style
    call figure(1000, 600)
    
    ! Plot omega_theta on primary scale
    call plot(eta_array, omega_theta_thin/1e6, label="\\omega_\\theta thin", linestyle="b-")
    call plot(eta_array, omega_theta_thick/1e6, label="\\omega_\\theta thick", linestyle="b--")
    
    ! Plot omega_phi on secondary scale (scaled by 100 for visibility)
    call plot(eta_array, omega_phi_thin/1e4, label="\\omega_\\phi thin (x10)", linestyle="r-")
    call plot(eta_array, omega_phi_thick/1e4, label="\\omega_\\phi thick (x10)", linestyle="r--")
    
    call xlabel("Pitch parameter \\eta")
    call ylabel("Frequency (MHz)")
    call title("Canonical Frequencies: Thin vs Thick Orbits (v = v_{thermal})")
    call legend()
    call savefig("all_frequencies.png")
    print *, 'Saved: all_frequencies.png'
    
    ! Summary statistics
    print *, ''
    print *, 'Summary Statistics:'
    print *, '==================='
    print *, 'Max relative difference in omega_theta:', &
             maxval(abs(omega_theta_thick - omega_theta_thin)/omega_theta_thin)*100, '%'
    print *, 'Max relative difference in omega_phi:', &
             maxval(abs(omega_phi_thick - omega_phi_thin)/omega_phi_thin)*100, '%'
    print *, 'Average q_eff (thin):', sum(q_eff_thin)/neta
    print *, 'Average q_eff (thick):', sum(q_eff_thick)/neta
    
    print *, ''
    print *, 'Plotting complete! Check the generated PNG files.'
    
end program plot_canonical_frequencies_demo