program plot_runtime_comparison
    ! Plot comparison between thin and thick orbit calculations
    ! using the new runtime configuration system
    
    use fortplot
    use runtime_config, only: set_use_thick_orbits, get_use_thick_orbits
    use orbit_interface, only: orbit_calculator_t, orbit_calculator_factory
    use field_interface, only: field_evaluator_t, field_evaluator_factory
    
    implicit none
    
    ! Parameters
    integer, parameter :: neta = 50
    real(8), parameter :: eta_min = 0.1d0, eta_max = 0.9d0
    real(8), parameter :: v_test = 1.0d6  ! m/s
    
    ! Arrays for plotting
    real(8) :: eta_array(neta)
    real(8) :: omega_theta_thin(neta), omega_phi_thin(neta)
    real(8) :: omega_theta_thick(neta), omega_phi_thick(neta)
    real(8) :: taub_thin(neta), taub_thick(neta)
    real(8) :: delphi_thin(neta), delphi_thick(neta)
    
    ! Local variables
    real(8) :: eta, deta, taub, delphi, Om_th, Om_ph, extraset(7)
    real(8) :: s_flux, theta_boozer, phi_boozer
    logical :: success
    integer :: i
    class(orbit_calculator_t), allocatable :: calculator
    
    print *, '======================================================='
    print *, 'Runtime Configuration Comparison: Thin vs Thick Orbits'
    print *, '======================================================='
    print *, 'v =', v_test, 'm/s'
    print *, 'Using runtime dispatch architecture...'
    
    deta = (eta_max - eta_min) / real(neta - 1, 8)
    
    ! Set up orbit parameters
    s_flux = 0.5d0
    theta_boozer = 0.0d0
    phi_boozer = 0.0d0
    
    ! Calculate thin orbit results
    print *, 'Calculating thin orbit results...'
    call set_use_thick_orbits(.false.)
    
    do i = 1, neta
        eta = eta_min + real(i - 1, 8) * deta
        eta_array(i) = eta
        
        ! Get thin orbit calculator
        calculator = orbit_calculator_factory(get_use_thick_orbits())
        
        ! Calculate bounce time
        call calculator%find_bounce(v_test, eta, s_flux, theta_boozer, phi_boozer, &
                                   taub, delphi, extraset, success)
        
        if (success) then
            taub_thin(i) = taub
            delphi_thin(i) = delphi
            
            ! Calculate frequencies
            call calculator%calculate_frequencies(v_test, eta, Om_th, Om_ph, success)
            if (success) then
                omega_theta_thin(i) = Om_th
                omega_phi_thin(i) = Om_ph
            else
                omega_theta_thin(i) = 0.0d0
                omega_phi_thin(i) = 0.0d0
            end if
        else
            taub_thin(i) = 0.0d0
            delphi_thin(i) = 0.0d0
            omega_theta_thin(i) = 0.0d0
            omega_phi_thin(i) = 0.0d0
        end if
        
        call calculator%cleanup()
        
        if (mod(i, 10) == 0) then
            print *, '  Progress:', i, '/', neta, '(thin orbit)'
        end if
    end do
    
    ! Calculate thick orbit results
    print *, 'Calculating thick orbit results...'
    call set_use_thick_orbits(.true.)
    
    do i = 1, neta
        eta = eta_array(i)
        
        ! Get thick orbit calculator
        calculator = orbit_calculator_factory(get_use_thick_orbits())
        
        ! Calculate bounce time
        call calculator%find_bounce(v_test, eta, s_flux, theta_boozer, phi_boozer, &
                                   taub, delphi, extraset, success)
        
        if (success) then
            taub_thick(i) = taub
            delphi_thick(i) = delphi
            
            ! Calculate frequencies
            call calculator%calculate_frequencies(v_test, eta, Om_th, Om_ph, success)
            if (success) then
                omega_theta_thick(i) = Om_th
                omega_phi_thick(i) = Om_ph
            else
                omega_theta_thick(i) = 0.0d0
                omega_phi_thick(i) = 0.0d0
            end if
        else
            taub_thick(i) = 0.0d0
            delphi_thick(i) = 0.0d0
            omega_theta_thick(i) = 0.0d0
            omega_phi_thick(i) = 0.0d0
        end if
        
        call calculator%cleanup()
        
        if (mod(i, 10) == 0) then
            print *, '  Progress:', i, '/', neta, '(thick orbit)'
        end if
    end do
    
    print *, 'Calculations complete. Creating plots...'
    
    ! Plot 1: Bounce time comparison
    call figure(800, 600)
    call plot(eta_array, taub_thin*1e6, label="Thin orbit", linestyle="b-")
    call plot(eta_array, taub_thick*1e6, label="Thick orbit", linestyle="r--")
    call xlabel("Pitch parameter \\eta")
    call ylabel("Bounce time (microseconds)")
    call title("Bounce Time vs Pitch Parameter (Runtime Dispatch)")
    call legend()
    call savefig("bounce_time_comparison.png")
    print *, 'Saved: bounce_time_comparison.png'
    
    ! Plot 2: Toroidal shift comparison
    call figure(800, 600)
    call plot(eta_array, delphi_thin, label="Thin orbit", linestyle="b-")
    call plot(eta_array, delphi_thick, label="Thick orbit", linestyle="r--")
    call xlabel("Pitch parameter \\eta")
    call ylabel("Toroidal shift \\Delta\\phi")
    call title("Toroidal Shift vs Pitch Parameter (Runtime Dispatch)")
    call legend()
    call savefig("toroidal_shift_comparison.png")
    print *, 'Saved: toroidal_shift_comparison.png'
    
    ! Plot 3: Poloidal frequency comparison
    call figure(800, 600)
    call plot(eta_array, omega_theta_thin/1e3, label="Thin orbit", linestyle="b-")
    call plot(eta_array, omega_theta_thick/1e3, label="Thick orbit", linestyle="r--")
    call xlabel("Pitch parameter \\eta")
    call ylabel("Poloidal frequency (kHz)")
    call title("Poloidal Frequency vs Pitch Parameter (Runtime Dispatch)")
    call legend()
    call savefig("poloidal_frequency_comparison.png")
    print *, 'Saved: poloidal_frequency_comparison.png'
    
    ! Plot 4: Toroidal frequency comparison
    call figure(800, 600)
    call plot(eta_array, omega_phi_thin/1e3, label="Thin orbit", linestyle="b-")
    call plot(eta_array, omega_phi_thick/1e3, label="Thick orbit", linestyle="r--")
    call xlabel("Pitch parameter \\eta")
    call ylabel("Toroidal frequency (kHz)")
    call title("Toroidal Frequency vs Pitch Parameter (Runtime Dispatch)")
    call legend()
    call savefig("toroidal_frequency_comparison.png")
    print *, 'Saved: toroidal_frequency_comparison.png'
    
    ! Plot 5: Relative differences
    call figure(800, 600)
    call plot(eta_array, 100*(taub_thick - taub_thin)/taub_thin, &
              label="\\Delta\\tau_b", linestyle="g-")
    call plot(eta_array, 100*(delphi_thick - delphi_thin)/abs(delphi_thin + 1e-10), &
              label="\\Delta\\phi", linestyle="m-")
    call xlabel("Pitch parameter \\eta")
    call ylabel("Relative difference (%)")
    call title("Orbit Width Effects (Runtime Dispatch)")
    call legend()
    call savefig("relative_differences.png")
    print *, 'Saved: relative_differences.png'
    
    ! Plot 6: All results on one plot
    call figure(1200, 800)
    
    ! Normalized plots for comparison
    call plot(eta_array, taub_thin/maxval(taub_thin), label="\\tau_b thin", linestyle="b-")
    call plot(eta_array, taub_thick/maxval(taub_thick), label="\\tau_b thick", linestyle="b--")
    call plot(eta_array, abs(delphi_thin)/maxval(abs(delphi_thin)), label="\\Delta\\phi thin", linestyle="r-")
    call plot(eta_array, abs(delphi_thick)/maxval(abs(delphi_thick)), label="\\Delta\\phi thick", linestyle="r--")
    
    call xlabel("Pitch parameter \\eta")
    call ylabel("Normalized values")
    call title("Runtime Dispatch Comparison: Thin vs Thick Orbits")
    call legend()
    call savefig("runtime_dispatch_comparison.png")
    print *, 'Saved: runtime_dispatch_comparison.png'
    
    ! Summary statistics
    print *, ''
    print *, 'Summary Statistics:'
    print *, '==================='
    print *, 'Average bounce time (thin):', sum(taub_thin)/neta * 1e6, 'microseconds'
    print *, 'Average bounce time (thick):', sum(taub_thick)/neta * 1e6, 'microseconds'
    print *, 'Average toroidal shift (thin):', sum(delphi_thin)/neta
    print *, 'Average toroidal shift (thick):', sum(delphi_thick)/neta
    print *, 'Max relative difference in bounce time:', &
             maxval(abs(taub_thick - taub_thin)/taub_thin)*100, '%'
    print *, 'Max relative difference in toroidal shift:', &
             maxval(abs(delphi_thick - delphi_thin)/abs(delphi_thin + 1e-10))*100, '%'
    
    print *, ''
    print *, 'Runtime dispatch plotting complete!'
    print *, 'Generated files:'
    print *, '  bounce_time_comparison.png'
    print *, '  toroidal_shift_comparison.png'
    print *, '  poloidal_frequency_comparison.png'
    print *, '  toroidal_frequency_comparison.png'
    print *, '  relative_differences.png'
    print *, '  runtime_dispatch_comparison.png'
    
end program plot_runtime_comparison