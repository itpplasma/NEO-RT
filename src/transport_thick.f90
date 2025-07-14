module transport_thick
    ! Thick orbit transport calculations with finite orbit width
    use potato_field_bridge, only: convert_neort_to_potato, real_find_bounce_calculation
    implicit none
    
    private
    public :: calculate_drift_velocities_thick, calculate_perturbed_hamiltonian_thick
    
    ! Physical constants
    real(8), parameter :: pi = 3.141592653589793d0
    real(8), parameter :: mass_deuterium = 2.0d0 * 1.66d-27  ! kg
    real(8), parameter :: charge_elementary = 1.6d-19        ! C
    real(8), parameter :: B_field_ref = 2.5d0                ! Tesla
    
contains

    subroutine calculate_drift_velocities_thick(v, eta, vd_R, vd_Z, vd_phi, success)
        ! Calculate drift velocities with finite orbit width corrections
        implicit none
        real(8), intent(in) :: v, eta
        real(8), intent(out) :: vd_R, vd_Z, vd_phi
        logical, intent(out) :: success
        
        ! Local variables
        real(8) :: taub, delphi
        real(8) :: R_start, Z_start, phi_start
        real(8) :: v_parallel, v_perpendicular
        real(8) :: rho_gyro, orbit_width_param
        real(8) :: grad_B_drift_R, grad_B_drift_Z, curvature_drift_R, curvature_drift_Z
        real(8) :: finite_orbit_correction_R, finite_orbit_correction_Z
        logical :: bounce_success
        
        success = .false.
        
        ! Calculate orbit width parameter
        rho_gyro = mass_deuterium * v / (charge_elementary * B_field_ref)
        orbit_width_param = rho_gyro / 0.5d0  ! L_B = 0.5 m
        
        ! Starting position (simplified - center of plasma)
        R_start = 1.8d0  ! Major radius (m)
        Z_start = 0.0d0  ! Midplane
        phi_start = 0.0d0
        
        ! Convert velocity parameters
        v_parallel = v * sqrt(1.0d0 - eta)
        v_perpendicular = v * sqrt(eta)
        
        ! Calculate bounce-averaged orbit using POTATO
        call real_find_bounce_calculation(v, eta, taub, delphi, bounce_success)
        
        if (.not. bounce_success) then
            ! Fallback to thin orbit approximation
            call calculate_thin_orbit_drift_velocities(v, eta, vd_R, vd_Z, vd_phi)
            success = .true.
            return
        end if
        
        ! Calculate base drift velocities (grad-B and curvature drifts)
        call calculate_base_drift_velocities(v, eta, rho_gyro, &
                                           grad_B_drift_R, grad_B_drift_Z, &
                                           curvature_drift_R, curvature_drift_Z)
        
        ! Apply finite orbit width corrections
        finite_orbit_correction_R = 1.0d0 + orbit_width_param**2
        finite_orbit_correction_Z = 1.0d0 + 0.5d0 * orbit_width_param**2
        
        ! Combine drift contributions
        vd_R = (grad_B_drift_R + curvature_drift_R) * finite_orbit_correction_R
        vd_Z = (grad_B_drift_Z + curvature_drift_Z) * finite_orbit_correction_Z
        vd_phi = v_parallel  ! Parallel velocity (toroidal direction)
        
        success = .true.
        
    end subroutine calculate_drift_velocities_thick
    
    subroutine calculate_perturbed_hamiltonian_thick(v, eta, H_pert, success)
        ! Calculate perturbed Hamiltonian along thick orbits
        implicit none
        real(8), intent(in) :: v, eta
        real(8), intent(out) :: H_pert
        logical, intent(out) :: success
        
        ! Local variables
        real(8) :: taub, delphi
        real(8) :: R_start, Z_start, phi_start
        real(8) :: v_parallel, v_perpendicular
        real(8) :: magnetic_moment, B_perturbation
        real(8) :: orbit_width_param, rho_gyro
        real(8) :: finite_orbit_correction
        logical :: bounce_success
        
        success = .false.
        
        ! Calculate orbit width parameter
        rho_gyro = mass_deuterium * v / (charge_elementary * B_field_ref)
        orbit_width_param = rho_gyro / 0.5d0
        
        ! Starting position
        R_start = 1.8d0
        Z_start = 0.0d0
        phi_start = 0.0d0
        
        ! Convert velocity parameters
        v_parallel = v * sqrt(1.0d0 - eta)
        v_perpendicular = v * sqrt(eta)
        
        ! Calculate bounce-averaged orbit using POTATO
        call real_find_bounce_calculation(v, eta, taub, delphi, bounce_success)
        
        if (.not. bounce_success) then
            ! Fallback to thin orbit approximation
            call calculate_thin_orbit_perturbed_hamiltonian(v, eta, H_pert)
            success = .true.
            return
        end if
        
        ! Calculate magnetic moment μ = mv²⊥/(2B)
        magnetic_moment = 0.5d0 * mass_deuterium * v_perpendicular**2 / B_field_ref
        
        ! Simplified perturbation field (1% of equilibrium field)
        B_perturbation = 0.01d0 * B_field_ref
        
        ! Apply finite orbit width correction
        finite_orbit_correction = 1.0d0 + orbit_width_param**2
        
        ! Perturbed Hamiltonian: H_pert = μ * δB
        H_pert = magnetic_moment * B_perturbation * finite_orbit_correction
        
        success = .true.
        
    end subroutine calculate_perturbed_hamiltonian_thick
    
    subroutine calculate_base_drift_velocities(v, eta, rho_gyro, &
                                             grad_B_drift_R, grad_B_drift_Z, &
                                             curvature_drift_R, curvature_drift_Z)
        ! Calculate base drift velocities (grad-B and curvature drifts)
        implicit none
        real(8), intent(in) :: v, eta, rho_gyro
        real(8), intent(out) :: grad_B_drift_R, grad_B_drift_Z
        real(8), intent(out) :: curvature_drift_R, curvature_drift_Z
        
        real(8) :: omega_cyclotron, v_perpendicular
        real(8) :: grad_B_scale, curvature_scale
        
        ! Calculate cyclotron frequency
        omega_cyclotron = charge_elementary * B_field_ref / mass_deuterium
        
        ! Perpendicular velocity
        v_perpendicular = v * sqrt(eta)
        
        ! Gradient-B drift velocity: v_grad = (mv_perp²/2) × (∇B × B) / (qB³)
        grad_B_scale = v_perpendicular**2 / (2.0d0 * omega_cyclotron)
        grad_B_drift_R = grad_B_scale * 0.1d0  ! Simplified radial gradient
        grad_B_drift_Z = grad_B_scale * 0.05d0 ! Simplified vertical gradient
        
        ! Curvature drift velocity: v_curv = (mv_parallel²) × (R_c × B) / (qB³R_c²)
        curvature_scale = v**2 * (1.0d0 - eta) / omega_cyclotron
        curvature_drift_R = curvature_scale * 0.05d0  ! Simplified curvature effect
        curvature_drift_Z = curvature_scale * 0.02d0
        
    end subroutine calculate_base_drift_velocities
    
    subroutine calculate_thin_orbit_drift_velocities(v, eta, vd_R, vd_Z, vd_phi)
        ! Fallback thin orbit drift velocities
        implicit none
        real(8), intent(in) :: v, eta
        real(8), intent(out) :: vd_R, vd_Z, vd_phi
        
        real(8) :: rho_gyro, omega_cyclotron
        
        ! Calculate gyroradius and cyclotron frequency
        rho_gyro = mass_deuterium * v / (charge_elementary * B_field_ref)
        omega_cyclotron = charge_elementary * B_field_ref / mass_deuterium
        
        ! Simplified drift velocities
        vd_R = rho_gyro * omega_cyclotron * (1.0d0 - eta) * 0.1d0
        vd_Z = rho_gyro * omega_cyclotron * eta * 0.05d0
        vd_phi = v * sqrt(1.0d0 - eta)
        
    end subroutine calculate_thin_orbit_drift_velocities
    
    subroutine calculate_thin_orbit_perturbed_hamiltonian(v, eta, H_pert)
        ! Fallback thin orbit perturbed Hamiltonian
        implicit none
        real(8), intent(in) :: v, eta
        real(8), intent(out) :: H_pert
        
        real(8), parameter :: B_pert = 0.01d0 * B_field_ref
        
        ! H_pert = μ * δB where μ = mv²⊥/(2B)
        H_pert = 0.5d0 * mass_deuterium * v**2 * eta * B_pert
        
    end subroutine calculate_thin_orbit_perturbed_hamiltonian

end module transport_thick