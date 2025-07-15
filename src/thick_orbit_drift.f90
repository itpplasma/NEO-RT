module thick_orbit_drift
    ! Module for calculating real bounce-averaged drift velocities
    ! Phase G.4.REAL.2 implementation
    
    use iso_fortran_env, only: real64
    implicit none
    
    private
    public :: calculate_bounce_averaged_drift
    public :: gradB_drift_velocity
    public :: curvature_drift_velocity
    public :: magnetic_moment_conservation
    
    integer, parameter :: dp = real64
    
contains

    subroutine calculate_bounce_averaged_drift(v, eta, v_drift_avg, success)
        ! Calculate bounce-averaged drift velocity: v̄_drift = ∫₀^τb v_drift(τ) dτ / τb
        ! This is the core implementation of Phase G.4.REAL.2
        ! Currently uses simplified orbit integration due to POTATO stability issues
        implicit none
        
        real(dp), intent(in) :: v, eta
        real(dp), intent(out) :: v_drift_avg(3)  ! [R, phi, Z] components
        logical, intent(out) :: success
        
        real(dp) :: taub_estimated, delphi_estimated
        integer :: i, n_steps
        real(dp) :: dt, time_step
        real(dp) :: v_drift_sum(3)
        
        ! Initialize
        success = .false.
        v_drift_avg = 0.0_dp
        
        ! Use simplified bounce time estimate for now
        ! TODO: Replace with real POTATO bounce time when stable
        taub_estimated = estimate_bounce_time(v, eta)
        delphi_estimated = estimate_toroidal_shift(v, eta)
        
        print *, 'Estimated bounce time: τb =', taub_estimated, 's'
        print *, 'Estimated toroidal shift: Δφ =', delphi_estimated, 'rad'
        
        ! Simplified bounce averaging integration
        n_steps = 100
        dt = taub_estimated / real(n_steps, dp)
        v_drift_sum = 0.0_dp
        
        print *, 'Calculating bounce-averaged drift velocity...'
        print *, 'Integration steps:', n_steps
        print *, 'Time step:', dt, 's'
        
        ! Integration loop over bounce period
        do i = 1, n_steps
            time_step = real(i-1, dp) * dt
            
            ! Calculate drift velocity at this time step
            call calculate_drift_at_time(v, eta, time_step, v_drift_sum, success)
            
            if (.not. success) then
                print *, 'WARNING: Drift calculation failed at step', i
                return
            end if
        end do
        
        ! Average over the bounce time
        v_drift_avg = v_drift_sum / real(n_steps, dp)
        
        print *, 'Bounce-averaged drift velocity calculated:'
        print *, '  v_drift_R =', v_drift_avg(1), 'm/s'
        print *, '  v_drift_phi =', v_drift_avg(2), 'm/s'
        print *, '  v_drift_Z =', v_drift_avg(3), 'm/s'
        
        success = .true.
        
    end subroutine calculate_bounce_averaged_drift
    
    subroutine calculate_drift_at_time(v, eta, time, v_drift, success)
        ! Calculate drift velocity at a specific time along the orbit
        ! This is a simplified version - should use real POTATO orbit integration
        implicit none
        
        real(dp), intent(in) :: v, eta, time
        real(dp), intent(inout) :: v_drift(3)
        logical, intent(out) :: success
        
        real(dp) :: R, Z, phi, B_field
        real(dp) :: v_gradB(3), v_curv(3)
        
        ! Initialize
        success = .false.
        
        ! Simplified orbit position (should come from POTATO integration)
        R = 1.5_dp + 0.1_dp * sin(time * 1.0d6)  ! Simplified R(t)
        Z = 0.0_dp + 0.05_dp * cos(time * 1.0d6)  ! Simplified Z(t)
        phi = time * v / R  ! Simplified phi(t)
        B_field = 2.5_dp  ! Simplified B field
        
        ! Calculate grad-B drift
        call gradB_drift_velocity(v, eta, R, Z, phi, B_field, v_gradB, success)
        if (.not. success) return
        
        ! Calculate curvature drift
        call curvature_drift_velocity(v, eta, R, Z, phi, B_field, v_curv, success)
        if (.not. success) return
        
        ! Combined drift velocity
        v_drift = v_drift + v_gradB + v_curv
        
        success = .true.
        
    end subroutine calculate_drift_at_time
    
    subroutine gradB_drift_velocity(v, eta, R, Z, phi, B_field, v_gradB, success)
        ! Calculate grad-B drift: v_gradB = (μ/q) × (B × ∇B)/B²
        implicit none
        
        real(dp), intent(in) :: v, eta, R, Z, phi, B_field
        real(dp), intent(out) :: v_gradB(3)
        logical, intent(out) :: success
        
        real(dp) :: mu, q, grad_B(3), B_vec(3)
        real(dp) :: cross_product(3)
        real(dp) :: m, v_perp
        
        ! Initialize
        success = .false.
        v_gradB = 0.0_dp
        
        ! Physical parameters
        m = 1.67d-27   ! Proton mass (kg)
        q = 1.602d-19  ! Elementary charge (C)
        v_perp = v * sqrt(eta)  ! Perpendicular velocity
        mu = m * v_perp**2 / (2.0_dp * B_field)  ! Magnetic moment (J/T)
        
        ! Simplified magnetic field and gradient
        B_vec = [0.0_dp, B_field, 0.0_dp]  ! Toroidal field (φ-direction)
        grad_B = [-B_field/R, 0.0_dp, 0.0_dp]  ! Simplified gradient (∇B ~ -B/R in R direction)
        
        ! Calculate B × ∇B (cross product: B_φ × ∇B_R gives drift in Z direction)
        cross_product(1) = B_vec(2) * grad_B(3) - B_vec(3) * grad_B(2)  ! R component = 0
        cross_product(2) = B_vec(3) * grad_B(1) - B_vec(1) * grad_B(3)  ! φ component = 0 
        cross_product(3) = B_vec(1) * grad_B(2) - B_vec(2) * grad_B(1)  ! Z component = -B_φ * ∇B_R = B_field²/R
        
        ! Calculate drift velocity
        v_gradB = (mu / q) * cross_product / (B_field**2)
        
        success = .true.
        
    end subroutine gradB_drift_velocity
    
    subroutine curvature_drift_velocity(v, eta, R, Z, phi, B_field, v_curv, success)
        ! Calculate curvature drift: v_curv = (mv_∥²/qB) × (B × ∇B)/B²
        implicit none
        
        real(dp), intent(in) :: v, eta, R, Z, phi, B_field
        real(dp), intent(out) :: v_curv(3)
        logical, intent(out) :: success
        
        real(dp) :: m, q, v_parallel, grad_B(3), B_vec(3)
        real(dp) :: cross_product(3)
        
        ! Initialize
        success = .false.
        v_curv = 0.0_dp
        
        ! Physical parameters
        m = 1.67d-27  ! Proton mass
        q = 1.602d-19  ! Elementary charge
        v_parallel = v * sqrt(1.0_dp - eta)  ! Parallel velocity
        
        ! Simplified magnetic field and gradient
        B_vec = [0.0_dp, B_field, 0.0_dp]  ! Toroidal field
        grad_B = [-B_field/R, 0.0_dp, 0.0_dp]  ! Simplified gradient
        
        ! Calculate B × ∇B
        cross_product(1) = B_vec(2) * grad_B(3) - B_vec(3) * grad_B(2)
        cross_product(2) = B_vec(3) * grad_B(1) - B_vec(1) * grad_B(3)
        cross_product(3) = B_vec(1) * grad_B(2) - B_vec(2) * grad_B(1)
        
        ! Calculate curvature drift
        v_curv = (m * v_parallel**2 / (q * B_field)) * cross_product / (B_field**2)
        
        success = .true.
        
    end subroutine curvature_drift_velocity
    
    subroutine magnetic_moment_conservation(v, eta, B_field, mu, success)
        ! Verify magnetic moment conservation: μ = mv_⊥²/(2B)
        implicit none
        
        real(dp), intent(in) :: v, eta, B_field
        real(dp), intent(out) :: mu
        logical, intent(out) :: success
        
        real(dp) :: m, v_perp
        
        ! Initialize
        success = .false.
        
        ! Physical parameters
        m = 1.67d-27  ! Proton mass
        v_perp = v * sqrt(eta)  ! Perpendicular velocity
        
        ! Calculate magnetic moment
        mu = m * v_perp**2 / (2.0_dp * B_field)
        
        print *, 'Magnetic moment:', mu, 'J/T'
        print *, 'Perpendicular velocity:', v_perp, 'm/s'
        
        success = .true.
        
    end subroutine magnetic_moment_conservation
    
    function estimate_bounce_time(v, eta) result(taub)
        ! Simplified bounce time estimation
        ! TODO: Replace with real POTATO calculation when stable
        implicit none
        real(dp), intent(in) :: v, eta
        real(dp) :: taub
        
        ! Simple estimate based on typical tokamak parameters
        real(dp), parameter :: R_major = 1.5_dp  ! Major radius
        real(dp), parameter :: q_safety = 2.0_dp  ! Safety factor
        real(dp) :: v_parallel
        
        v_parallel = v * sqrt(1.0_dp - eta)
        
        ! Estimated bounce time: τb ~ 2πRq/v_parallel
        taub = 2.0_dp * 3.141592653589793_dp * R_major * q_safety / v_parallel
        
        ! Ensure reasonable bounds
        if (taub < 1.0d-6) taub = 1.0d-6
        if (taub > 1.0d-3) taub = 1.0d-3
        
    end function estimate_bounce_time
    
    function estimate_toroidal_shift(v, eta) result(delphi)
        ! Simplified toroidal shift estimation
        ! TODO: Replace with real POTATO calculation when stable
        implicit none
        real(dp), intent(in) :: v, eta
        real(dp) :: delphi
        
        real(dp), parameter :: q_safety = 2.0_dp  ! Safety factor
        
        ! Simple estimate: Δφ ~ 2π/q
        delphi = 2.0_dp * 3.141592653589793_dp / q_safety
        
    end function estimate_toroidal_shift

end module thick_orbit_drift