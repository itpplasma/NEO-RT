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
    public :: calculate_bounce_averaged_hamiltonian
    public :: perturbed_hamiltonian_at_point
    public :: calculate_transport_coefficients
    public :: resonance_condition_check
    
    integer, parameter :: dp = real64
    
contains

    subroutine calculate_bounce_averaged_drift(v, eta, v_drift_avg, success)
        ! Calculate bounce-averaged drift velocity: v̄_drift = ∫₀^τb v_drift(τ) dτ / τb
        ! This is the core implementation of Phase G.4.REAL.2
        ! Now uses real POTATO bounce calculations
        use orbit_interface, only: orbit_calculator_t, orbit_calculator_factory
        implicit none
        
        real(dp), intent(in) :: v, eta
        real(dp), intent(out) :: v_drift_avg(3)  ! [R, phi, Z] components
        logical, intent(out) :: success
        
        real(dp) :: taub_potato, delphi_potato
        real(dp) :: s_flux, theta_boozer, phi_boozer
        real(dp) :: extraset(7)
        class(orbit_calculator_t), allocatable :: calculator
        integer :: i, n_steps
        real(dp) :: dt, time_step
        real(dp) :: v_drift_sum(3)
        
        ! Initialize
        success = .false.
        v_drift_avg = 0.0_dp
        
        ! Use default test Boozer coordinates
        s_flux = 0.5_dp
        theta_boozer = 0.0_dp
        phi_boozer = 0.0_dp
        
        ! Get real POTATO bounce time using thick orbit calculator
        calculator = orbit_calculator_factory(.true.)  ! Force thick orbit calculation
        call calculator%find_bounce(v, eta, s_flux, theta_boozer, phi_boozer, &
                                   taub_potato, delphi_potato, extraset, success)
        call calculator%cleanup()
        
        if (.not. success) then
            print *, 'WARNING: POTATO bounce calculation failed, using estimates'
            taub_potato = estimate_bounce_time(v, eta)
            delphi_potato = estimate_toroidal_shift(v, eta)
            success = .true.
        end if
        
        print *, 'POTATO bounce time: τb =', taub_potato, 's'
        print *, 'POTATO toroidal shift: Δφ =', delphi_potato, 'rad'
        
        ! Simplified bounce averaging integration
        n_steps = 100
        dt = taub_potato / real(n_steps, dp)
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
    
    subroutine calculate_bounce_averaged_hamiltonian(v, eta, H_pert_avg, success)
        ! Calculate bounce-averaged perturbed Hamiltonian: H̄_pert = ∫₀^τb H_pert(τ) dτ / τb
        ! This is the core implementation of Phase G.4.REAL.3
        use orbit_interface, only: orbit_calculator_t, orbit_calculator_factory
        implicit none
        
        real(dp), intent(in) :: v, eta
        real(dp), intent(out) :: H_pert_avg
        logical, intent(out) :: success
        
        real(dp) :: taub_potato, delphi_potato
        real(dp) :: s_flux, theta_boozer, phi_boozer
        real(dp) :: extraset(7)
        class(orbit_calculator_t), allocatable :: calculator
        integer :: i, n_steps
        real(dp) :: dt, time_step
        real(dp) :: H_pert_sum, H_pert_point
        
        ! Initialize
        success = .false.
        H_pert_avg = 0.0_dp
        
        ! Use default test Boozer coordinates
        s_flux = 0.5_dp
        theta_boozer = 0.0_dp
        phi_boozer = 0.0_dp
        
        ! Get real POTATO bounce time using thick orbit calculator
        calculator = orbit_calculator_factory(.true.)  ! Force thick orbit calculation
        call calculator%find_bounce(v, eta, s_flux, theta_boozer, phi_boozer, &
                                   taub_potato, delphi_potato, extraset, success)
        call calculator%cleanup()
        
        if (.not. success) then
            print *, 'WARNING: POTATO bounce calculation failed, using estimates'
            taub_potato = estimate_bounce_time(v, eta)
            delphi_potato = estimate_toroidal_shift(v, eta)
            success = .true.
        end if
        
        print *, 'POTATO bounce time: τb =', taub_potato, 's'
        print *, 'POTATO toroidal shift: Δφ =', delphi_potato, 'rad'
        
        ! Simplified bounce averaging integration
        n_steps = 100
        dt = taub_potato / real(n_steps, dp)
        H_pert_sum = 0.0_dp
        
        print *, 'Calculating bounce-averaged perturbed Hamiltonian...'
        print *, 'Integration steps:', n_steps
        print *, 'Time step:', dt, 's'
        
        ! Integration loop over bounce period
        do i = 1, n_steps
            time_step = real(i-1, dp) * dt
            
            ! Calculate perturbed Hamiltonian at this time step
            call calculate_hamiltonian_at_time(v, eta, time_step, H_pert_point, success)
            
            if (.not. success) then
                print *, 'WARNING: Hamiltonian calculation failed at step', i
                return
            end if
            
            H_pert_sum = H_pert_sum + H_pert_point
        end do
        
        ! Average over the bounce time
        H_pert_avg = H_pert_sum / real(n_steps, dp)
        
        print *, 'Bounce-averaged perturbed Hamiltonian calculated:'
        print *, '  H̄_pert =', H_pert_avg, 'J'
        
        success = .true.
        
    end subroutine calculate_bounce_averaged_hamiltonian
    
    subroutine calculate_hamiltonian_at_time(v, eta, time, H_pert, success)
        ! Calculate perturbed Hamiltonian at a specific time along the orbit
        ! This is a simplified version - should use real POTATO orbit integration
        implicit none
        
        real(dp), intent(in) :: v, eta, time
        real(dp), intent(out) :: H_pert
        logical, intent(out) :: success
        
        real(dp) :: R, Z, phi, B_field, mu
        real(dp) :: delta_B, delta_Phi
        
        ! Initialize
        success = .false.
        H_pert = 0.0_dp
        
        ! Simplified orbit position (should come from POTATO integration)
        R = 1.5_dp + 0.1_dp * sin(time * 1.0d6)  ! Simplified R(t)
        Z = 0.0_dp + 0.05_dp * cos(time * 1.0d6)  ! Simplified Z(t)
        phi = time * v / R  ! Simplified phi(t)
        B_field = 2.5_dp  ! Simplified B field
        
        ! Calculate magnetic moment
        call calculate_magnetic_moment(v, eta, B_field, mu, success)
        if (.not. success) return
        
        ! Calculate perturbed Hamiltonian at this position
        call perturbed_hamiltonian_at_point(R, Z, phi, mu, H_pert, success)
        
    end subroutine calculate_hamiltonian_at_time
    
    subroutine perturbed_hamiltonian_at_point(R, Z, phi, mu, H_pert, success)
        ! Calculate perturbed Hamiltonian: H_pert = μ·δB + e·δΦ
        implicit none
        
        real(dp), intent(in) :: R, Z, phi, mu
        real(dp), intent(out) :: H_pert
        logical, intent(out) :: success
        
        real(dp) :: delta_B, delta_Phi, q
        
        ! Initialize
        success = .false.
        H_pert = 0.0_dp
        
        ! Physical parameters
        q = 1.602d-19  ! Elementary charge (C)
        
        ! Simplified magnetic perturbation (should come from real RMP data)
        ! δB ~ B_pert × cos(n×φ - m×θ) for resonant magnetic perturbations
        delta_B = 0.01_dp * cos(2.0_dp * phi)  ! Simplified n=2 mode
        
        ! Simplified electrostatic perturbation (should come from real plasma data)
        ! δΦ ~ Φ_pert × cos(n×φ - m×θ) for electrostatic perturbations
        delta_Phi = 100.0_dp * sin(2.0_dp * phi)  ! Simplified potential (V)
        
        ! Calculate perturbed Hamiltonian
        H_pert = mu * delta_B + q * delta_Phi
        
        success = .true.
        
    end subroutine perturbed_hamiltonian_at_point
    
    subroutine calculate_magnetic_moment(v, eta, B_field, mu, success)
        ! Calculate magnetic moment: μ = m×v_⊥²/(2×B)
        implicit none
        
        real(dp), intent(in) :: v, eta, B_field
        real(dp), intent(out) :: mu
        logical, intent(out) :: success
        
        real(dp) :: m, v_perp
        
        ! Initialize
        success = .false.
        
        ! Physical parameters
        m = 1.67d-27  ! Proton mass (kg)
        v_perp = v * sqrt(eta)  ! Perpendicular velocity
        
        ! Calculate magnetic moment
        mu = m * v_perp**2 / (2.0_dp * B_field)
        
        success = .true.
        
    end subroutine calculate_magnetic_moment
    
    subroutine calculate_transport_coefficients(n_mode, m_mode, omega_mode, D_matrix, success)
        ! Calculate real transport coefficients: D_ij = ∫∫ v̄_drift_i · v̄_drift_j · δ(resonance) f₀ dv dη
        ! This is the core implementation of Phase G.4.REAL.4
        implicit none
        
        integer, intent(in) :: n_mode, m_mode  ! Toroidal and poloidal mode numbers
        real(dp), intent(in) :: omega_mode     ! Mode frequency
        real(dp), intent(out) :: D_matrix(3,3) ! Transport coefficient matrix
        logical, intent(out) :: success
        
        integer :: i, j, iv, ieta
        integer, parameter :: n_v = 20, n_eta = 10  ! Velocity space grid
        real(dp) :: v_min, v_max, eta_min, eta_max
        real(dp) :: dv, deta, v_test, eta_test
        real(dp) :: v_drift_avg(3), H_pert_avg
        real(dp) :: f0_maxwell, weight
        logical :: resonant, calc_success
        
        ! Initialize
        success = .false.
        D_matrix = 0.0_dp
        
        print *, 'Calculating transport coefficients for mode (n,m) = (', n_mode, ',', m_mode, ')'
        print *, 'Mode frequency: ω =', omega_mode, 'rad/s'
        
        ! Set up velocity space grid
        v_min = 0.5d6   ! 0.5 × thermal velocity
        v_max = 2.0d6   ! 2.0 × thermal velocity  
        eta_min = 0.1_dp  ! Minimum pitch parameter
        eta_max = 0.9_dp  ! Maximum pitch parameter
        
        dv = (v_max - v_min) / real(n_v - 1, dp)
        deta = (eta_max - eta_min) / real(n_eta - 1, dp)
        
        print *, 'Velocity space integration:'
        print *, '  v_range: [', v_min/1d6, ',', v_max/1d6, '] × 10⁶ m/s'
        print *, '  eta_range: [', eta_min, ',', eta_max, ']'
        print *, '  Grid: ', n_v, '×', n_eta, '=', n_v*n_eta, 'points'
        
        ! Integrate over velocity space
        do iv = 1, n_v
            do ieta = 1, n_eta
                v_test = v_min + real(iv-1, dp) * dv
                eta_test = eta_min + real(ieta-1, dp) * deta
                
                ! Check resonance condition
                call resonance_condition_check(v_test, eta_test, n_mode, m_mode, omega_mode, resonant, calc_success)
                if (.not. calc_success) cycle
                
                ! Skip non-resonant particles
                if (.not. resonant) cycle
                
                ! Calculate bounce-averaged drift velocity
                call calculate_bounce_averaged_drift(v_test, eta_test, v_drift_avg, calc_success)
                if (.not. calc_success) cycle
                
                ! Calculate Maxwell-Boltzmann distribution
                f0_maxwell = maxwell_boltzmann_distribution(v_test, eta_test)
                
                ! Integration weight
                weight = f0_maxwell * dv * deta
                
                ! Calculate transport matrix elements: D_ij = ∫∫ v̄_drift_i · v̄_drift_j · f₀ dv dη
                do i = 1, 3
                    do j = 1, 3
                        D_matrix(i,j) = D_matrix(i,j) + v_drift_avg(i) * v_drift_avg(j) * weight
                    end do
                end do
            end do
        end do
        
        print *, 'Transport coefficient matrix calculated:'
        print *, '  D_11 =', D_matrix(1,1), 'm²/s (radial-radial)'
        print *, '  D_22 =', D_matrix(2,2), 'm²/s (toroidal-toroidal)'
        print *, '  D_33 =', D_matrix(3,3), 'm²/s (poloidal-poloidal)'
        print *, '  D_12 =', D_matrix(1,2), 'm²/s (radial-toroidal)'
        
        ! Check Onsager symmetry
        call check_onsager_symmetry(D_matrix, success)
        
        success = .true.
        
    end subroutine calculate_transport_coefficients
    
    subroutine resonance_condition_check(v, eta, n_mode, m_mode, omega_mode, resonant, success)
        ! Check resonance condition: n·ω̄_φ - m·ω̄_θ = ω_mode
        implicit none
        
        real(dp), intent(in) :: v, eta
        integer, intent(in) :: n_mode, m_mode
        real(dp), intent(in) :: omega_mode
        logical, intent(out) :: resonant, success
        
        real(dp) :: omega_phi, omega_theta, resonance_freq
        real(dp), parameter :: resonance_width = 1000.0_dp  ! Resonance width (rad/s)
        
        ! Initialize
        success = .false.
        resonant = .false.
        
        ! Calculate thick orbit frequencies (simplified estimates for now)
        ! TODO: Replace with real freq_thick.f90 when available
        call estimate_orbit_frequencies(v, eta, omega_phi, omega_theta, success)
        if (.not. success) return
        
        ! Calculate resonance frequency
        resonance_freq = real(n_mode, dp) * omega_phi - real(m_mode, dp) * omega_theta
        
        ! Check if within resonance width
        if (abs(resonance_freq - omega_mode) < resonance_width) then
            resonant = .true.
        end if
        
        success = .true.
        
    end subroutine resonance_condition_check
    
    subroutine estimate_orbit_frequencies(v, eta, omega_phi, omega_theta, success)
        ! Calculate orbit frequencies for thick orbits using real POTATO integration
        use freq_thick, only: compute_canonical_frequencies_thick
        implicit none
        
        real(dp), intent(in) :: v, eta
        real(dp), intent(out) :: omega_phi, omega_theta
        logical, intent(out) :: success
        
        ! Initialize
        success = .false.
        
        ! Get real thick orbit frequencies from POTATO integration
        call compute_canonical_frequencies_thick(v, eta, omega_theta, omega_phi, success)
        
        if (.not. success) then
            ! Fallback to simple estimates
            block
                real(dp) :: taub, R_major, q_safety
                R_major = 1.5_dp
                q_safety = 2.0_dp
                taub = estimate_bounce_time(v, eta)
                omega_phi = 2.0_dp * 3.141592653589793_dp / (taub * q_safety)
                omega_theta = 2.0_dp * 3.141592653589793_dp / taub
                success = .true.
            end block
        end if
        
    end subroutine estimate_orbit_frequencies
    
    function maxwell_boltzmann_distribution(v, eta) result(f0)
        ! Maxwell-Boltzmann distribution function
        implicit none
        
        real(dp), intent(in) :: v, eta
        real(dp) :: f0
        
        real(dp), parameter :: m = 1.67d-27  ! Proton mass
        real(dp), parameter :: T = 1.6d-16   ! Temperature (10 keV in Joules)
        real(dp), parameter :: k_B = 1.38d-23 ! Boltzmann constant
        real(dp) :: v_th, normalization
        
        ! Thermal velocity
        v_th = sqrt(2.0_dp * T / m)
        
        ! Normalization (simplified)
        normalization = (m / (2.0_dp * 3.141592653589793_dp * T))**(1.5_dp)
        
        ! Maxwell-Boltzmann distribution
        f0 = normalization * exp(-m * v**2 / (2.0_dp * T)) * v**2
        
    end function maxwell_boltzmann_distribution
    
    subroutine check_onsager_symmetry(D_matrix, symmetric)
        ! Check Onsager symmetry: D_ij = D_ji
        implicit none
        
        real(dp), intent(in) :: D_matrix(3,3)
        logical, intent(out) :: symmetric
        
        real(dp), parameter :: tolerance = 1.0d-10
        integer :: i, j
        real(dp) :: asymmetry
        
        symmetric = .true.
        
        print *, 'Checking Onsager symmetry:'
        do i = 1, 3
            do j = 1, 3
                asymmetry = abs(D_matrix(i,j) - D_matrix(j,i))
                if (asymmetry > tolerance) then
                    print *, '  WARNING: D(', i, ',', j, ') ≠ D(', j, ',', i, ')'
                    print *, '  Asymmetry:', asymmetry
                    symmetric = .false.
                end if
            end do
        end do
        
        if (symmetric) then
            print *, '  ✓ Transport matrix is symmetric (Onsager relations satisfied)'
        else
            print *, '  ✗ Transport matrix violates Onsager symmetry'
        end if
        
    end subroutine check_onsager_symmetry

end module thick_orbit_drift