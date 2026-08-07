program test_gc_generated_invariant_conversions
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_physical_mu_symbolic, only: evaluate_neort_physical_mu
    use neort_buchholz_action_symbolic, only: evaluate_neort_buchholz_action
    use neort_buchholz_energy_symbolic, only: evaluate_neort_buchholz_energy
    use neort_buchholz_specific_energy_symbolic, only: &
        evaluate_neort_buchholz_specific_energy
    use neort_potato_velocity_symbolic, only: evaluate_neort_potato_velocity
    use neort_potato_mu_symbolic, only: evaluate_neort_potato_mu
    use neort_cylindrical_hamiltonian_symbolic, only: &
        evaluate_neort_cylindrical_hamiltonian
    use neort_cylindrical_canonical_symbolic, only: &
        evaluate_neort_cylindrical_canonical
    use neort_cylindrical_launch_symbolic, only: &
        evaluate_neort_cylindrical_launch
    use util_for_test, only: pass_test
    implicit none

    real(dp), parameter :: mass = 4.0_dp, c_light = 9.0_dp
    real(dp), parameter :: v_perp = 1.75_dp, bmod = 2.5_dp
    real(dp), parameter :: reference_velocity = 3.5_dp
    real(dp), parameter :: electrostatic_potential = 0.4_dp
    real(dp), parameter :: radius = 1.8_dp, p_parallel = -1.25_dp
    real(dp), parameter :: psi = 0.7_dp, bhat_phi = 0.65_dp
    real(dp), parameter :: charges(2) = [2.0_dp, -2.0_dp]
    real(dp) :: mu_phys, j_k, j_potato, expected_omega_c
    real(dp) :: mu_from_vperp, vperp_squared_from_mu
    real(dp) :: specific_energy_from_mu
    real(dp) :: j_k_from_mu, mu_from_j_k
    real(dp) :: omega_c, vperp_squared_from_j_k
    real(dp) :: specific_energy_from_j_k
    real(dp) :: specific_energy_from_j_k_separate
    real(dp) :: j_potato_from_vperp, vperp_squared_from_j_potato
    real(dp) :: j_potato_from_mu, mu_from_j_potato
    real(dp) :: hamiltonian, canonical_p_phi, psi_star
    real(dp) :: v_parallel_squared, p_parallel_from_p_phi
    real(dp) :: launch_energy_residual
    real(dp) :: expected_hamiltonian, expected_p_phi, expected_psi_star
    real(dp) :: expected_v_parallel_squared, expected_p_parallel
    real(dp) :: expected_residual, q_phi_energy
    real(dp) :: jk_by_charge(2), psi_star_by_charge(2)
    real(dp) :: q_phi_by_charge(2)
    real(dp) :: charge
    integer :: charge_index

    mu_phys = mass*v_perp**2/(2.0_dp*bmod)
    j_k = mass*c_light*mu_phys/abs(charges(1))
    j_potato = v_perp**2/(reference_velocity**2*bmod)
    expected_omega_c = abs(charges(1))*bmod/(mass*c_light)

    call evaluate_neort_physical_mu(v_perp, mass, bmod, mu_phys, &
        mu_from_vperp, vperp_squared_from_mu, specific_energy_from_mu)
    call require_close('dimensional mu', mu_from_vperp, mu_phys)
    call require_close('v_perp squared from dimensional mu', &
        vperp_squared_from_mu, v_perp**2)
    call require_close('specific energy from dimensional mu', &
        specific_energy_from_mu, mu_phys*bmod/mass)

    call evaluate_neort_potato_velocity(v_perp, j_potato, &
        reference_velocity, bmod, j_potato_from_vperp, &
        vperp_squared_from_j_potato)
    call require_close('POTATO J_perp from velocity', &
        j_potato_from_vperp, j_potato)
    call require_close('velocity squared from POTATO J_perp', &
        vperp_squared_from_j_potato, v_perp**2)

    call evaluate_neort_potato_mu(mu_phys, j_potato, mass, &
        reference_velocity, j_potato_from_mu, mu_from_j_potato)
    call require_close('POTATO J_perp from mu', j_potato_from_mu, &
        2.0_dp*mu_phys/(mass*reference_velocity**2))
    call require_close('mu from POTATO J_perp', mu_from_j_potato, mu_phys)

    do charge_index = 1, 2
        charge = charges(charge_index)
        j_k = mass*c_light*mu_phys/abs(charge)

        call evaluate_neort_buchholz_action(mu_phys, j_k, mass, charge, &
            c_light, j_k_from_mu, mu_from_j_k)
        call require_close('Buchholz J_K with mass factor', j_k_from_mu, &
            mass*c_light*mu_phys/abs(charge))
        call require_close('mu from Buchholz J_K', mu_from_j_k, mu_phys)

        call evaluate_neort_buchholz_energy(j_k, mass, charge, c_light, &
            bmod, omega_c, vperp_squared_from_j_k, &
            specific_energy_from_j_k)
        call require_close('Buchholz cyclotron frequency', omega_c, &
            abs(charge)*bmod/(mass*c_light))
        call require_close('velocity squared from Buchholz J_K', &
            vperp_squared_from_j_k, 2.0_dp*j_k*expected_omega_c/mass)
        call require_close('specific energy from Buchholz energy', &
            specific_energy_from_j_k, j_k*expected_omega_c/mass)

        call evaluate_neort_buchholz_specific_energy(j_k, &
            expected_omega_c, mass, specific_energy_from_j_k_separate)
        call require_close('separate Buchholz specific energy', &
            specific_energy_from_j_k_separate, j_k*expected_omega_c/mass)

        q_phi_energy = charge*electrostatic_potential
        expected_hamiltonian = p_parallel**2/(2.0_dp*mass) + &
            mu_phys*bmod + q_phi_energy
        expected_p_phi = charge/c_light*psi + &
            p_parallel*radius*bhat_phi
        expected_psi_star = c_light/charge*expected_p_phi

        call evaluate_neort_cylindrical_hamiltonian(mass, charge, mu_phys, &
            bmod, electrostatic_potential, p_parallel, hamiltonian)
        call require_close('Hamiltonian including q Phi', hamiltonian, &
            expected_hamiltonian)

        call evaluate_neort_cylindrical_canonical(charge, c_light, radius, &
            p_parallel, psi, bhat_phi, canonical_p_phi, psi_star)
        call require_close('canonical p_phi', canonical_p_phi, expected_p_phi)
        call require_close('psi_star', psi_star, expected_psi_star)

        expected_v_parallel_squared = 2.0_dp*(expected_hamiltonian - &
            mu_phys*bmod - q_phi_energy)/mass
        expected_p_parallel = (expected_p_phi - charge*psi/c_light)/ &
            (radius*bhat_phi)
        expected_residual = expected_p_parallel**2 - &
            mass**2*expected_v_parallel_squared

        call evaluate_neort_cylindrical_launch(mass, charge, c_light, &
            mu_phys, bmod, electrostatic_potential, radius, psi, &
            bhat_phi, hamiltonian, canonical_p_phi, v_parallel_squared, &
            p_parallel_from_p_phi, launch_energy_residual)
        call require_close('v_parallel squared from fixed H', &
            v_parallel_squared, expected_v_parallel_squared)
        call require_close('fixed-invariant p_parallel reconstruction', &
            p_parallel_from_p_phi, expected_p_parallel)
        call require_close('launch energy residual', launch_energy_residual, &
            expected_residual)

        jk_by_charge(charge_index) = j_k_from_mu
        psi_star_by_charge(charge_index) = psi_star
        q_phi_by_charge(charge_index) = q_phi_energy
    end do

    call require_close('J_K is even in charge', jk_by_charge(1), &
        jk_by_charge(2))
    call require_close('signed psi_star contribution', &
        psi_star_by_charge(1) - psi, &
        -(psi_star_by_charge(2) - psi))
    call require_close('signed q Phi reversal', q_phi_by_charge(1), &
        -q_phi_by_charge(2))

    write (*, '(A)') 'test_gc_generated_invariant_conversions OK'
    call pass_test

contains

    subroutine require_close(label, actual, expected)
        character(*), intent(in) :: label
        real(dp), intent(in) :: actual, expected
        real(dp) :: scale

        scale = max(1.0_dp, max(abs(actual), abs(expected)))
        if (abs(actual - expected) > 2.0e-12_dp*scale) then
            write (*, '(A,2(1X,ES24.16))') trim(label), actual, expected
            error stop 'generated invariant conversion oracle failed'
        end if
    end subroutine require_close

end program test_gc_generated_invariant_conversions
