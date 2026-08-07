program test_gc_full_fow_symbolic_kernel
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_axisymmetric_noether_symbolic, only: &
        evaluate_neort_axisymmetric_noether
    use neort_full_fow_action_symbolic, only: &
        evaluate_neort_action_normalization
    use neort_full_fow_eq17_symbolic, only: evaluate_neort_eq17_force
    use neort_full_fow_perturbation_symbolic, only: &
        evaluate_neort_perturbation_coefficient
    use neort_full_fow_normalization_symbolic, only: &
        evaluate_neort_full_fow_normalization
    use neort_full_fow_quadrature_map_symbolic, only: &
        evaluate_neort_full_fow_quadrature_map
    use neort_full_fow_resonance_symbolic, only: &
        evaluate_neort_resonance_weights
    use neort_polynomial_cell_enclosure_symbolic, only: &
        evaluate_neort_polynomial_cell_enclosure
    use util_for_test, only: pass_test
    implicit none

    real(dp), parameter :: mass = 2.0_dp, charge = -3.0_dp
    real(dp), parameter :: c_light = 5.0_dp, mu = 0.7_dp, bmod = 4.0_dp
    real(dp), parameter :: h = 9.0_dp, electrostatic_potential = 0.2_dp
    real(dp), parameter :: p_phi = -1.4_dp, p_parallel = 0.8_dp
    real(dp), parameter :: temperature = 1.5_dp
    real(dp), parameter :: n_mode = -4.0_dp, tau = 3.2_dp, g_prime = -0.7_dp
    real(dp), parameter :: a1 = 0.2_dp, a2 = -0.4_dp
    real(dp) :: j_k, omega_c, j_k_omega_c, j_k_candidate, j_k_max
    real(dp) :: d_j_k_d_h, d_j_k_d_potential, d_j_k_d_omega
    real(dp) :: eq4_coefficient, psi_star, d_psi_star_d_p_phi
    real(dp) :: q_phi_energy, e_on_shell, c_on_shell, eta_on_shell
    real(dp) :: perturbation_ratio, boozer_ratio
    real(dp) :: frequency_phase_derivative, frequency_root_weight
    real(dp) :: phase_root_weight
    real(dp) :: force_bracket
    real(dp) :: canonical_p_phi, d_p_phi_dt, d_lagrangian_d_phi_dot
    real(dp) :: canonical_p_phi_residual
    real(dp) :: normalization(13), quadrature_map(8), polynomial_enclosure(3)

    call evaluate_neort_action_normalization(mass, charge, c_light, mu, bmod, &
        h, electrostatic_potential, p_phi, 1.2_dp, j_k, omega_c, &
        j_k_omega_c, j_k_candidate, j_k_max, d_j_k_d_h, &
        d_j_k_d_potential, d_j_k_d_omega, eq4_coefficient, psi_star, &
        d_psi_star_d_p_phi)
    call require_close("J_K", j_k, mass*c_light*mu/abs(charge))
    call require_close("omega_c", omega_c, abs(charge)*bmod/(mass*c_light))
    call require_close("J_K omega_c", j_k_omega_c, mu*bmod)
    call require_close("positive J_K candidate", j_k_candidate, 8.0_dp)
    call require_close("positive-part J_K bound", j_k_max, 8.0_dp)
    call require_close("d candidate/dH", d_j_k_d_h, 1.0_dp/1.2_dp)
    call require_close("d candidate/dPhi", d_j_k_d_potential, 2.5_dp)
    call require_close("d candidate/domega", d_j_k_d_omega, -6.666666666666667_dp)
    call require_close("Eq4 coefficient", eq4_coefficient, 16.4_dp)
    call require_close("psi_star", psi_star, c_light/charge*p_phi)
    call require_close("d psi_star/dp_phi", d_psi_star_d_p_phi, c_light/charge)

    call evaluate_neort_perturbation_coefficient(mass, charge, mu, bmod, h, &
        electrostatic_potential, p_parallel, temperature, q_phi_energy, &
        e_on_shell, c_on_shell, eta_on_shell, perturbation_ratio, boozer_ratio)
    call require_close("q Phi energy", q_phi_energy, charge*electrostatic_potential)
    call require_close("on-shell E", e_on_shell, 2.96_dp)
    call require_close("on-shell C", c_on_shell, 3.12_dp)
    call require_close("C/E limiting expression", perturbation_ratio, &
        boozer_ratio)
    call require_close("eta", eta_on_shell, mu/e_on_shell)

    call evaluate_neort_resonance_weights(n_mode, tau, g_prime, &
        frequency_phase_derivative, frequency_root_weight, phase_root_weight)
    call require_close("F prime", frequency_phase_derivative, 0.875_dp)
    call require_close("root/phase weight identity", frequency_root_weight, &
        phase_root_weight)

    call evaluate_neort_eq17_force(h, charge, electrostatic_potential, &
        temperature, a1, a2, q_phi_energy, force_bracket)
    call require_close("Eq17 q Phi", q_phi_energy, charge*electrostatic_potential)
    call require_close("Eq17 force bracket", force_bracket, -2.36_dp)

    call evaluate_neort_axisymmetric_noether(charge, c_light, 1.1_dp, 0.4_dp, &
        p_parallel, 2.96_dp, 0.2_dp, canonical_p_phi, d_p_phi_dt, &
        d_lagrangian_d_phi_dot, canonical_p_phi_residual)
    call require_close("canonical p_phi", canonical_p_phi, -0.34_dp)
    call require_close("axisymmetric dp_phi/dt", d_p_phi_dt, 0.0_dp)
    call require_close("dL/dphi_dot", d_lagrangian_d_phi_dot, canonical_p_phi)
    call require_close("Noether residual", canonical_p_phi_residual, 0.0_dp)

    call evaluate_neort_full_fow_normalization(mass, charge, c_light, 2.0_dp, &
        4.0_dp, h, 10.0_dp, 8.0_dp, 2.0_dp, 0.5_dp, 6.0_dp, -2.0_dp, &
        1.0_dp, -0.5_dp, 4.0_dp, -6.0_dp, normalization(1), &
        normalization(2), normalization(3), normalization(4), &
        normalization(5), normalization(6), normalization(7), &
        normalization(8), normalization(9), normalization(10), &
        normalization(11), normalization(12), normalization(13))
    call require_close("Phi_eff", normalization(1), 5.0_dp/6.0_dp)
    call require_close("J_K scale", normalization(2), 20.0_dp/3.0_dp)
    call require_close("H hat", normalization(3), 4.5_dp)
    call require_close("J_K hat", normalization(4), 1.5_dp)
    call require_close("psi_star hat", normalization(5), 9.6_dp)
    call require_close("tau hat", normalization(7), 2.0_dp)
    call require_close("omega_b hat", normalization(8), 1.5_dp)
    call require_close("Hm real hat", normalization(12), 2.0_dp)
    call require_close("Hm imag hat", normalization(13), -3.0_dp)

    call evaluate_neort_full_fow_quadrature_map(mass, charge, c_light, 2.0_dp, &
        0.25_dp, 0.4_dp, 0.3_dp, 0.6_dp, -1.0_dp, 5.0_dp, &
        quadrature_map(1), quadrature_map(2), quadrature_map(3), &
        quadrature_map(4), quadrature_map(5), quadrature_map(6), &
        quadrature_map(7), quadrature_map(8))
    call require_close("mapped H hat", quadrature_map(1), -1.0_dp/6.0_dp)
    call require_close("mapped H physical", quadrature_map(2), -1.0_dp/3.0_dp)
    call require_close("energy Jacobian weight", quadrature_map(3), 32.0_dp/45.0_dp)
    call require_close("J_K hat maximum", quadrature_map(4), 0.75_dp)
    call require_close("mapped J_K physical", quadrature_map(6), 1.5_dp)
    call require_close("action Jacobian weight", quadrature_map(7), 0.45_dp)
    call require_close("paired normalized weight", quadrature_map(8), 0.32_dp)

    call evaluate_neort_polynomial_cell_enclosure(-10.0_dp, 1.0_dp, -2.0_dp, &
        0.5_dp, 0.0_dp, 1.0_dp, 0.5_dp, polynomial_enclosure(1), &
        polynomial_enclosure(2), polynomial_enclosure(3))
    call require_close("polynomial tail bound", polynomial_enclosure(1), &
        1.09375_dp)
    call require_close("negative-c0 lower bound", polynomial_enclosure(2), &
        -11.09375_dp)
    call require_close("negative-c0 upper bound", polynomial_enclosure(3), &
        -8.90625_dp)

    write (*, '(A)') 'test_gc_full_fow_symbolic_kernel OK'
    call pass_test

contains

    subroutine require_close(label, actual, expected)
        character(*), intent(in) :: label
        real(dp), intent(in) :: actual, expected
        real(dp) :: scale

        scale = max(1.0_dp, max(abs(actual), abs(expected)))
        if (abs(actual - expected) > 2.0e-12_dp*scale) then
            write (*, '(A,2(1X,ES24.16))') trim(label), actual, expected
            error stop 'full-FOW symbolic oracle failed'
        end if
    end subroutine require_close

end program test_gc_full_fow_symbolic_kernel
