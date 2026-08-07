module potato_symbolic_kernel_mod
    !! Cohesive semantic wrappers around the split Fortsym kernels.
    !!
    !! These procedures do no physics algebra.  They provide stable POTATO
    !! names and keep the generated module interfaces out of legacy callers;
    !! domain validation and status propagation remain in handwritten code.
    use, intrinsic :: iso_fortran_env, only : real64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use potato_gap_error_generated_mod, only : potato_gap_error_kernel_generated => &
        potato_gap_error_kernel
    use potato_hm_eq4_generated_mod, only : potato_hm_eq4_kernel_generated => &
        potato_hm_eq4_kernel
    use potato_jperp_domain_generated_mod, only : potato_jperp_domain_kernel_generated => &
        potato_jperp_domain_kernel
    use potato_limiting_generated_mod, only : potato_limiting_kernel_generated => &
        potato_limiting_kernel
    use potato_root_jacobian_generated_mod, only : potato_root_jacobian_kernel_generated => &
        potato_root_jacobian_kernel
    use potato_resonance_harmonic_generated_mod, only : &
        potato_resonance_harmonic_kernel_generated => &
        potato_resonance_harmonic_kernel
    use potato_resonance_extent_envelope_generated_mod, only : &
        potato_resonance_extent_envelope_kernel_generated => &
        potato_resonance_extent_envelope_kernel
    use potato_frequency_reduction_generated_mod, only : &
        potato_frequency_reduction_kernel_generated => &
        potato_frequency_reduction_kernel
    use potato_resonance_torque_generated_mod, only : potato_resonance_torque_kernel_generated => &
        potato_resonance_torque_kernel
    implicit none
    private

    public :: potato_jperp_kernel
    public :: potato_hm_kernel
    public :: potato_root_jacobian_kernel
    public :: potato_resonance_harmonic_kernel
    public :: potato_resonance_extent_envelope_kernel
    public :: potato_frequency_reduction_kernel
    public :: potato_resonance_torque_kernel
    public :: potato_gap_contribution_kernel
    public :: potato_limiting_kernel
    public :: potato_limiting_kernel_checked
    public :: potato_limiting_success, potato_limiting_invalid_reference, &
        potato_limiting_invalid_distance, potato_limiting_invalid_output

    integer, parameter :: potato_limiting_success=0
    integer, parameter :: potato_limiting_invalid_reference=1
    integer, parameter :: potato_limiting_invalid_distance=2
    integer, parameter :: potato_limiting_invalid_output=3

contains

    subroutine potato_jperp_kernel(energy_H,qPhi,magnetic_field_B, &
        qPhi_prime,magnetic_field_B_prime,Jperp_candidate,Jperp_positive_bound, &
        dJperp_dR)
        real(real64), intent(in) :: energy_H,qPhi,magnetic_field_B
        real(real64), intent(in) :: qPhi_prime,magnetic_field_B_prime
        real(real64), intent(out) :: Jperp_candidate,Jperp_positive_bound,dJperp_dR

        call potato_jperp_domain_kernel_generated(energy_H,qPhi, &
            magnetic_field_B,qPhi_prime,magnetic_field_B_prime,Jperp_candidate, &
            Jperp_positive_bound,dJperp_dR)
    end subroutine potato_jperp_kernel

    subroutine potato_hm_kernel(energy_H,qPhi,magnetic_field_B,Jperp, &
        deltaB_m_real,deltaB_m_imag,mode_phase,H_m_real,H_m_imag, &
        H_m_squared,hamiltonian_coefficient)
        real(real64), intent(in) :: energy_H,qPhi,magnetic_field_B,Jperp
        real(real64), intent(in) :: deltaB_m_real,deltaB_m_imag,mode_phase
        real(real64), intent(out) :: H_m_real,H_m_imag,H_m_squared
        real(real64), intent(out) :: hamiltonian_coefficient

        call potato_hm_eq4_kernel_generated(energy_H,qPhi,magnetic_field_B, &
            Jperp,deltaB_m_real,deltaB_m_imag,mode_phase,H_m_real,H_m_imag, &
            H_m_squared,hamiltonian_coefficient)
    end subroutine potato_hm_kernel

    subroutine potato_root_jacobian_kernel(dpsiast_dR,dresonance_dR, &
        root_jacobian)
        real(real64), intent(in) :: dpsiast_dR,dresonance_dR
        real(real64), intent(out) :: root_jacobian

        call potato_root_jacobian_kernel_generated(dpsiast_dR,dresonance_dR, &
            root_jacobian)
    end subroutine potato_root_jacobian_kernel

    subroutine potato_resonance_harmonic_kernel(mode_m,mode_n,delta_phi, &
        target_delphi,search_extent,resonance_g,no_root_margin)
        real(real64), intent(in) :: mode_m,mode_n,delta_phi
        real(real64), intent(out) :: target_delphi,search_extent
        real(real64), intent(out) :: resonance_g,no_root_margin

        call potato_resonance_harmonic_kernel_generated(mode_m,mode_n,delta_phi, &
            target_delphi,search_extent,resonance_g,no_root_margin)
    end subroutine potato_resonance_harmonic_kernel

    subroutine potato_resonance_extent_envelope_kernel(current_extent_envelope, &
        harmonic_extent,extent_envelope)
        real(real64), intent(in) :: current_extent_envelope,harmonic_extent
        real(real64), intent(out) :: extent_envelope

        call potato_resonance_extent_envelope_kernel_generated( &
            current_extent_envelope,harmonic_extent,extent_envelope)
    end subroutine potato_resonance_extent_envelope_kernel

    subroutine potato_frequency_reduction_kernel(mode_n, &
        resonance_q,resonance_q_prime,bounce_time,bounce_time_prime, &
        frequency_resonance,frequency_derivative,frequency_derivative_at_root, &
        frequency_root_remainder,frequency_delta_weight_from_n2, &
        collapsed_n_tau2_weight)
        real(real64), intent(in) :: mode_n,resonance_q
        real(real64), intent(in) :: resonance_q_prime,bounce_time
        real(real64), intent(in) :: bounce_time_prime
        real(real64), intent(out) :: frequency_resonance,frequency_derivative
        real(real64), intent(out) :: frequency_derivative_at_root
        real(real64), intent(out) :: frequency_root_remainder
        real(real64), intent(out) :: frequency_delta_weight_from_n2
        real(real64), intent(out) :: collapsed_n_tau2_weight

        call potato_frequency_reduction_kernel_generated(mode_n, &
            resonance_q,resonance_q_prime,bounce_time,bounce_time_prime, &
            frequency_resonance,frequency_derivative, &
            frequency_derivative_at_root,frequency_root_remainder, &
            frequency_delta_weight_from_n2,collapsed_n_tau2_weight)
    end subroutine potato_frequency_reduction_kernel

    subroutine potato_resonance_torque_kernel(dpsiast_dR,dresonance_dR,mode_n, &
        orbit_H_m_squared,maxwellian_weight,Phi_eff,bounce_time, &
        thermodynamic_force,reference_energy,delta_root_weight, &
        resonance_torque_weight)
        real(real64), intent(in) :: dpsiast_dR,dresonance_dR,mode_n
        real(real64), intent(in) :: orbit_H_m_squared
        real(real64), intent(in) :: maxwellian_weight,Phi_eff,bounce_time
        real(real64), intent(in) :: thermodynamic_force,reference_energy
        real(real64), intent(out) :: delta_root_weight
        real(real64), intent(out) :: resonance_torque_weight

        ! orbit_H_m_squared is the finite orbit-averaged provider output from
        ! pertham, not a second caller-supplied mode magnitude.  Its status
        ! and finiteness are checked before this kernel is reached.
        ! The root Jacobian is derived from the two physical derivatives inside
        ! the generated kernel; callers cannot spoof a precomputed magnitude.
        call potato_resonance_torque_kernel_generated(dpsiast_dR,dresonance_dR, &
            mode_n,orbit_H_m_squared,maxwellian_weight,Phi_eff,bounce_time, &
            thermodynamic_force,reference_energy,delta_root_weight, &
            resonance_torque_weight)
    end subroutine potato_resonance_torque_kernel

    subroutine potato_gap_contribution_kernel(integrand_envelope, &
        Jperp_gap_lo,Jperp_gap_hi,topology_gap_measure, &
        topology_contribution_error_bound)
        real(real64), intent(in) :: integrand_envelope,Jperp_gap_lo,Jperp_gap_hi
        real(real64), intent(out) :: topology_gap_measure
        real(real64), intent(out) :: topology_contribution_error_bound

        call potato_gap_error_kernel_generated(integrand_envelope,Jperp_gap_lo, &
            Jperp_gap_hi,topology_gap_measure, &
            topology_contribution_error_bound)
    end subroutine potato_gap_contribution_kernel

    subroutine potato_limiting_kernel(hessian_H_RR,hessian_H_Rp,hessian_H_pp, &
        regular_tau_value,cut_linear_slope,xpoint_cut_curvature,delta_R, &
        delta_R_ref, &
        hessian_determinant,lambda_local,C_tau,regular_action_offset, &
        regular_action_jacobian,xpoint_action_offset,xpoint_action_jacobian, &
        regular_tau_limit,separatrix_tau_scale,xpoint_tau_scale)
        real(real64), intent(in) :: hessian_H_RR,hessian_H_Rp,hessian_H_pp
        real(real64), intent(in) :: regular_tau_value,cut_linear_slope
        real(real64), intent(in) :: xpoint_cut_curvature,delta_R,delta_R_ref
        real(real64), intent(out) :: hessian_determinant,lambda_local,C_tau
        real(real64), intent(out) :: regular_action_offset
        real(real64), intent(out) :: regular_action_jacobian
        real(real64), intent(out) :: xpoint_action_offset
        real(real64), intent(out) :: xpoint_action_jacobian
        real(real64), intent(out) :: regular_tau_limit
        real(real64), intent(out) :: separatrix_tau_scale,xpoint_tau_scale

        integer :: ierr

        call potato_limiting_kernel_checked(hessian_H_RR,hessian_H_Rp, &
            hessian_H_pp,regular_tau_value,cut_linear_slope, &
            xpoint_cut_curvature,delta_R,delta_R_ref,hessian_determinant, &
            lambda_local,C_tau,regular_action_offset,regular_action_jacobian, &
            xpoint_action_offset,xpoint_action_jacobian,regular_tau_limit, &
            separatrix_tau_scale,xpoint_tau_scale,ierr)
        if(ierr /= potato_limiting_success) then
            error stop 'invalid POTATO limiting-expression inputs or output'
        endif
    end subroutine potato_limiting_kernel

    subroutine potato_limiting_kernel_checked(hessian_H_RR,hessian_H_Rp, &
        hessian_H_pp,regular_tau_value,cut_linear_slope,xpoint_cut_curvature, &
        delta_R,delta_R_ref,hessian_determinant,lambda_local,C_tau, &
        regular_action_offset,regular_action_jacobian,xpoint_action_offset, &
        xpoint_action_jacobian,regular_tau_limit,separatrix_tau_scale, &
        xpoint_tau_scale,ierr)
        real(real64), intent(in) :: hessian_H_RR,hessian_H_Rp,hessian_H_pp
        real(real64), intent(in) :: regular_tau_value,cut_linear_slope
        real(real64), intent(in) :: xpoint_cut_curvature,delta_R,delta_R_ref
        real(real64), intent(out) :: hessian_determinant,lambda_local,C_tau
        real(real64), intent(out) :: regular_action_offset
        real(real64), intent(out) :: regular_action_jacobian
        real(real64), intent(out) :: xpoint_action_offset
        real(real64), intent(out) :: xpoint_action_jacobian
        real(real64), intent(out) :: regular_tau_limit
        real(real64), intent(out) :: separatrix_tau_scale,xpoint_tau_scale
        integer, intent(out) :: ierr

        ierr=potato_limiting_success
        hessian_determinant=0.0_real64
        lambda_local=0.0_real64
        C_tau=0.0_real64
        regular_action_offset=0.0_real64
        regular_action_jacobian=0.0_real64
        xpoint_action_offset=0.0_real64
        xpoint_action_jacobian=0.0_real64
        regular_tau_limit=0.0_real64
        separatrix_tau_scale=0.0_real64
        xpoint_tau_scale=0.0_real64
        if(.not.ieee_is_finite(delta_R_ref) .or. delta_R_ref <= 0.0_real64) then
            ierr=potato_limiting_invalid_reference
            return
        endif
        if(.not.ieee_is_finite(delta_R) .or. abs(delta_R) <= 0.0_real64) then
            ierr=potato_limiting_invalid_distance
            return
        endif
        call potato_limiting_kernel_generated(hessian_H_RR,hessian_H_Rp, &
            hessian_H_pp,regular_tau_value,cut_linear_slope, &
            xpoint_cut_curvature,delta_R,delta_R_ref,hessian_determinant, &
            lambda_local,C_tau,regular_action_offset,regular_action_jacobian, &
            xpoint_action_offset,xpoint_action_jacobian,regular_tau_limit, &
            separatrix_tau_scale,xpoint_tau_scale)
        if(.not.ieee_is_finite(hessian_determinant) .or. &
           .not.ieee_is_finite(lambda_local) .or. lambda_local <= 0.0_real64 .or. &
           .not.ieee_is_finite(C_tau) .or. &
           .not.ieee_is_finite(regular_action_offset) .or. &
           .not.ieee_is_finite(regular_action_jacobian) .or. &
           .not.ieee_is_finite(xpoint_action_offset) .or. &
           .not.ieee_is_finite(xpoint_action_jacobian) .or. &
           .not.ieee_is_finite(regular_tau_limit) .or. &
           .not.ieee_is_finite(separatrix_tau_scale) .or. &
           .not.ieee_is_finite(xpoint_tau_scale)) then
            ierr=potato_limiting_invalid_output
        endif
    end subroutine potato_limiting_kernel_checked

end module potato_symbolic_kernel_mod
