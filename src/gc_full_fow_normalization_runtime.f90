module neort_gc_full_fow_normalization_runtime
    !! Runtime contract for the generated direct full-FOW normalizations.
    !!
    !! This module is deliberately an adapter, not a second symbolic layer.
    !! It validates the domain of every call, invokes the generated kernels,
    !! and validates their returned values.  The expressions for scales,
    !! actions, quadrature Jacobians, Eq. 17, and polynomial enclosures live
    !! only in the Fortsym-generated modules.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_full_fow_action_symbolic, only: &
        evaluate_neort_action_normalization
    use neort_full_fow_eq17_outer_symbolic, only: &
        evaluate_neort_eq17_outer_factor
    use neort_full_fow_eq17_symbolic, only: evaluate_neort_eq17_force
    use neort_full_fow_normalization_symbolic, only: &
        evaluate_neort_full_fow_normalization
    use neort_full_fow_quadrature_map_symbolic, only: &
        evaluate_neort_full_fow_quadrature_map
    use neort_full_fow_refinement_symbolic, only: &
        evaluate_neort_full_fow_refinement_scales
    use neort_polynomial_cell_enclosure_symbolic, only: &
        evaluate_neort_polynomial_cell_enclosure
    implicit none
    private

    integer, parameter, public :: GC_FULL_FOW_NORMALIZATION_SUCCESS = 0
    integer, parameter, public :: GC_FULL_FOW_NORMALIZATION_INVALID_INPUT = 1
    integer, parameter, public :: GC_FULL_FOW_NORMALIZATION_NONFINITE_INPUT = 2
    integer, parameter, public :: GC_FULL_FOW_NORMALIZATION_DOMAIN_ERROR = 3
    integer, parameter, public :: GC_FULL_FOW_NORMALIZATION_NOT_INITIALIZED = 4
    integer, parameter, public :: &
        GC_FULL_FOW_NORMALIZATION_CERTIFICATE_REQUIRED = 5
    integer, parameter, public :: &
        GC_FULL_FOW_NORMALIZATION_REFERENCE_MISMATCH = 6
    integer, parameter, public :: GC_FULL_FOW_NORMALIZATION_KERNEL_NONFINITE = 7
    integer, parameter, public :: &
        GC_FULL_FOW_NORMALIZATION_OUTPUT_INVALID = 8

    type, public :: gc_full_fow_reference_scales_t
        !! Physical inputs and the generated reference scales used by samples.
        logical :: initialized = .false.
        real(dp) :: mass = 0.0_dp
        real(dp) :: charge = 0.0_dp
        real(dp) :: c_light = 0.0_dp
        real(dp) :: e_ref = 0.0_dp
        real(dp) :: reference_velocity = 0.0_dp
        real(dp) :: phi_eff = 0.0_dp
        real(dp) :: jk_scale = 0.0_dp
    end type gc_full_fow_reference_scales_t

    type, public :: gc_full_fow_physical_sample_t
        !! One physical orbit sample before dimensionless normalization.
        real(dp) :: h_physical = 0.0_dp
        real(dp) :: jk_physical = 0.0_dp
        real(dp) :: psi_star_physical = 0.0_dp
        real(dp) :: dpsi_star_dx_physical = 0.0_dp
        real(dp) :: tau_b_physical = 0.0_dp
        real(dp) :: omega_b_physical = 0.0_dp
        real(dp) :: omega_phi_physical = 0.0_dp
        real(dp) :: domega_b_dx_physical = 0.0_dp
        real(dp) :: domega_phi_dx_physical = 0.0_dp
        complex(dp) :: hm_physical = (0.0_dp, 0.0_dp)
    end type gc_full_fow_physical_sample_t

    type, public :: gc_full_fow_normalized_sample_t
        !! Generated dimensionless sample, including the two reference scales.
        real(dp) :: phi_eff = 0.0_dp
        real(dp) :: jk_scale = 0.0_dp
        real(dp) :: h_hat = 0.0_dp
        real(dp) :: jk_hat = 0.0_dp
        real(dp) :: psi_star_hat = 0.0_dp
        real(dp) :: dpsi_star_hat_dx = 0.0_dp
        real(dp) :: tau_b_hat = 0.0_dp
        real(dp) :: omega_b_hat = 0.0_dp
        real(dp) :: omega_phi_hat = 0.0_dp
        real(dp) :: domega_b_hat_dx = 0.0_dp
        real(dp) :: domega_phi_hat_dx = 0.0_dp
        complex(dp) :: hm_hat = (0.0_dp, 0.0_dp)
    end type gc_full_fow_normalized_sample_t

    type, public :: gc_full_fow_phase_space_bound_certificate_t
        !! Certified lower bounds supplied by the equilibrium/domain layer.
        !! qphi_min is the energy q*Phi_min, not Phi_min itself.
        logical :: qphi_min_certified = .false.
        logical :: bmod_min_certified = .false.
        real(dp) :: qphi_min = 0.0_dp
        real(dp) :: bmod_min = 0.0_dp
        character(len=96) :: method = ''
        character(len=96) :: certificate_id = ''
    end type gc_full_fow_phase_space_bound_certificate_t

    type, public :: gc_full_fow_jk_envelope_t
        !! Positive-action envelope returned by the generated action kernel.
        logical :: certified = .false.
        real(dp) :: h_physical = 0.0_dp
        real(dp) :: qphi_min = 0.0_dp
        real(dp) :: bmod_min = 0.0_dp
        real(dp) :: omega_c_min = 0.0_dp
        real(dp) :: jk_candidate_physical = 0.0_dp
        real(dp) :: jk_max_physical = 0.0_dp
        real(dp) :: mass = 0.0_dp
        real(dp) :: charge = 0.0_dp
        real(dp) :: c_light = 0.0_dp
        character(len=96) :: method = ''
        character(len=96) :: certificate_id = ''
    end type gc_full_fow_jk_envelope_t

    type, public :: gc_full_fow_paired_quadrature_t
        !! One paired energy/action map and its two normalized Jacobians.
        real(dp) :: h_hat = 0.0_dp
        real(dp) :: h_physical = 0.0_dp
        real(dp) :: energy_normalized_weight = 0.0_dp
        real(dp) :: jk_hat_max = 0.0_dp
        real(dp) :: jk_hat = 0.0_dp
        real(dp) :: jk_physical = 0.0_dp
        real(dp) :: action_normalized_weight = 0.0_dp
        real(dp) :: paired_weight = 0.0_dp
    end type gc_full_fow_paired_quadrature_t

    type, public :: gc_full_fow_energy_quadrature_t
        !! Energy-only projection of the generated paired-domain map.
        real(dp) :: h_hat = 0.0_dp
        real(dp) :: h_physical = 0.0_dp
        real(dp) :: normalized_weight = 0.0_dp
    end type gc_full_fow_energy_quadrature_t

    type, public :: gc_full_fow_eq17_result_t
        !! Signed Eq. 17 outputs.  No sign is removed by this adapter.
        real(dp) :: qphi_energy = 0.0_dp
        real(dp) :: force_bracket = 0.0_dp
        !! Phi_eff is owned exactly once by outer_factor.
        real(dp) :: outer_factor = 0.0_dp
    end type gc_full_fow_eq17_result_t

    type, public :: gc_full_fow_degree5_enclosure_t
        !! Certified interval supplied by the generated degree-five kernel.
        real(dp) :: tail_absolute_bound = 0.0_dp
        real(dp) :: lower_bound = 0.0_dp
        real(dp) :: upper_bound = 0.0_dp
    end type gc_full_fow_degree5_enclosure_t

    public :: initialize_gc_full_fow_reference
    public :: evaluate_gc_full_fow_canonical_flux
    public :: normalize_gc_full_fow_sample
    public :: evaluate_gc_full_fow_jk_envelope
    public :: evaluate_gc_full_fow_scaled_magnitude
    public :: map_gc_full_fow_energy_quadrature
    public :: map_gc_full_fow_paired_quadrature
    public :: evaluate_gc_full_fow_eq17
    public :: evaluate_gc_full_fow_degree5_cell_enclosure

contains

    subroutine evaluate_gc_full_fow_canonical_flux(reference, p_phi, psi_star, &
            dpsi_star_dp_phi, status, message)
        !! Canonical conversion projected from the generated action kernel.
        !! The neutral action inputs do not enter either requested output.
        type(gc_full_fow_reference_scales_t), intent(in) :: reference
        real(dp), intent(in) :: p_phi
        real(dp), intent(out) :: psi_star, dpsi_star_dp_phi
        integer, intent(out) :: status
        character(len=*), intent(out) :: message

        real(dp) :: values(11)

        psi_star = 0.0_dp
        dpsi_star_dp_phi = 0.0_dp
        call set_result(status, message, GC_FULL_FOW_NORMALIZATION_INVALID_INPUT, &
            '')
        if (.not. reference%initialized) then
            call set_result(status, message, &
                GC_FULL_FOW_NORMALIZATION_NOT_INITIALIZED, &
                'reference scales are not initialized')
            return
        end if
        if (.not. valid_reference(reference)) then
            call set_result(status, message, &
                GC_FULL_FOW_NORMALIZATION_INVALID_INPUT, &
                'reference scale record is invalid')
            return
        end if
        if (.not. ieee_is_finite(p_phi)) then
            call set_result(status, message, &
                GC_FULL_FOW_NORMALIZATION_NONFINITE_INPUT, &
                'canonical momentum must be finite')
            return
        end if
        call evaluate_neort_action_normalization(reference%mass, &
            reference%charge, reference%c_light, 0.0_dp, 1.0_dp, 0.0_dp, &
            0.0_dp, p_phi, 1.0_dp, values(1), values(2), values(3), values(4), &
            values(5), values(6), values(7), values(8), values(9), values(10), &
            values(11))
        if (.not. all(ieee_is_finite(values))) then
            call set_result(status, message, &
                GC_FULL_FOW_NORMALIZATION_KERNEL_NONFINITE, &
                'generated canonical conversion returned a nonfinite value')
            return
        end if
        psi_star = values(10)
        dpsi_star_dp_phi = values(11)
        call set_result(status, message, GC_FULL_FOW_NORMALIZATION_SUCCESS, &
            'canonical flux evaluated')
    end subroutine evaluate_gc_full_fow_canonical_flux

    subroutine evaluate_gc_full_fow_scaled_magnitude(value, reference_scale, &
            scaled_magnitude, status, message)
        !! Scalar projection of the generated refinement-scale family.
        real(dp), intent(in) :: value, reference_scale
        real(dp), intent(out) :: scaled_magnitude
        integer, intent(out) :: status
        character(len=*), intent(out) :: message

        real(dp) :: scales(14)

        scaled_magnitude = 0.0_dp
        call set_result(status, message, GC_FULL_FOW_NORMALIZATION_INVALID_INPUT, &
            '')
        if (.not. all(ieee_is_finite([value, reference_scale]))) then
            call set_result(status, message, &
                GC_FULL_FOW_NORMALIZATION_NONFINITE_INPUT, &
                'refinement value and reference scale must be finite')
            return
        end if
        if (reference_scale <= 0.0_dp) then
            call set_result(status, message, &
                GC_FULL_FOW_NORMALIZATION_DOMAIN_ERROR, &
                'refinement reference scale must be positive')
            return
        end if
        call evaluate_neort_full_fow_refinement_scales(value, 0.0_dp, 0.0_dp, &
            0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
            0.0_dp, 0.0_dp, 0.0_dp, reference_scale, 1.0_dp, 1.0_dp, &
            1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, &
            1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, scales(1), scales(2), scales(3), &
            scales(4), scales(5), scales(6), scales(7), scales(8), scales(9), &
            scales(10), scales(11), scales(12), scales(13), scales(14))
        if (.not. all(ieee_is_finite(scales))) then
            call set_result(status, message, &
                GC_FULL_FOW_NORMALIZATION_KERNEL_NONFINITE, &
                'generated refinement scaling returned a nonfinite value')
            return
        end if
        scaled_magnitude = scales(1)
        call set_result(status, message, GC_FULL_FOW_NORMALIZATION_SUCCESS, &
            'refinement magnitude scaled')
    end subroutine evaluate_gc_full_fow_scaled_magnitude

    subroutine initialize_gc_full_fow_reference(mass, charge, c_light, e_ref, &
            reference_velocity, reference, status, message)
        real(dp), intent(in) :: mass, charge, c_light, e_ref
        real(dp), intent(in) :: reference_velocity
        type(gc_full_fow_reference_scales_t), intent(out) :: reference
        integer, intent(out) :: status
        character(len=*), intent(out) :: message

        real(dp) :: normalization(13)

        reference = gc_full_fow_reference_scales_t()
        call set_result(status, message, GC_FULL_FOW_NORMALIZATION_INVALID_INPUT, &
            '')
        if (.not. all(ieee_is_finite([mass, charge, c_light, e_ref, &
                reference_velocity]))) then
            call set_result(status, message, &
                GC_FULL_FOW_NORMALIZATION_NONFINITE_INPUT, &
                'reference inputs must be finite')
            return
        end if
        if (mass <= 0.0_dp .or. c_light <= 0.0_dp .or. e_ref <= 0.0_dp .or. &
                reference_velocity <= 0.0_dp) then
            call set_result(status, message, &
                GC_FULL_FOW_NORMALIZATION_DOMAIN_ERROR, &
                'mass, c_light, Eref, and reference velocity must be positive')
            return
        end if
        if (charge == 0.0_dp) then
            call set_result(status, message, &
                GC_FULL_FOW_NORMALIZATION_DOMAIN_ERROR, &
                'charge must be nonzero')
            return
        end if

        call call_normalization_kernel(mass, charge, c_light, e_ref, &
            reference_velocity, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
            0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, normalization)
        if (.not. all(ieee_is_finite(normalization))) then
            call set_result(status, message, &
                GC_FULL_FOW_NORMALIZATION_KERNEL_NONFINITE, &
                'generated normalization returned a nonfinite reference')
            return
        end if
        if (normalization(1) <= 0.0_dp .or. normalization(2) <= 0.0_dp) then
            call set_result(status, message, &
                GC_FULL_FOW_NORMALIZATION_OUTPUT_INVALID, &
                'generated reference scales are not positive')
            return
        end if

        reference%mass = mass
        reference%charge = charge
        reference%c_light = c_light
        reference%e_ref = e_ref
        reference%reference_velocity = reference_velocity
        reference%phi_eff = normalization(1)
        reference%jk_scale = normalization(2)
        reference%initialized = .true.
        call set_result(status, message, GC_FULL_FOW_NORMALIZATION_SUCCESS, &
            'reference scales initialized')
    end subroutine initialize_gc_full_fow_reference

    subroutine normalize_gc_full_fow_sample(reference, physical, normalized, &
            status, message)
        type(gc_full_fow_reference_scales_t), intent(in) :: reference
        type(gc_full_fow_physical_sample_t), intent(in) :: physical
        type(gc_full_fow_normalized_sample_t), intent(out) :: normalized
        integer, intent(out) :: status
        character(len=*), intent(out) :: message

        real(dp) :: values(13)

        normalized = gc_full_fow_normalized_sample_t()
        call set_result(status, message, GC_FULL_FOW_NORMALIZATION_INVALID_INPUT, &
            '')
        if (.not. reference%initialized) then
            call set_result(status, message, &
                GC_FULL_FOW_NORMALIZATION_NOT_INITIALIZED, &
                'reference scales are not initialized')
            return
        end if
        if (.not. valid_reference(reference)) then
            call set_result(status, message, &
                GC_FULL_FOW_NORMALIZATION_INVALID_INPUT, &
                'reference scale record is invalid')
            return
        end if
        if (.not. all(ieee_is_finite([physical%h_physical, &
                physical%jk_physical, physical%psi_star_physical, &
                physical%dpsi_star_dx_physical, physical%tau_b_physical, &
                physical%omega_b_physical, physical%omega_phi_physical, &
                physical%domega_b_dx_physical, &
                physical%domega_phi_dx_physical, real(physical%hm_physical, dp), &
                aimag(physical%hm_physical)]))) then
            call set_result(status, message, &
                GC_FULL_FOW_NORMALIZATION_NONFINITE_INPUT, &
                'physical sample must be finite')
            return
        end if
        if (physical%jk_physical < 0.0_dp .or. &
                physical%tau_b_physical <= 0.0_dp) then
            call set_result(status, message, &
                GC_FULL_FOW_NORMALIZATION_DOMAIN_ERROR, &
                'physical action must be nonnegative and bounce time positive')
            return
        end if

        call call_normalization_kernel(reference%mass, reference%charge, &
            reference%c_light, reference%e_ref, reference%reference_velocity, &
            physical%h_physical, physical%jk_physical, &
            physical%psi_star_physical, physical%dpsi_star_dx_physical, &
            physical%tau_b_physical, physical%omega_b_physical, &
            physical%omega_phi_physical, physical%domega_b_dx_physical, &
            physical%domega_phi_dx_physical, real(physical%hm_physical, dp), &
            aimag(physical%hm_physical), values)
        call validate_normalization_output(values, reference, status, message)
        if (status /= GC_FULL_FOW_NORMALIZATION_SUCCESS) return

        normalized%phi_eff = values(1)
        normalized%jk_scale = values(2)
        normalized%h_hat = values(3)
        normalized%jk_hat = values(4)
        normalized%psi_star_hat = values(5)
        normalized%dpsi_star_hat_dx = values(6)
        normalized%tau_b_hat = values(7)
        normalized%omega_b_hat = values(8)
        normalized%omega_phi_hat = values(9)
        normalized%domega_b_hat_dx = values(10)
        normalized%domega_phi_hat_dx = values(11)
        normalized%hm_hat = cmplx(values(12), values(13), kind=dp)
        call set_result(status, message, GC_FULL_FOW_NORMALIZATION_SUCCESS, &
            'sample normalized')
    end subroutine normalize_gc_full_fow_sample

    subroutine evaluate_gc_full_fow_jk_envelope(reference, h_physical, &
            certificate, envelope, status, message)
        type(gc_full_fow_reference_scales_t), intent(in) :: reference
        real(dp), intent(in) :: h_physical
        type(gc_full_fow_phase_space_bound_certificate_t), intent(in) :: &
            certificate
        type(gc_full_fow_jk_envelope_t), intent(out) :: envelope
        integer, intent(out) :: status
        character(len=*), intent(out) :: message

        real(dp) :: action_values(11)
        real(dp) :: electrostatic_potential_min

        envelope = gc_full_fow_jk_envelope_t()
        call set_result(status, message, GC_FULL_FOW_NORMALIZATION_INVALID_INPUT, &
            '')
        if (.not. reference%initialized) then
            call set_result(status, message, &
                GC_FULL_FOW_NORMALIZATION_NOT_INITIALIZED, &
                'reference scales are not initialized')
            return
        end if
        if (.not. valid_reference(reference)) then
            call set_result(status, message, &
                GC_FULL_FOW_NORMALIZATION_INVALID_INPUT, &
                'reference scale record is invalid')
            return
        end if
        if (.not. certificate%qphi_min_certified .or. &
                .not. certificate%bmod_min_certified) then
            call set_result(status, message, &
                GC_FULL_FOW_NORMALIZATION_CERTIFICATE_REQUIRED, &
                'both qPhi_min and B_min require certification')
            return
        end if
        if (len_trim(certificate%method) == 0) then
            call set_result(status, message, &
                GC_FULL_FOW_NORMALIZATION_CERTIFICATE_REQUIRED, &
                'phase-space bound method is required')
            return
        end if
        if (.not. ieee_is_finite(h_physical) .or. &
                .not. all(ieee_is_finite([certificate%qphi_min, &
                certificate%bmod_min]))) then
            call set_result(status, message, &
                GC_FULL_FOW_NORMALIZATION_NONFINITE_INPUT, &
                'Hamiltonian and phase-space bounds must be finite')
            return
        end if
        if (certificate%bmod_min <= 0.0_dp) then
            call set_result(status, message, &
                GC_FULL_FOW_NORMALIZATION_DOMAIN_ERROR, &
                'certified B_min must be positive')
            return
        end if
        electrostatic_potential_min = certificate%qphi_min/reference%charge
        if (.not. ieee_is_finite(electrostatic_potential_min)) then
            call set_result(status, message, &
                GC_FULL_FOW_NORMALIZATION_DOMAIN_ERROR, &
                'qPhi_min cannot be converted to a finite potential')
            return
        end if

        ! The generated action interface accepts Phi, while the certificate
        ! records q*Phi.  This is the sole unit-interface conversion here;
        ! the action/envelope expression remains generated.
        call evaluate_neort_action_normalization(reference%mass, &
            reference%charge, reference%c_light, 0.0_dp, &
            certificate%bmod_min, h_physical, &
            electrostatic_potential_min, 0.0_dp, 1.0_dp, &
            action_values(1), action_values(2), action_values(3), &
            action_values(4), action_values(5), action_values(6), &
            action_values(7), action_values(8), action_values(9), &
            action_values(10), action_values(11))
        if (.not. all(ieee_is_finite(action_values))) then
            call set_result(status, message, &
                GC_FULL_FOW_NORMALIZATION_KERNEL_NONFINITE, &
                'generated action envelope returned a nonfinite value')
            return
        end if
        if (action_values(2) <= 0.0_dp .or. action_values(5) < 0.0_dp) then
            call set_result(status, message, &
                GC_FULL_FOW_NORMALIZATION_OUTPUT_INVALID, &
                'generated action envelope is outside its domain')
            return
        end if

        envelope%h_physical = h_physical
        envelope%qphi_min = certificate%qphi_min
        envelope%bmod_min = certificate%bmod_min
        envelope%omega_c_min = action_values(2)
        envelope%jk_candidate_physical = action_values(4)
        envelope%jk_max_physical = action_values(5)
        envelope%mass = reference%mass
        envelope%charge = reference%charge
        envelope%c_light = reference%c_light
        envelope%method = certificate%method
        envelope%certificate_id = certificate%certificate_id
        envelope%certified = .true.
        call set_result(status, message, GC_FULL_FOW_NORMALIZATION_SUCCESS, &
            'J_K envelope evaluated')
    end subroutine evaluate_gc_full_fow_jk_envelope

    subroutine map_gc_full_fow_energy_quadrature(reference, qphi_min, &
            energy_unit_node, energy_unit_weight, mapped, status, message)
        !! Project the energy part of the generated paired map before the
        !! energy-dependent J_K envelope is known.  Zero action inputs are
        !! exact neutral inputs to that same generated kernel; no map algebra
        !! is duplicated here.
        type(gc_full_fow_reference_scales_t), intent(in) :: reference
        real(dp), intent(in) :: qphi_min, energy_unit_node, energy_unit_weight
        type(gc_full_fow_energy_quadrature_t), intent(out) :: mapped
        integer, intent(out) :: status
        character(len=*), intent(out) :: message

        real(dp) :: values(8)

        mapped = gc_full_fow_energy_quadrature_t()
        call set_result(status, message, GC_FULL_FOW_NORMALIZATION_INVALID_INPUT, &
            '')
        if (.not. reference%initialized) then
            call set_result(status, message, &
                GC_FULL_FOW_NORMALIZATION_NOT_INITIALIZED, &
                'reference scales are not initialized')
            return
        end if
        if (.not. valid_reference(reference)) then
            call set_result(status, message, &
                GC_FULL_FOW_NORMALIZATION_INVALID_INPUT, &
                'reference scale record is invalid')
            return
        end if
        if (.not. all(ieee_is_finite([qphi_min, energy_unit_node, &
                energy_unit_weight]))) then
            call set_result(status, message, &
                GC_FULL_FOW_NORMALIZATION_NONFINITE_INPUT, &
                'energy quadrature inputs must be finite')
            return
        end if
        if (energy_unit_node < 0.0_dp .or. energy_unit_node >= 1.0_dp .or. &
                energy_unit_weight < 0.0_dp) then
            call set_result(status, message, &
                GC_FULL_FOW_NORMALIZATION_DOMAIN_ERROR, &
                'energy quadrature node or weight is outside its unit domain')
            return
        end if

        call evaluate_neort_full_fow_quadrature_map(reference%mass, &
            reference%charge, reference%c_light, reference%e_ref, &
            energy_unit_node, energy_unit_weight, 0.0_dp, 0.0_dp, qphi_min, &
            0.0_dp, values(1), values(2), values(3), values(4), values(5), &
            values(6), values(7), values(8))
        if (.not. all(ieee_is_finite(values))) then
            call set_result(status, message, &
                GC_FULL_FOW_NORMALIZATION_KERNEL_NONFINITE, &
                'generated energy quadrature map returned a nonfinite value')
            return
        end if
        if (values(3) < 0.0_dp) then
            call set_result(status, message, &
                GC_FULL_FOW_NORMALIZATION_OUTPUT_INVALID, &
                'generated energy quadrature map returned a negative measure')
            return
        end if

        mapped%h_hat = values(1)
        mapped%h_physical = values(2)
        mapped%normalized_weight = values(3)
        call set_result(status, message, GC_FULL_FOW_NORMALIZATION_SUCCESS, &
            'energy quadrature mapped')
    end subroutine map_gc_full_fow_energy_quadrature

    subroutine map_gc_full_fow_paired_quadrature(reference, envelope, &
            energy_unit_node, energy_unit_weight, action_unit_node, &
            action_unit_weight, mapped, status, message)
        type(gc_full_fow_reference_scales_t), intent(in) :: reference
        type(gc_full_fow_jk_envelope_t), intent(in) :: envelope
        real(dp), intent(in) :: energy_unit_node, energy_unit_weight
        real(dp), intent(in) :: action_unit_node, action_unit_weight
        type(gc_full_fow_paired_quadrature_t), intent(out) :: mapped
        integer, intent(out) :: status
        character(len=*), intent(out) :: message

        real(dp) :: values(8)

        mapped = gc_full_fow_paired_quadrature_t()
        call set_result(status, message, GC_FULL_FOW_NORMALIZATION_INVALID_INPUT, &
            '')
        if (.not. reference%initialized) then
            call set_result(status, message, &
                GC_FULL_FOW_NORMALIZATION_NOT_INITIALIZED, &
                'reference scales are not initialized')
            return
        end if
        if (.not. valid_reference(reference)) then
            call set_result(status, message, &
                GC_FULL_FOW_NORMALIZATION_INVALID_INPUT, &
                'reference scale record is invalid')
            return
        end if
        if (.not. envelope%certified) then
            call set_result(status, message, &
                GC_FULL_FOW_NORMALIZATION_CERTIFICATE_REQUIRED, &
                'a certified J_K envelope is required')
            return
        end if
        if (.not. same_reference_identity(reference, envelope)) then
            call set_result(status, message, &
                GC_FULL_FOW_NORMALIZATION_REFERENCE_MISMATCH, &
                'J_K envelope belongs to another reference scale set')
            return
        end if
        if (.not. all(ieee_is_finite([energy_unit_node, energy_unit_weight, &
                action_unit_node, action_unit_weight, envelope%qphi_min, &
                envelope%jk_max_physical]))) then
            call set_result(status, message, &
                GC_FULL_FOW_NORMALIZATION_NONFINITE_INPUT, &
                'quadrature inputs must be finite')
            return
        end if
        if (energy_unit_node < 0.0_dp .or. energy_unit_node >= 1.0_dp .or. &
                action_unit_node < 0.0_dp .or. action_unit_node > 1.0_dp .or. &
                energy_unit_weight < 0.0_dp .or. action_unit_weight < 0.0_dp .or. &
                envelope%jk_max_physical < 0.0_dp) then
            call set_result(status, message, &
                GC_FULL_FOW_NORMALIZATION_DOMAIN_ERROR, &
                'quadrature nodes or weights are outside [0,1] domains')
            return
        end if

        call evaluate_neort_full_fow_quadrature_map(reference%mass, &
            reference%charge, reference%c_light, reference%e_ref, &
            energy_unit_node, energy_unit_weight, action_unit_node, &
            action_unit_weight, envelope%qphi_min, envelope%jk_max_physical, &
            values(1), values(2), values(3), values(4), values(5), values(6), &
            values(7), values(8))
        if (.not. all(ieee_is_finite(values))) then
            call set_result(status, message, &
                GC_FULL_FOW_NORMALIZATION_KERNEL_NONFINITE, &
                'generated quadrature map returned a nonfinite value')
            return
        end if
        if (values(4) < 0.0_dp .or. values(5) < 0.0_dp .or. &
                values(6) < 0.0_dp .or. values(3) < 0.0_dp .or. &
                values(7) < 0.0_dp .or. values(8) < 0.0_dp) then
            call set_result(status, message, &
                GC_FULL_FOW_NORMALIZATION_OUTPUT_INVALID, &
                'generated quadrature map returned a negative measure')
            return
        end if

        mapped%h_hat = values(1)
        mapped%h_physical = values(2)
        mapped%energy_normalized_weight = values(3)
        mapped%jk_hat_max = values(4)
        mapped%jk_hat = values(5)
        mapped%jk_physical = values(6)
        mapped%action_normalized_weight = values(7)
        mapped%paired_weight = values(8)
        call set_result(status, message, GC_FULL_FOW_NORMALIZATION_SUCCESS, &
            'paired quadrature mapped')
    end subroutine map_gc_full_fow_paired_quadrature

    subroutine evaluate_gc_full_fow_eq17(reference, h_physical, &
            electrostatic_potential, temperature, a1, a2, n0, residence, &
            result, status, message)
        type(gc_full_fow_reference_scales_t), intent(in) :: reference
        real(dp), intent(in) :: h_physical, electrostatic_potential
        real(dp), intent(in) :: temperature, a1, a2, n0, residence
        type(gc_full_fow_eq17_result_t), intent(out) :: result
        integer, intent(out) :: status
        character(len=*), intent(out) :: message

        result = gc_full_fow_eq17_result_t()
        call set_result(status, message, GC_FULL_FOW_NORMALIZATION_INVALID_INPUT, &
            '')
        if (.not. reference%initialized) then
            call set_result(status, message, &
                GC_FULL_FOW_NORMALIZATION_NOT_INITIALIZED, &
                'reference scales are not initialized')
            return
        end if
        if (.not. valid_reference(reference)) then
            call set_result(status, message, &
                GC_FULL_FOW_NORMALIZATION_INVALID_INPUT, &
                'reference scale record is invalid')
            return
        end if
        if (.not. all(ieee_is_finite([h_physical, electrostatic_potential, &
                temperature, a1, a2, n0, residence]))) then
            call set_result(status, message, &
                GC_FULL_FOW_NORMALIZATION_NONFINITE_INPUT, &
                'Eq. 17 inputs must be finite')
            return
        end if
        if (temperature <= 0.0_dp .or. n0 < 0.0_dp .or. residence < 0.0_dp) then
            call set_result(status, message, &
                GC_FULL_FOW_NORMALIZATION_DOMAIN_ERROR, &
                'temperature must be positive and density/residence nonnegative')
            return
        end if

        call evaluate_neort_eq17_force(h_physical, reference%charge, &
            electrostatic_potential, temperature, a1, a2, &
            result%qphi_energy, result%force_bracket)
        call evaluate_neort_eq17_outer_factor(reference%e_ref, &
            reference%phi_eff, n0, temperature, reference%charge, &
            electrostatic_potential, h_physical, residence, result%outer_factor)
        if (.not. all(ieee_is_finite([result%qphi_energy, &
                result%force_bracket, result%outer_factor]))) then
            call set_result(status, message, &
                GC_FULL_FOW_NORMALIZATION_KERNEL_NONFINITE, &
                'generated Eq. 17 kernels returned a nonfinite value')
            return
        end if
        call set_result(status, message, GC_FULL_FOW_NORMALIZATION_SUCCESS, &
            'Eq. 17 force and outer factor evaluated')
    end subroutine evaluate_gc_full_fow_eq17

    subroutine evaluate_gc_full_fow_degree5_cell_enclosure(coefficients, &
            cell_width, enclosure, status, message)
        real(dp), intent(in) :: coefficients(6), cell_width
        type(gc_full_fow_degree5_enclosure_t), intent(out) :: enclosure
        integer, intent(out) :: status
        character(len=*), intent(out) :: message

        real(dp) :: values(3)

        enclosure = gc_full_fow_degree5_enclosure_t()
        call set_result(status, message, GC_FULL_FOW_NORMALIZATION_INVALID_INPUT, &
            '')
        if (.not. all(ieee_is_finite(coefficients)) .or. &
                .not. ieee_is_finite(cell_width)) then
            call set_result(status, message, &
                GC_FULL_FOW_NORMALIZATION_NONFINITE_INPUT, &
                'polynomial coefficients and cell width must be finite')
            return
        end if
        if (cell_width < 0.0_dp) then
            call set_result(status, message, &
                GC_FULL_FOW_NORMALIZATION_DOMAIN_ERROR, &
                'cell width must be nonnegative')
            return
        end if

        call evaluate_neort_polynomial_cell_enclosure(coefficients(1), &
            coefficients(2), coefficients(3), coefficients(4), coefficients(5), &
            coefficients(6), cell_width, values(1), values(2), values(3))
        if (.not. all(ieee_is_finite(values))) then
            call set_result(status, message, &
                GC_FULL_FOW_NORMALIZATION_KERNEL_NONFINITE, &
                'generated polynomial enclosure returned a nonfinite value')
            return
        end if
        if (values(1) < 0.0_dp .or. values(2) > values(3)) then
            call set_result(status, message, &
                GC_FULL_FOW_NORMALIZATION_OUTPUT_INVALID, &
                'generated polynomial enclosure is not ordered')
            return
        end if

        enclosure%tail_absolute_bound = values(1)
        enclosure%lower_bound = values(2)
        enclosure%upper_bound = values(3)
        call set_result(status, message, GC_FULL_FOW_NORMALIZATION_SUCCESS, &
            'degree-five cell enclosed')
    end subroutine evaluate_gc_full_fow_degree5_cell_enclosure

    subroutine call_normalization_kernel(mass, charge, c_light, e_ref, &
            reference_velocity, h_physical, jk_physical, psi_star_physical, &
            dpsi_star_dx_physical, tau_b_physical, omega_b_physical, &
            omega_phi_physical, domega_b_dx_physical, domega_phi_dx_physical, &
            hm_real_physical, hm_imag_physical, values)
        real(dp), intent(in) :: mass, charge, c_light, e_ref
        real(dp), intent(in) :: reference_velocity, h_physical, jk_physical
        real(dp), intent(in) :: psi_star_physical, dpsi_star_dx_physical
        real(dp), intent(in) :: tau_b_physical, omega_b_physical
        real(dp), intent(in) :: omega_phi_physical, domega_b_dx_physical
        real(dp), intent(in) :: domega_phi_dx_physical
        real(dp), intent(in) :: hm_real_physical, hm_imag_physical
        real(dp), intent(out) :: values(13)

        call evaluate_neort_full_fow_normalization(mass, charge, c_light, &
            e_ref, reference_velocity, h_physical, jk_physical, &
            psi_star_physical, dpsi_star_dx_physical, tau_b_physical, &
            omega_b_physical, omega_phi_physical, domega_b_dx_physical, &
            domega_phi_dx_physical, hm_real_physical, hm_imag_physical, &
            values(1), values(2), values(3), values(4), values(5), values(6), &
            values(7), values(8), values(9), values(10), values(11), values(12), &
            values(13))
    end subroutine call_normalization_kernel

    subroutine validate_normalization_output(values, reference, status, message)
        real(dp), intent(in) :: values(13)
        type(gc_full_fow_reference_scales_t), intent(in) :: reference
        integer, intent(out) :: status
        character(len=*), intent(out) :: message

        call set_result(status, message, GC_FULL_FOW_NORMALIZATION_SUCCESS, '')
        if (.not. all(ieee_is_finite(values))) then
            call set_result(status, message, &
                GC_FULL_FOW_NORMALIZATION_KERNEL_NONFINITE, &
                'generated normalization returned a nonfinite sample')
            return
        end if
        if (values(1) <= 0.0_dp .or. values(2) <= 0.0_dp .or. &
                values(4) < 0.0_dp .or. values(7) <= 0.0_dp) then
            call set_result(status, message, &
                GC_FULL_FOW_NORMALIZATION_OUTPUT_INVALID, &
                'generated normalization returned an invalid scale or measure')
            return
        end if
        if (.not. close_value(values(1), reference%phi_eff) .or. &
                .not. close_value(values(2), reference%jk_scale)) then
            call set_result(status, message, &
                GC_FULL_FOW_NORMALIZATION_REFERENCE_MISMATCH, &
                'generated sample scales disagree with reference scales')
        end if
    end subroutine validate_normalization_output

    logical function valid_reference(reference)
        type(gc_full_fow_reference_scales_t), intent(in) :: reference

        valid_reference = .false.
        if (.not. reference%initialized) return
        if (.not. all(ieee_is_finite([reference%mass, reference%charge, &
                reference%c_light, reference%e_ref, &
                reference%reference_velocity, reference%phi_eff, &
                reference%jk_scale]))) return
        if (reference%mass <= 0.0_dp .or. reference%c_light <= 0.0_dp .or. &
                reference%e_ref <= 0.0_dp .or. &
                reference%reference_velocity <= 0.0_dp .or. &
                reference%phi_eff <= 0.0_dp .or. reference%jk_scale <= 0.0_dp) &
            return
        if (reference%charge == 0.0_dp) return
        valid_reference = .true.
    end function valid_reference

    logical function same_reference_identity(reference, envelope)
        type(gc_full_fow_reference_scales_t), intent(in) :: reference
        type(gc_full_fow_jk_envelope_t), intent(in) :: envelope

        same_reference_identity = .true.
        if (.not. close_value(reference%mass, envelope%mass)) &
            same_reference_identity = .false.
        if (.not. close_value(reference%charge, envelope%charge)) &
            same_reference_identity = .false.
        if (.not. close_value(reference%c_light, envelope%c_light)) &
            same_reference_identity = .false.
    end function same_reference_identity

    logical function close_value(first, second)
        real(dp), intent(in) :: first, second
        real(dp) :: scale

        scale = max(1.0_dp, max(abs(first), abs(second)))
        close_value = abs(first-second) <= 64.0_dp*epsilon(1.0_dp)*scale
    end function close_value

    subroutine set_result(status, message, new_status, new_message)
        integer, intent(out) :: status
        character(len=*), intent(out) :: message
        integer, intent(in) :: new_status
        character(len=*), intent(in) :: new_message

        status = new_status
        message = ''
        message = new_message
    end subroutine set_result

end module neort_gc_full_fow_normalization_runtime
