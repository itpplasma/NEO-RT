module neort_gc_nonlocal_transport_integral
    !! Outer Eq. 14--17 full-finite-orbit-width transport quadrature.
    !!
    !! The fixed-(H0,J_perp) root search is delegated to
    !! neort_gc_nonlocal_resonance_integral.  This module owns only the
    !! tensor-product H0,J_perp measure and the strict provider composition.
    !! A node is accepted only if the delegated kernel certifies every class.
    !! If one weighted node is unresolved, accepted contribution is reset to
    !! zero; diagnostic_contribution remains explicitly marked as diagnostic.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_gc_nonlocal_resonance_integral, only: &
        integrate_gc_nonlocal_resonance
    use neort_gc_nonlocal_resonance_types, only: &
        GC_NONLOCAL_CAPACITY, GC_NONLOCAL_CALLBACK_FAILURE, &
        GC_NONLOCAL_DERIVATIVE_MISSING, GC_NONLOCAL_FORCE_CONTRACT, &
        GC_NONLOCAL_INVALID_INPUT, GC_NONLOCAL_MAX_FORCE_VALUES, &
        GC_NONLOCAL_MAX_ROOT_ITERATIONS, GC_NONLOCAL_MAX_ROOTS, &
        GC_NONLOCAL_MAX_SCAN_INTERVALS, GC_NONLOCAL_NONFINITE, &
        GC_NONLOCAL_PARTIAL, GC_NONLOCAL_SAMPLE_VALID, GC_NONLOCAL_SUCCESS, &
        gc_nonlocal_orbit_sample_t, gc_nonlocal_resonance_options_t, &
        gc_nonlocal_resonance_result_t
    use neort_gc_nonlocal_transport_types, only: &
        GC_NONLOCAL_TRANSPORT_MAX_H0_NODES, &
        GC_NONLOCAL_TRANSPORT_MAX_JPERP_NODES, &
        GC_NONLOCAL_TRANSPORT_MAX_TOTAL_NODES, &
        GC_NONLOCAL_COMPLETE_CYCLE_FREQUENCIES, &
        gc_nonlocal_component_t, gc_nonlocal_transport_options_t, &
        gc_nonlocal_transport_provider_t, gc_nonlocal_transport_reference_t, &
        gc_nonlocal_transport_result_t

    implicit none
    private

    public :: integrate_gc_nonlocal_transport

contains

    subroutine integrate_gc_nonlocal_transport(provider, harmonic_m, harmonic_n, &
            options, result, status)
        class(gc_nonlocal_transport_provider_t), intent(inout) :: provider
        integer, intent(in) :: harmonic_m, harmonic_n
        type(gc_nonlocal_transport_options_t), intent(in) :: options
        type(gc_nonlocal_transport_result_t), intent(out) :: result
        integer, intent(out) :: status

        real(dp), allocatable :: h0_nodes(:), h0_weights(:)
        real(dp), allocatable :: jperp_nodes(:), jperp_weights(:)
        type(gc_nonlocal_component_t), allocatable :: components(:)
        type(gc_nonlocal_transport_reference_t) :: reference
        type(gc_nonlocal_resonance_result_t) :: node_result
        real(dp) :: accepted_sum(GC_NONLOCAL_MAX_FORCE_VALUES)
        real(dp) :: node_value(GC_NONLOCAL_MAX_FORCE_VALUES)
        real(dp) :: frequency_harmonic_factor
        integer :: h0_index, jperp_index, node_status, provider_status
        integer :: nforce, local_status, provider_callback_status
        logical :: failed

        result = gc_nonlocal_transport_result_t()
        result%harmonic_m = harmonic_m
        result%harmonic_n = harmonic_n
        accepted_sum = 0.0_dp
        failed = .false.
        status = validate_transport_request(harmonic_m, harmonic_n, options)
        if (status /= GC_NONLOCAL_SUCCESS) then
            result%status = status
            result%failure_status = status
            return
        end if

        call provider%get_section_reference(reference, provider_status)
        if (provider_status /= GC_NONLOCAL_SUCCESS) then
            result%status = provider_status
            result%failure_status = provider_status
            status = provider_status
            return
        end if
        local_status = validate_reference(reference)
        if (local_status /= GC_NONLOCAL_SUCCESS) then
            result%status = local_status
            result%failure_status = local_status
            status = local_status
            return
        end if

        call provider%get_quadrature(h0_nodes, h0_weights, jperp_nodes, &
            jperp_weights, provider_status)
        if (provider_status /= GC_NONLOCAL_SUCCESS) then
            result%status = provider_status
            result%failure_status = provider_status
            status = provider_status
            return
        end if
        local_status = validate_quadrature(h0_nodes, h0_weights, jperp_nodes, &
            jperp_weights, options)
        if (local_status /= GC_NONLOCAL_SUCCESS) then
            result%status = local_status
            result%failure_status = local_status
            status = local_status
            return
        end if

        nforce = options%resonance_options%force_count
        frequency_harmonic_factor = real(harmonic_n, dp)**2
        call initialize_result(result, harmonic_m, harmonic_n, nforce, h0_nodes, &
            h0_weights, jperp_nodes, jperp_weights)
        result%reference = reference

        do jperp_index = 1, result%njperp
            if (failed) exit
            do h0_index = 1, result%nh0
                result%weighted_nodes = result%weighted_nodes + 1
                call provider%get_components(reference, result%h0_nodes(h0_index), &
                    result%jperp_nodes(jperp_index), components, provider_status)
                if (provider_status /= GC_NONLOCAL_SUCCESS) then
                    call mark_node_failure(result, h0_index, jperp_index, &
                        provider_status)
                    failed = .true.
                    exit
                end if
                if (.not. allocated(components)) then
                    call mark_node_failure(result, h0_index, jperp_index, &
                        GC_NONLOCAL_INVALID_INPUT)
                    failed = .true.
                    exit
                end if
                if (size(components) < 1) then
                    call mark_node_failure(result, h0_index, jperp_index, &
                        GC_NONLOCAL_INVALID_INPUT)
                    failed = .true.
                    exit
                end if

                result%ncomponents = result%ncomponents + size(components)
                node_result = gc_nonlocal_resonance_result_t()
                provider_callback_status = GC_NONLOCAL_SUCCESS
                call integrate_gc_nonlocal_resonance(transport_evaluate, &
                    result%h0_nodes(h0_index), result%jperp_nodes(jperp_index), &
                    harmonic_m, harmonic_n, components, &
                    options%resonance_options, node_result, node_status)
                if (provider_callback_status /= GC_NONLOCAL_SUCCESS) then
                    node_status = provider_callback_status
                end if
                result%nroots = result%nroots + node_result%nroots
                if (node_result%nforce < nforce) then
                    node_value = 0.0_dp
                    local_status = GC_NONLOCAL_SUCCESS
                    if (node_status == GC_NONLOCAL_SUCCESS) then
                        local_status = GC_NONLOCAL_FORCE_CONTRACT
                    end if
                else
                    call scale_node_contribution(node_result, &
                        result%h0_weights(h0_index)* &
                        result%jperp_weights(jperp_index), nforce, node_value, &
                        local_status)
                end if
                if (local_status /= GC_NONLOCAL_SUCCESS) then
                    call mark_node_failure(result, h0_index, jperp_index, &
                        local_status)
                    failed = .true.
                    exit
                end if
                result%node_contribution(:, h0_index, jperp_index) = &
                    node_value(1:nforce)
                result%diagnostic_contribution(1:nforce) = &
                    result%diagnostic_contribution(1:nforce) + &
                    node_value(1:nforce)
                if (node_status /= GC_NONLOCAL_SUCCESS) then
                    call mark_node_failure(result, h0_index, jperp_index, &
                        node_status)
                    failed = .true.
                    exit
                end if
                if (.not. node_result%certified) then
                    call mark_node_failure(result, h0_index, jperp_index, &
                        GC_NONLOCAL_PARTIAL)
                    failed = .true.
                    exit
                end if

                result%certified_nodes = result%certified_nodes + 1
                result%node_status(h0_index, jperp_index) = GC_NONLOCAL_SUCCESS
                result%node_certified(h0_index, jperp_index) = .true.
                accepted_sum(1:nforce) = accepted_sum(1:nforce) + &
                    node_value(1:nforce)
            end do
        end do

        if (failed) then
            result%contribution = 0.0_dp
            result%certified = .false.
            status = result%status
            return
        end if

        if (.not. all(ieee_is_finite(accepted_sum(1:nforce)))) then
            result%contribution = 0.0_dp
            result%certified = .false.
            result%status = GC_NONLOCAL_NONFINITE
            result%failure_status = GC_NONLOCAL_NONFINITE
            status = GC_NONLOCAL_NONFINITE
            return
        end if
        result%contribution(1:nforce) = accepted_sum(1:nforce)
        result%certified = .true.
        result%status = GC_NONLOCAL_SUCCESS
        result%failure_status = GC_NONLOCAL_SUCCESS
        status = GC_NONLOCAL_SUCCESS

    contains

        subroutine transport_evaluate(h0, jperp, x, sigma, component_id, sample, &
                callback_status)
            real(dp), intent(in) :: h0, jperp, x
            integer, intent(in) :: sigma, component_id
            type(gc_nonlocal_orbit_sample_t), intent(out) :: sample
            integer, intent(out) :: callback_status

            real(dp) :: force_values(nforce)
            real(dp) :: outer_measure_factor
            integer :: local_callback_status

            sample = gc_nonlocal_orbit_sample_t()
            call provider%evaluate_orbit(reference, h0, jperp, x, sigma, &
                component_id, sample, local_callback_status)
            if (local_callback_status /= GC_NONLOCAL_SUCCESS) then
                provider_callback_status = local_callback_status
                callback_status = local_callback_status
                return
            end if

            ! A wall or unresolved sample is passed intact to the existing
            ! kernel, which records a diagnostic partial and refuses
            ! certification.  Profiles are only meaningful for a valid orbit.
            if (sample%status /= GC_NONLOCAL_SAMPLE_VALID) then
                callback_status = GC_NONLOCAL_SUCCESS
                return
            end if
            if (sample%nforce /= 0) then
                provider_callback_status = GC_NONLOCAL_FORCE_CONTRACT
                callback_status = GC_NONLOCAL_FORCE_CONTRACT
                return
            end if

            call provider%evaluate_profiles(reference, h0, jperp, x, sigma, &
                component_id, sample, nforce, force_values, local_callback_status)
            if (local_callback_status /= GC_NONLOCAL_SUCCESS) then
                provider_callback_status = local_callback_status
                callback_status = local_callback_status
                return
            end if
            if (.not. all(ieee_is_finite(force_values))) then
                provider_callback_status = GC_NONLOCAL_NONFINITE
                callback_status = GC_NONLOCAL_NONFINITE
                return
            end if

            call provider%evaluate_outer_measure_factor(reference, h0, jperp, x, &
                sigma, component_id, sample, outer_measure_factor, &
                local_callback_status)
            if (local_callback_status /= GC_NONLOCAL_SUCCESS) then
                provider_callback_status = local_callback_status
                callback_status = local_callback_status
                return
            end if
            if (.not. ieee_is_finite(outer_measure_factor)) then
                provider_callback_status = GC_NONLOCAL_NONFINITE
                callback_status = GC_NONLOCAL_NONFINITE
                return
            end if

            ! The patched resonance kernel multiplies its force slot by the
            ! frequency-form Eq. 17 root weight, including the sole
            ! |d psi_star/dx| and tau_b/|d r/dx| factors.  The m3^2 conversion
            ! from the phase form is applied here exactly once.  Composing
            ! only the residual signed Eq. 17 factor here keeps the root
            ! finder single-sourced without duplicating a canonical Jacobian.
            sample%thermodynamic_force(1:nforce) = frequency_harmonic_factor &
                *outer_measure_factor*force_values
            sample%nforce = nforce
            callback_status = GC_NONLOCAL_SUCCESS
        end subroutine transport_evaluate

    end subroutine integrate_gc_nonlocal_transport

    integer function validate_reference(reference) result(status)
        type(gc_nonlocal_transport_reference_t), intent(in) :: reference

        status = GC_NONLOCAL_INVALID_INPUT
        if (reference%section_id <= 0) return
        if (len_trim(reference%section_coordinate) == 0) return
        if (len_trim(reference%section_units) == 0) return
        if (.not. all(ieee_is_finite(reference%section_position))) return
        if (.not. ieee_is_finite(reference%section_flux)) return
        if (abs(reference%p_zeta_orientation) /= 1) return
        if (reference%frequency_semantics /= &
            GC_NONLOCAL_COMPLETE_CYCLE_FREQUENCIES) return
        if (.not. reference%hamiltonian_includes_phi) return
        if (.not. reference%fixed) return
        status = GC_NONLOCAL_SUCCESS
    end function validate_reference

    integer function validate_transport_request(harmonic_m, harmonic_n, options) &
            result(status)
        integer, intent(in) :: harmonic_m, harmonic_n
        type(gc_nonlocal_transport_options_t), intent(in) :: options

        status = GC_NONLOCAL_INVALID_INPUT
        if (harmonic_m == 0 .and. harmonic_n == 0) return
        if (harmonic_n == 0) return
        if (options%max_h0_nodes < 1) return
        if (options%max_h0_nodes > GC_NONLOCAL_TRANSPORT_MAX_H0_NODES) return
        if (options%max_jperp_nodes < 1) return
        if (options%max_jperp_nodes > GC_NONLOCAL_TRANSPORT_MAX_JPERP_NODES) return
        if (options%max_total_nodes < 1) return
        if (options%max_total_nodes > GC_NONLOCAL_TRANSPORT_MAX_TOTAL_NODES) return
        if (.not. valid_resonance_options(options%resonance_options)) return
        status = GC_NONLOCAL_SUCCESS
    end function validate_transport_request

    pure logical function valid_resonance_options(options)
        type(gc_nonlocal_resonance_options_t), intent(in) :: options

        valid_resonance_options = .false.
        if (options%scan_intervals < 1) return
        if (options%scan_intervals > GC_NONLOCAL_MAX_SCAN_INTERVALS) return
        if (options%max_root_iterations < 1) return
        if (options%max_root_iterations > GC_NONLOCAL_MAX_ROOT_ITERATIONS) return
        if (options%max_roots < 1) return
        if (options%max_roots > GC_NONLOCAL_MAX_ROOTS) return
        if (options%force_count < 1) return
        if (options%force_count > GC_NONLOCAL_MAX_FORCE_VALUES) return
        if (.not. ieee_is_finite(options%residual_tolerance)) return
        if (.not. ieee_is_finite(options%x_tolerance)) return
        if (.not. ieee_is_finite(options%derivative_tolerance)) return
        if (options%residual_tolerance <= 0.0_dp) return
        if (options%x_tolerance <= 0.0_dp) return
        if (options%derivative_tolerance <= 0.0_dp) return
        valid_resonance_options = .true.
    end function valid_resonance_options

    integer function validate_quadrature(h0_nodes, h0_weights, jperp_nodes, &
            jperp_weights, options) result(status)
        real(dp), allocatable, intent(in) :: h0_nodes(:), h0_weights(:)
        real(dp), allocatable, intent(in) :: jperp_nodes(:), jperp_weights(:)
        type(gc_nonlocal_transport_options_t), intent(in) :: options
        integer :: total_nodes

        status = GC_NONLOCAL_INVALID_INPUT
        if (.not. allocated(h0_nodes)) return
        if (.not. allocated(h0_weights)) return
        if (.not. allocated(jperp_nodes)) return
        if (.not. allocated(jperp_weights)) return
        if (size(h0_nodes) < 1) return
        if (size(jperp_nodes) < 1) return
        if (size(h0_nodes) /= size(h0_weights)) return
        if (size(jperp_nodes) /= size(jperp_weights)) return
        if (size(h0_nodes) > options%max_h0_nodes) return
        if (size(jperp_nodes) > options%max_jperp_nodes) return
        if (size(h0_nodes) > GC_NONLOCAL_TRANSPORT_MAX_H0_NODES) return
        if (size(jperp_nodes) > GC_NONLOCAL_TRANSPORT_MAX_JPERP_NODES) return
        if (size(h0_nodes) > huge(total_nodes)/size(jperp_nodes)) return
        total_nodes = size(h0_nodes)*size(jperp_nodes)
        if (total_nodes > options%max_total_nodes) return
        if (total_nodes > GC_NONLOCAL_TRANSPORT_MAX_TOTAL_NODES) return
        if (.not. all(ieee_is_finite(h0_nodes))) return
        if (.not. all(ieee_is_finite(h0_weights))) return
        if (.not. all(ieee_is_finite(jperp_nodes))) return
        if (.not. all(ieee_is_finite(jperp_weights))) return
        if (any(jperp_nodes < 0.0_dp)) return
        if (any(h0_weights <= 0.0_dp)) return
        if (any(jperp_weights <= 0.0_dp)) return
        status = GC_NONLOCAL_SUCCESS
    end function validate_quadrature

    subroutine initialize_result(result, harmonic_m, harmonic_n, nforce, h0_nodes, &
            h0_weights, jperp_nodes, jperp_weights)
        type(gc_nonlocal_transport_result_t), intent(inout) :: result
        integer, intent(in) :: harmonic_m, harmonic_n, nforce
        real(dp), intent(in) :: h0_nodes(:), h0_weights(:), jperp_nodes(:), &
            jperp_weights(:)

        result%harmonic_m = harmonic_m
        result%harmonic_n = harmonic_n
        result%nforce = nforce
        result%nh0 = size(h0_nodes)
        result%njperp = size(jperp_nodes)
        allocate(result%h0_nodes(result%nh0), result%h0_weights(result%nh0))
        allocate(result%jperp_nodes(result%njperp), &
            result%jperp_weights(result%njperp))
        allocate(result%node_status(result%nh0, result%njperp))
        allocate(result%node_certified(result%nh0, result%njperp))
        allocate(result%node_contribution(nforce, result%nh0, result%njperp))
        result%h0_nodes = h0_nodes
        result%h0_weights = h0_weights
        result%jperp_nodes = jperp_nodes
        result%jperp_weights = jperp_weights
        result%node_status = GC_NONLOCAL_INVALID_INPUT
        result%node_certified = .false.
        result%node_contribution = 0.0_dp
    end subroutine initialize_result

    subroutine scale_node_contribution(node_result, weight, nforce, node_value, &
            status)
        type(gc_nonlocal_resonance_result_t), intent(in) :: node_result
        real(dp), intent(in) :: weight
        integer, intent(in) :: nforce
        real(dp), intent(out) :: node_value(GC_NONLOCAL_MAX_FORCE_VALUES)
        integer, intent(out) :: status

        node_value = 0.0_dp
        status = GC_NONLOCAL_INVALID_INPUT
        if (nforce < 1 .or. nforce > GC_NONLOCAL_MAX_FORCE_VALUES) return
        if (node_result%nforce < nforce) return
        if (.not. ieee_is_finite(weight)) then
            status = GC_NONLOCAL_NONFINITE
            return
        end if
        if (.not. all(ieee_is_finite(node_result%contribution(1:nforce)))) then
            status = GC_NONLOCAL_NONFINITE
            return
        end if
        node_value(1:nforce) = weight*node_result%contribution(1:nforce)
        if (.not. all(ieee_is_finite(node_value(1:nforce)))) then
            status = GC_NONLOCAL_NONFINITE
            node_value = 0.0_dp
            return
        end if
        status = GC_NONLOCAL_SUCCESS
    end subroutine scale_node_contribution

    subroutine mark_node_failure(result, h0_index, jperp_index, status)
        type(gc_nonlocal_transport_result_t), intent(inout) :: result
        integer, intent(in) :: h0_index, jperp_index, status

        result%status = status
        result%failure_status = status
        result%unresolved_nodes = result%unresolved_nodes + 1
        result%node_status(h0_index, jperp_index) = status
        result%node_certified(h0_index, jperp_index) = .false.
        if (result%failed_h0_index == 0) then
            result%failed_h0_index = h0_index
            result%failed_jperp_index = jperp_index
        end if
    end subroutine mark_node_failure

end module neort_gc_nonlocal_transport_integral
