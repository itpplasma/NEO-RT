module neort_gc_nonlocal_transport_integral
    !! Physical paired (H0,J_K) outer integral for the full-orbit kernel.
    !!
    !! The provider supplies transformed nodes and their *paired* weights.
    !! This layer never constructs a Cartesian J_K domain and never calls a
    !! local eta evaluator.  A failed root, tangent, wall, or class
    !! certification fails the whole accepted result closed.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_gc_nonlocal_resonance_integral, only: &
        integrate_gc_nonlocal_resonance
    use neort_gc_nonlocal_resonance_types, only: &
        GC_NONLOCAL_CALLBACK_FAILURE, GC_NONLOCAL_CAPACITY, &
        GC_NONLOCAL_FORCE_CONTRACT, &
        GC_NONLOCAL_INVALID_INPUT, GC_NONLOCAL_MAX_FORCE_VALUES, &
        GC_NONLOCAL_MAX_ROOT_ITERATIONS, GC_NONLOCAL_MAX_ROOTS, &
        GC_NONLOCAL_MAX_SCAN_INTERVALS, GC_NONLOCAL_NONFINITE, &
        GC_NONLOCAL_PARTIAL, GC_NONLOCAL_SUCCESS, &
        gc_nonlocal_orbit_sample_t, gc_nonlocal_resonance_options_t, &
        gc_nonlocal_resonance_result_t
    use neort_gc_nonlocal_transport_types, only: &
        GC_NONLOCAL_CLASS_COUNT, GC_NONLOCAL_TRANSPORT_MAX_H0_NODES, &
        GC_NONLOCAL_TRANSPORT_MAX_JPERP_NODES, &
        GC_NONLOCAL_TRANSPORT_MAX_TOTAL_NODES, gc_nonlocal_component_t, &
        gc_nonlocal_transport_options_t, gc_nonlocal_transport_provider_t, &
        gc_nonlocal_transport_quadrature_t, &
        gc_nonlocal_transport_reference_t, gc_nonlocal_transport_result_t
    use neort_full_fow_torque_assembly_symbolic, only: &
        evaluate_neort_full_fow_torque_assembly

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

        type(gc_nonlocal_transport_result_t) :: base_result, refined_result
        type(gc_nonlocal_transport_options_t) :: refined_options

        result = gc_nonlocal_transport_result_t()
        call integrate_gc_nonlocal_transport_order(provider, harmonic_m, harmonic_n, &
            options, options%h0_order, options%jk_order, base_result, status)
        if (status /= GC_NONLOCAL_SUCCESS) then
            result = base_result
            return
        end if
        if (.not. options%require_converged) then
            result = base_result
            return
        end if
        refined_options = options
        refined_options%h0_order = 2*options%h0_order
        refined_options%jk_order = 2*options%jk_order
        if (refined_options%h0_order > refined_options%max_h0_nodes .or. &
                refined_options%jk_order > refined_options%max_jperp_nodes .or. &
                refined_options%h0_order*refined_options%jk_order > &
                refined_options%max_total_nodes) then
            result = base_result
            result%certified = .false.
            result%status = GC_NONLOCAL_CAPACITY
            result%failure_status = GC_NONLOCAL_CAPACITY
            status = GC_NONLOCAL_CAPACITY
            return
        end if
        call integrate_gc_nonlocal_transport_order(provider, harmonic_m, harmonic_n, &
            refined_options, refined_options%h0_order, refined_options%jk_order, &
            refined_result, status)
        if (status /= GC_NONLOCAL_SUCCESS) then
            result = refined_result
            return
        end if
        if (.not. result_converged(base_result, refined_result, options)) then
            result = refined_result
            result%certified = .false.
            result%status = GC_NONLOCAL_PARTIAL
            result%failure_status = GC_NONLOCAL_PARTIAL
            status = GC_NONLOCAL_PARTIAL
            return
        end if
        result = refined_result
        result%quadrature%converged = .true.
        status = GC_NONLOCAL_SUCCESS
    end subroutine integrate_gc_nonlocal_transport

    subroutine integrate_gc_nonlocal_transport_order(provider, harmonic_m, harmonic_n, &
            options, h0_order, jk_order, result, status)
        class(gc_nonlocal_transport_provider_t), intent(inout) :: provider
        integer, intent(in) :: harmonic_m, harmonic_n
        type(gc_nonlocal_transport_options_t), intent(in) :: options
        integer, intent(in) :: h0_order, jk_order
        type(gc_nonlocal_transport_result_t), intent(out) :: result
        integer, intent(out) :: status

        type(gc_nonlocal_transport_reference_t) :: reference
        type(gc_nonlocal_transport_quadrature_t) :: quadrature
        type(gc_nonlocal_component_t), allocatable :: components(:)
        type(gc_nonlocal_resonance_result_t) :: node_result
        real(dp) :: accepted_sum(GC_NONLOCAL_MAX_FORCE_VALUES)
        real(dp) :: node_value(GC_NONLOCAL_MAX_FORCE_VALUES)
        real(dp) :: class_node(GC_NONLOCAL_CLASS_COUNT, &
            GC_NONLOCAL_MAX_FORCE_VALUES)
        real(dp) :: node_weight
        integer :: node_index, node_status, provider_status
        integer :: local_status, nforce, class_kind, root_index
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

        call provider%set_harmonic(harmonic_m, harmonic_n, provider_status)
        if (provider_status /= GC_NONLOCAL_SUCCESS) then
            call fail_result(result, provider_status, status)
            return
        end if
        call provider%get_section_reference(reference, provider_status)
        if (provider_status /= GC_NONLOCAL_SUCCESS) then
            call fail_result(result, provider_status, status)
            return
        end if
        local_status = validate_reference(reference)
        if (local_status /= GC_NONLOCAL_SUCCESS) then
            call fail_result(result, local_status, status)
            return
        end if
        call provider%get_quadrature(h0_order, jk_order, quadrature, provider_status)
        if (provider_status /= GC_NONLOCAL_SUCCESS) then
            call fail_result(result, provider_status, status)
            return
        end if
        local_status = validate_quadrature(quadrature, options)
        if (local_status /= GC_NONLOCAL_SUCCESS) then
            call fail_result(result, local_status, status)
            return
        end if

        nforce = options%resonance_options%force_count
        call initialize_result(result, reference, quadrature, nforce)
        do node_index = 1, quadrature%n_nodes
            if (failed) exit
            result%weighted_nodes = result%weighted_nodes + 1
            if (quadrature%weight(node_index) == 0.0_dp) then
                ! Exact H=H_min/J_K,max=0 nodes have zero measure.  They are
                ! certified without invoking an orbit/class callback.
                result%node_status(node_index) = GC_NONLOCAL_SUCCESS
                result%node_certified(node_index) = .true.
                result%certified_nodes = result%certified_nodes + 1
                cycle
            end if
            call provider%get_components(reference, quadrature%h0(node_index), &
                quadrature%j_k(node_index), components, provider_status)
            if (provider_status /= GC_NONLOCAL_SUCCESS) then
                call mark_node_failure(result, node_index, provider_status)
                failed = .true.
                exit
            end if
            if (.not. allocated(components)) then
                call mark_node_failure(result, node_index, GC_NONLOCAL_INVALID_INPUT)
                failed = .true.
                exit
            end if
            if (size(components) == 0) then
                ! A certified over-enclosed action node can have no allowed
                ! Poincare component.  Its physical contribution is exactly
                ! zero; it is not a callback failure.
                result%node_status(node_index) = GC_NONLOCAL_SUCCESS
                result%node_certified(node_index) = .true.
                result%certified_nodes = result%certified_nodes + 1
                cycle
            end if
            result%ncomponents = result%ncomponents + size(components)
            node_result = gc_nonlocal_resonance_result_t()
            call integrate_gc_nonlocal_resonance(transport_evaluate, &
                quadrature%h0(node_index), quadrature%j_k(node_index), harmonic_m, &
                harmonic_n, components, options%resonance_options, node_result, &
                node_status)
            result%nroots = result%nroots + node_result%nroots
            if (node_result%nforce < nforce) then
                call mark_node_failure(result, node_index, GC_NONLOCAL_FORCE_CONTRACT)
                failed = .true.
                exit
            end if
            node_weight = quadrature%weight(node_index)
            call scale_node_contribution(node_result, node_weight, nforce, &
                node_value, local_status)
            if (local_status /= GC_NONLOCAL_SUCCESS) then
                call mark_node_failure(result, node_index, local_status)
                failed = .true.
                exit
            end if
            result%node_contribution(1:nforce, node_index) = &
                node_value(1:nforce)
            result%diagnostic_contribution(1:nforce) = &
                result%diagnostic_contribution(1:nforce) + node_value(1:nforce)
            if (node_status /= GC_NONLOCAL_SUCCESS .or. &
                    .not. node_result%certified) then
                local_status = node_status
                if (local_status == GC_NONLOCAL_SUCCESS) then
                    local_status = GC_NONLOCAL_PARTIAL
                end if
                call mark_node_failure(result, node_index, local_status)
                failed = .true.
                exit
            end if

            class_node = 0.0_dp
            do root_index = 1, node_result%nroots
                call provider%get_class_kind(reference, quadrature%h0(node_index), &
                    quadrature%j_k(node_index), node_result%root_x(root_index), &
                    node_result%root_sigma(root_index), &
                    node_result%root_component_id(root_index), class_kind, &
                    provider_status)
                if (provider_status /= GC_NONLOCAL_SUCCESS) then
                    call mark_node_failure(result, node_index, provider_status)
                    failed = .true.
                    exit
                end if
                if (class_kind < 1 .or. class_kind > GC_NONLOCAL_CLASS_COUNT) then
                    call mark_node_failure(result, node_index, GC_NONLOCAL_INVALID_INPUT)
                    failed = .true.
                    exit
                end if
                if (.not. allocated(node_result%root_class_kind)) then
                    call mark_node_failure(result, node_index, &
                        GC_NONLOCAL_INVALID_INPUT)
                    failed = .true.
                    exit
                end if
                node_result%root_class_kind(root_index) = class_kind
                class_node(class_kind, 1:nforce) = &
                    class_node(class_kind, 1:nforce) + &
                    node_result%root_force_contribution(1:nforce, root_index)
            end do
            if (failed) exit
            class_node = node_weight*class_node
            if (.not. class_sum_matches(class_node, node_value, nforce)) then
                call mark_node_failure(result, node_index, GC_NONLOCAL_INVALID_INPUT)
                failed = .true.
                exit
            end if
            result%class_contribution(:, 1:nforce) = &
                result%class_contribution(:, 1:nforce) + class_node(:, 1:nforce)
            result%certified_nodes = result%certified_nodes + 1
            result%node_status(node_index) = GC_NONLOCAL_SUCCESS
            result%node_certified(node_index) = .true.
            accepted_sum(1:nforce) = accepted_sum(1:nforce) + node_value(1:nforce)
        end do

        if (failed) then
            result%contribution = 0.0_dp
            result%class_contribution = 0.0_dp
            result%certified = .false.
            status = result%status
            return
        end if
        if (.not. all(ieee_is_finite(accepted_sum(1:nforce)))) then
            call fail_result(result, GC_NONLOCAL_NONFINITE, status)
            return
        end if
        if (.not. class_sum_matches(result%class_contribution, accepted_sum, &
                nforce)) then
            call fail_result(result, GC_NONLOCAL_INVALID_INPUT, status)
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

            real(dp) :: force_values(nforce), outer_measure_factor
            integer :: local_callback_status
            real(dp) :: torque_value
            integer :: force_index

            sample = gc_nonlocal_orbit_sample_t()
            call provider%evaluate_orbit(reference, h0, jperp, x, sigma, &
                component_id, sample, local_callback_status)
            if (local_callback_status /= GC_NONLOCAL_SUCCESS) then
                callback_status = local_callback_status
                return
            end if
            if (sample%status /= 0) then
                callback_status = GC_NONLOCAL_SUCCESS
                return
            end if
            if (sample%nforce /= 0) then
                callback_status = GC_NONLOCAL_FORCE_CONTRACT
                return
            end if
            call provider%evaluate_profiles(reference, h0, jperp, x, sigma, &
                component_id, sample, nforce, force_values, local_callback_status)
            if (local_callback_status /= GC_NONLOCAL_SUCCESS) then
                callback_status = local_callback_status
                return
            end if
            call provider%evaluate_outer_measure_factor(reference, h0, jperp, x, &
                sigma, component_id, sample, outer_measure_factor, &
                local_callback_status)
            if (local_callback_status /= GC_NONLOCAL_SUCCESS) then
                callback_status = local_callback_status
                return
            end if
            if (.not. all(ieee_is_finite(force_values))) then
                callback_status = GC_NONLOCAL_NONFINITE
                return
            end if
            if (.not. ieee_is_finite(outer_measure_factor)) then
                callback_status = GC_NONLOCAL_NONFINITE
                return
            end if
            ! Eq. 10 supplies n^2.  outer_measure_factor is the signed native
            ! W_outer prefactor; force_values retain the signed component too.
            do force_index = 1, nforce
                call evaluate_neort_full_fow_torque_assembly( &
                    real(harmonic_n, dp), outer_measure_factor, &
                    force_values(force_index), torque_value)
                sample%thermodynamic_force(force_index) = torque_value
            end do
            sample%nforce = nforce
            callback_status = GC_NONLOCAL_SUCCESS
        end subroutine transport_evaluate

    end subroutine integrate_gc_nonlocal_transport_order

    integer function validate_reference(reference) result(status)
        type(gc_nonlocal_transport_reference_t), intent(in) :: reference

        status = GC_NONLOCAL_INVALID_INPUT
        if (reference%section_id <= 0) return
        if (len_trim(reference%section_coordinate) == 0) return
        if (len_trim(reference%section_units) == 0) return
        if (.not. all(ieee_is_finite(reference%section_position))) return
        if (.not. ieee_is_finite(reference%section_flux)) return
        if (.not. all(ieee_is_finite(reference%force_scale))) return
        if (any(reference%force_scale <= 0.0_dp)) return
        if (abs(reference%p_zeta_orientation) /= 1) return
        if (reference%frequency_semantics /= &
                1) return
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
        if (options%h0_order < 2 .or. options%h0_order > &
                GC_NONLOCAL_TRANSPORT_MAX_H0_NODES) return
        if (options%jk_order < 2 .or. options%jk_order > &
                GC_NONLOCAL_TRANSPORT_MAX_JPERP_NODES) return
        if (options%max_h0_nodes < options%h0_order) return
        if (options%max_jperp_nodes < options%jk_order) return
        if (options%max_total_nodes < options%h0_order*options%jk_order) return
        if (options%max_total_nodes > GC_NONLOCAL_TRANSPORT_MAX_TOTAL_NODES) return
        if (.not. ieee_is_finite(options%quadrature_relative_tolerance)) return
        if (.not. ieee_is_finite(options%quadrature_absolute_tolerance)) return
        if (options%quadrature_relative_tolerance <= 0.0_dp) return
        if (options%quadrature_absolute_tolerance <= 0.0_dp) return
        if (.not. valid_resonance_options(options%resonance_options)) return
        status = GC_NONLOCAL_SUCCESS
    end function validate_transport_request

    pure logical function valid_resonance_options(options)
        type(gc_nonlocal_resonance_options_t), intent(in) :: options

        valid_resonance_options = .false.
        if (options%scan_intervals < 1 .or. options%scan_intervals > &
                GC_NONLOCAL_MAX_SCAN_INTERVALS) return
        if (options%max_root_iterations < 1 .or. options%max_root_iterations > &
                GC_NONLOCAL_MAX_ROOT_ITERATIONS) return
        if (options%max_roots < 1 .or. options%max_roots > GC_NONLOCAL_MAX_ROOTS) &
            return
        if (options%force_count < 1 .or. options%force_count > &
                GC_NONLOCAL_MAX_FORCE_VALUES) return
        if (.not. ieee_is_finite(options%residual_tolerance)) return
        if (.not. ieee_is_finite(options%x_tolerance)) return
        if (.not. ieee_is_finite(options%derivative_tolerance)) return
        if (options%residual_tolerance <= 0.0_dp) return
        if (options%x_tolerance <= 0.0_dp) return
        if (options%derivative_tolerance <= 0.0_dp) return
        valid_resonance_options = .true.
    end function valid_resonance_options

    integer function validate_quadrature(quadrature, options) result(status)
        type(gc_nonlocal_transport_quadrature_t), intent(in) :: quadrature
        type(gc_nonlocal_transport_options_t), intent(in) :: options

        status = GC_NONLOCAL_INVALID_INPUT
        if (.not. quadrature%paired_domain) return
        if (.not. quadrature%domain_certified) return
        if (quadrature%h0_order < options%h0_order) return
        if (quadrature%jk_order < options%jk_order) return
        if (quadrature%n_nodes < 4) return
        if (quadrature%n_nodes > options%max_total_nodes) return
        if (.not. ieee_is_finite(quadrature%h0_min)) return
        if (.not. ieee_is_finite(quadrature%h0_scale)) return
        if (quadrature%h0_scale <= 0.0_dp) return
        if (.not. allocated(quadrature%h0)) return
        if (.not. allocated(quadrature%j_k)) return
        if (.not. allocated(quadrature%weight)) return
        if (.not. allocated(quadrature%j_k_upper_bound)) return
        if (size(quadrature%h0) /= quadrature%n_nodes) return
        if (size(quadrature%j_k) /= quadrature%n_nodes) return
        if (size(quadrature%weight) /= quadrature%n_nodes) return
        if (size(quadrature%j_k_upper_bound) /= quadrature%n_nodes) return
        if (.not. all(ieee_is_finite(quadrature%h0))) return
        if (.not. all(ieee_is_finite(quadrature%j_k))) return
        if (.not. all(ieee_is_finite(quadrature%weight))) return
        if (.not. all(ieee_is_finite(quadrature%j_k_upper_bound))) return
        if (any(quadrature%weight <= 0.0_dp)) return
        if (any(quadrature%h0 < quadrature%h0_min)) return
        if (any(quadrature%j_k_upper_bound <= 0.0_dp)) return
        if (any(quadrature%j_k < 0.0_dp)) return
        if (any(quadrature%j_k > quadrature%j_k_upper_bound)) return
        status = GC_NONLOCAL_SUCCESS
    end function validate_quadrature

    subroutine initialize_result(result, reference, quadrature, nforce)
        type(gc_nonlocal_transport_result_t), intent(inout) :: result
        type(gc_nonlocal_transport_reference_t), intent(in) :: reference
        type(gc_nonlocal_transport_quadrature_t), intent(in) :: quadrature
        integer, intent(in) :: nforce

        result%reference = reference
        result%quadrature = quadrature
        result%nforce = nforce
        result%nh0 = quadrature%h0_order
        result%njperp = quadrature%jk_order
        allocate(result%h0_nodes(quadrature%n_nodes), &
            result%jperp_nodes(quadrature%n_nodes), &
            result%node_weights(quadrature%n_nodes), &
            result%node_status(quadrature%n_nodes), &
            result%node_certified(quadrature%n_nodes), &
            result%node_contribution(nforce, quadrature%n_nodes))
        result%h0_nodes = quadrature%h0
        result%jperp_nodes = quadrature%j_k
        result%node_weights = quadrature%weight
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
            node_value = 0.0_dp
            status = GC_NONLOCAL_NONFINITE
            return
        end if
        status = GC_NONLOCAL_SUCCESS
    end subroutine scale_node_contribution

    logical function class_sum_matches(class_values, total, nforce)
        real(dp), intent(in) :: class_values(:, :)
        real(dp), intent(in) :: total(:)
        integer, intent(in) :: nforce
        real(dp) :: class_total(GC_NONLOCAL_MAX_FORCE_VALUES)
        real(dp) :: scale

        class_total = 0.0_dp
        class_total(1:nforce) = sum(class_values(:, 1:nforce), dim=1)
        scale = max(tiny(1.0_dp), maxval(abs(total(1:nforce))))
        scale = max(scale, maxval(abs(class_total(1:nforce))))
        class_sum_matches = maxval(abs(class_total(1:nforce) - &
            total(1:nforce))) <= 2.0e-11_dp*scale
    end function class_sum_matches

    subroutine mark_node_failure(result, node_index, status)
        type(gc_nonlocal_transport_result_t), intent(inout) :: result
        integer, intent(in) :: node_index, status

        result%status = status
        result%failure_status = status
        result%unresolved_nodes = result%unresolved_nodes + 1
        result%node_status(node_index) = status
        result%node_certified(node_index) = .false.
        if (result%failed_node == 0) result%failed_node = node_index
    end subroutine mark_node_failure

    subroutine fail_result(result, failure_status, status)
        type(gc_nonlocal_transport_result_t), intent(inout) :: result
        integer, intent(in) :: failure_status
        integer, intent(out) :: status

        result%status = failure_status
        result%failure_status = failure_status
        result%certified = .false.
        result%contribution = 0.0_dp
        result%class_contribution = 0.0_dp
        status = failure_status
    end subroutine fail_result

    logical function result_converged(coarse, refined, options)
        type(gc_nonlocal_transport_result_t), intent(in) :: coarse, refined
        type(gc_nonlocal_transport_options_t), intent(in) :: options
        real(dp) :: scale
        real(dp) :: absolute_scale
        integer :: force_index, class_index

        result_converged = .false.
        if (.not. coarse%certified .or. .not. refined%certified) return
        if (.not. ieee_is_finite(options%quadrature_absolute_tolerance)) return
        if (.not. ieee_is_finite(options%quadrature_relative_tolerance)) return
        if (options%quadrature_absolute_tolerance <= 0.0_dp) return
        if (options%quadrature_relative_tolerance <= 0.0_dp) return
        if (.not. all(ieee_is_finite(coarse%contribution))) return
        if (.not. all(ieee_is_finite(refined%contribution))) return
        if (.not. all(ieee_is_finite(coarse%class_contribution))) return
        if (.not. all(ieee_is_finite(refined%class_contribution))) return
        do force_index = 1, GC_NONLOCAL_MAX_FORCE_VALUES
            scale = max(abs(coarse%contribution(force_index)), &
                abs(refined%contribution(force_index)))
            absolute_scale = options%quadrature_absolute_tolerance * &
                refined%reference%force_scale(force_index)
            if (abs(refined%contribution(force_index) - &
                    coarse%contribution(force_index)) > &
                    absolute_scale + &
                    options%quadrature_relative_tolerance*scale) return
        end do
        do class_index = 1, GC_NONLOCAL_CLASS_COUNT
            do force_index = 1, GC_NONLOCAL_MAX_FORCE_VALUES
                scale = max(abs(coarse%class_contribution(class_index, force_index)), &
                    abs(refined%class_contribution(class_index, force_index)))
                absolute_scale = options%quadrature_absolute_tolerance * &
                    refined%reference%force_scale(force_index)
                if (abs(refined%class_contribution(class_index, force_index) - &
                        coarse%class_contribution(class_index, force_index)) > &
                        absolute_scale + &
                        options%quadrature_relative_tolerance*scale) return
            end do
        end do
        result_converged = .true.
    end function result_converged

end module neort_gc_nonlocal_transport_integral
