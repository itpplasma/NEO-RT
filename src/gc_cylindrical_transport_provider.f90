module neort_gc_cylindrical_transport_provider
    !! Concrete clean-architecture adapter for cylindrical full-FOW transport.
    !!
    !! The outer transport contract varies (H0,Jperp), while both existing
    !! cylindrical seams are deliberately fixed-invariant.  A node factory
    !! therefore creates a fresh class adapter and a fresh callback context
    !! for every weighted node.  The class adapter is the authority for the
    !! disconnected Poincare classes; the nonlocal callback context is the
    !! authority for complete-cycle dynamics, phase-averaged H_m, and signed
    !! native forces.  Their class identities and intervals are cross-checked
    !! before a node is exposed to the resonance kernel.
    !!
    !! No local eta, electric-frequency correction, sign alignment, or
    !! normalization fit is performed here.  The callback context must report
    !! frequencies from the Hamiltonian which already contains Phi.  The only
    !! canonical conversion accepted by this adapter is
    !! psi_star=(c/q)*p_phi, applied once by the existing callback provider.
    !! W_outer is intentionally a bare signed factor: it excludes tau_b,
    !! n^2, and the sole canonical Jacobian, which belong to the transport and
    !! resonance layers respectively.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use neort_gc_cylindrical_class_adapter, only: &
        GC_CYL_CLASS_SUCCESS, gc_cylindrical_class_adapter_t, &
        gc_cylindrical_class_interval_t, gc_cylindrical_class_launch_t, &
        gc_cylindrical_class_result_t, enumerate_gc_cylindrical_classes, &
        launch_gc_cylindrical_class
    use neort_gc_cylindrical_nonlocal_provider, only: &
        GC_CYL_NONLOCAL_SUCCESS, &
        evaluate_gc_cylindrical_nonlocal_sample, &
        enumerate_gc_cylindrical_nonlocal_components, &
        gc_cylindrical_nonlocal_context_t
    use neort_gc_nonlocal_resonance_types, only: &
        GC_NONLOCAL_CALLBACK_FAILURE, GC_NONLOCAL_COMPONENT_IDENTITY, &
        GC_NONLOCAL_INVALID_INPUT, GC_NONLOCAL_MAX_FORCE_VALUES, &
        GC_NONLOCAL_NONFINITE, GC_NONLOCAL_SAMPLE_INVALID, &
        GC_NONLOCAL_SAMPLE_UNRESOLVED, GC_NONLOCAL_SAMPLE_VALID, &
        GC_NONLOCAL_SAMPLE_WALL, GC_NONLOCAL_SUCCESS, &
        gc_nonlocal_component_t, gc_nonlocal_orbit_sample_t
    use neort_gc_nonlocal_transport_types, only: &
        GC_NONLOCAL_COMPLETE_CYCLE_FREQUENCIES, &
        gc_nonlocal_transport_provider_t, &
        gc_nonlocal_transport_reference_t

    implicit none
    private

    integer, parameter, public :: GC_CYL_TRANSPORT_SUCCESS = &
        GC_NONLOCAL_SUCCESS
    integer, parameter, public :: GC_CYL_TRANSPORT_INVALID_INPUT = &
        GC_NONLOCAL_INVALID_INPUT
    integer, parameter, public :: GC_CYL_TRANSPORT_UNAVAILABLE = &
        GC_NONLOCAL_CALLBACK_FAILURE
    integer, parameter, public :: GC_CYL_TRANSPORT_CLASS_UNCERTIFIED = 201
    integer, parameter, public :: GC_CYL_TRANSPORT_NODE_MISMATCH = 202
    integer, parameter, public :: GC_CYL_TRANSPORT_REFERENCE_MISMATCH = 203
    integer, parameter, public :: GC_CYL_TRANSPORT_NORMALIZATION_ERROR = 204

    abstract interface
        subroutine gc_cylindrical_transport_node_factory_i(h0, jperp, &
                user_data, adapter, context, status)
            import :: dp, gc_cylindrical_class_adapter_t, &
                gc_cylindrical_nonlocal_context_t
            real(dp), intent(in) :: h0, jperp
            class(*), pointer, intent(inout) :: user_data
            type(gc_cylindrical_class_adapter_t), intent(out) :: adapter
            type(gc_cylindrical_nonlocal_context_t), intent(out) :: context
            integer, intent(out) :: status
        end subroutine gc_cylindrical_transport_node_factory_i

        subroutine gc_cylindrical_transport_outer_factor_i(reference, h0, &
                jperp, x, sigma, component_id, launch, sample, user_data, &
                outer_factor, status)
            import :: dp, gc_cylindrical_class_launch_t, &
                gc_nonlocal_orbit_sample_t, gc_nonlocal_transport_reference_t
            type(gc_nonlocal_transport_reference_t), intent(in) :: reference
            real(dp), intent(in) :: h0, jperp, x
            integer, intent(in) :: sigma, component_id
            type(gc_cylindrical_class_launch_t), intent(in) :: launch
            type(gc_nonlocal_orbit_sample_t), intent(in) :: sample
            class(*), pointer, intent(inout) :: user_data
            real(dp), intent(out) :: outer_factor
            integer, intent(out) :: status
        end subroutine gc_cylindrical_transport_outer_factor_i
    end interface

    type, public, extends(gc_nonlocal_transport_provider_t) :: &
            gc_cylindrical_transport_provider_t
        !! One configured harmonic and one explicit H0,Jperp quadrature.
        integer :: harmonic_m = 0
        integer :: harmonic_n = 0
        integer :: force_count = 0
        type(gc_nonlocal_transport_reference_t) :: reference
        character(len=64) :: section_reference_id = ''
        real(dp), allocatable :: h0_nodes(:), h0_weights(:)
        real(dp), allocatable :: jperp_nodes(:), jperp_weights(:)
        procedure(gc_cylindrical_transport_node_factory_i), pointer, nopass :: &
            node_factory => null()
        procedure(gc_cylindrical_transport_outer_factor_i), pointer, nopass :: &
            outer_factor_callback => null()
        class(*), pointer :: user_data => null()
        logical :: initialized = .false.
        logical :: node_ready = .false.
        real(dp) :: node_h0 = 0.0_dp
        real(dp) :: node_jperp = 0.0_dp
        type(gc_cylindrical_class_adapter_t) :: node_adapter
        type(gc_cylindrical_nonlocal_context_t) :: node_context
        type(gc_cylindrical_class_result_t) :: node_class_result
        type(gc_nonlocal_component_t), allocatable :: node_components(:)
        logical :: cache_valid = .false.
        real(dp) :: cache_h0 = 0.0_dp
        real(dp) :: cache_jperp = 0.0_dp
        real(dp) :: cache_x = 0.0_dp
        integer :: cache_sigma = 0
        integer :: cache_component_id = 0
        type(gc_nonlocal_orbit_sample_t) :: cache_sample
        type(gc_cylindrical_class_launch_t) :: cache_launch
        real(dp) :: cache_force(GC_NONLOCAL_MAX_FORCE_VALUES) = 0.0_dp
    contains
        procedure :: get_section_reference => &
            cylindrical_get_section_reference
        procedure :: get_quadrature => cylindrical_get_quadrature
        procedure :: get_components => cylindrical_get_components
        procedure :: evaluate_orbit => cylindrical_evaluate_orbit
        procedure :: evaluate_profiles => cylindrical_evaluate_profiles
        procedure :: evaluate_outer_measure_factor => &
            cylindrical_evaluate_outer_factor
    end type gc_cylindrical_transport_provider_t

    public :: gc_cylindrical_transport_node_factory_i
    public :: gc_cylindrical_transport_outer_factor_i
    public :: initialize_gc_cylindrical_transport_provider
    public :: clear_gc_cylindrical_transport_provider

contains

    subroutine initialize_gc_cylindrical_transport_provider(provider, h0_nodes, &
            h0_weights, jperp_nodes, jperp_weights, reference, harmonic_m, &
            harmonic_n, force_count, node_factory, outer_factor_callback, &
            status, user_data, section_reference_id)
        type(gc_cylindrical_transport_provider_t), intent(out) :: provider
        real(dp), intent(in) :: h0_nodes(:), h0_weights(:)
        real(dp), intent(in) :: jperp_nodes(:), jperp_weights(:)
        type(gc_nonlocal_transport_reference_t), intent(in) :: reference
        integer, intent(in) :: harmonic_m, harmonic_n, force_count
        procedure(gc_cylindrical_transport_node_factory_i), optional :: &
            node_factory
        procedure(gc_cylindrical_transport_outer_factor_i), optional :: &
            outer_factor_callback
        integer, intent(out) :: status
        class(*), target, intent(inout), optional :: user_data
        character(len=*), intent(in), optional :: section_reference_id

        provider%harmonic_m = 0
        provider%harmonic_n = 0
        provider%force_count = 0
        provider%reference = gc_nonlocal_transport_reference_t()
        provider%section_reference_id = ''
        if (allocated(provider%h0_nodes)) deallocate(provider%h0_nodes)
        if (allocated(provider%h0_weights)) deallocate(provider%h0_weights)
        if (allocated(provider%jperp_nodes)) deallocate(provider%jperp_nodes)
        if (allocated(provider%jperp_weights)) deallocate(provider%jperp_weights)
        nullify(provider%node_factory)
        nullify(provider%outer_factor_callback)
        nullify(provider%user_data)
        provider%initialized = .false.
        provider%node_ready = .false.
        provider%cache_valid = .false.
        if (allocated(provider%node_components)) &
            deallocate(provider%node_components)

        status = GC_CYL_TRANSPORT_INVALID_INPUT
        if (harmonic_m == 0 .and. harmonic_n == 0) return
        if (harmonic_n == 0) return
        if (force_count < 1 .or. force_count > GC_NONLOCAL_MAX_FORCE_VALUES) &
            return
        if (.not. valid_reference(reference)) return
        if (.not. present(section_reference_id)) then
            status = GC_CYL_TRANSPORT_REFERENCE_MISMATCH
            return
        end if
        if (len_trim(section_reference_id) == 0) then
            status = GC_CYL_TRANSPORT_REFERENCE_MISMATCH
            return
        end if
        if (.not. valid_quadrature_arrays(h0_nodes, h0_weights, jperp_nodes, &
            jperp_weights)) return
        if (.not. present(node_factory)) then
            status = GC_CYL_TRANSPORT_UNAVAILABLE
            return
        end if
        if (.not. present(outer_factor_callback)) then
            status = GC_CYL_TRANSPORT_UNAVAILABLE
            return
        end if

        allocate(provider%h0_nodes(size(h0_nodes)), &
            provider%h0_weights(size(h0_weights)), &
            provider%jperp_nodes(size(jperp_nodes)), &
            provider%jperp_weights(size(jperp_weights)))
        provider%h0_nodes = h0_nodes
        provider%h0_weights = h0_weights
        provider%jperp_nodes = jperp_nodes
        provider%jperp_weights = jperp_weights
        provider%reference = reference
        provider%section_reference_id = adjustl(section_reference_id)
        provider%harmonic_m = harmonic_m
        provider%harmonic_n = harmonic_n
        provider%force_count = force_count
        provider%node_factory => node_factory
        provider%outer_factor_callback => outer_factor_callback
        if (present(user_data)) provider%user_data => user_data
        provider%initialized = .true.
        status = GC_CYL_TRANSPORT_SUCCESS
    end subroutine initialize_gc_cylindrical_transport_provider

    subroutine clear_gc_cylindrical_transport_provider(provider)
        type(gc_cylindrical_transport_provider_t), intent(inout) :: provider

        if (allocated(provider%h0_nodes)) deallocate(provider%h0_nodes)
        if (allocated(provider%h0_weights)) deallocate(provider%h0_weights)
        if (allocated(provider%jperp_nodes)) deallocate(provider%jperp_nodes)
        if (allocated(provider%jperp_weights)) deallocate(provider%jperp_weights)
        if (allocated(provider%node_components)) &
            deallocate(provider%node_components)
        nullify(provider%node_factory)
        nullify(provider%outer_factor_callback)
        nullify(provider%user_data)
        provider%reference = gc_nonlocal_transport_reference_t()
        provider%section_reference_id = ''
        provider%harmonic_m = 0
        provider%harmonic_n = 0
        provider%force_count = 0
        provider%initialized = .false.
        provider%node_ready = .false.
        provider%cache_valid = .false.
    end subroutine clear_gc_cylindrical_transport_provider

    subroutine cylindrical_get_section_reference(provider, reference, status)
        class(gc_cylindrical_transport_provider_t), intent(inout) :: provider
        type(gc_nonlocal_transport_reference_t), intent(out) :: reference
        integer, intent(out) :: status

        reference = gc_nonlocal_transport_reference_t()
        status = GC_CYL_TRANSPORT_INVALID_INPUT
        if (.not. provider%initialized) return
        reference = provider%reference
        status = GC_CYL_TRANSPORT_SUCCESS
    end subroutine cylindrical_get_section_reference

    subroutine cylindrical_get_quadrature(provider, h0_nodes, h0_weights, &
            jperp_nodes, jperp_weights, status)
        class(gc_cylindrical_transport_provider_t), intent(inout) :: provider
        real(dp), allocatable, intent(out) :: h0_nodes(:), h0_weights(:)
        real(dp), allocatable, intent(out) :: jperp_nodes(:), jperp_weights(:)
        integer, intent(out) :: status

        if (allocated(h0_nodes)) deallocate(h0_nodes)
        if (allocated(h0_weights)) deallocate(h0_weights)
        if (allocated(jperp_nodes)) deallocate(jperp_nodes)
        if (allocated(jperp_weights)) deallocate(jperp_weights)
        status = GC_CYL_TRANSPORT_INVALID_INPUT
        if (.not. provider%initialized) return
        allocate(h0_nodes(size(provider%h0_nodes)), &
            h0_weights(size(provider%h0_weights)), &
            jperp_nodes(size(provider%jperp_nodes)), &
            jperp_weights(size(provider%jperp_weights)))
        h0_nodes = provider%h0_nodes
        h0_weights = provider%h0_weights
        jperp_nodes = provider%jperp_nodes
        jperp_weights = provider%jperp_weights
        status = GC_CYL_TRANSPORT_SUCCESS
    end subroutine cylindrical_get_quadrature

    subroutine cylindrical_get_components(provider, reference, h0, jperp, &
            components, status)
        class(gc_cylindrical_transport_provider_t), intent(inout) :: provider
        type(gc_nonlocal_transport_reference_t), intent(in) :: reference
        real(dp), intent(in) :: h0, jperp
        type(gc_nonlocal_component_t), allocatable, intent(out) :: components(:)
        integer, intent(out) :: status

        if (allocated(components)) deallocate(components)
        status = GC_CYL_TRANSPORT_INVALID_INPUT
        if (.not. provider%initialized) return
        if (.not. same_reference(reference, provider%reference)) then
            status = GC_CYL_TRANSPORT_REFERENCE_MISMATCH
            return
        end if
        call prepare_node(provider, h0, jperp, status)
        if (status /= GC_CYL_TRANSPORT_SUCCESS) return
        allocate(components(size(provider%node_components)))
        components = provider%node_components
        status = GC_CYL_TRANSPORT_SUCCESS
    end subroutine cylindrical_get_components

    subroutine cylindrical_evaluate_orbit(provider, reference, h0, jperp, x, &
            sigma, component_id, sample, status)
        class(gc_cylindrical_transport_provider_t), intent(inout) :: provider
        type(gc_nonlocal_transport_reference_t), intent(in) :: reference
        real(dp), intent(in) :: h0, jperp, x
        integer, intent(in) :: sigma, component_id
        type(gc_nonlocal_orbit_sample_t), intent(out) :: sample
        integer, intent(out) :: status

        type(gc_cylindrical_class_launch_t) :: launch
        integer :: class_status, local_status
        real(dp) :: expected_psi_star, expected_dpsi_star_dx

        sample = gc_nonlocal_orbit_sample_t()
        provider%cache_valid = .false.
        status = GC_CYL_TRANSPORT_INVALID_INPUT
        if (.not. provider%initialized) return
        if (.not. same_reference(reference, provider%reference)) then
            status = GC_CYL_TRANSPORT_REFERENCE_MISMATCH
            return
        end if
        call ensure_node(provider, h0, jperp, local_status)
        if (local_status /= GC_CYL_TRANSPORT_SUCCESS) then
            status = local_status
            return
        end if
        if (.not. point_belongs_to_class(provider%node_class_result%classes, x, &
            sigma, component_id)) then
            status = GC_NONLOCAL_COMPONENT_IDENTITY
            return
        end if

        call launch_gc_cylindrical_class(provider%node_adapter, x, sigma, &
            component_id, launch, class_status)
        if (class_status /= GC_CYL_CLASS_SUCCESS) then
            status = GC_CYL_TRANSPORT_CLASS_UNCERTIFIED
            return
        end if
        if (.not. launch%class_certified) then
            status = GC_CYL_TRANSPORT_CLASS_UNCERTIFIED
            return
        end if
        call evaluate_gc_cylindrical_nonlocal_sample(provider%node_context, x, &
            sigma, component_id, provider%harmonic_m, provider%harmonic_n, &
            provider%force_count, sample, local_status)
        if (local_status /= GC_CYL_NONLOCAL_SUCCESS) then
            status = GC_NONLOCAL_CALLBACK_FAILURE
            return
        end if
        if (sample%status == GC_NONLOCAL_SAMPLE_WALL .or. &
            sample%status == GC_NONLOCAL_SAMPLE_UNRESOLVED) then
            sample%nforce = 0
            sample%thermodynamic_force = 0.0_dp
            status = GC_CYL_TRANSPORT_SUCCESS
            return
        end if
        if (sample%status /= GC_NONLOCAL_SAMPLE_VALID) then
            status = GC_NONLOCAL_CALLBACK_FAILURE
            return
        end if
        if (sample%nforce /= provider%force_count) then
            status = GC_NONLOCAL_CALLBACK_FAILURE
            return
        end if
        if (.not. provider%node_context%charge_c_locked) then
            status = GC_CYL_TRANSPORT_NORMALIZATION_ERROR
            return
        end if
        ! The class adapter has already evaluated psi_star=(c/q)*p_phi.
        ! Its launch derivative is therefore the same normalized quantity;
        ! multiplying it by c/q here would apply the canonical map twice.
        expected_psi_star = launch%psi_star
        expected_dpsi_star_dx = launch%dpsi_star_drc
        if (.not. close_values(sample%psi_star, expected_psi_star, &
            2.0e-10_dp)) then
            status = GC_CYL_TRANSPORT_NORMALIZATION_ERROR
            return
        end if
        if (.not. close_values(sample%dpsi_star_dx, expected_dpsi_star_dx, &
            2.0e-10_dp)) then
            status = GC_CYL_TRANSPORT_NORMALIZATION_ERROR
            return
        end if
        if (.not. all(ieee_is_finite(sample%thermodynamic_force(1: &
            provider%force_count)))) then
            status = GC_NONLOCAL_NONFINITE
            return
        end if

        provider%cache_h0 = h0
        provider%cache_jperp = jperp
        provider%cache_x = x
        provider%cache_sigma = sigma
        provider%cache_component_id = component_id
        provider%cache_sample = sample
        provider%cache_launch = launch
        provider%cache_force = 0.0_dp
        provider%cache_force(1:provider%force_count) = &
            sample%thermodynamic_force(1:provider%force_count)
        provider%cache_valid = .true.
        sample%nforce = 0
        sample%thermodynamic_force = 0.0_dp
        status = GC_CYL_TRANSPORT_SUCCESS
    end subroutine cylindrical_evaluate_orbit

    subroutine cylindrical_evaluate_profiles(provider, reference, h0, jperp, &
            x, sigma, component_id, sample, force_count, force_values, status)
        class(gc_cylindrical_transport_provider_t), intent(inout) :: provider
        type(gc_nonlocal_transport_reference_t), intent(in) :: reference
        real(dp), intent(in) :: h0, jperp, x
        integer, intent(in) :: sigma, component_id, force_count
        type(gc_nonlocal_orbit_sample_t), intent(in) :: sample
        real(dp), intent(out) :: force_values(force_count)
        integer, intent(out) :: status

        force_values = 0.0_dp
        status = GC_CYL_TRANSPORT_INVALID_INPUT
        if (.not. provider%initialized) return
        if (force_count /= provider%force_count) return
        if (sample%status /= GC_NONLOCAL_SAMPLE_VALID) return
        if (.not. same_reference(reference, provider%reference)) then
            status = GC_CYL_TRANSPORT_REFERENCE_MISMATCH
            return
        end if
        if (.not. cache_matches(provider, h0, jperp, x, sigma, component_id, &
            sample)) then
            status = GC_CYL_TRANSPORT_NODE_MISMATCH
            return
        end if
        force_values = provider%cache_force(1:force_count)
        if (.not. all(ieee_is_finite(force_values))) then
            status = GC_NONLOCAL_NONFINITE
            return
        end if
        status = GC_CYL_TRANSPORT_SUCCESS
    end subroutine cylindrical_evaluate_profiles

    subroutine cylindrical_evaluate_outer_factor(provider, reference, h0, &
            jperp, x, sigma, component_id, sample, outer_factor, status)
        class(gc_cylindrical_transport_provider_t), intent(inout) :: provider
        type(gc_nonlocal_transport_reference_t), intent(in) :: reference
        real(dp), intent(in) :: h0, jperp, x
        integer, intent(in) :: sigma, component_id
        type(gc_nonlocal_orbit_sample_t), intent(in) :: sample
        real(dp), intent(out) :: outer_factor
        integer, intent(out) :: status

        outer_factor = 0.0_dp
        status = GC_CYL_TRANSPORT_INVALID_INPUT
        if (.not. provider%initialized) return
        if (.not. associated(provider%outer_factor_callback)) then
            status = GC_CYL_TRANSPORT_UNAVAILABLE
            return
        end if
        if (.not. same_reference(reference, provider%reference)) then
            status = GC_CYL_TRANSPORT_REFERENCE_MISMATCH
            return
        end if
        if (.not. cache_matches(provider, h0, jperp, x, sigma, component_id, &
            sample)) then
            status = GC_CYL_TRANSPORT_NODE_MISMATCH
            return
        end if
        call provider%outer_factor_callback(reference, h0, jperp, x, sigma, &
            component_id, provider%cache_launch, sample, provider%user_data, &
            outer_factor, status)
        if (status /= GC_CYL_TRANSPORT_SUCCESS) then
            status = GC_NONLOCAL_CALLBACK_FAILURE
            return
        end if
        if (.not. ieee_is_finite(outer_factor)) then
            outer_factor = 0.0_dp
            status = GC_NONLOCAL_NONFINITE
            return
        end if
        status = GC_CYL_TRANSPORT_SUCCESS
    end subroutine cylindrical_evaluate_outer_factor

    subroutine prepare_node(provider, h0, jperp, status)
        type(gc_cylindrical_transport_provider_t), intent(inout) :: provider
        real(dp), intent(in) :: h0, jperp
        integer, intent(out) :: status

        call ensure_node(provider, h0, jperp, status, force_refresh=.true.)
    end subroutine prepare_node

    subroutine ensure_node(provider, h0, jperp, status, force_refresh)
        type(gc_cylindrical_transport_provider_t), intent(inout) :: provider
        real(dp), intent(in) :: h0, jperp
        integer, intent(out) :: status
        logical, intent(in), optional :: force_refresh

        type(gc_cylindrical_class_adapter_t) :: candidate_adapter
        type(gc_cylindrical_nonlocal_context_t) :: candidate_context
        type(gc_cylindrical_class_result_t) :: candidate_classes
        type(gc_nonlocal_component_t), allocatable :: callback_components(:)
        integer :: factory_status, class_status, callback_status
        logical :: refresh

        refresh = .false.
        if (present(force_refresh)) refresh = force_refresh
        status = GC_CYL_TRANSPORT_INVALID_INPUT
        if (.not. provider%initialized) return
        if (.not. ieee_is_finite(h0)) return
        if (.not. ieee_is_finite(jperp)) return
        if (jperp < 0.0_dp) return
        if (.not. refresh) then
            if (same_real_exact(provider%node_h0, h0) .and. &
                same_real_exact(provider%node_jperp, jperp)) then
                status = GC_CYL_TRANSPORT_SUCCESS
                return
            end if
        end if
        if (.not. associated(provider%node_factory)) then
            status = GC_CYL_TRANSPORT_UNAVAILABLE
            return
        end if
        call provider%node_factory(h0, jperp, provider%user_data, &
            candidate_adapter, candidate_context, factory_status)
        if (factory_status /= GC_CYL_TRANSPORT_SUCCESS) then
            status = GC_NONLOCAL_CALLBACK_FAILURE
            return
        end if
        status = validate_node_objects(candidate_adapter, candidate_context, h0, &
            jperp, provider%reference, provider%section_reference_id)
        if (status /= GC_CYL_TRANSPORT_SUCCESS) return
        call enumerate_gc_cylindrical_classes(candidate_adapter, candidate_classes, &
            class_status)
        if (class_status /= GC_CYL_CLASS_SUCCESS) then
            status = GC_CYL_TRANSPORT_CLASS_UNCERTIFIED
            return
        end if
        if (.not. candidate_classes%class_complete) then
            status = GC_CYL_TRANSPORT_CLASS_UNCERTIFIED
            return
        end if
        status = validate_class_result(candidate_classes)
        if (status /= GC_CYL_TRANSPORT_SUCCESS) return
        call enumerate_gc_cylindrical_nonlocal_components(candidate_context, &
            callback_components, callback_status)
        if (callback_status /= GC_CYL_NONLOCAL_SUCCESS) then
            status = GC_NONLOCAL_CALLBACK_FAILURE
            return
        end if
        status = compare_components(candidate_classes%classes, callback_components)
        if (status /= GC_CYL_TRANSPORT_SUCCESS) return

        provider%node_adapter = candidate_adapter
        provider%node_context = candidate_context
        provider%node_class_result = candidate_classes
        if (allocated(provider%node_components)) &
            deallocate(provider%node_components)
        allocate(provider%node_components(size(callback_components)))
        provider%node_components = callback_components
        provider%node_h0 = h0
        provider%node_jperp = jperp
        provider%node_ready = .true.
        provider%cache_valid = .false.
        status = GC_CYL_TRANSPORT_SUCCESS
    end subroutine ensure_node

    integer function validate_node_objects(adapter, context, h0, jperp, reference, &
            section_reference_id) &
            result(status)
        type(gc_cylindrical_class_adapter_t), intent(in) :: adapter
        type(gc_cylindrical_nonlocal_context_t), intent(in) :: context
        real(dp), intent(in) :: h0, jperp
        type(gc_nonlocal_transport_reference_t), intent(in) :: reference
        character(len=*), intent(in) :: section_reference_id

        status = GC_CYL_TRANSPORT_INVALID_INPUT
        if (.not. adapter%initialized) return
        if (.not. context%initialized) return
        if (.not. same_real_exact(adapter%h0, h0)) return
        if (.not. same_real_exact(adapter%jperp, jperp)) return
        if (.not. same_real_exact(context%h0, h0)) return
        if (.not. same_real_exact(context%jperp, jperp)) return
        if (adapter%mass <= 0.0_dp) return
        if (abs(adapter%charge) <= tiny(adapter%charge)) return
        if (adapter%c_light <= 0.0_dp) return
        if (.not. context%charge_c_locked) then
            status = GC_CYL_TRANSPORT_NORMALIZATION_ERROR
            return
        end if
        if (.not. same_real_exact(adapter%charge, &
            context%particle_charge)) return
        if (.not. same_real_exact(adapter%c_light, context%c_light)) return
        if (.not. same_context_section(context, reference, section_reference_id)) then
            status = GC_CYL_TRANSPORT_REFERENCE_MISMATCH
            return
        end if
        if (associated(context%canonical_conversion_provider)) then
            status = GC_CYL_TRANSPORT_NORMALIZATION_ERROR
            return
        end if
        if (.not. associated(context%component_provider)) return
        if (.not. associated(context%orbit_provider)) return
        if (.not. associated(context%harmonic_provider)) return
        if (.not. associated(context%force_provider)) return
        status = GC_CYL_TRANSPORT_SUCCESS
    end function validate_node_objects

    integer function validate_class_result(result) result(status)
        type(gc_cylindrical_class_result_t), intent(in) :: result
        integer :: i

        status = GC_CYL_TRANSPORT_CLASS_UNCERTIFIED
        if (.not. result%class_complete) return
        if (result%nclasses < 1) return
        if (.not. allocated(result%classes)) return
        if (size(result%classes) /= result%nclasses) return
        do i = 1, result%nclasses
            if (.not. result%classes(i)%allowed_interval) return
            if (.not. result%classes(i)%topology_certified) return
            if (.not. result%classes(i)%orbit_return_certified) return
            if (result%classes(i)%component_id <= 0) return
            if (abs(result%classes(i)%sigma) /= 1) return
            if (result%classes(i)%rc_max <= result%classes(i)%rc_min) return
        end do
        status = GC_CYL_TRANSPORT_SUCCESS
    end function validate_class_result

    integer function compare_components(classes, components) result(status)
        type(gc_cylindrical_class_interval_t), intent(in) :: classes(:)
        type(gc_nonlocal_component_t), intent(in) :: components(:)
        logical :: found
        integer :: i, j

        status = GC_NONLOCAL_COMPONENT_IDENTITY
        if (size(classes) /= size(components)) return
        do i = 1, size(classes)
            found = .false.
            do j = 1, size(components)
                if (classes(i)%component_id /= components(j)%component_id) cycle
                if (classes(i)%sigma /= components(j)%sigma) cycle
                if (.not. close_values(classes(i)%rc_min, components(j)%x_min, &
                    2.0e-10_dp)) return
                if (.not. close_values(classes(i)%rc_max, components(j)%x_max, &
                    2.0e-10_dp)) return
                found = .true.
                exit
            end do
            if (.not. found) return
        end do
        status = GC_CYL_TRANSPORT_SUCCESS
    end function compare_components

    logical function point_belongs_to_class(classes, x, sigma, component_id)
        type(gc_cylindrical_class_interval_t), intent(in) :: classes(:)
        real(dp), intent(in) :: x
        integer, intent(in) :: sigma, component_id
        integer :: i

        point_belongs_to_class = .false.
        do i = 1, size(classes)
            if (classes(i)%sigma /= sigma) cycle
            if (classes(i)%component_id /= component_id) cycle
            if (x < classes(i)%rc_min) return
            if (x > classes(i)%rc_max) return
            point_belongs_to_class = .true.
            return
        end do
    end function point_belongs_to_class

    logical function cache_matches(provider, h0, jperp, x, sigma, component_id, &
            sample)
        type(gc_cylindrical_transport_provider_t), intent(in) :: provider
        real(dp), intent(in) :: h0, jperp, x
        integer, intent(in) :: sigma, component_id
        type(gc_nonlocal_orbit_sample_t), intent(in) :: sample

        cache_matches = provider%cache_valid
        if (.not. cache_matches) return
        if (.not. same_real_exact(provider%cache_h0, h0)) then
            cache_matches = .false.
            return
        end if
        if (.not. same_real_exact(provider%cache_jperp, jperp)) then
            cache_matches = .false.
            return
        end if
        if (.not. same_real_exact(provider%cache_x, x)) then
            cache_matches = .false.
            return
        end if
        if (provider%cache_sigma /= sigma) then
            cache_matches = .false.
            return
        end if
        if (provider%cache_component_id /= component_id) then
            cache_matches = .false.
            return
        end if
        cache_matches = same_sample(sample, provider%cache_sample)
    end function cache_matches

    logical function same_sample(first, second)
        type(gc_nonlocal_orbit_sample_t), intent(in) :: first, second

        same_sample = .true.
        if (first%status /= second%status) same_sample = .false.
        if (first%component_id /= second%component_id) same_sample = .false.
        if (first%sigma /= second%sigma) same_sample = .false.
        if (.not. close_values(first%psi_star, second%psi_star, 2.0e-10_dp)) &
            same_sample = .false.
        if (.not. close_values(first%dpsi_star_dx, second%dpsi_star_dx, &
            2.0e-10_dp)) same_sample = .false.
        if (.not. close_values(first%tau_b, second%tau_b, 2.0e-10_dp)) &
            same_sample = .false.
        if (.not. close_values(first%omega_b, second%omega_b, 2.0e-10_dp)) &
            same_sample = .false.
        if (.not. close_values(first%omega_phi, second%omega_phi, 2.0e-10_dp)) &
            same_sample = .false.
        if (.not. close_values(first%domega_b_dx, second%domega_b_dx, &
            2.0e-10_dp)) same_sample = .false.
        if (.not. close_values(first%domega_phi_dx, second%domega_phi_dx, &
            2.0e-10_dp)) same_sample = .false.
        if (abs(first%h_m - second%h_m) > 2.0e-10_dp) same_sample = .false.
    end function same_sample

    logical function valid_reference(reference)
        type(gc_nonlocal_transport_reference_t), intent(in) :: reference

        valid_reference = .false.
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
        valid_reference = .true.
    end function valid_reference

    logical function valid_quadrature_arrays(h0_nodes, h0_weights, jperp_nodes, &
            jperp_weights)
        real(dp), intent(in) :: h0_nodes(:), h0_weights(:), jperp_nodes(:), &
            jperp_weights(:)

        valid_quadrature_arrays = .false.
        if (size(h0_nodes) < 1 .or. size(jperp_nodes) < 1) return
        if (size(h0_nodes) /= size(h0_weights)) return
        if (size(jperp_nodes) /= size(jperp_weights)) return
        if (.not. all(ieee_is_finite(h0_nodes))) return
        if (.not. all(ieee_is_finite(h0_weights))) return
        if (.not. all(ieee_is_finite(jperp_nodes))) return
        if (.not. all(ieee_is_finite(jperp_weights))) return
        if (any(jperp_nodes < 0.0_dp)) return
        valid_quadrature_arrays = .true.
    end function valid_quadrature_arrays

    logical function same_reference(first, second)
        type(gc_nonlocal_transport_reference_t), intent(in) :: first, second
        integer :: i

        same_reference = first%section_id == second%section_id
        if (.not. same_reference) return
        if (trim(first%section_coordinate) /= trim(second%section_coordinate)) &
            same_reference = .false.
        if (.not. same_reference) return
        if (trim(first%section_units) /= trim(second%section_units)) &
            same_reference = .false.
        if (.not. same_reference) return
        if (first%p_zeta_orientation /= second%p_zeta_orientation) &
            same_reference = .false.
        if (.not. same_reference) return
        if (first%frequency_semantics /= second%frequency_semantics) &
            same_reference = .false.
        if (.not. same_reference) return
        if (first%hamiltonian_includes_phi .neqv. &
            second%hamiltonian_includes_phi) same_reference = .false.
        if (first%fixed .neqv. second%fixed) same_reference = .false.
        if (.not. same_reference) return
        if (.not. close_values(first%section_flux, second%section_flux, &
            2.0e-12_dp)) same_reference = .false.
        do i = 1, 3
            if (.not. close_values(first%section_position(i), &
                second%section_position(i), 2.0e-12_dp)) then
                same_reference = .false.
                return
            end if
        end do
    end function same_reference

    logical function same_context_section(context, reference, section_reference_id)
        type(gc_cylindrical_nonlocal_context_t), intent(in) :: context
        type(gc_nonlocal_transport_reference_t), intent(in) :: reference
        character(len=*), intent(in) :: section_reference_id

        integer :: i

        same_context_section = context%section%locked
        if (.not. same_context_section) return
        if (trim(context%section%coordinate) /= &
            trim(reference%section_coordinate)) then
            same_context_section = .false.
            return
        end if
        if (trim(context%section%reference_id) /= trim(section_reference_id)) then
            same_context_section = .false.
            return
        end if
        do i = 1, 3
            if (.not. close_values(context%section%reference(i), &
                reference%section_position(i), 2.0e-12_dp)) then
                same_context_section = .false.
                return
            end if
        end do
    end function same_context_section

    logical function close_values(value, expected, tolerance)
        real(dp), intent(in) :: value, expected, tolerance

        close_values = abs(value - expected) <= tolerance*max(1.0_dp, &
            max(abs(value), abs(expected)))
    end function close_values

    pure logical function same_real_exact(first, second)
        real(dp), intent(in) :: first, second
        integer(int64) :: first_bits, second_bits

        first_bits = transfer(first, first_bits)
        second_bits = transfer(second, second_bits)
        same_real_exact = first_bits == second_bits
    end function same_real_exact

end module neort_gc_cylindrical_transport_provider
