module neort_gc_nonlocal_transport_types
    !! Contracts for the physical full-orbit transport integral.
    !!
    !! Buchholz et al. Eq. 14--18 are represented explicitly here.  The
    !! phase-space node is a *paired* (H0,J_K) point: J_K is allowed to have
    !! a different upper bound at every H0 node.  In particular, the
    !! quadrature is not a Cartesian tensor of independently bounded arrays.
    !! H0 is the Hamiltonian including q*Phi and J_K is the positive
    !! positive action J_K=m*c*mu_phys/abs(q), so that
    !! e*J_K/(m*c*Eref)=p_perp_hat**2/B, with the
    !! charge magnitude used in that latter normalization.  The signed q is
    !! retained in psi_star=(c/q)*p_phi and in the Maxwellian derivative.
    !!
    !! The root kernel owns |d psi_star/dx|*|H_m|^2*tau_b/|dF/dx|, where
    !! F=m2*omega_b+m3*omega_phi.  This is algebraically identical to the
    !! Eq. 17 phase-root weight |m3|*|d psi_star/dx|*|H_m|^2*tau_b**2/
    !! |d(Delta_phi_b+2*pi*m2/m3)/dx|.  The transport layer owns m3**2 and
    !! the physical signed Eq. 10 prefactor exactly once.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_gc_nonlocal_resonance_types, only: &
        GC_NONLOCAL_CAPACITY, GC_NONLOCAL_CALLBACK_FAILURE, &
        GC_NONLOCAL_COMPONENT_IDENTITY, GC_NONLOCAL_DERIVATIVE_MISSING, &
        GC_NONLOCAL_FORCE_CONTRACT, GC_NONLOCAL_INVALID_INPUT, &
        GC_NONLOCAL_MAX_FORCE_VALUES, GC_NONLOCAL_NONFINITE, &
        GC_NONLOCAL_PARTIAL, GC_NONLOCAL_ROOT_NOT_CONVERGED, &
        GC_NONLOCAL_SAMPLE_INVALID, GC_NONLOCAL_SAMPLE_UNRESOLVED, &
        GC_NONLOCAL_SAMPLE_VALID, GC_NONLOCAL_SAMPLE_WALL, &
        GC_NONLOCAL_SINGULAR_RESONANCE, GC_NONLOCAL_SUCCESS, &
        gc_nonlocal_component_t, gc_nonlocal_orbit_sample_t, &
        gc_nonlocal_resonance_options_t

    implicit none
    private

    integer, parameter, public :: GC_NONLOCAL_TRANSPORT_MAX_H0_NODES = 256
    integer, parameter, public :: GC_NONLOCAL_TRANSPORT_MAX_JPERP_NODES = 256
    integer, parameter, public :: GC_NONLOCAL_TRANSPORT_MAX_TOTAL_NODES = 65536
    integer, parameter, public :: GC_NONLOCAL_COMPLETE_CYCLE_FREQUENCIES = 1

    integer, parameter, public :: GC_NONLOCAL_CLASS_COPASSING = 1
    integer, parameter, public :: GC_NONLOCAL_CLASS_COUNTERPASSING = 2
    integer, parameter, public :: GC_NONLOCAL_CLASS_TRAPPED = 3
    integer, parameter, public :: GC_NONLOCAL_CLASS_COUNT = 3

    public :: GC_NONLOCAL_CAPACITY, GC_NONLOCAL_CALLBACK_FAILURE
    public :: GC_NONLOCAL_COMPONENT_IDENTITY, GC_NONLOCAL_DERIVATIVE_MISSING
    public :: GC_NONLOCAL_FORCE_CONTRACT, GC_NONLOCAL_INVALID_INPUT
    public :: GC_NONLOCAL_MAX_FORCE_VALUES, GC_NONLOCAL_NONFINITE
    public :: GC_NONLOCAL_PARTIAL, GC_NONLOCAL_ROOT_NOT_CONVERGED
    public :: GC_NONLOCAL_SAMPLE_INVALID, GC_NONLOCAL_SAMPLE_UNRESOLVED
    public :: GC_NONLOCAL_SAMPLE_VALID, GC_NONLOCAL_SAMPLE_WALL
    public :: GC_NONLOCAL_SINGULAR_RESONANCE, GC_NONLOCAL_SUCCESS
    public :: dp, gc_nonlocal_component_t, gc_nonlocal_orbit_sample_t
    public :: gc_nonlocal_resonance_options_t

    type, public :: gc_nonlocal_transport_reference_t
        integer :: section_id = 0
        character(len=32) :: section_coordinate = ''
        character(len=16) :: section_units = ''
        real(dp) :: section_position(3) = 0.0_dp
        real(dp) :: section_flux = 0.0_dp
        integer :: p_zeta_orientation = 0
        integer :: frequency_semantics = 0
        !! Normalization in Buchholz Eq. 17--18.  Physical providers pass
        !! dimensional H0,J_K to their dynamics and expose normalized
        !! sample quantities to the root kernel.
        real(dp) :: energy_scale = 1.0_dp
        real(dp) :: velocity_scale = 1.0_dp
        real(dp) :: psi_star_scale = 1.0_dp
        !! Dimensional ownership for quadrature_absolute_tolerance.  The
        !! tolerance itself is dimensionless; each force slot supplies its
        !! physical scale here, so no mixed-unit absolute comparison is
        !! possible at the convergence gate.
        real(dp) :: force_scale(GC_NONLOCAL_MAX_FORCE_VALUES) = 1.0_dp
        logical :: hamiltonian_includes_phi = .false.
        logical :: fixed = .false.
    end type gc_nonlocal_transport_reference_t

    type, public :: gc_nonlocal_transport_quadrature_t
        !! One physical pair per entry.  jperp_upper_bound(i) is the actual
        !! paired-domain bound used to map the J_K node at H0(i).
        integer :: h0_order = 0
        integer :: jk_order = 0
        integer :: n_nodes = 0
        logical :: paired_domain = .false.
        logical :: domain_certified = .false.
        logical :: upper_bound_is_envelope = .false.
        character(len=64) :: domain_certificate = ''
        logical :: converged = .false.
        real(dp) :: h0_min = 0.0_dp
        real(dp) :: h0_scale = 0.0_dp
        real(dp), allocatable :: h0(:)
        real(dp), allocatable :: j_k(:)
        real(dp), allocatable :: weight(:)
        real(dp), allocatable :: j_k_upper_bound(:)
    end type gc_nonlocal_transport_quadrature_t

    type, public :: gc_nonlocal_transport_options_t
        !! Physical outer quadrature orders, not a node-count shortcut.
        integer :: h0_order = 16
        integer :: jk_order = 16
        integer :: max_h0_nodes = 256
        integer :: max_jperp_nodes = 256
        integer :: max_total_nodes = 65536
        real(dp) :: quadrature_relative_tolerance = 1.0e-7_dp
        real(dp) :: quadrature_absolute_tolerance = 1.0e-12_dp
        logical :: require_converged = .true.
        type(gc_nonlocal_resonance_options_t) :: resonance_options = &
            gc_nonlocal_resonance_options_t()
    end type gc_nonlocal_transport_options_t

    type, public :: gc_nonlocal_transport_observed_evidence_t
        !! Provider-observed counters.  Each axis is conserved independently:
        !! invariant_successes+invariant_failures, return_successes+
        !! return_no_return+return_radial_domain+return_wall_loss+
        !! return_errors, and wall_not_hit+wall_loss+wall_errors must each
        !! equal physical_return_attempts.
        integer :: physical_return_attempts = 0
        integer :: invariant_successes = 0
        integer :: invariant_failures = 0
        integer :: return_successes = 0
        integer :: return_no_return = 0
        integer :: return_radial_domain = 0
        integer :: return_wall_loss = 0
        integer :: return_errors = 0
        integer :: wall_not_hit = 0
        integer :: wall_loss = 0
        integer :: wall_errors = 0
        integer :: topology_certification_attempts = 0
        integer :: topology_certification_successes = 0
        integer :: harmonic_average_successes = 0
        integer :: residence_average_successes = 0
        logical :: return_accounting_complete = .false.
    end type gc_nonlocal_transport_observed_evidence_t

    type, public :: gc_nonlocal_transport_result_t
        integer :: status = GC_NONLOCAL_INVALID_INPUT
        integer :: harmonic_m = 0
        integer :: harmonic_n = 0
        integer :: nforce = 0
        integer :: nh0 = 0
        integer :: njperp = 0
        integer :: weighted_nodes = 0
        integer :: certified_nodes = 0
        integer :: unresolved_nodes = 0
        integer :: ncomponents = 0
        integer :: nroots = 0
        integer :: failed_node = 0
        integer :: failure_status = GC_NONLOCAL_INVALID_INPUT
        type(gc_nonlocal_transport_reference_t) :: reference
        type(gc_nonlocal_transport_quadrature_t) :: quadrature
        real(dp) :: contribution(GC_NONLOCAL_MAX_FORCE_VALUES) = 0.0_dp
        real(dp) :: diagnostic_contribution(GC_NONLOCAL_MAX_FORCE_VALUES) = 0.0_dp
        real(dp) :: class_contribution(GC_NONLOCAL_CLASS_COUNT, &
            GC_NONLOCAL_MAX_FORCE_VALUES) = 0.0_dp
        logical :: certified = .false.
        real(dp), allocatable :: h0_nodes(:)
        real(dp), allocatable :: jperp_nodes(:)
        real(dp), allocatable :: node_weights(:)
        integer, allocatable :: node_status(:)
        logical, allocatable :: node_certified(:)
        real(dp), allocatable :: node_contribution(:, :)
    end type gc_nonlocal_transport_result_t

    type, abstract, public :: gc_nonlocal_transport_provider_t
    contains
        procedure(gc_nonlocal_transport_reference_i), deferred :: &
            get_section_reference
        procedure(gc_nonlocal_transport_quadrature_i), deferred :: &
            get_quadrature
        procedure(gc_nonlocal_transport_components_i), deferred :: &
            get_components
        procedure(gc_nonlocal_transport_orbit_i), deferred :: evaluate_orbit
        procedure(gc_nonlocal_transport_profiles_i), deferred :: &
            evaluate_profiles
        procedure(gc_nonlocal_transport_outer_factor_i), deferred :: &
            evaluate_outer_measure_factor
        procedure :: set_harmonic => default_set_harmonic
        procedure :: get_class_kind => default_get_class_kind
        procedure :: begin_execution => default_begin_execution
        procedure :: get_execution_evidence => default_get_execution_evidence
    end type gc_nonlocal_transport_provider_t

    abstract interface
        subroutine gc_nonlocal_transport_reference_i(provider, reference, status)
            import :: gc_nonlocal_transport_provider_t, &
                gc_nonlocal_transport_reference_t
            class(gc_nonlocal_transport_provider_t), intent(inout) :: provider
            type(gc_nonlocal_transport_reference_t), intent(out) :: reference
            integer, intent(out) :: status
        end subroutine gc_nonlocal_transport_reference_i

        subroutine gc_nonlocal_transport_quadrature_i(provider, h0_order, jk_order, &
                quadrature, status)
            import :: gc_nonlocal_transport_provider_t, &
                gc_nonlocal_transport_quadrature_t
            class(gc_nonlocal_transport_provider_t), intent(inout) :: provider
            integer, intent(in) :: h0_order, jk_order
            type(gc_nonlocal_transport_quadrature_t), intent(out) :: quadrature
            integer, intent(out) :: status
        end subroutine gc_nonlocal_transport_quadrature_i

        subroutine gc_nonlocal_transport_components_i(provider, reference, h0, &
                jperp, components, status)
            import :: dp, gc_nonlocal_component_t, &
                gc_nonlocal_transport_provider_t, &
                gc_nonlocal_transport_reference_t
            class(gc_nonlocal_transport_provider_t), intent(inout) :: provider
            type(gc_nonlocal_transport_reference_t), intent(in) :: reference
            real(dp), intent(in) :: h0, jperp
            type(gc_nonlocal_component_t), allocatable, intent(out) :: &
                components(:)
            integer, intent(out) :: status
        end subroutine gc_nonlocal_transport_components_i

        subroutine gc_nonlocal_transport_orbit_i(provider, reference, h0, jperp, &
                x, sigma, component_id, sample, status)
            import :: dp, gc_nonlocal_orbit_sample_t, &
                gc_nonlocal_transport_provider_t, &
                gc_nonlocal_transport_reference_t
            class(gc_nonlocal_transport_provider_t), intent(inout) :: provider
            type(gc_nonlocal_transport_reference_t), intent(in) :: reference
            real(dp), intent(in) :: h0, jperp, x
            integer, intent(in) :: sigma, component_id
            type(gc_nonlocal_orbit_sample_t), intent(out) :: sample
            integer, intent(out) :: status
        end subroutine gc_nonlocal_transport_orbit_i

        subroutine gc_nonlocal_transport_profiles_i(provider, reference, h0, &
                jperp, x, sigma, component_id, sample, force_count, &
                force_values, status)
            import :: dp, gc_nonlocal_orbit_sample_t, &
                gc_nonlocal_transport_provider_t, &
                gc_nonlocal_transport_reference_t
            class(gc_nonlocal_transport_provider_t), intent(inout) :: provider
            type(gc_nonlocal_transport_reference_t), intent(in) :: reference
            real(dp), intent(in) :: h0, jperp, x
            integer, intent(in) :: sigma, component_id, force_count
            type(gc_nonlocal_orbit_sample_t), intent(in) :: sample
            real(dp), intent(out) :: force_values(force_count)
            integer, intent(out) :: status
        end subroutine gc_nonlocal_transport_profiles_i

        subroutine gc_nonlocal_transport_outer_factor_i(provider, reference, h0, &
                jperp, x, sigma, component_id, sample, outer_factor, status)
            import :: dp, gc_nonlocal_orbit_sample_t, &
                gc_nonlocal_transport_provider_t, &
                gc_nonlocal_transport_reference_t
            class(gc_nonlocal_transport_provider_t), intent(inout) :: provider
            type(gc_nonlocal_transport_reference_t), intent(in) :: reference
            real(dp), intent(in) :: h0, jperp, x
            integer, intent(in) :: sigma, component_id
            type(gc_nonlocal_orbit_sample_t), intent(in) :: sample
            real(dp), intent(out) :: outer_factor
            integer, intent(out) :: status
        end subroutine gc_nonlocal_transport_outer_factor_i

        subroutine gc_nonlocal_transport_set_harmonic_i(provider, harmonic_m, &
                harmonic_n, status)
            import :: gc_nonlocal_transport_provider_t
            class(gc_nonlocal_transport_provider_t), intent(inout) :: provider
            integer, intent(in) :: harmonic_m, harmonic_n
            integer, intent(out) :: status
        end subroutine gc_nonlocal_transport_set_harmonic_i

        subroutine gc_nonlocal_transport_class_kind_i(provider, reference, h0, &
                jperp, x, sigma, component_id, class_kind, status)
            import :: dp, gc_nonlocal_transport_provider_t, &
                gc_nonlocal_transport_reference_t
            class(gc_nonlocal_transport_provider_t), intent(inout) :: provider
            type(gc_nonlocal_transport_reference_t), intent(in) :: reference
            real(dp), intent(in) :: h0, jperp, x
            integer, intent(in) :: sigma, component_id
            integer, intent(out) :: class_kind, status
        end subroutine gc_nonlocal_transport_class_kind_i

        subroutine gc_nonlocal_transport_begin_execution_i(provider, status)
            import :: gc_nonlocal_transport_provider_t
            class(gc_nonlocal_transport_provider_t), intent(inout) :: provider
            integer, intent(out) :: status
        end subroutine gc_nonlocal_transport_begin_execution_i

        subroutine gc_nonlocal_transport_evidence_i(provider, evidence, status)
            import :: gc_nonlocal_transport_observed_evidence_t, &
                gc_nonlocal_transport_provider_t
            class(gc_nonlocal_transport_provider_t), intent(inout) :: provider
            type(gc_nonlocal_transport_observed_evidence_t), intent(out) :: evidence
            integer, intent(out) :: status
        end subroutine gc_nonlocal_transport_evidence_i
    end interface

    public :: gc_nonlocal_transport_reference_i
    public :: gc_nonlocal_transport_components_i
    public :: gc_nonlocal_transport_orbit_i
    public :: gc_nonlocal_transport_profiles_i
    public :: gc_nonlocal_transport_outer_factor_i
    public :: gc_nonlocal_transport_set_harmonic_i
    public :: gc_nonlocal_transport_class_kind_i
    public :: gc_nonlocal_transport_begin_execution_i
    public :: gc_nonlocal_transport_evidence_i

contains

    subroutine default_set_harmonic(provider, harmonic_m, harmonic_n, status)
        class(gc_nonlocal_transport_provider_t), intent(inout) :: provider
        integer, intent(in) :: harmonic_m, harmonic_n
        integer, intent(out) :: status

        associate (unused_provider => provider, unused_m => harmonic_m, &
                unused_n => harmonic_n)
        end associate
        status = GC_NONLOCAL_SUCCESS
    end subroutine default_set_harmonic

    subroutine default_get_class_kind(provider, reference, h0, jperp, x, sigma, &
            component_id, class_kind, status)
        class(gc_nonlocal_transport_provider_t), intent(inout) :: provider
        type(gc_nonlocal_transport_reference_t), intent(in) :: reference
        real(dp), intent(in) :: h0, jperp, x
        integer, intent(in) :: sigma, component_id
        integer, intent(out) :: class_kind, status

        associate (unused_provider => provider, unused_reference => reference, &
                unused_h0 => h0, unused_jperp => jperp, unused_x => x, &
                unused_component => component_id)
        end associate
        class_kind = GC_NONLOCAL_CLASS_COPASSING
        if (sigma < 0) class_kind = GC_NONLOCAL_CLASS_COUNTERPASSING
        status = GC_NONLOCAL_SUCCESS
    end subroutine default_get_class_kind

    subroutine default_begin_execution(provider, status)
        class(gc_nonlocal_transport_provider_t), intent(inout) :: provider
        integer, intent(out) :: status

        associate (unused_provider => provider)
        end associate
        status = GC_NONLOCAL_SUCCESS
    end subroutine default_begin_execution

    subroutine default_get_execution_evidence(provider, evidence, status)
        class(gc_nonlocal_transport_provider_t), intent(inout) :: provider
        type(gc_nonlocal_transport_observed_evidence_t), intent(out) :: evidence
        integer, intent(out) :: status

        associate (unused_provider => provider)
        end associate
        evidence = gc_nonlocal_transport_observed_evidence_t()
        status = GC_NONLOCAL_SUCCESS
    end subroutine default_get_execution_evidence

end module neort_gc_nonlocal_transport_types
