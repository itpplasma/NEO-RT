module neort_gc_nonlocal_transport_types
    !! Contracts for the outer full-finite-orbit-width transport integral.
    !!
    !! The provider is deliberately an abstract type.  It owns the physical
    !! cylindrical implementation and must expose all quantities needed by
    !! Buchholz et al. Eq. 14--17 without converting to a local Boozer eta
    !! coordinate:
    !!
    !!   r = m2*omega_b + m3*omega_phi,
    !!   g = Delta_phi_b + 2*pi*m2/m3,
    !!
    !!   sum_i,j w_H(i) w_J(j) sum_sigma int dx
    !!       m3^2 W_outer |d psi_star/dx| |H_m|^2 tau_b delta(r)
    !!
    !! is identical at a simple root to the phase-residual form
    !!
    !!   sum_i,j w_H(i) w_J(j) sum_sigma int dx
    !!       |m3| W_outer |d psi_star/dx| |H_m|^2 tau_b^2 delta(g).
    !!
    !! The inner resonance kernel owns the frequency-residual factor
    !! tau_b/|d r/dx|.  The outer layer supplies m3^2 exactly once.  The
    !! phase-equivalent |m3| form is an oracle, not a second multiplicative
    !! factor in the production path.  Eq. 17 requires m3 /= 0.
    !!
    !! The orbit evaluator and the existing fixed-node resonance kernel own
    !! the sole |d psi_star/dx| factor.  The outer-measure callback returns
    !! only W_outer: the explicitly signed Eq. 17 constants, physical
    !! distribution factors, and any orbit-volume selector not already in
    !! A_*, |H_m|^2, tau_b, or the harmonic conversion.  It must not return
    !! |d p_phi/dx|, tau_b, m3^2, or any equivalent canonical Jacobian.  Eq.
    !! 14's |d p_phi/dR_c| has already become the Eq. 17 |d psi_star/dx|
    !! factor through Eqs. 15--16 and may not be multiplied in a second time.
    !!
    !! The profile callback returns signed native thermodynamic forces.  The
    !! orbit callback returns native signed frequencies and complex
    !! perturbation amplitude.  No sign alignment or normalization is
    !! performed here.
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

    integer, parameter, public :: GC_NONLOCAL_TRANSPORT_MAX_H0_NODES = 1024
    integer, parameter, public :: GC_NONLOCAL_TRANSPORT_MAX_JPERP_NODES = 1024
    integer, parameter, public :: GC_NONLOCAL_TRANSPORT_MAX_TOTAL_NODES = 262144
    integer, parameter, public :: &
        GC_NONLOCAL_COMPLETE_CYCLE_FREQUENCIES = 1

    ! The outer layer uses the status vocabulary of the fixed-(H0,J_perp)
    ! resonance kernel.  Re-exporting it makes a provider adapter's contract
    ! self-contained while preserving exact status propagation.
    public :: GC_NONLOCAL_CAPACITY, GC_NONLOCAL_CALLBACK_FAILURE
    public :: GC_NONLOCAL_COMPONENT_IDENTITY, GC_NONLOCAL_DERIVATIVE_MISSING
    public :: GC_NONLOCAL_FORCE_CONTRACT, GC_NONLOCAL_INVALID_INPUT
    public :: GC_NONLOCAL_MAX_FORCE_VALUES, GC_NONLOCAL_NONFINITE
    public :: GC_NONLOCAL_PARTIAL, GC_NONLOCAL_ROOT_NOT_CONVERGED
    public :: GC_NONLOCAL_SAMPLE_INVALID, GC_NONLOCAL_SAMPLE_UNRESOLVED
    public :: GC_NONLOCAL_SAMPLE_VALID, GC_NONLOCAL_SAMPLE_WALL
    public :: GC_NONLOCAL_SINGULAR_RESONANCE, GC_NONLOCAL_SUCCESS
    public :: dp
    public :: gc_nonlocal_component_t, gc_nonlocal_orbit_sample_t
    public :: gc_nonlocal_resonance_options_t

    type, public :: gc_nonlocal_transport_reference_t
        !! Fixed Poincare section and frequency convention for one run.
        !!
        !! section_position is in section_units and is held fixed for all
        !! H0,J_perp nodes and every disconnected class.  The section
        !! coordinate x parameterizes the orbit label P_zeta on that fixed
        !! cut (possibly through a regularized transform of Rc); it is not an
        !! instantaneous geometric angle.
        integer :: section_id = 0
        character(len=32) :: section_coordinate = ''
        character(len=16) :: section_units = ''
        real(dp) :: section_position(3) = 0.0_dp
        real(dp) :: section_flux = 0.0_dp
        integer :: p_zeta_orientation = 0
        integer :: frequency_semantics = 0
        logical :: hamiltonian_includes_phi = .false.
        logical :: fixed = .false.
    end type gc_nonlocal_transport_reference_t

    type, public :: gc_nonlocal_transport_options_t
        !! Bounds for the tensor-product outer quadrature and inner kernel.
        integer :: max_h0_nodes = 64
        integer :: max_jperp_nodes = 64
        integer :: max_total_nodes = 4096
        type(gc_nonlocal_resonance_options_t) :: resonance_options = &
            gc_nonlocal_resonance_options_t()
    end type gc_nonlocal_transport_options_t

    type, public :: gc_nonlocal_transport_result_t
        !! Accepted values are present only when certified is true.
        !!
        !! diagnostic_contribution and node_contribution can contain a
        !! partial sum after a failure.  They are intentionally separate from
        !! contribution, which is reset to zero on every uncertified return.
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
        integer :: failed_h0_index = 0
        integer :: failed_jperp_index = 0
        integer :: failure_status = GC_NONLOCAL_INVALID_INPUT
        type(gc_nonlocal_transport_reference_t) :: reference
        real(dp) :: contribution(GC_NONLOCAL_MAX_FORCE_VALUES) = 0.0_dp
        real(dp) :: diagnostic_contribution(GC_NONLOCAL_MAX_FORCE_VALUES) = 0.0_dp
        logical :: certified = .false.
        real(dp), allocatable :: h0_nodes(:)
        real(dp), allocatable :: h0_weights(:)
        real(dp), allocatable :: jperp_nodes(:)
        real(dp), allocatable :: jperp_weights(:)
        integer, allocatable :: node_status(:, :)
        logical, allocatable :: node_certified(:, :)
        real(dp), allocatable :: node_contribution(:, :, :)
    end type gc_nonlocal_transport_result_t

    type, abstract, public :: gc_nonlocal_transport_provider_t
        !! Physical provider for one complete Eq. 14--17 calculation.
        !!
        !! A concrete provider must return a new complete class list for each
        !! fixed H0,J_perp pair.  It must not merge opposite-sigma branches,
        !! fit signs, or silently replace an unresolved orbit by a nearby one.
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
    end type gc_nonlocal_transport_provider_t

    abstract interface
        subroutine gc_nonlocal_transport_reference_i(provider, reference, status)
            import :: gc_nonlocal_transport_provider_t, &
                gc_nonlocal_transport_reference_t
            class(gc_nonlocal_transport_provider_t), intent(inout) :: provider
            type(gc_nonlocal_transport_reference_t), intent(out) :: reference
            integer, intent(out) :: status
        end subroutine gc_nonlocal_transport_reference_i

        subroutine gc_nonlocal_transport_quadrature_i(provider, h0_nodes, &
                h0_weights, jperp_nodes, jperp_weights, status)
            import :: dp, gc_nonlocal_transport_provider_t
            class(gc_nonlocal_transport_provider_t), intent(inout) :: provider
            real(dp), allocatable, intent(out) :: h0_nodes(:), h0_weights(:)
            real(dp), allocatable, intent(out) :: jperp_nodes(:), jperp_weights(:)
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
    end interface

    abstract interface
        subroutine gc_nonlocal_transport_local_reference_i(h0, jperp, &
                force_count, values, status)
            import :: dp
            real(dp), intent(in) :: h0, jperp
            integer, intent(in) :: force_count
            real(dp), intent(out) :: values(force_count)
            integer, intent(out) :: status
        end subroutine gc_nonlocal_transport_local_reference_i
    end interface

    public :: gc_nonlocal_transport_quadrature_i
    public :: gc_nonlocal_transport_reference_i
    public :: gc_nonlocal_transport_components_i
    public :: gc_nonlocal_transport_orbit_i
    public :: gc_nonlocal_transport_profiles_i
    public :: gc_nonlocal_transport_outer_factor_i
    public :: gc_nonlocal_transport_local_reference_i

end module neort_gc_nonlocal_transport_types
