module neort_gc_cylindrical_nonlocal_provider
    !! Adapter contract between cylindrical full-orbit dynamics and the
    !! nonlocal Eq. 14--17 phase-space integral.
    !!
    !! The class coordinate is an explicitly named physical Poincare-section
    !! coordinate.  It defaults to R_c.  H0 and J_perp are passed unchanged
    !! through every callback.  In particular, this module contains no eta
    !! conversion and does not turn an energy-allowed interval into a
    !! nonlocal invariant measure.
    !!
    !! The callbacks are deliberately strict.  The existing cylindrical
    !! orbit integrator supplies the fixed-invariant return map once a caller
    !! supplies its physical launch map.  It does not, by itself, supply the
    !! Poincare component enumeration, P_phi(R_c), dP_phi/dR_c, drift-surface
    !! forces, or the phase-averaged perturbation.  Until those callbacks are
    !! supplied, the provider returns GC_CYL_NONLOCAL_UNAVAILABLE.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_gc_callback_context, only: gc_callback_context_t
    use neort_gc_nonlocal_resonance_types, only: &
        GC_NONLOCAL_MAX_FORCE_VALUES, &
        GC_NONLOCAL_SAMPLE_INVALID, &
        GC_NONLOCAL_SAMPLE_UNRESOLVED, GC_NONLOCAL_SAMPLE_VALID, &
        GC_NONLOCAL_SAMPLE_WALL, GC_NONLOCAL_SUCCESS, &
        gc_nonlocal_component_t, gc_nonlocal_orbit_sample_t, &
        gc_nonlocal_resonance_options_t, gc_nonlocal_resonance_result_t
    use neort_gc_nonlocal_transport_types, only: &
        GC_NONLOCAL_CLASS_COUNTERPASSING, GC_NONLOCAL_CLASS_COPASSING, &
        GC_NONLOCAL_CLASS_TRAPPED
    use neort_full_fow_action_symbolic, only: &
        evaluate_neort_action_normalization

    implicit none
    private

    integer, parameter, public :: GC_CYL_NONLOCAL_SUCCESS = 0
    integer, parameter, public :: GC_CYL_NONLOCAL_INVALID_INPUT = 100
    integer, parameter, public :: GC_CYL_NONLOCAL_UNAVAILABLE = 101
    integer, parameter, public :: GC_CYL_NONLOCAL_CALLBACK_FAILURE = 102
    integer, parameter, public :: GC_CYL_NONLOCAL_COMPONENT_IDENTITY = 103
    integer, parameter, public :: GC_CYL_NONLOCAL_ORBIT_ERROR = 104
    integer, parameter, public :: GC_CYL_NONLOCAL_DERIVATIVE_MISSING = 105
    integer, parameter, public :: GC_CYL_NONLOCAL_NONFINITE = 106
    integer, parameter, public :: GC_CYL_NONLOCAL_FORCE_ERROR = 107
    integer, parameter, public :: GC_CYL_NONLOCAL_SECTION_ERROR = 108
    integer, parameter, public :: GC_CYL_NONLOCAL_NO_COMPONENTS = 109
    integer, parameter, public :: GC_CYL_NONLOCAL_REFERENCE_MISMATCH = 110
    integer, parameter, public :: GC_CYL_NONLOCAL_MAPPING_UNAVAILABLE = 111
    integer, parameter, public :: GC_CYL_NONLOCAL_INTEGRAL_UNAVAILABLE = 112
    integer, parameter, public :: GC_CYL_NONLOCAL_NORMALIZATION_UNAVAILABLE = 113

    integer, parameter, public :: GC_CYL_NONLOCAL_PHASE_REPRESENTATION = 1
    integer, parameter, public :: GC_CYL_NONLOCAL_FREQUENCY_REPRESENTATION = 2

    integer, parameter, public :: GC_CYL_NONLOCAL_ORBIT_VALID = 0
    integer, parameter, public :: GC_CYL_NONLOCAL_ORBIT_WALL = 1
    integer, parameter, public :: GC_CYL_NONLOCAL_ORBIT_UNRESOLVED = 2
    integer, parameter, public :: GC_CYL_NONLOCAL_ORBIT_ERROR_STATUS = 3

    integer, parameter, public :: GC_CYL_NONLOCAL_WALL_CLEAR = 0
    integer, parameter, public :: GC_CYL_NONLOCAL_WALL_HIT = 1
    integer, parameter, public :: GC_CYL_NONLOCAL_WALL_UNRESOLVED = 2

    character(len=*), parameter, public :: &
        GC_CYL_NONLOCAL_DEFAULT_SECTION = 'R_c'

    type, public :: gc_cylindrical_nonlocal_section_t
        !! Locked physical reference for the Poincare section.
        !!
        !! reference is ordered (R,Z,phi), matching the cylindrical backend.
        !! It is metadata, not an additional dynamical coordinate.
        character(len=32) :: coordinate = GC_CYL_NONLOCAL_DEFAULT_SECTION
        real(dp) :: reference(3) = 0.0_dp
        character(len=64) :: reference_id = ''
        logical :: locked = .false.
        integer :: required_return_crossings = 1
        integer :: return_crossings = 0
    end type gc_cylindrical_nonlocal_section_t

    type, public :: gc_cylindrical_nonlocal_orbit_t
        !! Physical output of the fixed-(H0,J_perp,R_c,sigma) orbit callback.
        !!
        !! p_phi is the dimensional canonical toroidal momentum P_phi.  The
        !! provider never copies it directly into psi_star: it applies the
        !! locked psi_star=(c/q)*p_phi conversion, or one certified callback,
        !! exactly once.
        integer :: status = GC_CYL_NONLOCAL_ORBIT_ERROR_STATUS
        integer :: component_id = 0
        integer :: sigma = 0
        integer :: winding = 0
        integer :: class_kind = 0
        integer :: parallel_sign_changes = 0
        integer :: section_return_crossings = 0
        integer :: intersection_orientations(2) = 0
        real(dp) :: intersection_times(2) = 0.0_dp
        real(dp) :: intersection_rates(2) = 0.0_dp
        logical :: intersection_multiplicity_certified = .false.
        logical :: winding_available = .false.
        logical :: section_return_available = .false.
        logical :: complete_cycle_return = .false.
        logical :: class_behavior_certified = .false.
        logical :: p_phi_mapping_certified = .false.
        integer :: mapping_orientation = 0
        logical :: electric_potential_included = .false.
        real(dp) :: delta_phi_b = 0.0_dp
        real(dp) :: ddelta_phi_b_dx = 0.0_dp
        logical :: phase_advance_available = .false.
        logical :: period_derivative_available = .false.
        logical :: wall_checked = .false.
        integer :: wall_status = GC_CYL_NONLOCAL_WALL_UNRESOLVED
        type(gc_cylindrical_nonlocal_section_t) :: section
        real(dp) :: p_phi = 0.0_dp
        real(dp) :: dp_phi_dx = 0.0_dp
        real(dp) :: tau_b = 0.0_dp
        real(dp) :: dtau_b_dx = 0.0_dp
        real(dp) :: omega_b = 0.0_dp
        real(dp) :: omega_phi = 0.0_dp
        real(dp) :: domega_b_dx = 0.0_dp
        real(dp) :: domega_phi_dx = 0.0_dp
        logical :: derivatives_available = .false.
    end type gc_cylindrical_nonlocal_orbit_t

    type, public :: gc_cylindrical_nonlocal_evaluation_t
        !! Provider output including topology and wall identity which the
        !! existing nonlocal sample type intentionally does not carry.
        type(gc_nonlocal_orbit_sample_t) :: sample
        integer :: winding = 0
        integer :: class_kind = 0
        integer :: parallel_sign_changes = 0
        integer :: section_return_crossings = 0
        integer :: intersection_orientations(2) = 0
        real(dp) :: intersection_times(2) = 0.0_dp
        real(dp) :: intersection_rates(2) = 0.0_dp
        logical :: intersection_multiplicity_certified = .false.
        integer :: wall_status = GC_CYL_NONLOCAL_WALL_UNRESOLVED
        logical :: winding_available = .false.
        logical :: section_return_available = .false.
        logical :: complete_cycle_return = .false.
        logical :: class_behavior_certified = .false.
        logical :: p_phi_mapping_certified = .false.
        integer :: mapping_orientation = 0
        logical :: electric_potential_included = .false.
        real(dp) :: p_phi = 0.0_dp
        real(dp) :: dp_phi_dx = 0.0_dp
        real(dp) :: dtau_b_dx = 0.0_dp
        real(dp) :: delta_phi_b = 0.0_dp
        real(dp) :: ddelta_phi_b_dx = 0.0_dp
        logical :: phase_advance_available = .false.
        logical :: period_derivative_available = .false.
        integer :: frequency_representation = &
            GC_CYL_NONLOCAL_FREQUENCY_REPRESENTATION
        logical :: canonical_normalization_certified = .false.
        character(len=64) :: psi_star_units = ''
        character(len=64) :: psi_star_definition = ''
        logical :: wall_checked = .false.
        type(gc_cylindrical_nonlocal_section_t) :: section
    end type gc_cylindrical_nonlocal_evaluation_t

    abstract interface
        subroutine gc_cylindrical_nonlocal_component_provider_i(h0, jperp, &
                user_data, components, status)
            import :: dp, gc_callback_context_t, gc_nonlocal_component_t
            real(dp), intent(in) :: h0, jperp
            class(gc_callback_context_t), pointer, intent(inout) :: user_data
            type(gc_nonlocal_component_t), allocatable, intent(out) :: components(:)
            integer, intent(out) :: status
        end subroutine gc_cylindrical_nonlocal_component_provider_i

        subroutine gc_cylindrical_nonlocal_orbit_provider_i(h0, jperp, x, &
                sigma, component_id, user_data, orbit, status)
            import :: dp, gc_callback_context_t, &
                gc_cylindrical_nonlocal_orbit_t
            real(dp), intent(in) :: h0, jperp, x
            integer, intent(in) :: sigma, component_id
            class(gc_callback_context_t), pointer, intent(inout) :: user_data
            type(gc_cylindrical_nonlocal_orbit_t), intent(out) :: orbit
            integer, intent(out) :: status
        end subroutine gc_cylindrical_nonlocal_orbit_provider_i

        subroutine gc_cylindrical_nonlocal_harmonic_provider_i(h0, jperp, x, &
                sigma, component_id, harmonic_m, harmonic_n, orbit, user_data, &
                h_m, status)
            import :: dp, gc_callback_context_t, &
                gc_cylindrical_nonlocal_orbit_t
            real(dp), intent(in) :: h0, jperp, x
            integer, intent(in) :: sigma, component_id, harmonic_m, harmonic_n
            type(gc_cylindrical_nonlocal_orbit_t), intent(in) :: orbit
            class(gc_callback_context_t), pointer, intent(inout) :: user_data
            complex(dp), intent(out) :: h_m
            integer, intent(out) :: status
        end subroutine gc_cylindrical_nonlocal_harmonic_provider_i

        subroutine gc_cylindrical_nonlocal_force_provider_i(h0, jperp, x, &
                sigma, component_id, orbit, user_data, force, status)
            import :: dp, gc_callback_context_t, &
                gc_cylindrical_nonlocal_orbit_t
            real(dp), intent(in) :: h0, jperp, x
            integer, intent(in) :: sigma, component_id
            type(gc_cylindrical_nonlocal_orbit_t), intent(in) :: orbit
            class(gc_callback_context_t), pointer, intent(inout) :: user_data
            real(dp), intent(out) :: force(:)
            integer, intent(out) :: status
        end subroutine gc_cylindrical_nonlocal_force_provider_i

        subroutine gc_cylindrical_nonlocal_canonical_conversion_i(p_phi, &
                dp_phi_dx, user_data, psi_star, dpsi_star_dx, units, certified, &
                status)
            import :: dp, gc_callback_context_t
            real(dp), intent(in) :: p_phi, dp_phi_dx
            class(gc_callback_context_t), pointer, intent(inout) :: user_data
            real(dp), intent(out) :: psi_star, dpsi_star_dx
            character(len=*), intent(out) :: units
            logical, intent(out) :: certified
            integer, intent(out) :: status
        end subroutine gc_cylindrical_nonlocal_canonical_conversion_i
    end interface

    type, public :: gc_cylindrical_nonlocal_context_t
        !! Immutable physical invariants and callback-owned cylindrical seams.
        real(dp) :: h0 = 0.0_dp
        real(dp) :: jperp = 0.0_dp
        real(dp) :: particle_charge = 0.0_dp
        real(dp) :: c_light = 0.0_dp
        logical :: charge_c_locked = .false.
        type(gc_cylindrical_nonlocal_section_t) :: section
        logical :: initialized = .false.
        logical :: components_enumerated = .false.
        class(gc_callback_context_t), pointer :: user_data => null()
        procedure(gc_cylindrical_nonlocal_component_provider_i), pointer, nopass :: &
            component_provider => null()
        procedure(gc_cylindrical_nonlocal_orbit_provider_i), pointer, nopass :: &
            orbit_provider => null()
        procedure(gc_cylindrical_nonlocal_harmonic_provider_i), pointer, nopass :: &
            harmonic_provider => null()
        procedure(gc_cylindrical_nonlocal_force_provider_i), pointer, nopass :: &
            force_provider => null()
        procedure(gc_cylindrical_nonlocal_canonical_conversion_i), pointer, &
            nopass :: canonical_conversion_provider => null()
        type(gc_nonlocal_component_t), allocatable :: components(:)
    end type gc_cylindrical_nonlocal_context_t

    public :: clear_gc_cylindrical_nonlocal_provider
    public :: initialize_gc_cylindrical_nonlocal_provider
    public :: enumerate_gc_cylindrical_nonlocal_components
    public :: evaluate_gc_cylindrical_nonlocal
    public :: evaluate_gc_cylindrical_nonlocal_sample
    public :: integrate_gc_cylindrical_nonlocal_resonance
    public :: gc_cylindrical_nonlocal_resonance_residuals
    public :: gc_cylindrical_nonlocal_resonance_weights
    public :: gc_cylindrical_nonlocal_component_provider_i
    public :: gc_cylindrical_nonlocal_orbit_provider_i
    public :: gc_cylindrical_nonlocal_harmonic_provider_i
    public :: gc_cylindrical_nonlocal_force_provider_i
    public :: gc_cylindrical_nonlocal_canonical_conversion_i

contains

    subroutine initialize_gc_cylindrical_nonlocal_provider(h0, jperp, context, &
            status, particle_charge, c_light, component_provider, orbit_provider, &
            harmonic_provider, force_provider, canonical_conversion_provider, &
            section_coordinate, section_reference, section_reference_id, user_data, &
            required_return_crossings)
        real(dp), intent(in) :: h0, jperp
        type(gc_cylindrical_nonlocal_context_t), intent(out) :: context
        integer, intent(out) :: status
        real(dp), intent(in), optional :: particle_charge, c_light
        procedure(gc_cylindrical_nonlocal_component_provider_i), optional :: &
            component_provider
        procedure(gc_cylindrical_nonlocal_orbit_provider_i), optional :: &
            orbit_provider
        procedure(gc_cylindrical_nonlocal_harmonic_provider_i), optional :: &
            harmonic_provider
        procedure(gc_cylindrical_nonlocal_force_provider_i), optional :: &
            force_provider
        procedure(gc_cylindrical_nonlocal_canonical_conversion_i), optional :: &
            canonical_conversion_provider
        character(len=*), intent(in), optional :: section_coordinate
        real(dp), intent(in), optional :: section_reference(3)
        character(len=*), intent(in), optional :: section_reference_id
        class(gc_callback_context_t), target, intent(inout), optional :: user_data
        integer, intent(in), optional :: required_return_crossings

        context%h0 = 0.0_dp
        context%jperp = 0.0_dp
        context%particle_charge = 0.0_dp
        context%c_light = 0.0_dp
        context%charge_c_locked = .false.
        context%section = gc_cylindrical_nonlocal_section_t()
        context%initialized = .false.
        context%components_enumerated = .false.
        nullify(context%user_data)
        nullify(context%component_provider)
        nullify(context%orbit_provider)
        nullify(context%harmonic_provider)
        nullify(context%force_provider)
        nullify(context%canonical_conversion_provider)
        if (allocated(context%components)) deallocate(context%components)
        if (present(required_return_crossings)) then
            if (required_return_crossings < 1 .or. required_return_crossings > 2) then
                status = GC_CYL_NONLOCAL_SECTION_ERROR
                return
            end if
            context%section%required_return_crossings = required_return_crossings
            context%section%return_crossings = required_return_crossings
        end if

        status = GC_CYL_NONLOCAL_INVALID_INPUT
        if (.not. ieee_is_finite(h0)) return
        if (.not. ieee_is_finite(jperp)) return
        if (jperp < 0.0_dp) return
        if (present(particle_charge)) then
            if (.not. present(c_light)) then
                status = GC_CYL_NONLOCAL_NORMALIZATION_UNAVAILABLE
                return
            end if
            if (.not. ieee_is_finite(particle_charge)) return
            if (.not. ieee_is_finite(c_light)) return
            if (abs(particle_charge) <= tiny(particle_charge) .or. &
                c_light <= 0.0_dp) then
                status = GC_CYL_NONLOCAL_NORMALIZATION_UNAVAILABLE
                return
            end if
            context%particle_charge = particle_charge
            context%c_light = c_light
            context%charge_c_locked = .true.
        else
            if (present(c_light)) then
                status = GC_CYL_NONLOCAL_NORMALIZATION_UNAVAILABLE
                return
            end if
        end if
        if (present(section_coordinate)) then
            if (.not. valid_section_coordinate(section_coordinate)) then
                status = GC_CYL_NONLOCAL_SECTION_ERROR
                return
            end if
            context%section%coordinate = adjustl(section_coordinate)
        end if
        if (present(section_reference)) then
            if (.not. all(ieee_is_finite(section_reference))) return
            if (.not. present(section_reference_id)) then
                status = GC_CYL_NONLOCAL_SECTION_ERROR
                return
            end if
            if (len_trim(section_reference_id) == 0) then
                status = GC_CYL_NONLOCAL_SECTION_ERROR
                return
            end if
            context%section%reference = section_reference
            context%section%reference_id = adjustl(section_reference_id)
            context%section%locked = .true.
        else
            if (present(section_reference_id)) then
                status = GC_CYL_NONLOCAL_SECTION_ERROR
                return
            end if
        end if

        context%h0 = h0
        context%jperp = jperp
        if (present(user_data)) context%user_data => user_data
        if (present(component_provider)) then
            context%component_provider => component_provider
        end if
        if (present(orbit_provider)) then
            context%orbit_provider => orbit_provider
        end if
        if (present(harmonic_provider)) then
            context%harmonic_provider => harmonic_provider
        end if
        if (present(force_provider)) then
            context%force_provider => force_provider
        end if
        if (present(canonical_conversion_provider)) then
            context%canonical_conversion_provider => canonical_conversion_provider
        end if
        context%initialized = .true.
        status = GC_CYL_NONLOCAL_SUCCESS
    end subroutine initialize_gc_cylindrical_nonlocal_provider

    subroutine clear_gc_cylindrical_nonlocal_provider(context)
        type(gc_cylindrical_nonlocal_context_t), intent(inout) :: context

        if (allocated(context%components)) deallocate(context%components)
        nullify(context%user_data)
        nullify(context%component_provider)
        nullify(context%orbit_provider)
        nullify(context%harmonic_provider)
        nullify(context%force_provider)
        nullify(context%canonical_conversion_provider)
        context%h0 = 0.0_dp
        context%jperp = 0.0_dp
        context%particle_charge = 0.0_dp
        context%c_light = 0.0_dp
        context%charge_c_locked = .false.
        context%section = gc_cylindrical_nonlocal_section_t()
        context%initialized = .false.
        context%components_enumerated = .false.
    end subroutine clear_gc_cylindrical_nonlocal_provider

    subroutine enumerate_gc_cylindrical_nonlocal_components(context, components, &
            status)
        type(gc_cylindrical_nonlocal_context_t), intent(inout) :: context
        type(gc_nonlocal_component_t), allocatable, intent(out) :: components(:)
        integer, intent(out) :: status

        type(gc_nonlocal_component_t), allocatable :: proposed(:)
        integer :: callback_status

        if (allocated(components)) deallocate(components)
        status = GC_CYL_NONLOCAL_INVALID_INPUT
        if (.not. context%initialized) return
        if (.not. associated(context%component_provider)) then
            status = GC_CYL_NONLOCAL_UNAVAILABLE
            return
        end if

        call context%component_provider(context%h0, context%jperp, &
            context%user_data, proposed, callback_status)
        if (callback_status /= GC_CYL_NONLOCAL_SUCCESS) then
            status = GC_CYL_NONLOCAL_CALLBACK_FAILURE
            return
        end if
        if (.not. allocated(proposed)) then
            status = GC_CYL_NONLOCAL_UNAVAILABLE
            return
        end if
        status = validate_components(proposed)
        if (status /= GC_CYL_NONLOCAL_SUCCESS) return

        allocate(context%components(size(proposed)))
        context%components = proposed
        context%components_enumerated = .true.
        allocate(components(size(proposed)))
        components = proposed
        status = GC_CYL_NONLOCAL_SUCCESS
    end subroutine enumerate_gc_cylindrical_nonlocal_components

    subroutine evaluate_gc_cylindrical_nonlocal(context, x, sigma, component_id, &
            harmonic_m, harmonic_n, nforce, evaluation, status)
        type(gc_cylindrical_nonlocal_context_t), intent(inout) :: context
        real(dp), intent(in) :: x
        integer, intent(in) :: sigma, component_id, harmonic_m, harmonic_n, nforce
        type(gc_cylindrical_nonlocal_evaluation_t), intent(out) :: evaluation
        integer, intent(out) :: status

        type(gc_nonlocal_component_t) :: component
        type(gc_cylindrical_nonlocal_orbit_t) :: orbit
        real(dp) :: force(GC_NONLOCAL_MAX_FORCE_VALUES)
        complex(dp) :: h_m
        integer :: component_index, callback_status, local_status

        evaluation = gc_cylindrical_nonlocal_evaluation_t()
        evaluation%sample = gc_nonlocal_orbit_sample_t()
        h_m = cmplx(0.0_dp, 0.0_dp, kind=dp)
        force = 0.0_dp
        status = GC_CYL_NONLOCAL_INVALID_INPUT
        if (.not. context%initialized) return
        if (.not. ieee_is_finite(x)) return
        if (abs(sigma) /= 1) return
        if (nforce < 1) return
        if (nforce > size(force)) return
        if (.not. context%components_enumerated) then
            status = GC_CYL_NONLOCAL_UNAVAILABLE
            return
        end if
        call locate_component(context%components, x, sigma, component_id, &
            component_index, component, local_status)
        if (local_status /= GC_CYL_NONLOCAL_SUCCESS) then
            status = local_status
            return
        end if
        associate (unused_component_index => component_index)
        end associate
        call initialize_evaluation_identity(evaluation, component, context%section)

        if (.not. associated(context%orbit_provider)) then
            status = GC_CYL_NONLOCAL_UNAVAILABLE
            return
        end if
        call context%orbit_provider(context%h0, context%jperp, x, sigma, &
            component_id, context%user_data, orbit, callback_status)
        if (callback_status /= GC_CYL_NONLOCAL_SUCCESS) then
            status = GC_CYL_NONLOCAL_CALLBACK_FAILURE
            return
        end if
        call copy_orbit_metadata(orbit, evaluation)
        if (orbit%component_id /= component_id) then
            status = GC_CYL_NONLOCAL_COMPONENT_IDENTITY
            return
        end if
        if (orbit%sigma /= sigma) then
            status = GC_CYL_NONLOCAL_COMPONENT_IDENTITY
            return
        end if
        call validate_section_metadata(orbit%section, context%section, local_status)
        if (local_status /= GC_CYL_NONLOCAL_SUCCESS) then
            status = local_status
            return
        end if
        if (orbit%status == GC_CYL_NONLOCAL_ORBIT_WALL) then
            if (.not. orbit%wall_checked) then
                status = GC_CYL_NONLOCAL_UNAVAILABLE
                return
            end if
            if (orbit%wall_status /= GC_CYL_NONLOCAL_WALL_HIT) then
                status = GC_CYL_NONLOCAL_ORBIT_ERROR
                return
            end if
            evaluation%sample%status = GC_NONLOCAL_SAMPLE_WALL
            status = GC_CYL_NONLOCAL_SUCCESS
            return
        end if
        if (orbit%status == GC_CYL_NONLOCAL_ORBIT_UNRESOLVED) then
            evaluation%sample%status = GC_NONLOCAL_SAMPLE_UNRESOLVED
            status = GC_CYL_NONLOCAL_SUCCESS
            return
        end if
        if (orbit%status /= GC_CYL_NONLOCAL_ORBIT_VALID) then
            status = GC_CYL_NONLOCAL_ORBIT_ERROR
            return
        end if
        call validate_valid_orbit(orbit, context%section%required_return_crossings, &
            local_status)
        if (local_status /= GC_CYL_NONLOCAL_SUCCESS) then
            status = local_status
            return
        end if

        evaluation%sample%status = GC_NONLOCAL_SAMPLE_VALID
        evaluation%sample%component_id = component_id
        evaluation%sample%sigma = sigma
        evaluation%sample%class_kind = orbit%class_kind
        evaluation%sample%class_behavior_certified = &
            orbit%class_behavior_certified
        evaluation%p_phi = orbit%p_phi
        evaluation%dp_phi_dx = orbit%dp_phi_dx
        evaluation%sample%tau_b = orbit%tau_b
        evaluation%sample%omega_b = orbit%omega_b
        evaluation%sample%omega_phi = orbit%omega_phi
        evaluation%sample%domega_b_dx = orbit%domega_b_dx
        evaluation%sample%domega_phi_dx = orbit%domega_phi_dx
        evaluation%sample%nforce = nforce
        evaluation%sample%derivatives_available = orbit%derivatives_available

        call convert_canonical_coordinate(context, orbit%p_phi, orbit%dp_phi_dx, &
            evaluation%sample%psi_star, evaluation%sample%dpsi_star_dx, &
            evaluation%psi_star_units, evaluation%psi_star_definition, &
            evaluation%canonical_normalization_certified, local_status)
        if (local_status /= GC_CYL_NONLOCAL_SUCCESS) then
            evaluation%sample%status = GC_NONLOCAL_SAMPLE_INVALID
            status = local_status
            return
        end if

        if (.not. associated(context%harmonic_provider)) then
            evaluation%sample%status = GC_NONLOCAL_SAMPLE_INVALID
            status = GC_CYL_NONLOCAL_UNAVAILABLE
            return
        end if
        call context%harmonic_provider(context%h0, context%jperp, x, sigma, &
            component_id, harmonic_m, harmonic_n, orbit, context%user_data, h_m, &
            callback_status)
        if (callback_status /= GC_CYL_NONLOCAL_SUCCESS) then
            evaluation%sample%status = GC_NONLOCAL_SAMPLE_INVALID
            status = GC_CYL_NONLOCAL_CALLBACK_FAILURE
            return
        end if
        if (.not. all(ieee_is_finite([real(h_m, dp), aimag(h_m)]))) then
            evaluation%sample%status = GC_NONLOCAL_SAMPLE_INVALID
            status = GC_CYL_NONLOCAL_NONFINITE
            return
        end if
        evaluation%sample%h_m = h_m

        if (.not. associated(context%force_provider)) then
            evaluation%sample%status = GC_NONLOCAL_SAMPLE_INVALID
            status = GC_CYL_NONLOCAL_UNAVAILABLE
            return
        end if
        call context%force_provider(context%h0, context%jperp, x, sigma, &
            component_id, orbit, context%user_data, force(1:nforce), &
            callback_status)
        if (callback_status /= GC_CYL_NONLOCAL_SUCCESS) then
            evaluation%sample%status = GC_NONLOCAL_SAMPLE_INVALID
            status = GC_CYL_NONLOCAL_FORCE_ERROR
            return
        end if
        if (.not. all(ieee_is_finite(force(1:nforce)))) then
            evaluation%sample%status = GC_NONLOCAL_SAMPLE_INVALID
            status = GC_CYL_NONLOCAL_NONFINITE
            return
        end if
        evaluation%sample%thermodynamic_force(1:nforce) = force(1:nforce)
        status = GC_CYL_NONLOCAL_SUCCESS
    end subroutine evaluate_gc_cylindrical_nonlocal

    subroutine evaluate_gc_cylindrical_nonlocal_sample(context, x, sigma, &
            component_id, harmonic_m, harmonic_n, nforce, sample, status)
        type(gc_cylindrical_nonlocal_context_t), intent(inout) :: context
        real(dp), intent(in) :: x
        integer, intent(in) :: sigma, component_id, harmonic_m, harmonic_n, nforce
        type(gc_nonlocal_orbit_sample_t), intent(out) :: sample
        integer, intent(out) :: status

        type(gc_cylindrical_nonlocal_evaluation_t) :: evaluation

        sample = gc_nonlocal_orbit_sample_t()
        call evaluate_gc_cylindrical_nonlocal(context, x, sigma, component_id, &
            harmonic_m, harmonic_n, nforce, evaluation, status)
        sample = evaluation%sample
    end subroutine evaluate_gc_cylindrical_nonlocal_sample

    subroutine integrate_gc_cylindrical_nonlocal_resonance(context, harmonic_m, &
            harmonic_n, options, result, status)
        !! Deliberately keep this bounded provider tranche fail closed.  The
        !! repository's corrected inner kernel uses the frequency-delta
        !! tau_b/|d(m*omega_b+n*omega_phi)/dx| weight; this seam does not
        !! duplicate that factor or claim that physical callbacks are present.
        !! A production integration adapter belongs after component, orbit,
        !! phase-average, force, and outer-measure callbacks are supplied.
        type(gc_cylindrical_nonlocal_context_t), intent(inout) :: context
        integer, intent(in) :: harmonic_m, harmonic_n
        type(gc_nonlocal_resonance_options_t), intent(in) :: options
        type(gc_nonlocal_resonance_result_t), intent(out) :: result
        integer, intent(out) :: status

        result = gc_nonlocal_resonance_result_t()
        result%h0 = context%h0
        result%jperp = context%jperp
        result%harmonic_m = harmonic_m
        result%harmonic_n = harmonic_n
        associate (unused_context => context, unused_options => options)
        end associate
        status = GC_CYL_NONLOCAL_INTEGRAL_UNAVAILABLE
        result%status = status
        result%certified = .false.
    end subroutine integrate_gc_cylindrical_nonlocal_resonance

    subroutine gc_cylindrical_nonlocal_resonance_residuals(evaluation, harmonic_m, &
            harmonic_n, phase_residual, phase_derivative, frequency_residual, &
            frequency_derivative, status)
        !! Return both equivalent resonance representations.
        !!
        !! g = Delta_phi_b + 2*pi*m/n
        !! F = m*omega_b + n*omega_phi = n*g/tau_b.
        !! The derivative identity is used at a resonance root; away from a
        !! root the derivative of tau_b contributes to F'.
        type(gc_cylindrical_nonlocal_evaluation_t), intent(in) :: evaluation
        integer, intent(in) :: harmonic_m, harmonic_n
        real(dp), intent(out) :: phase_residual, phase_derivative
        real(dp), intent(out) :: frequency_residual, frequency_derivative
        integer, intent(out) :: status

        real(dp), parameter :: two_pi = 6.28318530717958647692528676656_dp

        phase_residual = 0.0_dp
        phase_derivative = 0.0_dp
        frequency_residual = 0.0_dp
        frequency_derivative = 0.0_dp
        status = GC_CYL_NONLOCAL_INVALID_INPUT
        if (harmonic_n == 0) return
        if (evaluation%sample%status /= GC_NONLOCAL_SAMPLE_VALID) return
        if (.not. evaluation%phase_advance_available) return
        if (evaluation%frequency_representation &
            /= GC_CYL_NONLOCAL_FREQUENCY_REPRESENTATION) return
        if (evaluation%sample%tau_b <= 0.0_dp) return
        phase_residual = evaluation%delta_phi_b + two_pi*real(harmonic_m, dp) &
            /real(harmonic_n, dp)
        phase_derivative = evaluation%ddelta_phi_b_dx
        frequency_residual = real(harmonic_m, dp)*evaluation%sample%omega_b &
            +real(harmonic_n, dp)*evaluation%sample%omega_phi
        frequency_derivative = real(harmonic_m, dp) &
            *evaluation%sample%domega_b_dx &
            +real(harmonic_n, dp)*evaluation%sample%domega_phi_dx
        if (.not. all(ieee_is_finite([phase_residual, phase_derivative, &
            frequency_residual, frequency_derivative]))) then
            phase_residual = 0.0_dp
            phase_derivative = 0.0_dp
            frequency_residual = 0.0_dp
            frequency_derivative = 0.0_dp
            status = GC_CYL_NONLOCAL_NONFINITE
            return
        end if
        status = GC_CYL_NONLOCAL_SUCCESS
    end subroutine gc_cylindrical_nonlocal_resonance_residuals

    subroutine gc_cylindrical_nonlocal_resonance_weights(evaluation, harmonic_m, &
            harmonic_n, phase_weight, frequency_inner_weight, frequency_prefactor, &
            status)
        !! Convert the phase-delta Eq. 17 weight to the frequency-delta form.
        !!
        !! At a simple root, Eq. 17 has
        !!   |n|*tau_b**2/|g'|,
        !! while the equivalent frequency form has
        !!   n**2 * [tau_b/|F'|].
        !! The returned frequency_prefactor is n**2 and the returned inner
        !! weight contains one, not two, powers of tau_b.
        type(gc_cylindrical_nonlocal_evaluation_t), intent(in) :: evaluation
        integer, intent(in) :: harmonic_m, harmonic_n
        real(dp), intent(out) :: phase_weight, frequency_inner_weight
        integer, intent(out) :: frequency_prefactor
        integer, intent(out) :: status

        real(dp) :: phase_residual, phase_derivative
        real(dp) :: frequency_residual, frequency_derivative
        real(dp), parameter :: root_tolerance = 1.0e-10_dp

        phase_weight = 0.0_dp
        frequency_inner_weight = 0.0_dp
        frequency_prefactor = 0
        call gc_cylindrical_nonlocal_resonance_residuals(evaluation, harmonic_m, &
            harmonic_n, phase_residual, phase_derivative, frequency_residual, &
            frequency_derivative, status)
        if (status /= GC_CYL_NONLOCAL_SUCCESS) return
        if (abs(phase_residual) > root_tolerance) then
            status = GC_CYL_NONLOCAL_SECTION_ERROR
            return
        end if
        if (abs(phase_derivative) <= tiny(phase_derivative)) then
            status = GC_CYL_NONLOCAL_DERIVATIVE_MISSING
            return
        end if
        if (abs(frequency_derivative) <= tiny(frequency_derivative)) then
            status = GC_CYL_NONLOCAL_DERIVATIVE_MISSING
            return
        end if
        phase_weight = abs(real(harmonic_n, dp))*evaluation%sample%tau_b**2 &
            /abs(phase_derivative)
        frequency_inner_weight = evaluation%sample%tau_b/abs(frequency_derivative)
        frequency_prefactor = harmonic_n*harmonic_n
        status = GC_CYL_NONLOCAL_SUCCESS
    end subroutine gc_cylindrical_nonlocal_resonance_weights

    subroutine initialize_evaluation_identity(evaluation, component, section)
        type(gc_cylindrical_nonlocal_evaluation_t), intent(inout) :: evaluation
        type(gc_nonlocal_component_t), intent(in) :: component
        type(gc_cylindrical_nonlocal_section_t), intent(in) :: section

        evaluation%sample%component_id = component%component_id
        evaluation%sample%sigma = component%sigma
        evaluation%wall_status = GC_CYL_NONLOCAL_WALL_UNRESOLVED
        evaluation%section = section
    end subroutine initialize_evaluation_identity

    subroutine copy_orbit_metadata(orbit, evaluation)
        type(gc_cylindrical_nonlocal_orbit_t), intent(in) :: orbit
        type(gc_cylindrical_nonlocal_evaluation_t), intent(inout) :: evaluation

        evaluation%winding = orbit%winding
        evaluation%class_kind = orbit%class_kind
        evaluation%parallel_sign_changes = orbit%parallel_sign_changes
        evaluation%section_return_crossings = orbit%section_return_crossings
        evaluation%intersection_orientations = orbit%intersection_orientations
        evaluation%intersection_times = orbit%intersection_times
        evaluation%intersection_rates = orbit%intersection_rates
        evaluation%intersection_multiplicity_certified = &
            orbit%intersection_multiplicity_certified
        evaluation%wall_status = orbit%wall_status
        evaluation%winding_available = orbit%winding_available
        evaluation%section_return_available = orbit%section_return_available
        evaluation%complete_cycle_return = orbit%complete_cycle_return
        evaluation%class_behavior_certified = orbit%class_behavior_certified
        evaluation%p_phi_mapping_certified = orbit%p_phi_mapping_certified
        evaluation%mapping_orientation = orbit%mapping_orientation
        evaluation%electric_potential_included = orbit%electric_potential_included
        evaluation%p_phi = orbit%p_phi
        evaluation%dp_phi_dx = orbit%dp_phi_dx
        evaluation%dtau_b_dx = orbit%dtau_b_dx
        evaluation%delta_phi_b = orbit%delta_phi_b
        evaluation%ddelta_phi_b_dx = orbit%ddelta_phi_b_dx
        evaluation%phase_advance_available = orbit%phase_advance_available
        evaluation%period_derivative_available = orbit%period_derivative_available
        evaluation%wall_checked = orbit%wall_checked
        evaluation%section = orbit%section
    end subroutine copy_orbit_metadata

    subroutine convert_canonical_coordinate(context, p_phi, dp_phi_dx, psi_star, &
            dpsi_star_dx, units, definition, certified, status)
        type(gc_cylindrical_nonlocal_context_t), intent(inout) :: context
        real(dp), intent(in) :: p_phi, dp_phi_dx
        real(dp), intent(out) :: psi_star, dpsi_star_dx
        character(len=*), intent(out) :: units, definition
        logical, intent(out) :: certified
        integer, intent(out) :: status

        real(dp) :: action_values(11)
        integer :: callback_status

        psi_star = 0.0_dp
        dpsi_star_dx = 0.0_dp
        units = ''
        definition = ''
        certified = .false.
        status = GC_CYL_NONLOCAL_NORMALIZATION_UNAVAILABLE
        if (.not. all(ieee_is_finite([p_phi, dp_phi_dx]))) then
            status = GC_CYL_NONLOCAL_NONFINITE
            return
        end if
        if (associated(context%canonical_conversion_provider)) then
            call context%canonical_conversion_provider(p_phi, dp_phi_dx, &
                context%user_data, psi_star, dpsi_star_dx, units, certified, &
                callback_status)
            if (callback_status /= GC_CYL_NONLOCAL_SUCCESS) then
                status = GC_CYL_NONLOCAL_CALLBACK_FAILURE
                return
            end if
            if (.not. certified) return
            if (len_trim(units) == 0) return
            definition = 'certified canonical conversion callback'
        else
            if (.not. context%charge_c_locked) return
            call evaluate_neort_action_normalization(1.0_dp, &
                context%particle_charge, context%c_light, 0.0_dp, 1.0_dp, &
                0.0_dp, 0.0_dp, p_phi, 1.0_dp, action_values(1), &
                action_values(2), action_values(3), action_values(4), &
                action_values(5), action_values(6), action_values(7), &
                action_values(8), action_values(9), action_values(10), &
                action_values(11))
            if (.not. all(ieee_is_finite(action_values))) return
            psi_star = action_values(10)
            dpsi_star_dx = action_values(11)*dp_phi_dx
            units = 'psi_star=(c/q)*p_phi'
            definition = 'Fortsym psi_star=(c/q)*p_phi'
            certified = .true.
        end if
        if (.not. all(ieee_is_finite([psi_star, dpsi_star_dx]))) then
            psi_star = 0.0_dp
            dpsi_star_dx = 0.0_dp
            certified = .false.
            status = GC_CYL_NONLOCAL_NONFINITE
            return
        end if
        status = GC_CYL_NONLOCAL_SUCCESS
    end subroutine convert_canonical_coordinate

    subroutine validate_valid_orbit(orbit, required_return_crossings, status)
        type(gc_cylindrical_nonlocal_orbit_t), intent(in) :: orbit
        integer, intent(in) :: required_return_crossings
        integer, intent(out) :: status

        real(dp), parameter :: two_pi = 6.28318530717958647692528676656_dp
        real(dp), parameter :: consistency_tolerance = 1.0e-10_dp
        real(dp) :: values(10), expected(4)

        status = GC_CYL_NONLOCAL_INVALID_INPUT
        if (.not. orbit%wall_checked) then
            status = GC_CYL_NONLOCAL_UNAVAILABLE
            return
        end if
        if (orbit%wall_status /= GC_CYL_NONLOCAL_WALL_CLEAR) then
            status = GC_CYL_NONLOCAL_ORBIT_ERROR
            return
        end if
        if (.not. orbit%winding_available) then
            status = GC_CYL_NONLOCAL_UNAVAILABLE
            return
        end if
        if (.not. orbit%section_return_available) then
            status = GC_CYL_NONLOCAL_UNAVAILABLE
            return
        end if
        if (.not. orbit%phase_advance_available) then
            status = GC_CYL_NONLOCAL_UNAVAILABLE
            return
        end if
        if (.not. orbit%period_derivative_available) then
            status = GC_CYL_NONLOCAL_UNAVAILABLE
            return
        end if
        if (.not. orbit%complete_cycle_return) then
            status = GC_CYL_NONLOCAL_UNAVAILABLE
            return
        end if
        ! class_kind, winding, and parallel-sign history are descriptive
        ! diagnostics only.  They never certify or split a cut component.
        ! Completeness is established by the exhaustive certified cut atlas
        ! and its boundary/intersection evidence above.
        if (.not. orbit%p_phi_mapping_certified) then
            status = GC_CYL_NONLOCAL_MAPPING_UNAVAILABLE
            return
        end if
        if (abs(orbit%mapping_orientation) /= 1) then
            status = GC_CYL_NONLOCAL_MAPPING_UNAVAILABLE
            return
        end if
        if (.not. orbit%electric_potential_included) then
            status = GC_CYL_NONLOCAL_UNAVAILABLE
            return
        end if
        if (required_return_crossings < 1 .or. required_return_crossings > 2) then
            status = GC_CYL_NONLOCAL_SECTION_ERROR
            return
        end if
        if (orbit%section_return_crossings /= required_return_crossings) then
            status = GC_CYL_NONLOCAL_SECTION_ERROR
            return
        end if
        if (required_return_crossings == 2) then
            if (.not. orbit%intersection_multiplicity_certified) then
                status = GC_CYL_NONLOCAL_SECTION_ERROR
                return
            end if
            if (abs(orbit%intersection_orientations(1)) /= 1 .or. &
                    abs(orbit%intersection_orientations(2)) /= 1 .or. &
                    orbit%intersection_orientations(1) == &
                    orbit%intersection_orientations(2) .or. &
                    orbit%intersection_times(2) <= orbit%intersection_times(1)) then
                status = GC_CYL_NONLOCAL_SECTION_ERROR
                return
            end if
            if (.not. all(ieee_is_finite(orbit%intersection_rates))) then
                status = GC_CYL_NONLOCAL_SECTION_ERROR
                return
            end if
            if (any(abs(orbit%intersection_rates) <= tiny(1.0_dp))) then
                status = GC_CYL_NONLOCAL_SECTION_ERROR
                return
            end if
        end if
        values = [orbit%p_phi, orbit%dp_phi_dx, orbit%tau_b, orbit%omega_b, &
            orbit%omega_phi, orbit%domega_b_dx, orbit%domega_phi_dx, &
            orbit%delta_phi_b, orbit%ddelta_phi_b_dx, orbit%dtau_b_dx]
        if (.not. all(ieee_is_finite(values))) then
            status = GC_CYL_NONLOCAL_NONFINITE
            return
        end if
        if (orbit%tau_b <= 0.0_dp) return
        if (.not. orbit%derivatives_available) then
            status = GC_CYL_NONLOCAL_DERIVATIVE_MISSING
            return
        end if
        expected = [two_pi/orbit%tau_b, orbit%delta_phi_b/orbit%tau_b, &
            -two_pi*orbit%dtau_b_dx/orbit%tau_b**2, &
            (orbit%ddelta_phi_b_dx*orbit%tau_b &
            -orbit%delta_phi_b*orbit%dtau_b_dx)/orbit%tau_b**2]
        if (.not. close_values(orbit%omega_b, expected(1), &
            consistency_tolerance)) then
            status = GC_CYL_NONLOCAL_ORBIT_ERROR
            return
        end if
        if (.not. close_values(orbit%omega_phi, expected(2), &
            consistency_tolerance)) then
            status = GC_CYL_NONLOCAL_ORBIT_ERROR
            return
        end if
        if (.not. close_values(orbit%domega_b_dx, expected(3), &
            consistency_tolerance)) then
            status = GC_CYL_NONLOCAL_ORBIT_ERROR
            return
        end if
        if (.not. close_values(orbit%domega_phi_dx, expected(4), &
            consistency_tolerance)) then
            status = GC_CYL_NONLOCAL_ORBIT_ERROR
            return
        end if
        status = GC_CYL_NONLOCAL_SUCCESS
    end subroutine validate_valid_orbit

    pure logical function close_values(value, reference, relative_tolerance)
        real(dp), intent(in) :: value, reference, relative_tolerance

        close_values = abs(value - reference) <= relative_tolerance*max(1.0_dp, &
            max(abs(value), abs(reference)))
    end function close_values

    subroutine validate_section_metadata(orbit_section, context_section, status)
        type(gc_cylindrical_nonlocal_section_t), intent(in) :: orbit_section
        type(gc_cylindrical_nonlocal_section_t), intent(in) :: context_section
        integer, intent(out) :: status

        real(dp), parameter :: reference_tolerance = 1.0e-12_dp
        real(dp) :: scale
        integer :: i

        status = GC_CYL_NONLOCAL_UNAVAILABLE
        if (.not. context_section%locked) return
        if (.not. orbit_section%locked) return
        if (trim(orbit_section%coordinate) /= trim(context_section%coordinate)) then
            status = GC_CYL_NONLOCAL_REFERENCE_MISMATCH
            return
        end if
        if (trim(orbit_section%reference_id) /= &
            trim(context_section%reference_id)) then
            status = GC_CYL_NONLOCAL_REFERENCE_MISMATCH
            return
        end if
        do i = 1, 3
            scale = max(1.0_dp, abs(context_section%reference(i)))
            if (abs(orbit_section%reference(i) - context_section%reference(i)) &
                > reference_tolerance*scale) then
                status = GC_CYL_NONLOCAL_REFERENCE_MISMATCH
                return
            end if
        end do
        if (orbit_section%required_return_crossings /= &
                context_section%required_return_crossings) then
            status = GC_CYL_NONLOCAL_SECTION_ERROR
            return
        end if
        status = GC_CYL_NONLOCAL_SUCCESS
    end subroutine validate_section_metadata

    subroutine locate_component(components, x, sigma, component_id, index, &
            component, status)
        type(gc_nonlocal_component_t), intent(in) :: components(:)
        real(dp), intent(in) :: x
        integer, intent(in) :: sigma, component_id
        integer, intent(out) :: index
        type(gc_nonlocal_component_t), intent(out) :: component
        integer, intent(out) :: status

        integer :: i

        index = 0
        component = gc_nonlocal_component_t()
        status = GC_CYL_NONLOCAL_COMPONENT_IDENTITY
        do i = 1, size(components)
            if (components(i)%component_id /= component_id) cycle
            if (components(i)%sigma /= sigma) cycle
            index = i
            component = components(i)
            if (x < component%x_min) then
                status = GC_CYL_NONLOCAL_SECTION_ERROR
                return
            end if
            if (x > component%x_max) then
                status = GC_CYL_NONLOCAL_SECTION_ERROR
                return
            end if
            status = GC_CYL_NONLOCAL_SUCCESS
            return
        end do
    end subroutine locate_component

    pure integer function validate_components(components) result(status)
        type(gc_nonlocal_component_t), intent(in) :: components(:)

        integer :: i, j

        status = GC_CYL_NONLOCAL_INVALID_INPUT
        if (size(components) < 1) then
            status = GC_CYL_NONLOCAL_NO_COMPONENTS
            return
        end if
        do i = 1, size(components)
            if (components(i)%component_id <= 0) return
            if (abs(components(i)%sigma) /= 1) return
            if (.not. ieee_is_finite(components(i)%x_min)) return
            if (.not. ieee_is_finite(components(i)%x_max)) return
            if (components(i)%x_max <= components(i)%x_min) return
            do j = 1, i - 1
                if (same_identity(components(i), components(j))) then
                    status = GC_CYL_NONLOCAL_COMPONENT_IDENTITY
                    return
                end if
                if (components(i)%sigma /= components(j)%sigma) cycle
                if (intervals_overlap(components(i), components(j))) return
            end do
        end do
        status = GC_CYL_NONLOCAL_SUCCESS
    end function validate_components

    pure logical function same_identity(first, second)
        type(gc_nonlocal_component_t), intent(in) :: first, second

        same_identity = first%component_id == second%component_id
        if (same_identity) same_identity = first%sigma == second%sigma
    end function same_identity

    pure logical function intervals_overlap(first, second)
        type(gc_nonlocal_component_t), intent(in) :: first, second

        intervals_overlap = min(first%x_max, second%x_max) &
            > max(first%x_min, second%x_min)
    end function intervals_overlap

    pure logical function valid_section_coordinate(name)
        character(len=*), intent(in) :: name
        character(len=64) :: lowered

        lowered = lower_ascii(adjustl(trim(name)))
        valid_section_coordinate = len_trim(lowered) > 0
        if (.not. valid_section_coordinate) return
        if (index(lowered, 'eta') > 0) valid_section_coordinate = .false.
        if (index(lowered, 'pitch') > 0) valid_section_coordinate = .false.
    end function valid_section_coordinate

    pure function lower_ascii(value) result(lowered)
        character(len=*), intent(in) :: value
        character(len=len(value)) :: lowered
        integer :: i, code

        lowered = value
        do i = 1, len(value)
            code = iachar(lowered(i:i))
            if (code >= iachar('A') .and. code <= iachar('Z')) then
                lowered(i:i) = achar(code + iachar('a') - iachar('A'))
            end if
        end do
    end function lower_ascii

end module neort_gc_cylindrical_nonlocal_provider
