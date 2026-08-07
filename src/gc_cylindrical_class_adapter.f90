module neort_gc_cylindrical_class_adapter
    !! Fixed-(H0,J_K) Poincare-cut class and launch adapter.
    !!
    !! The class coordinate is the physical cut parameter Rc.  H0 and the
    !! compatibility field named jperp are the action variables of Buchholz
    !! et al. (2022).  In this adapter jperp means J_K, with
    !! J_K=m*c*mu_phys/abs(q); it is not the repository's magnetic-moment
    !! storage convention and it is not POTATO's normalized invariant.
    !! The adapter evaluates
    !!
    !!   mu_phys = abs(q) J_K / (m c),
    !!   v_parallel^2 = 2 (H0 - mu_phys B - q Phi) / m
    !!                = 2 (H0 - J_K omega_c - q Phi) / m,
    !!   psi_star = psi + (m c/q) v_parallel (R B_phi/B),
    !!
    !! where B_phi in the second expression is the covariant cylindrical
    !! component R B_hatphi.  Every allowed Rc interval is retained.  A
    !! separate orbit-return/fixed-point splitter is deliberately required
    !! before the intervals can be called complete topological classes.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_gc_callback_context, only: gc_callback_context_t
    use neort_gc_cylindrical_model, only: &
        GC_CYL_SUCCESS, gc_cylindrical_field_sample_t, &
        gc_cylindrical_field_t, gc_cylindrical_potential_t, &
        gc_cylindrical_state_t, gc_cylindrical_allowed_component_t, &
        canonical_toroidal_momentum_from_state
    use neort_gc_perpendicular_invariant, only: &
        gc_buchholz_jk_from_mu_phys, gc_mu_phys_from_buchholz_jk
    use neort_gc_cylindrical_topology, only: &
        find_gc_cylindrical_allowed_regions, &
        gc_cylindrical_allowed_region_set_t

    implicit none
    private

    integer, parameter, public :: GC_CYL_CLASS_SUCCESS = 0
    integer, parameter, public :: GC_CYL_CLASS_INVALID_INPUT = 1
    integer, parameter, public :: GC_CYL_CLASS_FIELD_ERROR = 2
    integer, parameter, public :: GC_CYL_CLASS_POTENTIAL_ERROR = 3
    integer, parameter, public :: GC_CYL_CLASS_CUT_ERROR = 4
    integer, parameter, public :: GC_CYL_CLASS_INTERIOR_INVALID = 5
    integer, parameter, public :: GC_CYL_CLASS_NO_ALLOWED_INTERVAL = 6
    integer, parameter, public :: GC_CYL_CLASS_SPLITTER_UNAVAILABLE = 7
    integer, parameter, public :: GC_CYL_CLASS_SPLITTER_FAILURE = 8
    integer, parameter, public :: GC_CYL_CLASS_NOT_ENUMERATED = 9
    integer, parameter, public :: GC_CYL_CLASS_NOT_ALLOWED = 10
    integer, parameter, public :: GC_CYL_CLASS_COMPONENT_ERROR = 11
    integer, parameter, public :: GC_CYL_CLASS_LAUNCH_ERROR = 12

    character(len=*), parameter, public :: GC_CYL_CLASS_PSI_STAR_UNITS = &
        'psi_star=(c/q)*p_phi'
    character(len=*), parameter, public :: GC_CYL_CLASS_H0_UNITS = 'energy'
    character(len=*), parameter, public :: GC_CYL_CLASS_JK_UNITS = &
        'J_K=mass*c*mu_phys/abs(charge)'
    ! Keep the historical public name for source compatibility.  Its value
    ! is explicitly Buchholz J_K, never 2*m*mu_phys and never POTATO J_perp.
    character(len=*), parameter, public :: GC_CYL_CLASS_JPERP_UNITS = &
        GC_CYL_CLASS_JK_UNITS
    character(len=*), parameter, public :: GC_CYL_CLASS_RC_UNITS = 'cut-coordinate'
    character(len=*), parameter, public :: GC_CYL_CLASS_VPARALLEL_UNITS = 'velocity'

    type, public :: gc_cylindrical_class_options_t
        integer :: scan_points = 257
        integer :: validation_points = 1025
        real(dp) :: allowed_tolerance = 1.0e-11_dp
        real(dp) :: tangency_tolerance = 1.0e-10_dp
    end type gc_cylindrical_class_options_t

    type, public :: gc_cylindrical_class_interval_t
        integer :: component_id = 0
        integer :: sigma = 0
        real(dp) :: rc_min = 0.0_dp
        real(dp) :: rc_max = 0.0_dp
        real(dp) :: psi_star_min = 0.0_dp
        real(dp) :: psi_star_max = 0.0_dp
        real(dp) :: canonical_measure = 0.0_dp
        logical :: lower_root = .false.
        logical :: upper_root = .false.
        logical :: lower_tangent = .false.
        logical :: upper_tangent = .false.
        logical :: allowed_interval = .true.
        logical :: topology_certified = .false.
        logical :: orbit_return_certified = .false.
        logical :: root_isolation_certified = .false.
        character(len=24) :: lower_boundary_kind = 'unresolved'
        character(len=24) :: upper_boundary_kind = 'unresolved'
        character(len=32) :: limiting_chart = 'unresolved'
    end type gc_cylindrical_class_interval_t

    type, public :: gc_cylindrical_class_point_t
        real(dp) :: rc = 0.0_dp
        integer :: sigma = 0
        real(dp) :: position(3) = 0.0_dp
        real(dp) :: dposition_drc(3) = 0.0_dp
        type(gc_cylindrical_field_sample_t) :: field
        real(dp) :: potential = 0.0_dp
        real(dp) :: potential_gradient(3) = 0.0_dp
        real(dp) :: omega_c = 0.0_dp
        real(dp) :: domega_c_drc = 0.0_dp
        real(dp) :: vparallel_squared = 0.0_dp
        real(dp) :: dvparallel_squared_drc = 0.0_dp
        real(dp) :: vparallel = 0.0_dp
        real(dp) :: dvparallel_drc = 0.0_dp
        real(dp) :: psi_star = 0.0_dp
        real(dp) :: dpsi_star_drc = 0.0_dp
        logical :: allowed = .false.
        logical :: at_turning_point = .false.
        logical :: derivative_available = .false.
        logical :: derivative_singular = .false.
    end type gc_cylindrical_class_point_t

    type, public :: gc_cylindrical_class_launch_t
        integer :: status = GC_CYL_CLASS_LAUNCH_ERROR
        integer :: component_id = 0
        integer :: sigma = 0
        real(dp) :: rc = 0.0_dp
        real(dp) :: h0 = 0.0_dp
        real(dp) :: jperp = 0.0_dp
        real(dp) :: p_phi = 0.0_dp
        real(dp) :: energy_residual = 0.0_dp
        real(dp) :: jperp_residual = 0.0_dp
        real(dp) :: canonical_residual = 0.0_dp
        real(dp) :: psi_star = 0.0_dp
        real(dp) :: dpsi_star_drc = 0.0_dp
        real(dp) :: vparallel_squared = 0.0_dp
        real(dp) :: vparallel = 0.0_dp
        real(dp) :: omega_c = 0.0_dp
        real(dp) :: position(3) = 0.0_dp
        type(gc_cylindrical_state_t) :: state
        logical :: class_certified = .false.
        logical :: endpoint_tangent = .false.
        logical :: derivative_available = .false.
        character(len=32) :: h0_units = GC_CYL_CLASS_H0_UNITS
        character(len=64) :: jperp_units = GC_CYL_CLASS_JPERP_UNITS
        character(len=32) :: rc_units = GC_CYL_CLASS_RC_UNITS
        character(len=32) :: vparallel_units = GC_CYL_CLASS_VPARALLEL_UNITS
        character(len=32) :: psi_star_units = GC_CYL_CLASS_PSI_STAR_UNITS
    end type gc_cylindrical_class_launch_t

    type, public :: gc_cylindrical_class_result_t
        integer :: status = GC_CYL_CLASS_INVALID_INPUT
        integer :: nallowed_intervals = 0
        integer :: nclasses = 0
        logical :: splitter_present = .false.
        logical :: class_complete = .false.
        type(gc_cylindrical_class_interval_t), allocatable :: &
            allowed_intervals(:)
        type(gc_cylindrical_class_interval_t), allocatable :: classes(:)
    end type gc_cylindrical_class_result_t

    abstract interface
        subroutine gc_cylindrical_allowed_region_provider_i(h0, jperp, sigma, &
                user_data, regions, status)
            import :: dp, gc_callback_context_t, &
                gc_cylindrical_allowed_region_set_t
            real(dp), intent(in) :: h0, jperp
            integer, intent(in) :: sigma
            class(gc_callback_context_t), pointer, intent(inout) :: user_data
            type(gc_cylindrical_allowed_region_set_t), intent(out) :: regions
            integer, intent(out) :: status
        end subroutine gc_cylindrical_allowed_region_provider_i

        subroutine gc_cylindrical_allowed_region_verifier_i(h0, jperp, sigma, &
                rc_min, rc_max, regions, user_data, certificate_id, status)
            import :: dp, gc_callback_context_t, &
                gc_cylindrical_allowed_region_set_t
            real(dp), intent(in) :: h0, jperp, rc_min, rc_max
            integer, intent(in) :: sigma
            type(gc_cylindrical_allowed_region_set_t), intent(in) :: regions
            class(gc_callback_context_t), pointer, intent(inout) :: user_data
            integer, intent(out) :: certificate_id
            integer, intent(out) :: status
        end subroutine gc_cylindrical_allowed_region_verifier_i

        subroutine gc_cylindrical_class_cut_map_i(rc, user_data, position, &
                dposition_drc, status)
            import :: dp, gc_callback_context_t
            real(dp), intent(in) :: rc
            class(gc_callback_context_t), pointer, intent(inout) :: user_data
            real(dp), intent(out) :: position(3)
            real(dp), intent(out) :: dposition_drc(3)
            integer, intent(out) :: status
        end subroutine gc_cylindrical_class_cut_map_i

        subroutine gc_cylindrical_class_splitter_i(h0, jperp, sigma, candidate, &
                user_data, split_classes, certified, status)
            import :: dp, gc_callback_context_t, &
                gc_cylindrical_class_interval_t
            real(dp), intent(in) :: h0, jperp
            integer, intent(in) :: sigma
            type(gc_cylindrical_class_interval_t), intent(in) :: candidate
            class(gc_callback_context_t), pointer, intent(inout) :: user_data
            type(gc_cylindrical_class_interval_t), allocatable, intent(out) :: &
                split_classes(:)
            logical, intent(out) :: certified
            integer, intent(out) :: status
        end subroutine gc_cylindrical_class_splitter_i
    end interface

    type, public :: gc_cylindrical_class_adapter_t
        real(dp) :: h0 = 0.0_dp
        ! Compatibility component name: this value is Buchholz J_K.
        real(dp) :: jperp = 0.0_dp
        real(dp) :: mass = 0.0_dp
        real(dp) :: charge = 0.0_dp
        real(dp) :: c_light = 0.0_dp
        real(dp) :: rc_min = 0.0_dp
        real(dp) :: rc_max = 0.0_dp
        type(gc_cylindrical_class_options_t) :: options
        class(gc_cylindrical_field_t), pointer :: field => null()
        class(gc_cylindrical_potential_t), pointer :: potential => null()
        procedure(gc_cylindrical_class_cut_map_i), pointer, nopass :: &
            cut_map => null()
        procedure(gc_cylindrical_class_splitter_i), pointer, nopass :: &
            splitter => null()
        procedure(gc_cylindrical_allowed_region_provider_i), pointer, nopass :: &
            allowed_region_provider => null()
        procedure(gc_cylindrical_allowed_region_verifier_i), pointer, nopass :: &
            allowed_region_verifier => null()
        integer :: allowed_region_certificate_id = 0
        ! Borrowed procedure associations.  The provider and verifier targets
        ! must remain valid until this adapter is cleared; they are never
        ! copied into owned executable state.
        ! Borrowed target; the caller must keep it alive until this adapter is
        ! cleared and must not copy the adapter past that target's lifetime.
        class(gc_callback_context_t), pointer :: user_data => null()
        ! Provider data is only input to the typed verifier; its flags and
        ! labels are never treated as a topology proof.
        logical :: initialized = .false.
        logical :: classes_enumerated = .false.
        logical :: class_complete = .false.
        type(gc_cylindrical_class_interval_t), allocatable :: intervals(:)
    end type gc_cylindrical_class_adapter_t

    public :: gc_cylindrical_class_cut_map_i
    public :: gc_cylindrical_class_splitter_i
    public :: gc_cylindrical_allowed_region_provider_i
    public :: gc_cylindrical_allowed_region_verifier_i
    public :: initialize_gc_cylindrical_class_adapter
    public :: clear_gc_cylindrical_class_adapter
    public :: evaluate_gc_cylindrical_class_point
    public :: enumerate_gc_cylindrical_classes
    public :: launch_gc_cylindrical_class

contains

    subroutine initialize_gc_cylindrical_class_adapter(field_model, potential_model, &
            h0, jperp, mass, charge, c_light, rc_min, rc_max, cut_map, adapter, &
            status, options, splitter, user_data, allowed_region_provider, &
            allowed_region_verifier, allowed_region_certificate_id)
        class(gc_cylindrical_field_t), target, intent(in) :: field_model
        class(gc_cylindrical_potential_t), target, intent(in) :: potential_model
        real(dp), intent(in) :: h0, jperp, mass, charge, c_light, rc_min, rc_max
        procedure(gc_cylindrical_class_cut_map_i) :: cut_map
        type(gc_cylindrical_class_adapter_t), intent(out) :: adapter
        integer, intent(out) :: status
        type(gc_cylindrical_class_options_t), intent(in), optional :: options
        procedure(gc_cylindrical_class_splitter_i), optional :: splitter
        ! Stored as a borrowed pointer; user_data must outlive adapter and all
        ! callback use, and the adapter must be cleared before it is invalidated.
        class(gc_callback_context_t), target, intent(inout), optional :: user_data
        procedure(gc_cylindrical_allowed_region_provider_i), optional :: &
            allowed_region_provider
        procedure(gc_cylindrical_allowed_region_verifier_i), optional :: &
            allowed_region_verifier
        integer, intent(in), optional :: allowed_region_certificate_id

        ! Optional callback procedures and user_data are borrowed for the
        ! lifetime of adapter.  Callers must clear the adapter before any of
        ! those associations or their defining scopes cease to exist.

        adapter%h0 = 0.0_dp
        adapter%jperp = 0.0_dp
        adapter%mass = 0.0_dp
        adapter%charge = 0.0_dp
        adapter%c_light = 0.0_dp
        adapter%rc_min = 0.0_dp
        adapter%rc_max = 0.0_dp
        adapter%options = gc_cylindrical_class_options_t()
        nullify(adapter%field)
        nullify(adapter%potential)
        nullify(adapter%cut_map)
        nullify(adapter%splitter)
        nullify(adapter%allowed_region_provider)
        nullify(adapter%allowed_region_verifier)
        nullify(adapter%user_data)
        adapter%allowed_region_certificate_id = 0
        if (allocated(adapter%intervals)) deallocate(adapter%intervals)
        adapter%initialized = .false.
        adapter%classes_enumerated = .false.
        adapter%class_complete = .false.

        status = GC_CYL_CLASS_INVALID_INPUT
        if (.not. all(ieee_is_finite([h0, jperp, mass, charge, c_light, &
            rc_min, rc_max]))) return
        if (mass <= 0.0_dp) return
        if (jperp < 0.0_dp) return
        if (abs(charge) <= tiny(charge)) return
        if (c_light <= 0.0_dp) return
        if (rc_max <= rc_min) return
        if (present(allowed_region_provider)) then
            if (.not. present(allowed_region_verifier)) return
            if (.not. present(allowed_region_certificate_id)) return
            if (allowed_region_certificate_id <= 0) return
        else
            if (present(allowed_region_verifier)) return
            if (present(allowed_region_certificate_id)) return
        end if
        if (present(options)) then
            if (options%scan_points < 2) return
            if (options%validation_points < 2) return
            if (options%allowed_tolerance <= 0.0_dp) return
            if (options%tangency_tolerance <= 0.0_dp) return
            adapter%options = options
        end if

        adapter%h0 = h0
        adapter%jperp = jperp
        adapter%mass = mass
        adapter%charge = charge
        adapter%c_light = c_light
        adapter%rc_min = rc_min
        adapter%rc_max = rc_max
        adapter%field => field_model
        adapter%potential => potential_model
        adapter%cut_map => cut_map
        if (present(splitter)) adapter%splitter => splitter
        if (present(allowed_region_provider)) then
            adapter%allowed_region_provider => allowed_region_provider
            adapter%allowed_region_verifier => allowed_region_verifier
            adapter%allowed_region_certificate_id = allowed_region_certificate_id
        end if
        if (present(user_data)) adapter%user_data => user_data
        adapter%initialized = .true.
        status = GC_CYL_CLASS_SUCCESS
    end subroutine initialize_gc_cylindrical_class_adapter

    subroutine clear_gc_cylindrical_class_adapter(adapter)
        type(gc_cylindrical_class_adapter_t), intent(inout) :: adapter

        if (allocated(adapter%intervals)) deallocate(adapter%intervals)
        nullify(adapter%field)
        nullify(adapter%potential)
        nullify(adapter%cut_map)
        nullify(adapter%splitter)
        nullify(adapter%allowed_region_provider)
        nullify(adapter%allowed_region_verifier)
        nullify(adapter%user_data)
        adapter%allowed_region_certificate_id = 0
        adapter%initialized = .false.
        adapter%classes_enumerated = .false.
        adapter%class_complete = .false.
    end subroutine clear_gc_cylindrical_class_adapter

    subroutine evaluate_gc_cylindrical_class_point(adapter, rc, sigma, point, &
            status)
        type(gc_cylindrical_class_adapter_t), intent(inout) :: adapter
        real(dp), intent(in) :: rc
        integer, intent(in) :: sigma
        type(gc_cylindrical_class_point_t), intent(out) :: point
        integer, intent(out) :: status

        real(dp) :: tangent_physical(3), db_drc(3), bphi_covariant
        real(dp) :: dbphi_covariant_drc, dbmod_drc, ratio, dratio_drc
        real(dp) :: kinetic_energy, v2_scale, velocity_scale, mu_phys
        real(dp) :: charge_magnitude, coefficient
        integer :: local_status

        point = gc_cylindrical_class_point_t()
        point%rc = rc
        point%sigma = sigma
        status = GC_CYL_CLASS_INVALID_INPUT
        if (.not. adapter%initialized) return
        if (.not. ieee_is_finite(rc)) return
        if (abs(sigma) /= 1) return
        if (.not. associated(adapter%field)) return
        if (.not. associated(adapter%potential)) return
        if (.not. associated(adapter%cut_map)) return

        call adapter%cut_map(rc, adapter%user_data, point%position, &
            point%dposition_drc, local_status)
        if (local_status /= GC_CYL_SUCCESS) then
            status = GC_CYL_CLASS_CUT_ERROR
            return
        end if
        if (.not. all(ieee_is_finite([point%position, point%dposition_drc]))) then
            status = GC_CYL_CLASS_CUT_ERROR
            return
        end if
        if (point%position(1) <= 0.0_dp) then
            status = GC_CYL_CLASS_CUT_ERROR
            return
        end if

        call adapter%field%evaluate(point%position, point%field, local_status)
        if (local_status /= GC_CYL_SUCCESS) then
            status = GC_CYL_CLASS_FIELD_ERROR
            return
        end if
        call adapter%potential%evaluate(point%position, point%field, &
            point%potential, point%potential_gradient, local_status)
        if (local_status /= GC_CYL_SUCCESS) then
            status = GC_CYL_CLASS_POTENTIAL_ERROR
            return
        end if
        if (.not. all(ieee_is_finite([point%potential, &
            point%potential_gradient]))) then
            status = GC_CYL_CLASS_POTENTIAL_ERROR
            return
        end if

        ! dposition_drc is in (R,Z,phi) coordinates.  The field sample uses
        ! the orthonormal (R,phi,Z) basis and db_dq is differentiated in
        ! physical distance, so the phi tangent carries the factor R.
        tangent_physical = [point%dposition_drc(1), &
            point%position(1)*point%dposition_drc(3), point%dposition_drc(2)]
        db_drc = point%field%db_dq(:, 1)*tangent_physical(1) &
            +point%field%db_dq(:, 2)*tangent_physical(2) &
            +point%field%db_dq(:, 3)*tangent_physical(3)
        dbmod_drc = dot_product(point%field%bhat, db_drc)
        if (.not. all(ieee_is_finite([tangent_physical, db_drc, dbmod_drc]))) then
            status = GC_CYL_CLASS_FIELD_ERROR
            return
        end if

        point%omega_c = abs(adapter%charge)*point%field%bmod &
            /(adapter%mass*adapter%c_light)
        point%domega_c_drc = abs(adapter%charge)*dbmod_drc &
            /(adapter%mass*adapter%c_light)
        mu_phys = gc_mu_phys_from_buchholz_jk(adapter%jperp, adapter%mass, &
            adapter%charge, adapter%c_light)
        kinetic_energy = adapter%h0 - mu_phys*point%field%bmod &
            -adapter%charge*point%potential
        point%vparallel_squared = 2.0_dp*kinetic_energy/adapter%mass
        point%dvparallel_squared_drc = -2.0_dp/adapter%mass &
            *(mu_phys*dbmod_drc &
            +adapter%charge*dot_product(point%potential_gradient, &
            tangent_physical))
        v2_scale = max(1.0_dp, abs(2.0_dp*adapter%h0/adapter%mass))
        velocity_scale = sqrt(max(v2_scale, tiny(v2_scale)))
        point%at_turning_point = abs(point%vparallel_squared) <= &
            adapter%options%allowed_tolerance*v2_scale
        point%allowed = point%vparallel_squared >= &
            -adapter%options%allowed_tolerance*v2_scale

        point%psi_star = point%field%psi
        point%dpsi_star_drc = dot_product(point%field%grad_psi, tangent_physical)
        if (point%vparallel_squared > 0.0_dp) then
            point%vparallel = real(sigma, dp)*sqrt(point%vparallel_squared)
            if (sqrt(point%vparallel_squared) > &
                100.0_dp*sqrt(tiny(point%vparallel_squared))) then
                point%dvparallel_drc = real(sigma, dp) &
                    *point%dvparallel_squared_drc &
                    /(2.0_dp*sqrt(point%vparallel_squared))
                point%derivative_available = .true.
            else
                point%dvparallel_drc = 0.0_dp
                point%derivative_singular = .true.
            end if
        else
            point%vparallel = 0.0_dp
            point%dvparallel_drc = 0.0_dp
            if (point%vparallel_squared < 0.0_dp) then
                point%derivative_singular = .true.
            else
                if (abs(point%dvparallel_squared_drc) > &
                    adapter%options%tangency_tolerance*v2_scale) then
                    point%derivative_singular = .true.
                else
                    point%derivative_available = .true.
                end if
            end if
        end if

        bphi_covariant = point%position(1)*point%field%b(2)
        dbphi_covariant_drc = tangent_physical(1)*point%field%b(2) &
            +point%position(1)*db_drc(2)
        ratio = bphi_covariant/point%field%bmod
        dratio_drc = dbphi_covariant_drc/point%field%bmod &
            -bphi_covariant*dbmod_drc/(point%field%bmod**2)
        charge_magnitude = abs(adapter%charge)
        coefficient = adapter%mass*adapter%c_light/adapter%charge
        if (point%vparallel_squared > 0.0_dp) then
            point%psi_star = point%field%psi + coefficient*point%vparallel*ratio
            if (point%derivative_available) then
                point%dpsi_star_drc = &
                    dot_product(point%field%grad_psi, tangent_physical) &
                    +coefficient*(point%dvparallel_drc*ratio &
                    +point%vparallel*dratio_drc)
            else
                point%dpsi_star_drc = &
                    dot_product(point%field%grad_psi, tangent_physical)
            end if
        else
            ! The value at a turning point is exact.  A simple turning root
            ! has a square-root derivative; retain a finite endpoint value
            ! but explicitly mark that derivative unavailable.
            point%psi_star = point%field%psi
            point%dpsi_star_drc = &
                dot_product(point%field%grad_psi, tangent_physical)
        end if
        if (.not. all(ieee_is_finite([point%omega_c, point%domega_c_drc, &
            point%vparallel_squared, point%dvparallel_squared_drc, &
            point%vparallel, point%psi_star, point%dpsi_star_drc]))) then
            status = GC_CYL_CLASS_FIELD_ERROR
            return
        end if
        if (charge_magnitude <= tiny(charge_magnitude)) then
            status = GC_CYL_CLASS_INVALID_INPUT
            return
        end if
        if (velocity_scale <= 0.0_dp) then
            status = GC_CYL_CLASS_INVALID_INPUT
            return
        end if
        status = GC_CYL_CLASS_SUCCESS
    end subroutine evaluate_gc_cylindrical_class_point

    subroutine enumerate_gc_cylindrical_classes(adapter, result, status)
        type(gc_cylindrical_class_adapter_t), intent(inout) :: adapter
        type(gc_cylindrical_class_result_t), intent(out) :: result
        integer, intent(out) :: status

        type(gc_cylindrical_allowed_region_set_t) :: regions
        type(gc_cylindrical_class_interval_t), allocatable :: raw(:)
        type(gc_cylindrical_class_interval_t), allocatable :: split(:)
        type(gc_cylindrical_class_interval_t) :: interval
        integer :: sigma, i, j, local_status, nraw, nclass
        logical :: certified
        integer :: certificate_id

        result = gc_cylindrical_class_result_t()
        adapter%classes_enumerated = .false.
        adapter%class_complete = .false.
        if (allocated(adapter%intervals)) deallocate(adapter%intervals)
        status = validate_adapter(adapter)
        if (status /= GC_CYL_CLASS_SUCCESS) then
            result%status = status
            return
        end if
        call validate_cut_samples(adapter, local_status)
        if (local_status /= GC_CYL_CLASS_SUCCESS) then
            result%status = local_status
            status = local_status
            return
        end if

        nraw = 0
        do sigma = -1, 1, 2
            regions = gc_cylindrical_allowed_region_set_t()
            local_status = GC_CYL_CLASS_INTERIOR_INVALID
            certificate_id = 0
            if (associated(adapter%allowed_region_provider)) then
                call adapter%allowed_region_provider(adapter%h0, adapter%jperp, &
                    sigma, adapter%user_data, regions, local_status)
                if (local_status /= GC_CYL_SUCCESS) then
                    status = GC_CYL_CLASS_INTERIOR_INVALID
                    result%status = status
                    return
                end if
                local_status = GC_CYL_CLASS_SPLITTER_FAILURE
                certificate_id = 0
                call adapter%allowed_region_verifier(adapter%h0, adapter%jperp, &
                    sigma, adapter%rc_min, adapter%rc_max, regions, &
                    adapter%user_data, certificate_id, local_status)
                if (local_status /= GC_CYL_SUCCESS) then
                    result%status = GC_CYL_CLASS_SPLITTER_FAILURE
                    status = GC_CYL_CLASS_SPLITTER_FAILURE
                    return
                end if
                if (certificate_id /= adapter%allowed_region_certificate_id) then
                    result%status = GC_CYL_CLASS_SPLITTER_FAILURE
                    status = GC_CYL_CLASS_SPLITTER_FAILURE
                    return
                end if
                call validate_certified_regions(adapter, sigma, regions, &
                    local_status)
                if (local_status /= GC_CYL_CLASS_SUCCESS) then
                    result%status = local_status
                    status = local_status
                    return
                end if
            else
                call find_gc_cylindrical_allowed_regions(evaluate_sigma, &
                    adapter%rc_min, adapter%rc_max, adapter%options%scan_points, &
                    sigma, regions, local_status)
                if (local_status /= GC_CYL_SUCCESS) then
                    status = GC_CYL_CLASS_INTERIOR_INVALID
                    result%status = status
                    return
                end if
                if (.not. regions%topology_certified) then
                    ! The sampler returns useful diagnostics, but finite values
                    ! and finite-difference slopes cannot exclude a missed even
                    ! root, an X point, or a separatrix.  Only a generated
                    ! interval/root-isolation certificate may cross this gate.
                    status = GC_CYL_CLASS_SPLITTER_FAILURE
                    result%status = status
                    return
                end if
            end if
            do i = 1, regions%ncomponents
                call interval_from_topology(adapter, regions%components(i), &
                    interval, local_status, &
                    associated(adapter%allowed_region_provider))
                if (local_status /= GC_CYL_CLASS_SUCCESS) then
                    result%status = local_status
                    status = local_status
                    return
                end if
                call append_interval(raw, nraw, interval, local_status)
                if (local_status /= GC_CYL_CLASS_SUCCESS) then
                    result%status = local_status
                    status = local_status
                    return
                end if
            end do
        end do
        if (nraw == 0) then
            result%status = GC_CYL_CLASS_NO_ALLOWED_INTERVAL
            status = GC_CYL_CLASS_NO_ALLOWED_INTERVAL
            return
        end if

        allocate(result%allowed_intervals(nraw))
        result%allowed_intervals = raw(1:nraw)
        result%nallowed_intervals = nraw
        result%splitter_present = associated(adapter%splitter)

        if (.not. associated(adapter%splitter)) then
            allocate(result%classes(nraw))
            result%classes = raw(1:nraw)
            result%nclasses = nraw
            result%class_complete = .false.
            adapter%class_complete = .false.
            adapter%classes_enumerated = .true.
            call copy_intervals(adapter%intervals, raw, nraw)
            result%status = GC_CYL_CLASS_SPLITTER_UNAVAILABLE
            status = GC_CYL_CLASS_SPLITTER_UNAVAILABLE
            return
        end if

        nclass = 0
        do i = 1, nraw
            certified = .false.
            if (allocated(split)) deallocate(split)
            call adapter%splitter(adapter%h0, adapter%jperp, raw(i)%sigma, &
                raw(i), adapter%user_data, split, certified, local_status)
            if (local_status /= GC_CYL_CLASS_SUCCESS) then
                result%status = GC_CYL_CLASS_SPLITTER_FAILURE
                status = GC_CYL_CLASS_SPLITTER_FAILURE
                return
            end if
            if (.not. certified) then
                result%status = GC_CYL_CLASS_SPLITTER_FAILURE
                status = GC_CYL_CLASS_SPLITTER_FAILURE
                return
            end if
            call validate_split(adapter, raw(i), split, certified, local_status)
            if (local_status /= GC_CYL_CLASS_SUCCESS) then
                result%status = local_status
                status = local_status
                return
            end if
            do j = 1, size(split)
                call append_interval(result%classes, nclass, split(j), local_status)
                if (local_status /= GC_CYL_CLASS_SUCCESS) then
                    result%status = local_status
                    status = local_status
                    return
                end if
            end do
        end do
        result%nclasses = nclass
        result%class_complete = .true.
        result%status = GC_CYL_CLASS_SUCCESS
        adapter%class_complete = .true.
        adapter%classes_enumerated = .true.
        call copy_intervals(adapter%intervals, result%classes, nclass)
        status = GC_CYL_CLASS_SUCCESS

    contains

        subroutine evaluate_sigma(x, value, derivative, psi_star, &
                dpsi_star_dx, callback_status)
            real(dp), intent(in) :: x
            real(dp), intent(out) :: value, derivative, psi_star
            real(dp), intent(out) :: dpsi_star_dx
            integer, intent(out) :: callback_status

            type(gc_cylindrical_class_point_t) :: point
            integer :: point_status

            call evaluate_gc_cylindrical_class_point(adapter, x, sigma, point, &
                point_status)
            if (point_status /= GC_CYL_CLASS_SUCCESS) then
                value = 0.0_dp
                derivative = 0.0_dp
                psi_star = 0.0_dp
                dpsi_star_dx = 0.0_dp
                callback_status = GC_CYL_CLASS_INTERIOR_INVALID
                return
            end if
            value = point%vparallel_squared
            derivative = point%dvparallel_squared_drc
            psi_star = point%psi_star
            dpsi_star_dx = point%dpsi_star_drc
            callback_status = GC_CYL_SUCCESS
        end subroutine evaluate_sigma

    end subroutine enumerate_gc_cylindrical_classes

    subroutine validate_certified_regions(adapter, sigma, regions, status)
        type(gc_cylindrical_class_adapter_t), intent(in) :: adapter
        integer, intent(in) :: sigma
        type(gc_cylindrical_allowed_region_set_t), intent(in) :: regions
        integer, intent(out) :: status

        integer :: i, j
        real(dp) :: bound_tolerance, measure_sum, measure_tolerance
        real(dp) :: canonical_scale, canonical_tolerance
        real(dp) :: previous_end
        logical :: root_found
        integer :: root_index, root_matches

        status = GC_CYL_CLASS_SPLITTER_FAILURE
        ! topology_certified and certificate_method are provider metadata, not
        ! proofs.  The typed verifier and its adapter-owned certificate ID
        ! have already established certification before this structural check.
        if (regions%nroots < 0) return
        if (regions%ncomponents < 0) return
        if (.not. ieee_is_finite(regions%total_canonical_measure)) return
        if (regions%total_canonical_measure < 0.0_dp) return

        bound_tolerance = 100.0_dp*adapter%options%allowed_tolerance*max(1.0_dp, &
            max(abs(adapter%rc_min), abs(adapter%rc_max)))
        if (.not. ieee_is_finite(bound_tolerance)) return

        if (regions%nroots > 0) then
            if (.not. allocated(regions%roots)) return
            if (.not. allocated(regions%root_canonical)) return
            if (size(regions%roots) /= regions%nroots) return
            if (size(regions%root_canonical) /= regions%nroots) return
        else
            if (allocated(regions%roots)) then
                if (size(regions%roots) /= 0) return
            end if
            if (allocated(regions%root_canonical)) then
                if (size(regions%root_canonical) /= 0) return
            end if
        end if
        if (regions%nroots > 0) then
            if (.not. all(ieee_is_finite(regions%roots))) return
            if (.not. all(ieee_is_finite(regions%root_canonical))) return
            do i = 1, regions%nroots
                if (regions%roots(i) < adapter%rc_min - bound_tolerance) return
                if (regions%roots(i) > adapter%rc_max + bound_tolerance) return
                if (i > 1) then
                    if (regions%roots(i) <= regions%roots(i - 1)) return
                end if
            end do
        end if

        if (regions%ncomponents > 0) then
            if (.not. allocated(regions%components)) return
            if (size(regions%components) /= regions%ncomponents) return
        else
            if (allocated(regions%components)) then
                if (size(regions%components) /= 0) return
            end if
            if (regions%total_canonical_measure > bound_tolerance) return
            status = GC_CYL_CLASS_SUCCESS
            return
        end if

        canonical_scale = max(1.0_dp, abs(regions%total_canonical_measure))
        if (regions%nroots > 0) then
            do i = 1, regions%nroots
                canonical_scale = max(canonical_scale, &
                    abs(regions%root_canonical(i)))
            end do
        end if
        do i = 1, regions%ncomponents
            if (regions%components(i)%component_id <= 0) return
            if (regions%components(i)%sigma /= sigma) return
            if (.not. all(ieee_is_finite([regions%components(i)%x_begin, &
                regions%components(i)%x_end, &
                regions%components(i)%canonical_begin, &
                regions%components(i)%canonical_end, &
                regions%components(i)%canonical_measure]))) return
            if (regions%components(i)%x_end <= regions%components(i)%x_begin) return
            if (regions%components(i)%canonical_measure <= 0.0_dp) return
            canonical_scale = max(canonical_scale, &
                abs(regions%components(i)%canonical_begin))
            canonical_scale = max(canonical_scale, &
                abs(regions%components(i)%canonical_end))
            canonical_scale = max(canonical_scale, &
                abs(regions%components(i)%canonical_measure))
            do j = 1, i - 1
                if (regions%components(j)%component_id == &
                    regions%components(i)%component_id) return
            end do
        end do
        canonical_tolerance = 100.0_dp*adapter%options%allowed_tolerance &
            *max(1.0_dp, canonical_scale)
        if (.not. ieee_is_finite(canonical_tolerance)) return

        previous_end = adapter%rc_min
        measure_sum = 0.0_dp
        do i = 1, regions%ncomponents
            if (regions%components(i)%x_begin < adapter%rc_min - &
                bound_tolerance) return
            if (regions%components(i)%x_end > adapter%rc_max + &
                bound_tolerance) return
            if (i == 1) then
                if (regions%components(i)%x_begin > adapter%rc_min + &
                        bound_tolerance) then
                    if (.not. regions%components(i)%lower_root) return
                end if
            else
                if (regions%components(i)%x_begin < previous_end - &
                        bound_tolerance) return
                if (regions%components(i)%x_begin > previous_end + &
                        bound_tolerance) then
                    if (.not. regions%components(i - 1)%upper_root) return
                    if (.not. regions%components(i)%lower_root) return
                end if
            end if
            if (regions%components(i)%lower_root) then
                root_found = .false.
                root_index = 0
                root_matches = 0
                do j = 1, regions%nroots
                    if (abs(regions%roots(j) - regions%components(i)%x_begin) &
                        <= bound_tolerance) then
                        root_found = .true.
                        root_matches = root_matches + 1
                        root_index = j
                    end if
                end do
                if (.not. root_found) return
                if (root_matches /= 1) return
                if (abs(regions%root_canonical(root_index) - &
                        regions%components(i)%canonical_begin) > &
                    canonical_tolerance) return
            end if
            if (regions%components(i)%upper_root) then
                root_found = .false.
                root_index = 0
                root_matches = 0
                do j = 1, regions%nroots
                    if (abs(regions%roots(j) - regions%components(i)%x_end) &
                        <= bound_tolerance) then
                        root_found = .true.
                        root_matches = root_matches + 1
                        root_index = j
                    end if
                end do
                if (.not. root_found) return
                if (root_matches /= 1) return
                if (abs(regions%root_canonical(root_index) - &
                        regions%components(i)%canonical_end) > &
                    canonical_tolerance) return
            end if
            measure_sum = measure_sum + regions%components(i)%canonical_measure
            previous_end = regions%components(i)%x_end
        end do
        if (previous_end < adapter%rc_max - bound_tolerance) then
            if (.not. regions%components(regions%ncomponents)%upper_root) return
        end if
        if (previous_end > adapter%rc_max + bound_tolerance) return
        if (measure_sum <= 0.0_dp) return
        measure_tolerance = 1000.0_dp*epsilon(measure_sum)*max(1.0_dp, &
            max(abs(measure_sum), abs(regions%total_canonical_measure)))
        if (.not. ieee_is_finite(measure_tolerance)) return
        if (abs(measure_sum - regions%total_canonical_measure) > &
            measure_tolerance) return
        status = GC_CYL_CLASS_SUCCESS
    end subroutine validate_certified_regions

    subroutine launch_gc_cylindrical_class(adapter, rc, sigma, component_id, &
            launch, status)
        type(gc_cylindrical_class_adapter_t), intent(inout) :: adapter
        real(dp), intent(in) :: rc
        integer, intent(in) :: sigma, component_id
        type(gc_cylindrical_class_launch_t), intent(out) :: launch
        integer, intent(out) :: status

        type(gc_cylindrical_class_point_t) :: point
        integer :: i, point_status
        logical :: found
        real(dp) :: tolerance, mu_phys, energy, canonical_momentum

        launch = gc_cylindrical_class_launch_t()
        launch%rc = rc
        launch%sigma = sigma
        launch%component_id = component_id
        launch%h0 = adapter%h0
        launch%jperp = adapter%jperp
        status = GC_CYL_CLASS_LAUNCH_ERROR
        if (.not. adapter%initialized) return
        if (.not. adapter%classes_enumerated) then
            status = GC_CYL_CLASS_NOT_ENUMERATED
            launch%status = status
            return
        end if
        if (.not. adapter%class_complete) then
            status = GC_CYL_CLASS_SPLITTER_UNAVAILABLE
            launch%status = status
            return
        end if
        if (abs(sigma) /= 1) return
        found = .false.
        tolerance = 100.0_dp*adapter%options%allowed_tolerance &
            *max(1.0_dp, max(abs(adapter%rc_min), abs(adapter%rc_max)))
        if (allocated(adapter%intervals)) then
            do i = 1, size(adapter%intervals)
                if (adapter%intervals(i)%sigma /= sigma) cycle
                if (adapter%intervals(i)%component_id /= component_id) cycle
                if (rc < adapter%intervals(i)%rc_min - tolerance) cycle
                if (rc > adapter%intervals(i)%rc_max + tolerance) cycle
                found = .true.
                exit
            end do
        end if
        if (.not. found) then
            status = GC_CYL_CLASS_COMPONENT_ERROR
            launch%status = status
            return
        end if

        call evaluate_gc_cylindrical_class_point(adapter, rc, sigma, point, &
            point_status)
        if (point_status /= GC_CYL_CLASS_SUCCESS) then
            status = GC_CYL_CLASS_LAUNCH_ERROR
            launch%status = status
            return
        end if
        if (.not. point%allowed) then
            status = GC_CYL_CLASS_NOT_ALLOWED
            launch%status = status
            return
        end if

        launch%position = point%position
        launch%vparallel_squared = max(0.0_dp, point%vparallel_squared)
        launch%vparallel = point%vparallel
        launch%omega_c = point%omega_c
        launch%psi_star = point%psi_star
        launch%dpsi_star_drc = point%dpsi_star_drc
        launch%derivative_available = point%derivative_available
        launch%endpoint_tangent = point%at_turning_point .and. &
            abs(point%dvparallel_squared_drc) <= &
            adapter%options%tangency_tolerance &
            *max(1.0_dp, abs(2.0_dp*adapter%h0/adapter%mass))
        launch%state%R = point%position(1)
        launch%state%Z = point%position(2)
        launch%state%phi = point%position(3)
        launch%state%p_parallel = adapter%mass*point%vparallel
        mu_phys = gc_mu_phys_from_buchholz_jk(adapter%jperp, adapter%mass, &
            adapter%charge, adapter%c_light)
        launch%state%mu = mu_phys
        energy = launch%state%p_parallel**2/(2.0_dp*adapter%mass) &
            +mu_phys*point%field%bmod + adapter%charge*point%potential
        launch%p_phi = adapter%charge/adapter%c_light*point%psi_star
        canonical_momentum = canonical_toroidal_momentum_from_state( &
            launch%state, point%field, adapter%charge, adapter%c_light)
        launch%energy_residual = energy - adapter%h0
        launch%jperp_residual = gc_buchholz_jk_from_mu_phys(mu_phys, &
            adapter%mass, adapter%charge, adapter%c_light) - adapter%jperp
        launch%canonical_residual = canonical_momentum - launch%p_phi
        if (.not. all(ieee_is_finite([launch%energy_residual, &
            launch%jperp_residual, launch%canonical_residual]))) then
            status = GC_CYL_CLASS_LAUNCH_ERROR
            launch%status = status
            return
        end if
        launch%class_certified = adapter%class_complete
        launch%status = GC_CYL_CLASS_SUCCESS
        status = GC_CYL_CLASS_SUCCESS
    end subroutine launch_gc_cylindrical_class

    integer function validate_adapter(adapter)
        type(gc_cylindrical_class_adapter_t), intent(in) :: adapter

        validate_adapter = GC_CYL_CLASS_INVALID_INPUT
        if (.not. adapter%initialized) return
        if (.not. associated(adapter%field)) return
        if (.not. associated(adapter%potential)) return
        if (.not. associated(adapter%cut_map)) return
        if (.not. all(ieee_is_finite([adapter%h0, adapter%jperp, adapter%mass, &
            adapter%charge, adapter%c_light, adapter%rc_min, &
            adapter%rc_max]))) return
        if (adapter%mass <= 0.0_dp) return
        if (adapter%jperp < 0.0_dp) return
        if (abs(adapter%charge) <= tiny(adapter%charge)) return
        if (adapter%c_light <= 0.0_dp) return
        if (adapter%rc_max <= adapter%rc_min) return
        if (associated(adapter%allowed_region_provider)) then
            if (.not. associated(adapter%allowed_region_verifier)) return
            if (adapter%allowed_region_certificate_id <= 0) return
        else
            if (associated(adapter%allowed_region_verifier)) return
            if (adapter%allowed_region_certificate_id /= 0) return
        end if
        validate_adapter = GC_CYL_CLASS_SUCCESS
    end function validate_adapter

    subroutine validate_cut_samples(adapter, status)
        type(gc_cylindrical_class_adapter_t), intent(inout) :: adapter
        integer, intent(out) :: status

        type(gc_cylindrical_class_point_t) :: point
        integer :: i, n, local_status
        real(dp) :: rc

        n = max(adapter%options%validation_points, &
            adapter%options%scan_points)
        status = GC_CYL_CLASS_SUCCESS
        do i = 0, n
            rc = adapter%rc_min + (adapter%rc_max - adapter%rc_min) &
                *real(i, dp)/real(n, dp)
            call evaluate_gc_cylindrical_class_point(adapter, rc, 1, point, &
                local_status)
            if (local_status /= GC_CYL_CLASS_SUCCESS) then
                status = GC_CYL_CLASS_INTERIOR_INVALID
                return
            end if
            call evaluate_gc_cylindrical_class_point(adapter, rc, -1, point, &
                local_status)
            if (local_status /= GC_CYL_CLASS_SUCCESS) then
                status = GC_CYL_CLASS_INTERIOR_INVALID
                return
            end if
        end do
    end subroutine validate_cut_samples

    subroutine interval_from_topology(adapter, component, interval, status, &
            certified)
        type(gc_cylindrical_class_adapter_t), intent(inout) :: adapter
        type(gc_cylindrical_allowed_component_t), intent(in) :: component
        type(gc_cylindrical_class_interval_t), intent(out) :: interval
        integer, intent(out) :: status
        logical, intent(in) :: certified

        type(gc_cylindrical_class_point_t) :: lower_point, upper_point
        real(dp) :: scale
        integer :: lower_status, upper_status

        interval = gc_cylindrical_class_interval_t()
        interval%component_id = component%component_id
        interval%sigma = component%sigma
        interval%rc_min = component%x_begin
        interval%rc_max = component%x_end
        interval%psi_star_min = component%canonical_begin
        interval%psi_star_max = component%canonical_end
        interval%canonical_measure = component%canonical_measure
        interval%lower_root = component%lower_root
        interval%upper_root = component%upper_root
        ! Certification reaches this helper only after the typed provider's
        ! interval/root contract has passed validation.  The legacy finite scan
        ! never reaches it because its topology certificate is false.
        interval%topology_certified = certified
        interval%root_isolation_certified = certified
        interval%lower_boundary_kind = 'unresolved'
        interval%upper_boundary_kind = 'unresolved'
        interval%limiting_chart = 'unresolved'
        if (certified) then
            if (component%lower_root) then
                interval%lower_boundary_kind = 'certified-root'
            else
                interval%lower_boundary_kind = 'certified-bound'
            end if
            if (component%upper_root) then
                interval%upper_boundary_kind = 'certified-root'
            else
                interval%upper_boundary_kind = 'certified-bound'
            end if
            interval%limiting_chart = 'certified-provider'
        end if
        scale = max(1.0_dp, abs(2.0_dp*adapter%h0/adapter%mass))
        call evaluate_gc_cylindrical_class_point(adapter, interval%rc_min, &
            interval%sigma, lower_point, lower_status)
        call evaluate_gc_cylindrical_class_point(adapter, interval%rc_max, &
            interval%sigma, upper_point, upper_status)
        if (lower_status /= GC_CYL_CLASS_SUCCESS) then
            status = GC_CYL_CLASS_INTERIOR_INVALID
            return
        end if
        if (upper_status /= GC_CYL_CLASS_SUCCESS) then
            status = GC_CYL_CLASS_INTERIOR_INVALID
            return
        end if
        interval%lower_tangent = interval%lower_root
        if (interval%lower_tangent) then
            interval%lower_tangent = abs(lower_point%dvparallel_squared_drc) &
                <= adapter%options%tangency_tolerance*scale
        end if
        interval%upper_tangent = interval%upper_root
        if (interval%upper_tangent) then
            interval%upper_tangent = abs(upper_point%dvparallel_squared_drc) &
                <= adapter%options%tangency_tolerance*scale
        end if
        if (certified) then
            status = GC_CYL_CLASS_SUCCESS
        else
            status = GC_CYL_CLASS_SPLITTER_FAILURE
        end if
    end subroutine interval_from_topology

    subroutine validate_split(adapter, candidate, split, certified, status)
        type(gc_cylindrical_class_adapter_t), intent(inout) :: adapter
        type(gc_cylindrical_class_interval_t), intent(in) :: candidate
        type(gc_cylindrical_class_interval_t), allocatable, intent(inout) :: split(:)
        logical, intent(in) :: certified
        integer, intent(out) :: status

        type(gc_cylindrical_class_point_t) :: point
        real(dp) :: tolerance, canonical_tolerance, measure_tolerance
        real(dp) :: canonical_scale, previous_end, measure_sum
        integer :: i, point_status, j

        status = GC_CYL_CLASS_SPLITTER_FAILURE
        if (.not. certified) return
        if (.not. allocated(split)) return
        if (size(split) < 1) return
        if (candidate%component_id <= 0) return
        if (abs(candidate%sigma) /= 1) return
        if (.not. candidate%allowed_interval) return
        if (.not. all(ieee_is_finite([candidate%rc_min, candidate%rc_max, &
            candidate%psi_star_min, candidate%psi_star_max, &
            candidate%canonical_measure]))) return
        if (candidate%rc_max <= candidate%rc_min) return
        if (candidate%canonical_measure <= 0.0_dp) return
        tolerance = 100.0_dp*adapter%options%allowed_tolerance &
            *max(1.0_dp, max(abs(candidate%rc_min), abs(candidate%rc_max)))
        if (.not. ieee_is_finite(tolerance)) return
        canonical_scale = max(1.0_dp, abs(candidate%psi_star_min))
        canonical_scale = max(canonical_scale, abs(candidate%psi_star_max))
        canonical_scale = max(canonical_scale, abs(candidate%canonical_measure))
        do i = 1, size(split)
            if (.not. all(ieee_is_finite([split(i)%rc_min, split(i)%rc_max, &
                split(i)%psi_star_min, split(i)%psi_star_max, &
                split(i)%canonical_measure]))) return
            if (.not. split(i)%allowed_interval) return
            if (split(i)%canonical_measure <= 0.0_dp) return
            canonical_scale = max(canonical_scale, &
                abs(split(i)%psi_star_min))
            canonical_scale = max(canonical_scale, &
                abs(split(i)%psi_star_max))
            canonical_scale = max(canonical_scale, &
                abs(split(i)%canonical_measure))
        end do
        canonical_tolerance = 100.0_dp*adapter%options%allowed_tolerance &
            *max(1.0_dp, canonical_scale)
        if (.not. ieee_is_finite(canonical_tolerance)) return
        previous_end = candidate%rc_min
        measure_sum = 0.0_dp
        do i = 1, size(split)
            if (split(i)%sigma /= candidate%sigma) return
            if (split(i)%component_id <= 0) return
            if (split(i)%rc_max <= split(i)%rc_min) return
            if (split(i)%rc_min < candidate%rc_min - tolerance) return
            if (split(i)%rc_max > candidate%rc_max + tolerance) return
            if (i == 1) then
                if (abs(split(i)%rc_min - candidate%rc_min) > tolerance) return
            else
                if (abs(split(i)%rc_min - previous_end) > tolerance) return
                if (abs(split(i)%psi_star_min - split(i - 1)%psi_star_max) &
                        > canonical_tolerance) return
            end if
            call evaluate_gc_cylindrical_class_point(adapter, split(i)%rc_min, &
                split(i)%sigma, point, point_status)
            if (point_status /= GC_CYL_CLASS_SUCCESS) return
            if (.not. point%allowed) return
            call evaluate_gc_cylindrical_class_point(adapter, split(i)%rc_max, &
                split(i)%sigma, point, point_status)
            if (point_status /= GC_CYL_CLASS_SUCCESS) return
            if (.not. point%allowed) return
            previous_end = split(i)%rc_max
            measure_sum = measure_sum + split(i)%canonical_measure
            do j = 1, i - 1
                if (split(j)%component_id == split(i)%component_id) return
            end do
        end do
        if (abs(split(1)%psi_star_min - candidate%psi_star_min) > &
                canonical_tolerance) return
        if (abs(split(size(split))%psi_star_max - candidate%psi_star_max) &
                > canonical_tolerance) return
        if (abs(previous_end - candidate%rc_max) > tolerance) return
        if (.not. ieee_is_finite(measure_sum)) return
        measure_tolerance = 100.0_dp*adapter%options%allowed_tolerance &
            *max(1.0_dp, max(abs(candidate%canonical_measure), &
            abs(measure_sum)))
        if (.not. ieee_is_finite(measure_tolerance)) return
        if (abs(measure_sum - candidate%canonical_measure) > measure_tolerance) return
        if (split(1)%lower_root .neqv. candidate%lower_root) return
        if (split(size(split))%upper_root .neqv. candidate%upper_root) return
        do i = 1, size(split)
            split(i)%topology_certified = .true.
            split(i)%root_isolation_certified = candidate%root_isolation_certified
            split(i)%orbit_return_certified = .true.
        end do
        status = GC_CYL_CLASS_SUCCESS
    end subroutine validate_split

    subroutine append_interval(intervals, count, interval, status)
        type(gc_cylindrical_class_interval_t), allocatable, intent(inout) :: &
            intervals(:)
        integer, intent(inout) :: count
        type(gc_cylindrical_class_interval_t), intent(in) :: interval
        integer, intent(out) :: status

        type(gc_cylindrical_class_interval_t), allocatable :: enlarged(:)

        status = GC_CYL_CLASS_INVALID_INPUT
        if (count < 0) return
        if (.not. allocated(intervals)) then
            allocate(intervals(1))
            intervals(1) = interval
            count = 1
            status = GC_CYL_CLASS_SUCCESS
            return
        end if
        if (count /= size(intervals)) return
        allocate(enlarged(count + 1))
        if (count > 0) enlarged(1:count) = intervals
        enlarged(count + 1) = interval
        call move_alloc(enlarged, intervals)
        count = count + 1
        status = GC_CYL_CLASS_SUCCESS
    end subroutine append_interval

    subroutine copy_intervals(destination, source, count)
        type(gc_cylindrical_class_interval_t), allocatable, intent(out) :: &
            destination(:)
        type(gc_cylindrical_class_interval_t), allocatable, intent(in) :: &
            source(:)
        integer, intent(in) :: count

        if (allocated(destination)) deallocate(destination)
        if (count < 1) return
        allocate(destination(count))
        destination = source(1:count)
    end subroutine copy_intervals

end module neort_gc_cylindrical_class_adapter
