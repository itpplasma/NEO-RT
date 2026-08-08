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
        gc_cylindrical_allowed_region_set_t, &
        gc_cylindrical_root_coordinate_map_t
    use neort_gc_certified_interval_roots, only: &
        GC_INTERVAL_ROOT_EXTREMAL, GC_INTERVAL_ROOT_SIMPLE, &
        GC_INTERVAL_ROOT_TANGENT, GC_INTERVAL_ROOT_TRANSVERSE, &
        gc_interval_root_box_t, gc_interval_t

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
        type(gc_interval_t) :: rc_min_enclosure
        type(gc_interval_t) :: rc_max_enclosure
        type(gc_interval_t) :: psi_star_min_enclosure
        type(gc_interval_t) :: psi_star_max_enclosure
        type(gc_interval_t) :: canonical_measure_enclosure
        type(gc_interval_root_box_t) :: lower_root_certificate
        type(gc_interval_root_box_t) :: upper_root_certificate
        type(gc_cylindrical_root_coordinate_map_t) :: &
            lower_root_coordinate_map
        type(gc_cylindrical_root_coordinate_map_t) :: &
            upper_root_coordinate_map
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
                call interval_from_topology(regions, regions%components(i), &
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
        real(dp) :: measure_lower_sum, measure_upper_sum
        real(dp) :: canonical_scale, canonical_tolerance
        real(dp) :: previous_end
        integer :: root_index

        status = GC_CYL_CLASS_SPLITTER_FAILURE
        ! topology_certified and certificate_method are provider metadata, not
        ! proofs.  The typed verifier and its adapter-owned certificate ID
        ! have already established certification before this structural check.
        if (regions%nroots < 0) return
        if (regions%ncomponents < 0) return
        if (.not. ieee_is_finite(regions%total_canonical_measure)) return
        if (regions%total_canonical_measure < 0.0_dp) return
        if (.not. valid_interval(regions%total_canonical_measure_enclosure)) return
        if (regions%total_canonical_measure_enclosure%lo < 0.0_dp) return
        if (.not. interval_contains( &
                regions%total_canonical_measure_enclosure, &
                regions%total_canonical_measure)) return

        bound_tolerance = 100.0_dp*adapter%options%allowed_tolerance*max(1.0_dp, &
            max(abs(adapter%rc_min), abs(adapter%rc_max)))
        if (.not. ieee_is_finite(bound_tolerance)) return

        if (regions%nroots > 0) then
            if (.not. allocated(regions%roots)) return
            if (.not. allocated(regions%root_canonical)) return
            if (.not. allocated(regions%root_boxes)) return
            if (.not. allocated(regions%root_coordinate_maps)) return
            if (.not. allocated(regions%root_canonical_enclosures)) return
            if (size(regions%roots) /= regions%nroots) return
            if (size(regions%root_canonical) /= regions%nroots) return
            if (size(regions%root_boxes) /= regions%nroots) return
            if (size(regions%root_coordinate_maps) /= regions%nroots) return
            if (size(regions%root_canonical_enclosures) /= regions%nroots) return
        else
            if (allocated(regions%roots)) then
                if (size(regions%roots) /= 0) return
            end if
            if (allocated(regions%root_canonical)) then
                if (size(regions%root_canonical) /= 0) return
            end if
            if (allocated(regions%root_boxes)) then
                if (size(regions%root_boxes) /= 0) return
            end if
            if (allocated(regions%root_coordinate_maps)) then
                if (size(regions%root_coordinate_maps) /= 0) return
            end if
            if (allocated(regions%root_canonical_enclosures)) then
                if (size(regions%root_canonical_enclosures) /= 0) return
            end if
        end if
        if (regions%nroots > 0) then
            if (.not. all(ieee_is_finite(regions%roots))) return
            if (.not. all(ieee_is_finite(regions%root_canonical))) return
            do i = 1, regions%nroots
                if (.not. valid_certified_root_box(regions%root_boxes(i), &
                        adapter%rc_min, adapter%rc_max)) return
                if (.not. valid_root_coordinate_map( &
                        regions%root_coordinate_maps(i), &
                        regions%root_boxes(i), regions%roots(i))) return
                if (.not. interval_contains(gc_interval_t( &
                        regions%root_boxes(i)%lo, regions%root_boxes(i)%hi), &
                        regions%roots(i))) return
                if (.not. valid_interval( &
                        regions%root_canonical_enclosures(i))) return
                if (.not. interval_contains( &
                        regions%root_canonical_enclosures(i), &
                        regions%root_canonical(i))) return
                if (i > 1) then
                    if (regions%roots(i) <= regions%roots(i - 1)) return
                    if (regions%root_boxes(i)%lo <= &
                            regions%root_boxes(i - 1)%hi) return
                    if (.not. consistent_root_coordinate_maps( &
                            regions%root_coordinate_maps(i - 1), &
                            regions%root_coordinate_maps(i))) return
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
            if (.not. interval_contains( &
                    regions%total_canonical_measure_enclosure, 0.0_dp)) return
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
                regions%components(i)%canonical_measure, &
                regions%components(i)%canonical_measure_lower, &
                regions%components(i)%canonical_measure_upper]))) return
            if (regions%components(i)%x_end <= regions%components(i)%x_begin) return
            if (regions%components(i)%canonical_measure <= 0.0_dp) return
            if (regions%components(i)%canonical_measure_lower <= 0.0_dp) return
            if (regions%components(i)%canonical_measure_upper < &
                    regions%components(i)%canonical_measure_lower) return
            if (regions%components(i)%canonical_measure < &
                    regions%components(i)%canonical_measure_lower) return
            if (regions%components(i)%canonical_measure > &
                    regions%components(i)%canonical_measure_upper) return
            if (regions%components(i)%lower_root) then
                root_index = regions%components(i)%lower_root_index
                if (root_index < 1 .or. root_index > regions%nroots) return
                if (.not. interval_contains(gc_interval_t( &
                        regions%root_boxes(root_index)%lo, &
                        regions%root_boxes(root_index)%hi), &
                        regions%components(i)%x_begin)) return
                if (.not. interval_contains( &
                        regions%root_canonical_enclosures(root_index), &
                        regions%components(i)%canonical_begin)) return
            else
                if (regions%components(i)%lower_root_index /= 0) return
            end if
            if (regions%components(i)%upper_root) then
                root_index = regions%components(i)%upper_root_index
                if (root_index < 1 .or. root_index > regions%nroots) return
                if (.not. interval_contains(gc_interval_t( &
                        regions%root_boxes(root_index)%lo, &
                        regions%root_boxes(root_index)%hi), &
                        regions%components(i)%x_end)) return
                if (.not. interval_contains( &
                        regions%root_canonical_enclosures(root_index), &
                        regions%components(i)%canonical_end)) return
            else
                if (regions%components(i)%upper_root_index /= 0) return
            end if
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
        measure_lower_sum = 0.0_dp
        measure_upper_sum = 0.0_dp
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
                root_index = regions%components(i)%lower_root_index
                if (abs(regions%root_canonical(root_index) - &
                        regions%components(i)%canonical_begin) > &
                        canonical_tolerance) return
            end if
            if (regions%components(i)%upper_root) then
                root_index = regions%components(i)%upper_root_index
                if (abs(regions%root_canonical(root_index) - &
                        regions%components(i)%canonical_end) > &
                        canonical_tolerance) return
            end if
            measure_sum = measure_sum + regions%components(i)%canonical_measure
            measure_lower_sum = measure_lower_sum &
                +regions%components(i)%canonical_measure_lower
            measure_upper_sum = measure_upper_sum &
                +regions%components(i)%canonical_measure_upper
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
        if (regions%total_canonical_measure_enclosure%lo > &
                measure_lower_sum) return
        if (regions%total_canonical_measure_enclosure%hi < &
                measure_upper_sum) return
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
        if (adapter%intervals(i)%lower_root) then
            if (rc < adapter%intervals(i)%lower_root_certificate%hi) then
                status = GC_CYL_CLASS_NOT_ALLOWED
                launch%status = status
                return
            end if
        end if
        if (adapter%intervals(i)%upper_root) then
            if (rc > adapter%intervals(i)%upper_root_certificate%lo) then
                status = GC_CYL_CLASS_NOT_ALLOWED
                launch%status = status
                return
            end if
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
        launch%endpoint_tangent = .false.
        if (adapter%intervals(i)%lower_tangent) then
            launch%endpoint_tangent = &
                adapter%intervals(i)%lower_root_certificate%lo == rc
            launch%endpoint_tangent = launch%endpoint_tangent .and. &
                adapter%intervals(i)%lower_root_certificate%hi == rc
        end if
        if (adapter%intervals(i)%upper_tangent) then
            if (adapter%intervals(i)%upper_root_certificate%lo == rc .and. &
                    adapter%intervals(i)%upper_root_certificate%hi == rc) then
                launch%endpoint_tangent = .true.
            end if
        end if
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

    pure logical function valid_interval(value)
        type(gc_interval_t), intent(in) :: value

        valid_interval = ieee_is_finite(value%lo)
        valid_interval = valid_interval .and. ieee_is_finite(value%hi)
        valid_interval = valid_interval .and. value%lo <= value%hi
    end function valid_interval

    pure logical function interval_contains(interval, value)
        type(gc_interval_t), intent(in) :: interval
        real(dp), intent(in) :: value

        interval_contains = valid_interval(interval)
        interval_contains = interval_contains .and. ieee_is_finite(value)
        interval_contains = interval_contains .and. value >= interval%lo
        interval_contains = interval_contains .and. value <= interval%hi
    end function interval_contains

    pure logical function valid_certified_root_box(root, domain_lo, domain_hi)
        type(gc_interval_root_box_t), intent(in) :: root
        real(dp), intent(in) :: domain_lo, domain_hi

        logical :: derivative_excludes_zero, second_derivative_excludes_zero

        valid_certified_root_box = .false.
        if (.not. all(ieee_is_finite([root%lo, root%hi, &
                root%derivative_enclosure%lo, &
                root%derivative_enclosure%hi, &
                root%second_derivative_enclosure%lo, &
                root%second_derivative_enclosure%hi]))) return
        if (root%lo > root%hi) return
        if (root%lo < domain_lo .or. root%hi > domain_hi) return
        if (root%enclosure_certificate_id <= 0) return
        if (.not. root%classification_certified) return
        if (.not. root%bracket_certified .and. &
                .not. root%interval_newton_certified) return
        if (.not. valid_interval(root%derivative_enclosure)) return
        if (.not. valid_interval(root%second_derivative_enclosure)) return
        derivative_excludes_zero = root%derivative_enclosure%hi < 0.0_dp
        derivative_excludes_zero = derivative_excludes_zero .or. &
            root%derivative_enclosure%lo > 0.0_dp
        second_derivative_excludes_zero = &
            root%second_derivative_enclosure%hi < 0.0_dp
        second_derivative_excludes_zero = &
            second_derivative_excludes_zero .or. &
            root%second_derivative_enclosure%lo > 0.0_dp

        select case (root%kind)
            case (GC_INTERVAL_ROOT_SIMPLE)
                if (root%multiplicity_lower /= 1) return
                if (root%multiplicity_upper /= 1) return
                if (.not. root%derivative_excludes_zero) return
                if (.not. derivative_excludes_zero) return
                if (.not. root%transversality_certified) return
                if (root%transversality_kind /= &
                        GC_INTERVAL_ROOT_TRANSVERSE) return
            case (GC_INTERVAL_ROOT_TANGENT)
                if (root%multiplicity_lower /= 2) return
                if (root%multiplicity_upper /= 2) return
                if (.not. root%stationary_certified) return
                if (root%stationary_certificate_id <= 0) return
                if (.not. second_derivative_excludes_zero) return
                if (root%transversality_certified) return
                if (root%transversality_kind /= &
                        GC_INTERVAL_ROOT_EXTREMAL) return
                if (root%lo /= root%hi) return
            case default
                return
        end select
        valid_certified_root_box = .true.
    end function valid_certified_root_box

    pure logical function valid_root_coordinate_map(mapping, class_root, &
            class_representative)
        type(gc_cylindrical_root_coordinate_map_t), intent(in) :: mapping
        type(gc_interval_root_box_t), intent(in) :: class_root
        real(dp), intent(in) :: class_representative

        integer :: class_derivative_sign, source_derivative_sign
        integer :: class_second_sign, source_second_sign

        valid_root_coordinate_map = .false.
        if (.not. valid_interval(mapping%source_domain_enclosure)) return
        if (.not. valid_interval(mapping%mapped_class_enclosure)) return
        if (.not. valid_certified_root_box(mapping%source_root_certificate, &
                mapping%source_domain_enclosure%lo, &
                mapping%source_domain_enclosure%hi)) return
        if (mapping%map_certificate_id <= 0) return
        if (abs(mapping%monotonicity_sign) /= 1) return
        if (.not. mapping%strict_monotonicity_certified) return
        if (.not. mapping%mapping_enclosure_certified) return
        if (len_trim(mapping%source_coordinate) == 0) return
        if (len_trim(mapping%source_units) == 0) return
        if (len_trim(mapping%class_coordinate) == 0) return
        if (len_trim(mapping%class_units) == 0) return
        if (mapping%mapped_class_enclosure%lo /= class_root%lo) return
        if (mapping%mapped_class_enclosure%hi /= class_root%hi) return
        if (.not. interval_contains(mapping%mapped_class_enclosure, &
                class_representative)) return
        if (mapping%source_root_certificate%kind /= class_root%kind) return
        if (mapping%source_root_certificate%multiplicity_lower /= &
                class_root%multiplicity_lower) return
        if (mapping%source_root_certificate%multiplicity_upper /= &
                class_root%multiplicity_upper) return
        if (mapping%source_root_certificate%transversality_kind /= &
                class_root%transversality_kind) return

        select case (class_root%kind)
            case (GC_INTERVAL_ROOT_SIMPLE)
                class_derivative_sign = interval_strict_sign( &
                    class_root%derivative_enclosure)
                source_derivative_sign = interval_strict_sign( &
                    mapping%source_root_certificate%derivative_enclosure)
                if (class_derivative_sign == 0) return
                if (source_derivative_sign == 0) return
                if (class_derivative_sign /= source_derivative_sign &
                        *mapping%monotonicity_sign) return
            case (GC_INTERVAL_ROOT_TANGENT)
                class_second_sign = interval_strict_sign( &
                    class_root%second_derivative_enclosure)
                source_second_sign = interval_strict_sign( &
                    mapping%source_root_certificate%second_derivative_enclosure)
                if (class_second_sign == 0) return
                if (source_second_sign == 0) return
                if (class_second_sign /= source_second_sign) return
            case default
                return
        end select
        valid_root_coordinate_map = .true.
    end function valid_root_coordinate_map

    pure logical function consistent_root_coordinate_maps(previous, current)
        type(gc_cylindrical_root_coordinate_map_t), intent(in) :: previous
        type(gc_cylindrical_root_coordinate_map_t), intent(in) :: current

        consistent_root_coordinate_maps = .false.
        if (current%map_certificate_id /= previous%map_certificate_id) return
        if (current%monotonicity_sign /= previous%monotonicity_sign) return
        if (current%source_domain_enclosure%lo /= &
                previous%source_domain_enclosure%lo) return
        if (current%source_domain_enclosure%hi /= &
                previous%source_domain_enclosure%hi) return
        if (trim(current%source_coordinate) /= &
                trim(previous%source_coordinate)) return
        if (trim(current%source_units) /= trim(previous%source_units)) return
        if (trim(current%class_coordinate) /= &
                trim(previous%class_coordinate)) return
        if (trim(current%class_units) /= trim(previous%class_units)) return
        if (current%monotonicity_sign > 0) then
            if (current%source_root_certificate%lo <= &
                    previous%source_root_certificate%hi) return
        else
            if (current%source_root_certificate%hi >= &
                    previous%source_root_certificate%lo) return
        end if
        consistent_root_coordinate_maps = .true.
    end function consistent_root_coordinate_maps

    pure integer function interval_strict_sign(value)
        type(gc_interval_t), intent(in) :: value

        interval_strict_sign = 0
        if (value%lo > 0.0_dp) interval_strict_sign = 1
        if (value%hi < 0.0_dp) interval_strict_sign = -1
    end function interval_strict_sign

    subroutine interval_from_topology(regions, component, interval, status, &
            certified)
        type(gc_cylindrical_allowed_region_set_t), intent(in) :: regions
        type(gc_cylindrical_allowed_component_t), intent(in) :: component
        type(gc_cylindrical_class_interval_t), intent(out) :: interval
        integer, intent(out) :: status
        logical, intent(in) :: certified

        integer :: root_index

        interval = gc_cylindrical_class_interval_t()
        interval%component_id = component%component_id
        interval%sigma = component%sigma
        interval%rc_min = component%x_begin
        interval%rc_max = component%x_end
        interval%psi_star_min = component%canonical_begin
        interval%psi_star_max = component%canonical_end
        interval%canonical_measure = component%canonical_measure
        interval%rc_min_enclosure = gc_interval_t(component%x_begin, &
            component%x_begin)
        interval%rc_max_enclosure = gc_interval_t(component%x_end, &
            component%x_end)
        interval%psi_star_min_enclosure = gc_interval_t( &
            component%canonical_begin, component%canonical_begin)
        interval%psi_star_max_enclosure = gc_interval_t( &
            component%canonical_end, component%canonical_end)
        interval%canonical_measure_enclosure = gc_interval_t( &
            component%canonical_measure_lower, &
            component%canonical_measure_upper)
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
                root_index = component%lower_root_index
                interval%lower_root_certificate = &
                    regions%root_boxes(root_index)
                interval%lower_root_coordinate_map = &
                    regions%root_coordinate_maps(root_index)
                interval%rc_min_enclosure = gc_interval_t( &
                    regions%root_boxes(root_index)%lo, &
                    regions%root_boxes(root_index)%hi)
                interval%psi_star_min_enclosure = &
                    regions%root_canonical_enclosures(root_index)
                interval%lower_tangent = regions%root_boxes(root_index)%kind &
                    == GC_INTERVAL_ROOT_TANGENT
                if (interval%lower_tangent) then
                    interval%lower_boundary_kind = 'certified-tangent-root'
                else
                    interval%lower_boundary_kind = 'certified-simple-root'
                end if
            else
                interval%lower_boundary_kind = 'certified-bound'
            end if
            if (component%upper_root) then
                root_index = component%upper_root_index
                interval%upper_root_certificate = &
                    regions%root_boxes(root_index)
                interval%upper_root_coordinate_map = &
                    regions%root_coordinate_maps(root_index)
                interval%rc_max_enclosure = gc_interval_t( &
                    regions%root_boxes(root_index)%lo, &
                    regions%root_boxes(root_index)%hi)
                interval%psi_star_max_enclosure = &
                    regions%root_canonical_enclosures(root_index)
                interval%upper_tangent = regions%root_boxes(root_index)%kind &
                    == GC_INTERVAL_ROOT_TANGENT
                if (interval%upper_tangent) then
                    interval%upper_boundary_kind = 'certified-tangent-root'
                else
                    interval%upper_boundary_kind = 'certified-simple-root'
                end if
            else
                interval%upper_boundary_kind = 'certified-bound'
            end if
            interval%limiting_chart = 'fortsym-certified-root-box'
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
        if (candidate%lower_root) then
            if (.not. same_root_coordinate_map( &
                    split(1)%lower_root_coordinate_map, &
                    candidate%lower_root_coordinate_map)) return
        end if
        if (candidate%upper_root) then
            if (.not. same_root_coordinate_map( &
                    split(size(split))%upper_root_coordinate_map, &
                    candidate%upper_root_coordinate_map)) return
        end if
        do i = 2, size(split)
            if (split(i)%lower_root) return
        end do
        do i = 1, size(split) - 1
            if (split(i)%upper_root) return
        end do
        do i = 1, size(split)
            split(i)%topology_certified = .true.
            split(i)%root_isolation_certified = candidate%root_isolation_certified
            split(i)%orbit_return_certified = .true.
        end do
        status = GC_CYL_CLASS_SUCCESS
    end subroutine validate_split

    pure logical function same_root_coordinate_map(left, right)
        type(gc_cylindrical_root_coordinate_map_t), intent(in) :: left, right

        same_root_coordinate_map = same_root_box( &
            left%source_root_certificate, right%source_root_certificate)
        same_root_coordinate_map = same_root_coordinate_map .and. &
            same_interval(left%source_domain_enclosure, &
                right%source_domain_enclosure)
        same_root_coordinate_map = same_root_coordinate_map .and. &
            same_interval(left%mapped_class_enclosure, &
                right%mapped_class_enclosure)
        same_root_coordinate_map = same_root_coordinate_map .and. &
            left%map_certificate_id == right%map_certificate_id
        same_root_coordinate_map = same_root_coordinate_map .and. &
            left%monotonicity_sign == right%monotonicity_sign
        same_root_coordinate_map = same_root_coordinate_map .and. &
            (left%strict_monotonicity_certified .eqv. &
                right%strict_monotonicity_certified)
        same_root_coordinate_map = same_root_coordinate_map .and. &
            (left%mapping_enclosure_certified .eqv. &
                right%mapping_enclosure_certified)
        same_root_coordinate_map = same_root_coordinate_map .and. &
            trim(left%source_coordinate) == trim(right%source_coordinate)
        same_root_coordinate_map = same_root_coordinate_map .and. &
            trim(left%source_units) == trim(right%source_units)
        same_root_coordinate_map = same_root_coordinate_map .and. &
            trim(left%class_coordinate) == trim(right%class_coordinate)
        same_root_coordinate_map = same_root_coordinate_map .and. &
            trim(left%class_units) == trim(right%class_units)
    end function same_root_coordinate_map

    pure logical function same_root_box(left, right)
        type(gc_interval_root_box_t), intent(in) :: left, right

        same_root_box = left%lo == right%lo .and. left%hi == right%hi
        same_root_box = same_root_box .and. left%kind == right%kind
        same_root_box = same_root_box .and. left%cut_id == right%cut_id
        same_root_box = same_root_box .and. &
            left%boundary_role == right%boundary_role
        same_root_box = same_root_box .and. &
            left%multiplicity_lower == right%multiplicity_lower
        same_root_box = same_root_box .and. &
            left%multiplicity_upper == right%multiplicity_upper
        same_root_box = same_root_box .and. &
            (left%derivative_excludes_zero .eqv. &
                right%derivative_excludes_zero)
        same_root_box = same_root_box .and. &
            (left%stationary_certified .eqv. right%stationary_certified)
        same_root_box = same_root_box .and. &
            (left%bracket_certified .eqv. right%bracket_certified)
        same_root_box = same_root_box .and. &
            (left%classification_certified .eqv. &
                right%classification_certified)
        same_root_box = same_root_box .and. &
            (left%transversality_certified .eqv. &
                right%transversality_certified)
        same_root_box = same_root_box .and. &
            left%transversality_kind == right%transversality_kind
        same_root_box = same_root_box .and. &
            (left%interval_newton_certified .eqv. &
                right%interval_newton_certified)
        same_root_box = same_root_box .and. &
            (left%left_endpoint_root .eqv. right%left_endpoint_root)
        same_root_box = same_root_box .and. &
            (left%right_endpoint_root .eqv. right%right_endpoint_root)
        same_root_box = same_root_box .and. &
            left%enclosure_certificate_id == right%enclosure_certificate_id
        same_root_box = same_root_box .and. &
            left%stationary_certificate_id == right%stationary_certificate_id
        same_root_box = same_root_box .and. &
            same_interval(left%derivative_enclosure, &
                right%derivative_enclosure)
        same_root_box = same_root_box .and. &
            same_interval(left%second_derivative_enclosure, &
                right%second_derivative_enclosure)
    end function same_root_box

    pure logical function same_interval(left, right)
        type(gc_interval_t), intent(in) :: left, right

        same_interval = left%lo == right%lo .and. left%hi == right%hi
    end function same_interval

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
