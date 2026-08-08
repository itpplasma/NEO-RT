module neort_gc_eqdsk_certified_allowed_provider
    !! Certified fixed-(H0,J_K) allowed components on the complete physical-R cut.
    !!
    !! This module owns no physics formula.  It delegates every value and jet
    !! to the generated EQDSK interval evaluator, and uses the generic
    !! interval/root engine for coverage.  A component is returned only when
    !! its canonical measure is closed by generated interval data.  Turning
    !! and tangent limits for which that data is not available fail closed.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_next_after
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_gc_callback_context, only: gc_callback_context_t
    use neort_gc_cylindrical_model, only: &
        gc_cylindrical_allowed_component_t
    use neort_gc_cylindrical_topology, only: &
        gc_cylindrical_allowed_region_set_t, &
        gc_cylindrical_root_coordinate_map_t
    use neort_gc_certified_interval_roots, only: &
        GC_INTERVAL_ROOT_SUCCESS, gc_interval_callback_result_t, &
        gc_interval_enclosure_verifier_i, gc_interval_root_box_t, &
        gc_interval_root_callback_i, gc_interval_root_options_t, &
        gc_interval_root_result_t, gc_interval_stationary_verifier_i, &
        gc_interval_t, isolate_gc_interval_roots
    use neort_gc_eqdsk_allowed_region_cut_box, only: &
        EQDSK_CUT_BOX_SUCCESS, eqdsk_allowed_region_cut_provenance_t, &
        eqdsk_potential_profile_nodes_t, &
        evaluate_eqdsk_allowed_region_cut_box
    use neort_gc_eqdsk_allowed_region_interval, only: &
        eqdsk_allowed_interval_result_t
    use neort_gc_eqdsk_composite_cut_atlas, only: &
        eqdsk_composite_cut_atlas_t
    use neort_gc_eqdsk_cut_graph_atlas, only: &
        eqdsk_cut_graph_atlas_t
    use neort_gc_eqdsk_composite_r_ownership, only: &
        EQDSK_R_OWNERSHIP_SUCCESS, eqdsk_composite_r_partition_t, &
        validate_eqdsk_composite_r_partition
    use neort_gc_outward_interval, only: &
        gc_outward_interval, gc_outward_interval_t, &
        gc_outward_interval_is_valid
    use neort_generated_certificate_registry, only: certificate_index

    implicit none
    private

    abstract interface
        subroutine gc_eqdsk_stationary_root_certificate_i(x_lo, x_hi, h0, &
                j_k, sigma, branch, stationary_point, certificate_id, status)
            import dp
            real(dp), intent(in) :: x_lo, x_hi, h0, j_k
            integer, intent(in) :: sigma, branch
            real(dp), intent(out) :: stationary_point
            integer, intent(out) :: certificate_id, status
        end subroutine gc_eqdsk_stationary_root_certificate_i
    end interface

    integer, parameter, public :: GC_EQDSK_ALLOWED_PROVIDER_SUCCESS = 0
    integer, parameter, public :: GC_EQDSK_ALLOWED_PROVIDER_INVALID_INPUT = 1
    integer, parameter, public :: GC_EQDSK_ALLOWED_PROVIDER_ROOT_FAILURE = 2
    integer, parameter, public :: GC_EQDSK_ALLOWED_PROVIDER_MEASURE_FAILURE = 3
    integer, parameter, public :: GC_EQDSK_ALLOWED_PROVIDER_CERTIFICATE_FAILURE = 4
    integer, parameter, public :: GC_EQDSK_ALLOWED_PROVIDER_CERTIFICATE_ID = 130031

    type, public, extends(gc_callback_context_t) :: &
            gc_eqdsk_certified_allowed_provider_context_t
        type(eqdsk_cut_graph_atlas_t) :: graph(3)
        type(eqdsk_composite_r_partition_t) :: partition
        type(eqdsk_potential_profile_nodes_t) :: profile
        real(dp) :: field_scale = 0.0_dp
        real(dp) :: raw_psi_sep = 0.0_dp
        real(dp) :: mass = 0.0_dp
        real(dp) :: charge = 0.0_dp
        real(dp) :: c_light = 0.0_dp
        integer :: root_initial_partition = 2
        integer :: root_max_depth = 24
        integer :: root_max_boxes = 8192
        integer :: root_max_roots = 1024
        real(dp) :: root_x_tolerance = 1.0e-10_dp
        integer :: measure_max_depth = 20
        !! Zero means that no generated fixed-(H0,J_K) stationary
        !! certificate is available.  A nonzero value is accepted only when
        !! it is the distinct registered stationary-kernel ID and the point
        !! query below closes exact f=f'=0 evidence.
        integer :: stationary_certificate_id = 0
        procedure(gc_eqdsk_stationary_root_certificate_i), pointer, nopass :: &
            stationary_root_certificate => null()
        logical :: initialized = .false.
    end type gc_eqdsk_certified_allowed_provider_context_t

    public :: initialize_gc_eqdsk_certified_allowed_provider
    public :: build_gc_eqdsk_certified_allowed_regions
    public :: verify_gc_eqdsk_certified_allowed_regions
    public :: gc_eqdsk_stationary_root_certificate_i
    public :: verify_gc_eqdsk_fixed_invariant_stationary_certificate

contains

    subroutine verify_gc_eqdsk_fixed_invariant_stationary_certificate( &
            x_lo, x_hi, candidate, point_value, expected_enclosure_id, &
            expected_stationary_id, status)
        !! Validate the exact point half of the stationary-root contract.
        !!
        !! The interval candidate only isolates a stationary point and
        !! encloses f.  It is not a tangent-root certificate.  Promotion is
        !! possible only when the generated fixed-invariant point query
        !! returns the same stationary identity and exact zero f and f'.
        !! This adapter deliberately does not turn an interval containing
        !! zero into an exact zero.
        real(dp), intent(in) :: x_lo, x_hi
        type(gc_interval_callback_result_t), intent(in) :: candidate
        type(gc_interval_callback_result_t), intent(in) :: point_value
        integer, intent(in) :: expected_enclosure_id, expected_stationary_id
        integer, intent(out) :: status
        integer :: registered_stationary_id

        status = GC_EQDSK_ALLOWED_PROVIDER_CERTIFICATE_FAILURE
        if (.not. all(ieee_is_finite([x_lo, x_hi]))) return
        if (x_hi < x_lo .or. expected_enclosure_id <= 0 .or. &
                expected_stationary_id <= 0) return
        registered_stationary_id = certificate_index( &
            'eqdsk_allowed_energy_stationary')
        if (registered_stationary_id <= 0 .or. &
                expected_stationary_id /= registered_stationary_id) return
        if (candidate%query_lo /= x_lo .or. candidate%query_hi /= x_hi) return
        if (candidate%status /= 0 .or. point_value%status /= 0) return
        if (candidate%enclosure_certificate_id /= expected_enclosure_id .or. &
                candidate%stationary_certificate_id /= expected_stationary_id) return
        if (.not. valid_interval(candidate%f) .or. &
                .not. valid_interval(candidate%df) .or. &
                .not. valid_interval(candidate%d2f)) return
        if (.not. ieee_is_finite(candidate%stationary_point) .or. &
                candidate%stationary_point < x_lo .or. &
                candidate%stationary_point > x_hi) return
        if (point_value%query_lo /= candidate%stationary_point .or. &
                point_value%query_hi /= candidate%stationary_point) return
        if (point_value%enclosure_certificate_id /= expected_enclosure_id .or. &
                point_value%stationary_certificate_id /= expected_stationary_id .or. &
                point_value%stationary_point /= candidate%stationary_point) return
        if (.not. valid_interval(point_value%f) .or. &
                .not. valid_interval(point_value%df) .or. &
                .not. valid_interval(point_value%d2f)) return
        if (.not. exact_zero(point_value%f) .or. &
                .not. exact_zero(point_value%df) .or. &
                .not. excludes_zero(candidate%d2f) .or. &
                .not. excludes_zero(point_value%d2f)) return
        status = GC_EQDSK_ALLOWED_PROVIDER_SUCCESS
    end subroutine verify_gc_eqdsk_fixed_invariant_stationary_certificate

    subroutine initialize_gc_eqdsk_certified_allowed_provider(context, atlas, &
            partition, profile, field_scale, raw_psi_sep, mass, charge, &
            c_light, status, stationary_root_certificate)
        type(gc_eqdsk_certified_allowed_provider_context_t), intent(out) :: context
        type(eqdsk_composite_cut_atlas_t), intent(in) :: atlas
        type(eqdsk_composite_r_partition_t), intent(in) :: partition
        type(eqdsk_potential_profile_nodes_t), intent(in) :: profile
        real(dp), intent(in) :: field_scale, raw_psi_sep, mass, charge, c_light
        integer, intent(out) :: status
        procedure(gc_eqdsk_stationary_root_certificate_i), optional :: &
            stationary_root_certificate

        integer :: branch, local_status

        context = gc_eqdsk_certified_allowed_provider_context_t()
        status = GC_EQDSK_ALLOWED_PROVIDER_INVALID_INPUT
        call validate_eqdsk_composite_r_partition(partition, local_status)
        if (local_status /= EQDSK_R_OWNERSHIP_SUCCESS) return
        if (.not. atlas%geometric_completeness_certified .or. &
                .not. atlas%branch_connectivity_certified .or. &
                .not. atlas%surface_intersection_pair_certified) return
        if (.not. profile%structurally_validated) return
        if (.not. allocated(profile%psi) .or. &
                .not. allocated(profile%phi) .or. &
                .not. allocated(profile%omega)) return
        if (.not. all(ieee_is_finite([field_scale, raw_psi_sep, mass, &
                charge, c_light]))) return
        if (field_scale <= 0.0_dp .or. raw_psi_sep <= 0.0_dp .or. &
                mass <= 0.0_dp .or. c_light <= 0.0_dp .or. &
                abs(charge) <= tiny(charge)) return

        context%graph(1) = atlas%inboard_graph
        context%graph(2) = atlas%axis_graph
        context%graph(3) = atlas%outboard_graph
        context%partition = partition
        context%profile = profile
        context%field_scale = field_scale
        context%raw_psi_sep = raw_psi_sep
        context%mass = mass
        context%charge = charge
        context%c_light = c_light
        if (present(stationary_root_certificate)) then
            context%stationary_root_certificate => stationary_root_certificate
        end if
        do branch = 1, 3
            if (context%graph(branch)%certificate_id <= 0) return
            if (.not. context%graph(branch)%global_completeness_certified) return
        end do
        context%initialized = .true.
        status = GC_EQDSK_ALLOWED_PROVIDER_SUCCESS
    end subroutine initialize_gc_eqdsk_certified_allowed_provider

    subroutine build_gc_eqdsk_certified_allowed_regions(context, h0, j_k, sigma, &
            regions, status)
        type(gc_eqdsk_certified_allowed_provider_context_t), intent(inout) :: context
        real(dp), intent(in) :: h0, j_k
        integer, intent(in) :: sigma
        type(gc_cylindrical_allowed_region_set_t), intent(out) :: regions
        integer, intent(out) :: status

        type(gc_interval_root_result_t) :: root_results(3)
        type(gc_interval_root_options_t) :: root_options
        integer :: branch, local_status, root_count, component_count
        integer :: energy_certificate_id
        integer :: root_offset
        integer :: active_branch, active_sigma
        real(dp) :: active_h0, active_j_k, active_lo, active_hi
        type(gc_cylindrical_allowed_component_t), allocatable :: components(:)
        type(gc_interval_root_box_t), allocatable :: roots(:)
        type(gc_interval_t), allocatable :: root_canonical_enclosures(:)
        real(dp), allocatable :: root_coordinates(:), root_canonical(:)
        type(gc_cylindrical_root_coordinate_map_t), allocatable :: root_maps(:)

        regions = gc_cylindrical_allowed_region_set_t()
        status = GC_EQDSK_ALLOWED_PROVIDER_INVALID_INPUT
        if (.not. context%initialized) return
        if (.not. all(ieee_is_finite([h0, j_k]))) return
        if (j_k < 0.0_dp .or. abs(sigma) /= 1) return

        energy_certificate_id = certificate_index('eqdsk_allowed_energy')
        if (energy_certificate_id <= 0) then
            status = GC_EQDSK_ALLOWED_PROVIDER_CERTIFICATE_FAILURE
            return
        end if
        root_options = gc_interval_root_options_t()
        root_options%initial_partition = context%root_initial_partition
        root_options%max_depth = context%root_max_depth
        root_options%max_boxes = context%root_max_boxes
        root_options%max_roots = context%root_max_roots
        root_options%expected_enclosure_certificate_id = energy_certificate_id
        context%stationary_certificate_id = certificate_index( &
            'eqdsk_allowed_energy_stationary')
        if (context%stationary_certificate_id < 0) then
            status = GC_EQDSK_ALLOWED_PROVIDER_CERTIFICATE_FAILURE
            return
        end if
        ! The default remains zero until a real generated fixed-(H0,J_K)
        ! certificate is registered.  The generic engine then fails tangent
        ! candidates closed before any interval extremum can be promoted.
        root_options%expected_stationary_certificate_id = &
            context%stationary_certificate_id
        root_options%x_tolerance = context%root_x_tolerance

        root_count = 0
        do branch = 1, 3
            active_branch = branch
            active_sigma = sigma
            active_h0 = h0
            active_j_k = j_k
            call branch_domain(context, branch, active_lo, active_hi)
            call isolate_gc_interval_roots(evaluate_energy, verify_energy, &
                verify_energy_stationary, active_lo, active_hi, root_options, &
                root_results(branch))
            if (root_results(branch)%status /= GC_INTERVAL_ROOT_SUCCESS .or. &
                    .not. root_results(branch)%coverage_certified) then
                write (*, '(a,6(1x,i0),2(1x,es24.16))') &
                    'allowed root diagnostic=', branch, root_results(branch)%status, &
                    root_results(branch)%nroots, root_results(branch)%boxes_visited, &
                    root_results(branch)%boxes_discarded, &
                    root_results(branch)%unresolved_boxes, active_lo, active_hi
                status = GC_EQDSK_ALLOWED_PROVIDER_ROOT_FAILURE
                return
            end if
            root_count = root_count + root_results(branch)%nroots
        end do

        if (root_count > 0) then
            allocate(roots(root_count), root_coordinates(root_count), &
                root_canonical(root_count), root_canonical_enclosures(root_count), &
                root_maps(root_count))
        else
            allocate(roots(0), root_coordinates(0), root_canonical(0), &
                root_canonical_enclosures(0), root_maps(0))
        end if

        root_offset = 0
        do branch = 1, 3
            call copy_branch_roots(context, h0, j_k, sigma, branch, &
                root_results(branch), root_offset, roots, root_coordinates, &
                root_canonical, root_canonical_enclosures, root_maps, local_status)
            if (local_status /= GC_EQDSK_ALLOWED_PROVIDER_SUCCESS) then
                status = local_status
                return
            end if
        end do

        component_count = 0
        allocate(components(0))
        root_offset = 0
        do branch = 1, 3
            call build_branch_components(context, h0, j_k, sigma, branch, &
                root_results(branch), root_offset, components, component_count, &
                local_status)
            if (local_status /= GC_EQDSK_ALLOWED_PROVIDER_SUCCESS) then
                status = local_status
                return
            end if
            root_offset = root_offset + root_results(branch)%nroots
        end do
        if (component_count == 0) then
            status = GC_EQDSK_ALLOWED_PROVIDER_MEASURE_FAILURE
            return
        end if

        regions%nroots = root_count
        regions%ncomponents = component_count
        call move_alloc(root_coordinates, regions%roots)
        call move_alloc(root_canonical, regions%root_canonical)
        call move_alloc(roots, regions%root_boxes)
        call move_alloc(root_maps, regions%root_coordinate_maps)
        call move_alloc(root_canonical_enclosures, regions%root_canonical_enclosures)
        call move_alloc(components, regions%components)
        regions%total_canonical_measure = 0.0_dp
        regions%total_canonical_measure_enclosure = gc_interval_t(0.0_dp, 0.0_dp)
        do branch = 1, regions%ncomponents
            regions%total_canonical_measure = regions%total_canonical_measure + &
                regions%components(branch)%canonical_measure
            regions%total_canonical_measure_enclosure = add_interval( &
                regions%total_canonical_measure_enclosure, &
                regions%components(branch)%canonical_measure_enclosure)
        end do
        if (.not. ieee_is_finite(regions%total_canonical_measure) .or. &
                .not. valid_interval(regions%total_canonical_measure_enclosure)) then
            status = GC_EQDSK_ALLOWED_PROVIDER_MEASURE_FAILURE
            return
        end if
        regions%topology_certified = .true.
        regions%certificate_method = 'fortsym-energy-interval-root-and-measure'
        status = GC_EQDSK_ALLOWED_PROVIDER_SUCCESS

    contains

        subroutine branch_domain(ctx, branch_id, lo, hi)
            type(gc_eqdsk_certified_allowed_provider_context_t), intent(in) :: ctx
            integer, intent(in) :: branch_id
            real(dp), intent(out) :: lo, hi

            lo = ctx%partition%branches(branch_id)%r_lo
            hi = ctx%partition%branches(branch_id)%r_hi
            if (branch_id < 3) hi = ieee_next_after(hi, -huge(hi))
        end subroutine branch_domain

        subroutine evaluate_energy(x_lo, x_hi, value)
            real(dp), intent(in) :: x_lo, x_hi
            type(gc_interval_callback_result_t), intent(out) :: value
            type(eqdsk_allowed_interval_result_t) :: evaluated
            type(eqdsk_allowed_region_cut_provenance_t) :: provenance
            integer :: local_status

            value = gc_interval_callback_result_t()
            value%query_lo = x_lo
            value%query_hi = x_hi
            value%cut_id = active_branch
            value%enclosure_certificate_id = &
                root_options%expected_enclosure_certificate_id
            value%stationary_certificate_id = &
                root_options%expected_stationary_certificate_id
            if (value%stationary_certificate_id > 0) then
                ! A finite sentinel is required by the generic callback
                ! contract.  It is replaced only by the real certificate
                ! source below; a midpoint is never promoted to evidence.
                value%stationary_point = huge(1.0_dp)
            end if
            call evaluate_eqdsk_allowed_region_cut_box(context%graph(active_branch), &
                gc_outward_interval(x_lo, x_hi), context%field_scale, &
                context%raw_psi_sep, context%profile, active_h0, active_j_k, &
                context%mass, context%charge, context%c_light, active_sigma, &
                evaluated, provenance, local_status)
            if (local_status /= EQDSK_CUT_BOX_SUCCESS) then
                value%status = local_status
                return
            end if
            value%f = to_interval(evaluated%energy_margin)
            value%df = to_interval(evaluated%denergy_margin_dR)
            value%d2f = to_interval(evaluated%d2energy_margin_dR2)
            value%status = 0
            value%certified = provenance%certified
            if (value%stationary_certificate_id > 0 .and. &
                    contains_zero(value%df) .and. &
                    associated(context%stationary_root_certificate)) then
                call context%stationary_root_certificate(x_lo, x_hi, active_h0, &
                    active_j_k, active_sigma, active_branch, &
                    value%stationary_point, value%stationary_certificate_id, &
                    local_status)
                if (local_status /= 0 .or. &
                        value%stationary_certificate_id /= &
                            root_options%expected_stationary_certificate_id .or. &
                        .not. ieee_is_finite(value%stationary_point) .or. &
                        value%stationary_point < x_lo .or. &
                        value%stationary_point > x_hi) then
                    value%stationary_point = huge(1.0_dp)
                    value%stationary_certificate_id = &
                        root_options%expected_stationary_certificate_id
                end if
            end if
        end subroutine evaluate_energy

        subroutine verify_energy(x_lo, x_hi, value, expected_id, local_status)
            real(dp), intent(in) :: x_lo, x_hi
            type(gc_interval_callback_result_t), intent(in) :: value
            integer, intent(in) :: expected_id
            integer, intent(out) :: local_status
            type(gc_interval_callback_result_t) :: expected

            call evaluate_energy(x_lo, x_hi, expected)
            local_status = 1
            if (expected_id /= root_options%expected_enclosure_certificate_id) return
            if (value%query_lo /= x_lo .or. value%query_hi /= x_hi) return
            if (value%cut_id /= active_branch .or. &
                    value%enclosure_certificate_id /= expected_id .or. &
                    value%stationary_certificate_id /= &
                        root_options%expected_stationary_certificate_id) return
            if (.not. encloses(value%f, expected%f) .or. &
                    .not. encloses(value%df, expected%df) .or. &
                    .not. encloses(value%d2f, expected%d2f)) return
            local_status = 0
        end subroutine verify_energy

        subroutine verify_energy_stationary(x_lo, x_hi, point, value, &
                expected_enclosure_id, expected_stationary_id, local_status)
            real(dp), intent(in) :: x_lo, x_hi, point
            type(gc_interval_callback_result_t), intent(out) :: value
            integer, intent(in) :: expected_enclosure_id, expected_stationary_id
            integer, intent(out) :: local_status
            type(gc_interval_callback_result_t) :: candidate

            value = gc_interval_callback_result_t()
            local_status = 1
            if (expected_enclosure_id /= &
                    root_options%expected_enclosure_certificate_id) return
            if (expected_stationary_id <= 0 .or. &
                    expected_stationary_id /= &
                    root_options%expected_stationary_certificate_id) return
            if (expected_stationary_id /= &
                    certificate_index('eqdsk_allowed_energy_stationary')) return
            if (point < x_lo .or. point > x_hi) return
            call evaluate_energy(x_lo, x_hi, candidate)
            call evaluate_energy(point, point, value)
            call verify_gc_eqdsk_fixed_invariant_stationary_certificate( &
                x_lo, x_hi, candidate, value, expected_enclosure_id, &
                expected_stationary_id, local_status)
        end subroutine verify_energy_stationary

        subroutine copy_branch_roots(ctx, energy, jk, sign, branch_id, root_result, &
                offset, root_boxes, coordinates, canonical, canonical_boxes, maps, &
                local_status)
            type(gc_eqdsk_certified_allowed_provider_context_t), intent(inout) :: ctx
            real(dp), intent(in) :: energy, jk
            integer, intent(in) :: sign, branch_id, offset
            type(gc_interval_root_result_t), intent(in) :: root_result
            type(gc_interval_root_box_t), intent(inout) :: root_boxes(:)
            real(dp), intent(inout) :: coordinates(:), canonical(:)
            type(gc_interval_t), intent(inout) :: canonical_boxes(:)
            type(gc_cylindrical_root_coordinate_map_t), intent(inout) :: maps(:)
            integer, intent(out) :: local_status

            integer :: i, global_index
            type(gc_interval_t) :: endpoint

            local_status = GC_EQDSK_ALLOWED_PROVIDER_SUCCESS
            do i = 1, root_result%nroots
                global_index = offset+i
                root_boxes(global_index) = root_result%roots(i)
                coordinates(global_index) = midpoint(gc_interval_t( &
                    root_result%roots(i)%lo, root_result%roots(i)%hi))
                call root_canonical_enclosure(ctx, energy, jk, sign, branch_id, &
                    root_result%roots(i), endpoint, local_status)
                if (local_status /= GC_EQDSK_ALLOWED_PROVIDER_SUCCESS) return
                canonical_boxes(global_index) = endpoint
                canonical(global_index) = midpoint(endpoint)
                maps(global_index)%source_root_certificate = root_result%roots(i)
                maps(global_index)%source_domain_enclosure = gc_interval_t( &
                    ctx%partition%branches(branch_id)%source_r_lo, &
                    ctx%partition%branches(branch_id)%source_r_hi)
                maps(global_index)%mapped_class_enclosure = gc_interval_t( &
                    root_result%roots(i)%lo, root_result%roots(i)%hi)
                maps(global_index)%map_certificate_id = &
                    ctx%partition%certificate_id
                maps(global_index)%monotonicity_sign = 1
                maps(global_index)%strict_monotonicity_certified = .true.
                maps(global_index)%mapping_enclosure_certified = .true.
                maps(global_index)%source_coordinate = 'physical-R'
                maps(global_index)%source_units = 'cm'
                maps(global_index)%class_coordinate = 'physical-R'
                maps(global_index)%class_units = 'cm'
            end do
        end subroutine copy_branch_roots

        subroutine build_branch_components(ctx, energy, jk, sign, branch_id, &
                root_result, offset, components, count, local_status)
            type(gc_eqdsk_certified_allowed_provider_context_t), intent(inout) :: ctx
            real(dp), intent(in) :: energy, jk
            integer, intent(in) :: sign, branch_id, offset
            type(gc_interval_root_result_t), intent(in) :: root_result
            type(gc_cylindrical_allowed_component_t), allocatable, intent(inout) :: components(:)
            integer, intent(inout) :: count
            integer, intent(out) :: local_status

            integer :: i, n, global_lower, global_upper
            real(dp) :: domain_lo, domain_hi, left, right, probe
            type(gc_interval_t) :: begin_canonical, end_canonical, measure
            logical :: lower_root, upper_root
            type(gc_cylindrical_allowed_component_t) :: component

            local_status = GC_EQDSK_ALLOWED_PROVIDER_SUCCESS
            call branch_domain(ctx, branch_id, domain_lo, domain_hi)
            n = root_result%nroots
            left = domain_lo
            do i = 0, n
                lower_root = i > 0
                upper_root = i < n
                if (lower_root) then
                    left = root_result%roots(i)%hi
                else
                    left = domain_lo
                end if
                if (upper_root) then
                    right = root_result%roots(i+1)%lo
                else
                    right = domain_hi
                end if
                if (right <= left) cycle
                probe = left + 0.5_dp*(right-left)
                if (.not. ieee_is_finite(probe)) then
                    local_status = GC_EQDSK_ALLOWED_PROVIDER_MEASURE_FAILURE
                    return
                end if
                if (.not. interval_is_allowed(ctx, energy, jk, sign, branch_id, &
                        probe, probe)) cycle
                ! A generated canonical jet cannot yet certify the one-sided
                ! limit at a turning endpoint.  Do not replace it by a sample
                ! or a fitted square-root expression.
                if (lower_root .or. upper_root) then
                    local_status = GC_EQDSK_ALLOWED_PROVIDER_MEASURE_FAILURE
                    return
                end if
                call endpoint_canonical(ctx, energy, jk, sign, branch_id, left, &
                    begin_canonical, local_status)
                if (local_status /= GC_EQDSK_ALLOWED_PROVIDER_SUCCESS) return
                call endpoint_canonical(ctx, energy, jk, sign, branch_id, right, &
                    end_canonical, local_status)
                if (local_status /= GC_EQDSK_ALLOWED_PROVIDER_SUCCESS) return
                call certify_measure(ctx, energy, jk, sign, branch_id, left, right, &
                    measure, local_status)
                if (local_status /= GC_EQDSK_ALLOWED_PROVIDER_SUCCESS) return
                component = gc_cylindrical_allowed_component_t()
                count = count+1
                component%component_id = count
                component%sigma = sign
                component%x_begin = left
                component%x_end = right
                component%canonical_begin = midpoint(begin_canonical)
                component%canonical_end = midpoint(end_canonical)
                component%canonical_measure = midpoint(measure)
                component%canonical_measure_lower = measure%lo
                component%canonical_measure_upper = measure%hi
                component%x_begin_enclosure = gc_interval_t(left, left)
                component%x_end_enclosure = gc_interval_t(right, right)
                component%canonical_begin_enclosure = begin_canonical
                component%canonical_end_enclosure = end_canonical
                component%canonical_measure_enclosure = measure
                component%canonical_measure_certificate_id = &
                    certificate_index('eqdsk_canonical_cut')
                component%endpoint_evidence_certified = .true.
                component%canonical_measure_certified = .true.
                component%lower_root = lower_root
                component%upper_root = upper_root
                if (lower_root) then
                    global_lower = offset+i
                    component%lower_root_index = global_lower
                end if
                if (upper_root) then
                    global_upper = offset+i+1
                    component%upper_root_index = global_upper
                end if
                call append_component(components, component)
            end do
        end subroutine build_branch_components

        subroutine root_canonical_enclosure(ctx, energy, jk, sign, branch_id, root, &
                value, local_status)
            type(gc_eqdsk_certified_allowed_provider_context_t), intent(inout) :: ctx
            real(dp), intent(in) :: energy, jk
            integer, intent(in) :: sign, branch_id
            type(gc_interval_root_box_t), intent(in) :: root
            type(gc_interval_t), intent(out) :: value
            integer, intent(out) :: local_status
            type(eqdsk_allowed_interval_result_t) :: evaluated
            type(eqdsk_allowed_region_cut_provenance_t) :: provenance
            integer :: status_box

            value = gc_interval_t()
            call evaluate_eqdsk_allowed_region_cut_box(ctx%graph(branch_id), &
                gc_outward_interval(root%lo, root%hi), ctx%field_scale, &
                ctx%raw_psi_sep, ctx%profile, energy, jk, ctx%mass, ctx%charge, &
                ctx%c_light, sign, evaluated, provenance, status_box)
            local_status = GC_EQDSK_ALLOWED_PROVIDER_MEASURE_FAILURE
            if (status_box /= EQDSK_CUT_BOX_SUCCESS) return
            value = to_interval(evaluated%psi_physical)
            if (.not. valid_interval(value)) return
            local_status = GC_EQDSK_ALLOWED_PROVIDER_SUCCESS
        end subroutine root_canonical_enclosure

        subroutine endpoint_canonical(ctx, energy, jk, sign, branch_id, x, value, local_status)
            type(gc_eqdsk_certified_allowed_provider_context_t), intent(inout) :: ctx
            real(dp), intent(in) :: energy, jk, x
            integer, intent(in) :: sign, branch_id
            type(gc_interval_t), intent(out) :: value
            integer, intent(out) :: local_status
            type(eqdsk_allowed_interval_result_t) :: evaluated
            type(eqdsk_allowed_region_cut_provenance_t) :: provenance
            integer :: status_box

            value = gc_interval_t()
            call evaluate_eqdsk_allowed_region_cut_box(ctx%graph(branch_id), &
                gc_outward_interval(x, x), ctx%field_scale, ctx%raw_psi_sep, &
                ctx%profile, energy, jk, ctx%mass, ctx%charge, ctx%c_light, sign, &
                evaluated, provenance, status_box)
            local_status = GC_EQDSK_ALLOWED_PROVIDER_MEASURE_FAILURE
            if (status_box /= EQDSK_CUT_BOX_SUCCESS) return
            if (.not. evaluated%canonical_chart_certified) return
            value = to_interval(evaluated%psi_star)
            if (.not. valid_interval(value)) return
            local_status = GC_EQDSK_ALLOWED_PROVIDER_SUCCESS
        end subroutine endpoint_canonical

        subroutine certify_measure(ctx, energy, jk, sign, branch_id, left, right, &
                value, local_status)
            type(gc_eqdsk_certified_allowed_provider_context_t), intent(inout) :: ctx
            real(dp), intent(in) :: energy, jk, left, right
            integer, intent(in) :: sign, branch_id
            type(gc_interval_t), intent(out) :: value
            integer, intent(out) :: local_status

            type(gc_interval_t) :: left_value, right_value
            type(eqdsk_allowed_interval_result_t) :: evaluated
            type(eqdsk_allowed_region_cut_provenance_t) :: provenance
            integer :: status_box

            value = gc_interval_t()
            call evaluate_eqdsk_allowed_region_cut_box(ctx%graph(branch_id), &
                gc_outward_interval(left, right), ctx%field_scale, &
                ctx%raw_psi_sep, ctx%profile, energy, jk, ctx%mass, ctx%charge, &
                ctx%c_light, sign, evaluated, provenance, status_box)
            local_status = GC_EQDSK_ALLOWED_PROVIDER_MEASURE_FAILURE
            if (status_box /= EQDSK_CUT_BOX_SUCCESS) return
            if (.not. evaluated%canonical_chart_certified) return
            if (.not. excludes_zero(to_interval(evaluated%dpsi_star_dR))) return
            call endpoint_canonical(ctx, energy, jk, sign, branch_id, left, &
                left_value, local_status)
            if (local_status /= GC_EQDSK_ALLOWED_PROVIDER_SUCCESS) return
            call endpoint_canonical(ctx, energy, jk, sign, branch_id, right, &
                right_value, local_status)
            if (local_status /= GC_EQDSK_ALLOWED_PROVIDER_SUCCESS) return
            value = abs_difference(right_value, left_value)
            if (.not. valid_interval(value) .or. value%hi <= 0.0_dp) then
                local_status = GC_EQDSK_ALLOWED_PROVIDER_MEASURE_FAILURE
                return
            end if
            local_status = GC_EQDSK_ALLOWED_PROVIDER_SUCCESS
        end subroutine certify_measure

        logical function interval_is_allowed(ctx, energy, jk, sign, branch_id, lo, hi)
            type(gc_eqdsk_certified_allowed_provider_context_t), intent(inout) :: ctx
            real(dp), intent(in) :: energy, jk, lo, hi
            integer, intent(in) :: sign, branch_id
            type(eqdsk_allowed_interval_result_t) :: evaluated
            type(eqdsk_allowed_region_cut_provenance_t) :: provenance
            integer :: status_box

            call evaluate_eqdsk_allowed_region_cut_box(ctx%graph(branch_id), &
                gc_outward_interval(lo, hi), ctx%field_scale, ctx%raw_psi_sep, &
                ctx%profile, energy, jk, ctx%mass, ctx%charge, ctx%c_light, sign, &
                evaluated, provenance, status_box)
            interval_is_allowed = status_box == EQDSK_CUT_BOX_SUCCESS .and. &
                evaluated%energy_margin%lo > 0.0_dp
        end function interval_is_allowed

        subroutine append_component(array, item)
            type(gc_cylindrical_allowed_component_t), allocatable, intent(inout) :: array(:)
            type(gc_cylindrical_allowed_component_t), intent(in) :: item
            type(gc_cylindrical_allowed_component_t), allocatable :: extended(:)
            integer :: n

            n = size(array)
            allocate(extended(n+1))
            if (n > 0) extended(:n) = array
            extended(n+1) = item
            call move_alloc(extended, array)
        end subroutine append_component

    end subroutine build_gc_eqdsk_certified_allowed_regions

    subroutine verify_gc_eqdsk_certified_allowed_regions(context, h0, j_k, sigma, &
            rc_min, rc_max, regions, certificate_id, status)
        type(gc_eqdsk_certified_allowed_provider_context_t), intent(inout) :: context
        real(dp), intent(in) :: h0, j_k, rc_min, rc_max
        integer, intent(in) :: sigma
        type(gc_cylindrical_allowed_region_set_t), intent(in) :: regions
        integer, intent(out) :: certificate_id, status

        type(gc_cylindrical_allowed_region_set_t) :: expected
        integer :: i

        certificate_id = 0
        status = GC_EQDSK_ALLOWED_PROVIDER_CERTIFICATE_FAILURE
        if (.not. context%initialized) return
        if (.not. all(ieee_is_finite([h0, j_k, rc_min, rc_max]))) return
        if (rc_min /= context%partition%physical_r_lo .or. &
                rc_max /= context%partition%physical_r_hi) return
        call build_gc_eqdsk_certified_allowed_regions(context, h0, j_k, sigma, &
            expected, status)
        if (status /= GC_EQDSK_ALLOWED_PROVIDER_SUCCESS) then
            certificate_id = 0
            return
        end if
        if (regions%nroots /= expected%nroots .or. &
                regions%ncomponents /= expected%ncomponents) return
        if (regions%nroots > 0) then
            if (.not. all(regions%roots == expected%roots)) return
            if (.not. all(regions%root_canonical == expected%root_canonical)) return
        end if
        do i = 1, regions%ncomponents
            if (regions%components(i)%component_id /= &
                    expected%components(i)%component_id) return
            if (regions%components(i)%sigma /= expected%components(i)%sigma) return
            if (regions%components(i)%x_begin /= expected%components(i)%x_begin .or. &
                    regions%components(i)%x_end /= expected%components(i)%x_end) return
            if (regions%components(i)%canonical_measure /= &
                    expected%components(i)%canonical_measure) return
            if (regions%components(i)%canonical_measure_certificate_id /= &
                    expected%components(i)%canonical_measure_certificate_id) return
        end do
        certificate_id = GC_EQDSK_ALLOWED_PROVIDER_CERTIFICATE_ID
        status = GC_EQDSK_ALLOWED_PROVIDER_SUCCESS
    end subroutine verify_gc_eqdsk_certified_allowed_regions

    pure function to_interval(value) result(interval)
        type(gc_outward_interval_t), intent(in) :: value
        type(gc_interval_t) :: interval

        interval = gc_interval_t(value%lo, value%hi)
    end function to_interval

    pure function add_interval(left, right) result(value)
        type(gc_interval_t), intent(in) :: left, right
        type(gc_interval_t) :: value

        value = gc_interval_t(ieee_next_after(left%lo + right%lo, -huge(1.0_dp)), &
            ieee_next_after(left%hi + right%hi, huge(1.0_dp)))
    end function add_interval

    pure function abs_difference(right, left) result(value)
        type(gc_interval_t), intent(in) :: right, left
        type(gc_interval_t) :: value
        real(dp) :: lo, hi

        lo = min(abs(right%lo-left%lo), abs(right%hi-left%hi))
        hi = max(abs(right%lo-left%hi), abs(right%hi-left%lo))
        value = gc_interval_t(ieee_next_after(max(0.0_dp, lo), -huge(1.0_dp)), &
            ieee_next_after(hi, huge(1.0_dp)))
    end function abs_difference

    pure real(dp) function midpoint(value)
        type(gc_interval_t), intent(in) :: value
        midpoint = 0.5_dp*(value%lo+value%hi)
    end function midpoint

    pure logical function valid_interval(value)
        type(gc_interval_t), intent(in) :: value
        valid_interval = ieee_is_finite(value%lo) .and. &
            ieee_is_finite(value%hi) .and. value%lo <= value%hi
    end function valid_interval

    pure logical function excludes_zero(value)
        type(gc_interval_t), intent(in) :: value
        excludes_zero = valid_interval(value) .and. &
            (value%hi < 0.0_dp .or. value%lo > 0.0_dp)
    end function excludes_zero

    pure logical function contains_zero(value)
        type(gc_interval_t), intent(in) :: value
        contains_zero = valid_interval(value) .and. value%lo <= 0.0_dp .and. &
            value%hi >= 0.0_dp
    end function contains_zero

    pure logical function exact_zero(value)
        type(gc_interval_t), intent(in) :: value
        exact_zero = valid_interval(value) .and. value%lo == 0.0_dp .and. &
            value%hi == 0.0_dp
    end function exact_zero

    pure logical function encloses(value, expected)
        type(gc_interval_t), intent(in) :: value, expected
        encloses = valid_interval(value) .and. valid_interval(expected) .and. &
            value%lo <= expected%lo .and. value%hi >= expected%hi
    end function encloses

end module neort_gc_eqdsk_certified_allowed_provider
