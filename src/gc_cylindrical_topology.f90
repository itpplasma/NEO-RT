module neort_gc_cylindrical_topology
    !! Allowed-region and canonical-measure data for a real-space cut.
    !!
    !! The scan parameter is deliberately not called eta.  It may be R, a
    !! section coordinate, or any caller-owned class parameter.  Every
    !! positive component of v_parallel^2 is retained separately with its
    !! sigma label; no global topology state is used.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_gc_cylindrical_model, only: GC_CYL_INVALID_INPUT, &
        GC_CYL_SUCCESS, gc_cylindrical_allowed_component_t, &
        gc_cylindrical_allowed_value_i
    use neort_gc_certified_interval_roots, only: gc_interval_root_box_t, &
        gc_interval_t

    implicit none
    private

    type, public :: gc_cylindrical_root_coordinate_map_t
        !! A root is isolated in the provider's source coordinate before it
        !! is mapped into the class coordinate consumed by the adapter.  Keep
        !! both certificates: a mapped interval is not itself a root proof.
        type(gc_interval_root_box_t) :: source_root_certificate
        type(gc_interval_t) :: source_domain_enclosure
        type(gc_interval_t) :: mapped_class_enclosure
        integer :: map_certificate_id = 0
        integer :: monotonicity_sign = 0
        logical :: strict_monotonicity_certified = .false.
        logical :: mapping_enclosure_certified = .false.
        character(len=64) :: source_coordinate = ''
        character(len=64) :: source_units = ''
        character(len=64) :: class_coordinate = ''
        character(len=64) :: class_units = ''
    end type gc_cylindrical_root_coordinate_map_t

    type, public :: gc_cylindrical_allowed_region_set_t
        integer :: nroots = 0
        integer :: ncomponents = 0
        real(dp), allocatable :: roots(:)
        real(dp), allocatable :: root_canonical(:)
        type(gc_interval_root_box_t), allocatable :: root_boxes(:)
        type(gc_cylindrical_root_coordinate_map_t), allocatable :: &
            root_coordinate_maps(:)
        type(gc_interval_t), allocatable :: root_canonical_enclosures(:)
        type(gc_cylindrical_allowed_component_t), allocatable :: components(:)
        real(dp) :: total_canonical_measure = 0.0_dp
        type(gc_interval_t) :: total_canonical_measure_enclosure
        !! The adaptive sampler below is diagnostic only.  It is not an
        !! interval/root-isolation proof and must never be consumed as one.
        logical :: topology_certified = .false.
        character(len=64) :: certificate_method = ''
    end type gc_cylindrical_allowed_region_set_t

    public :: find_gc_cylindrical_allowed_regions
    public :: canonical_measure_density
    public :: canonical_measure_from_curve

    integer, parameter :: maximum_adaptive_samples = 32768
    integer, parameter :: maximum_adaptive_depth = 12

contains

    subroutine find_gc_cylindrical_allowed_regions(evaluate, x_min, x_max, &
            scan_points, sigma, regions, status)
        procedure(gc_cylindrical_allowed_value_i) :: evaluate
        real(dp), intent(in) :: x_min, x_max
        integer, intent(in) :: scan_points, sigma
        type(gc_cylindrical_allowed_region_set_t), intent(out) :: regions
        integer, intent(out) :: status

        real(dp), allocatable :: x(:), value(:), derivative(:), roots(:), root_psi(:)
        integer :: nroots, ncomponents, root_capacity, n_samples

        regions = gc_cylindrical_allowed_region_set_t()
        status = GC_CYL_INVALID_INPUT
        if (x_max <= x_min .or. scan_points < 2) return
        if (abs(sigma) /= 1) return
        root_capacity = max(8, 4*scan_points)
        allocate(x(maximum_adaptive_samples), value(maximum_adaptive_samples), &
            derivative(maximum_adaptive_samples), &
            roots(root_capacity), root_psi(root_capacity))
        call adaptive_sample_allowed_function(evaluate, x_min, x_max, scan_points, &
            x, value, derivative, n_samples, status)
        if (status /= GC_CYL_SUCCESS) return
        call collect_all_roots(evaluate, x(1:n_samples), value(1:n_samples), &
            derivative(1:n_samples), roots, root_psi, &
            nroots, status)
        if (status /= GC_CYL_SUCCESS) return
        call sort_and_deduplicate_roots(roots, root_psi, nroots)
        call build_components(evaluate, x_min, x_max, roots, root_psi, nroots, &
            sigma, regions%components, ncomponents, status)
        if (status /= GC_CYL_SUCCESS) return
        regions%nroots = nroots
        regions%ncomponents = ncomponents
        call copy_root_data(roots, root_psi, nroots, regions)
        call sum_component_measure(regions)
        ! A finite sample set can miss an even root, a stationary canonical
        ! momentum X point, or a separatrix between samples.  Keep the scan
        ! available as a diagnostic decomposition, but do not promote it to
        ! a physical topology certificate.  A future Fortsym interval/root
        ! isolation provider must set this flag and its method explicitly.
        regions%topology_certified = .false.
        regions%certificate_method = 'unresolved-finite-scan'
    end subroutine find_gc_cylindrical_allowed_regions

    subroutine adaptive_sample_allowed_function(evaluate, x_min, x_max, scan_points, &
            x, value, derivative, n_samples, status)
        procedure(gc_cylindrical_allowed_value_i) :: evaluate
        real(dp), intent(in) :: x_min, x_max
        integer, intent(in) :: scan_points
        real(dp), intent(out) :: x(:), value(:), derivative(:)
        integer, intent(out) :: n_samples, status

        real(dp) :: step, left, right, f_left, f_right, d_left, d_right
        real(dp) :: psi_star, dpsi_star_dx
        integer :: i, n

        status = GC_CYL_SUCCESS
        n_samples = 0
        n = max(2, scan_points)
        step = (x_max - x_min)/real(n - 1, dp)
        call evaluate(x_min, f_left, d_left, psi_star, dpsi_star_dx, status)
        if (status /= GC_CYL_SUCCESS) return
        call append_sample(x_min, f_left, d_left, x, value, derivative, &
            n_samples, status)
        if (status /= GC_CYL_SUCCESS) return
        do i = 1, n - 1
            left = x_min + real(i - 1, dp)*step
            right = x_min + real(i, dp)*step
            call evaluate(right, f_right, d_right, psi_star, dpsi_star_dx, status)
            if (status /= GC_CYL_SUCCESS) return
            call refine_allowed_interval(evaluate, left, f_left, d_left, right, &
                f_right, d_right, 0, x, value, derivative, n_samples, status)
            if (status /= GC_CYL_SUCCESS) return
            f_left = f_right
            d_left = d_right
        end do
    end subroutine adaptive_sample_allowed_function

    recursive subroutine refine_allowed_interval(evaluate, left, f_left, d_left, &
            right, f_right, d_right, depth, x, value, derivative, n_samples, status)
        procedure(gc_cylindrical_allowed_value_i) :: evaluate
        real(dp), intent(in) :: left, f_left, d_left, right, f_right, d_right
        integer, intent(in) :: depth
        real(dp), intent(inout) :: x(:), value(:), derivative(:)
        integer, intent(inout) :: n_samples, status

        real(dp) :: midpoint, f_mid, d_mid, psi_star, dpsi_star_dx
        logical :: refine

        if (status /= GC_CYL_SUCCESS) return
        midpoint = 0.5_dp*(left + right)
        call evaluate(midpoint, f_mid, d_mid, psi_star, dpsi_star_dx, status)
        if (status /= GC_CYL_SUCCESS) return
        if (.not. all(ieee_is_finite([left, right, f_left, f_mid, f_right, &
                d_left, d_mid, d_right, psi_star, dpsi_star_dx]))) then
            status = GC_CYL_INVALID_INPUT
            return
        end if
        refine = allowed_interval_needs_refinement(f_left, f_mid, f_right, &
            d_left, d_mid, d_right)
        if (refine .and. depth < maximum_adaptive_depth .and. &
                n_samples + 4 < size(x)) then
            call refine_allowed_interval(evaluate, left, f_left, d_left, midpoint, &
                f_mid, d_mid, depth + 1, x, value, derivative, n_samples, status)
            call refine_allowed_interval(evaluate, midpoint, f_mid, d_mid, right, &
                f_right, d_right, depth + 1, x, value, derivative, n_samples, status)
        else
            call append_sample(midpoint, f_mid, d_mid, x, value, derivative, &
                n_samples, status)
            call append_sample(right, f_right, d_right, x, value, derivative, &
                n_samples, status)
        end if
    end subroutine refine_allowed_interval

    logical function allowed_interval_needs_refinement(f_left, f_mid, f_right, &
            d_left, d_mid, d_right)
        real(dp), intent(in) :: f_left, f_mid, f_right, d_left, d_mid, d_right
        real(dp) :: curvature, scale, derivative_scale

        allowed_interval_needs_refinement = .false.
        if (f_left*f_mid <= 0.0_dp .or. f_mid*f_right <= 0.0_dp) then
            allowed_interval_needs_refinement = .true.
            return
        end if
        if (d_left*d_mid <= 0.0_dp .or. d_mid*d_right <= 0.0_dp) then
            allowed_interval_needs_refinement = .true.
            return
        end if
        scale = max(max(abs(f_left), abs(f_mid)), abs(f_right))
        curvature = abs(f_mid - 0.5_dp*(f_left + f_right))
        derivative_scale = max(max(abs(d_left), abs(d_mid)), abs(d_right))
        if (curvature > 0.125_dp*max(scale, tiny(scale)) .and. &
                derivative_scale > tiny(derivative_scale)) then
            allowed_interval_needs_refinement = .true.
        end if
    end function allowed_interval_needs_refinement

    subroutine append_sample(x_value, f_value, derivative_value, x, value, &
            derivative, n_samples, status)
        real(dp), intent(in) :: x_value, f_value, derivative_value
        real(dp), intent(inout) :: x(:), value(:), derivative(:)
        integer, intent(inout) :: n_samples, status

        status = GC_CYL_SUCCESS
        if (n_samples > 0) then
            if (x_value <= x(n_samples)) return
        end if
        if (n_samples >= size(x)) then
            status = GC_CYL_INVALID_INPUT
            return
        end if
        n_samples = n_samples + 1
        x(n_samples) = x_value
        value(n_samples) = f_value
        derivative(n_samples) = derivative_value
    end subroutine append_sample

    subroutine collect_all_roots(evaluate, x, value, derivative, roots, root_psi, &
            nroots, status)
        procedure(gc_cylindrical_allowed_value_i) :: evaluate
        real(dp), intent(in) :: x(:), value(:), derivative(:)
        real(dp), intent(out) :: roots(:), root_psi(:)
        integer, intent(out) :: nroots, status

        ! The value test must be tighter than the coordinate spacing created
        ! by adaptive refinement.  A looser tolerance turns a tangential
        ! root into a cloud of nearby roots before the sorted deduplication
        ! pass can see the stationary-point identity.
        real(dp), parameter :: root_tolerance = 1.0e-13_dp
        real(dp), parameter :: stationary_tolerance = 1.0e-12_dp
        real(dp) :: root, psi_star, root_value, scale, derivative_scale
        integer :: i, local_status

        nroots = 0
        status = GC_CYL_SUCCESS
        scale = max(1.0_dp, maxval(abs(value)))
        derivative_scale = max(1.0_dp, maxval(abs(derivative)))
        do i = 1, size(x) - 1
            ! Interior near-zero samples are not roots: on a tangency they
            ! form a whole cloud around the stationary root.  Exact sampled
            ! zeros are retained; all other simple roots come from a sign
            ! bracket and even roots from the derivative bracket below.
            if (value(i) == 0.0_dp) then
                call append_root(evaluate, x(i), roots, root_psi, nroots, status)
                if (status /= GC_CYL_SUCCESS) return
            end if
            if (value(i)*value(i + 1) < 0.0_dp) then
                call bisect_root(evaluate, x(i), x(i + 1), value(i), &
                    value(i + 1), root, psi_star, local_status)
                if (local_status /= GC_CYL_SUCCESS) then
                    status = local_status
                    return
                end if
                call append_root_value(root, psi_star, roots, root_psi, nroots, status)
                if (status /= GC_CYL_SUCCESS) return
            end if
            if (derivative(i)*derivative(i + 1) < 0.0_dp) then
                call bisect_stationary(evaluate, x(i), x(i + 1), derivative(i), &
                    derivative(i + 1), root, root_value, psi_star, local_status, &
                    stationary_tolerance*derivative_scale)
                if (local_status /= GC_CYL_SUCCESS) then
                    status = local_status
                    return
                end if
                if (abs(root_value) <= root_tolerance*scale) then
                    call append_root_value(root, psi_star, roots, root_psi, &
                        nroots, status)
                    if (status /= GC_CYL_SUCCESS) return
                end if
            end if
        end do
        if (value(size(value)) == 0.0_dp) then
            call append_root(evaluate, x(size(x)), roots, root_psi, nroots, status)
        end if
    end subroutine collect_all_roots

    subroutine bisect_root(evaluate, left, right, f_left, f_right, root, &
            psi_star, status)
        procedure(gc_cylindrical_allowed_value_i) :: evaluate
        real(dp), intent(in) :: left, right, f_left, f_right
        real(dp), intent(out) :: root, psi_star
        integer, intent(out) :: status

        real(dp), parameter :: coordinate_tolerance = 1.0e-12_dp
        real(dp) :: a, b, fa, fb, midpoint, fm, derivative, dpsi_star_dx
        integer :: iteration

        a = left
        b = right
        fa = f_left
        fb = f_right
        psi_star = 0.0_dp
        status = GC_CYL_SUCCESS
        do iteration = 1, 80
            midpoint = 0.5_dp*(a + b)
            call evaluate(midpoint, fm, derivative, psi_star, dpsi_star_dx, status)
            if (status /= GC_CYL_SUCCESS) return
            if (abs(fm) <= coordinate_tolerance .or. &
                abs(b - a) <= coordinate_tolerance*max(1.0_dp, abs(midpoint))) then
                root = midpoint
                return
            end if
            if (fa*fm <= 0.0_dp) then
                b = midpoint
                fb = fm
            else
                a = midpoint
                fa = fm
            end if
        end do
        root = 0.5_dp*(a + b)
        call evaluate(root, fm, derivative, psi_star, dpsi_star_dx, status)
        associate (unused_fb => fb)
        end associate
    end subroutine bisect_root

    subroutine append_root(evaluate, root, roots, root_psi, nroots, status)
        procedure(gc_cylindrical_allowed_value_i) :: evaluate
        real(dp), intent(in) :: root
        real(dp), intent(inout) :: roots(:), root_psi(:)
        integer, intent(inout) :: nroots
        integer, intent(out) :: status

        real(dp) :: value, derivative, psi_star, dpsi_star_dx

        call evaluate(root, value, derivative, psi_star, dpsi_star_dx, status)
        if (status /= GC_CYL_SUCCESS) return
        call append_root_value(root, psi_star, roots, root_psi, nroots, status)
    end subroutine append_root

    subroutine append_root_value(root, psi_star, roots, root_psi, nroots, status)
        real(dp), intent(in) :: root, psi_star
        real(dp), intent(inout) :: roots(:), root_psi(:)
        integer, intent(inout) :: nroots
        integer, intent(out) :: status

        real(dp), parameter :: duplicate_tolerance = 1.0e-9_dp

        status = GC_CYL_SUCCESS
        if (nroots > 0) then
            if (abs(root - roots(nroots)) <= duplicate_tolerance) return
        end if
        if (nroots >= size(roots)) then
            status = GC_CYL_INVALID_INPUT
            return
        end if
        nroots = nroots + 1
        roots(nroots) = root
        root_psi(nroots) = psi_star
    end subroutine append_root_value

    subroutine sort_and_deduplicate_roots(roots, root_psi, nroots)
        real(dp), intent(inout) :: roots(:), root_psi(:)
        integer, intent(inout) :: nroots

        real(dp), parameter :: duplicate_tolerance = 1.0e-8_dp
        real(dp) :: root, psi
        integer :: i, j, first

        do i = 1, nroots - 1
            first = i
            do j = i + 1, nroots
                if (roots(j) < roots(first)) first = j
            end do
            if (first /= i) then
                root = roots(i)
                roots(i) = roots(first)
                roots(first) = root
                psi = root_psi(i)
                root_psi(i) = root_psi(first)
                root_psi(first) = psi
            end if
        end do

        if (nroots < 2) return
        j = 1
        do i = 2, nroots
            if (abs(roots(i) - roots(j)) <= duplicate_tolerance) cycle
            j = j + 1
            roots(j) = roots(i)
            root_psi(j) = root_psi(i)
        end do
        nroots = j
    end subroutine sort_and_deduplicate_roots

    subroutine bisect_stationary(evaluate, left, right, dleft, dright, root, &
            value, psi_star, status, derivative_tolerance)
        procedure(gc_cylindrical_allowed_value_i) :: evaluate
        real(dp), intent(in) :: left, right, dleft, dright, derivative_tolerance
        real(dp), intent(out) :: root, value, psi_star
        integer, intent(out) :: status

        real(dp), parameter :: coordinate_tolerance = 1.0e-12_dp
        real(dp) :: a, b, fa, fb, midpoint, fm, derivative, dpsi_star_dx
        integer :: iteration

        a = left
        b = right
        fa = dleft
        fb = dright
        status = GC_CYL_SUCCESS
        value = 0.0_dp
        psi_star = 0.0_dp
        do iteration = 1, 80
            midpoint = 0.5_dp*(a + b)
            call evaluate(midpoint, value, derivative, psi_star, dpsi_star_dx, &
                status)
            if (status /= GC_CYL_SUCCESS) return
            if (abs(derivative) <= derivative_tolerance .or. &
                abs(b - a) <= coordinate_tolerance*max(1.0_dp, abs(midpoint))) then
                root = midpoint
                return
            end if
            if (fa*derivative <= 0.0_dp) then
                b = midpoint
                fb = derivative
            else
                a = midpoint
                fa = derivative
            end if
        end do
        root = 0.5_dp*(a + b)
        call evaluate(root, value, derivative, psi_star, dpsi_star_dx, status)
        associate (unused_fb => fb)
        end associate
    end subroutine bisect_stationary

    subroutine build_components(evaluate, x_min, x_max, roots, root_psi, nroots, &
            sigma, components, ncomponents, status)
        procedure(gc_cylindrical_allowed_value_i) :: evaluate
        real(dp), intent(in) :: x_min, x_max, roots(:), root_psi(:)
        integer, intent(in) :: nroots, sigma
        type(gc_cylindrical_allowed_component_t), allocatable, intent(out) :: components(:)
        integer, intent(out) :: ncomponents, status

        real(dp), allocatable :: edges(:)
        integer :: i, count
        real(dp) :: midpoint, value, derivative, psi_star, dpsi_star_dx

        allocate(edges(nroots + 2))
        edges(1) = x_min
        if (nroots > 0) edges(2:nroots + 1) = roots(1:nroots)
        edges(nroots + 2) = x_max
        count = 0
        do i = 1, size(edges) - 1
            midpoint = 0.5_dp*(edges(i) + edges(i + 1))
            call evaluate(midpoint, value, derivative, psi_star, dpsi_star_dx, status)
            if (status /= GC_CYL_SUCCESS) return
            if (value > 0.0_dp) count = count + 1
        end do
        ncomponents = count
        allocate(components(ncomponents))
        count = 0
        do i = 1, size(edges) - 1
            midpoint = 0.5_dp*(edges(i) + edges(i + 1))
            call evaluate(midpoint, value, derivative, psi_star, dpsi_star_dx, status)
            if (status /= GC_CYL_SUCCESS) return
            if (value <= 0.0_dp) cycle
            count = count + 1
            call fill_component(evaluate, edges(i), edges(i + 1), i - 1, &
                nroots, root_psi, sigma, count, components(count), status)
            if (status /= GC_CYL_SUCCESS) return
        end do
    end subroutine build_components

    subroutine fill_component(evaluate, x_begin, x_end, edge_index, nroots, &
            root_psi, sigma, component_id, component, status)
        procedure(gc_cylindrical_allowed_value_i) :: evaluate
        real(dp), intent(in) :: x_begin, x_end, root_psi(:)
        integer, intent(in) :: edge_index, nroots, sigma, component_id
        type(gc_cylindrical_allowed_component_t), intent(out) :: component
        integer, intent(out) :: status

        real(dp) :: value, derivative, psi_begin, psi_end, dpsi_begin, dpsi_end

        call evaluate(x_begin, value, derivative, psi_begin, dpsi_begin, status)
        if (status /= GC_CYL_SUCCESS) return
        call evaluate(x_end, value, derivative, psi_end, dpsi_end, status)
        if (status /= GC_CYL_SUCCESS) return
        component%component_id = component_id
        component%sigma = sigma
        component%x_begin = x_begin
        component%x_end = x_end
        component%canonical_begin = psi_begin
        component%canonical_end = psi_end
        component%lower_root = edge_index > 0
        component%upper_root = edge_index < nroots + 1
        if (component%lower_root) component%lower_root_index = edge_index
        if (component%upper_root) component%upper_root_index = edge_index + 1
        call integrate_canonical_measure(evaluate, x_begin, x_end, &
            component%canonical_measure, status)
        if (status /= GC_CYL_SUCCESS) return
        component%canonical_measure_lower = component%canonical_measure
        component%canonical_measure_upper = component%canonical_measure
    end subroutine fill_component

    subroutine integrate_canonical_measure(evaluate, x_begin, x_end, measure, &
            status)
        procedure(gc_cylindrical_allowed_value_i) :: evaluate
        real(dp), intent(in) :: x_begin, x_end
        real(dp), intent(out) :: measure
        integer, intent(out) :: status

        integer, parameter :: quadrature_points = 64
        real(dp) :: previous_x, current_x, previous_derivative, current_derivative
        real(dp) :: value, derivative, psi_star, dpsi_star_dx, step
        integer :: i, local_status

        measure = 0.0_dp
        status = GC_CYL_SUCCESS
        step = (x_end - x_begin)/real(quadrature_points, dp)
        previous_x = x_begin
        call evaluate(previous_x, value, derivative, psi_star, &
            previous_derivative, local_status)
        if (local_status /= GC_CYL_SUCCESS) then
            status = local_status
            return
        end if
        do i = 1, quadrature_points
            current_x = x_begin + real(i, dp)*step
            call evaluate(current_x, value, derivative, psi_star, &
                dpsi_star_dx, local_status)
            if (local_status /= GC_CYL_SUCCESS) then
                status = local_status
                return
            end if
            measure = measure &
                +0.5_dp*(abs(previous_derivative) + abs(dpsi_star_dx))*step
            previous_x = current_x
            previous_derivative = dpsi_star_dx
        end do
    end subroutine integrate_canonical_measure

    subroutine copy_root_data(roots, root_psi, nroots, regions)
        real(dp), intent(in) :: roots(:), root_psi(:)
        integer, intent(in) :: nroots
        type(gc_cylindrical_allowed_region_set_t), intent(inout) :: regions

        allocate(regions%roots(nroots), regions%root_canonical(nroots))
        if (nroots == 0) return
        regions%roots = roots(1:nroots)
        regions%root_canonical = root_psi(1:nroots)
    end subroutine copy_root_data

    subroutine sum_component_measure(regions)
        type(gc_cylindrical_allowed_region_set_t), intent(inout) :: regions

        regions%total_canonical_measure = 0.0_dp
        regions%total_canonical_measure_enclosure = gc_interval_t()
        if (regions%ncomponents == 0) return
        regions%total_canonical_measure = sum(regions%components%canonical_measure)
        regions%total_canonical_measure_enclosure = gc_interval_t( &
            regions%total_canonical_measure, regions%total_canonical_measure)
    end subroutine sum_component_measure

    pure function canonical_measure_density(dpsi_star_dx) result(value)
        real(dp), intent(in) :: dpsi_star_dx
        real(dp) :: value

        value = abs(dpsi_star_dx)
    end function canonical_measure_density

    subroutine canonical_measure_from_curve(x, psi_star, measure, status)
        real(dp), intent(in) :: x(:), psi_star(:)
        real(dp), intent(out) :: measure
        integer, intent(out) :: status

        integer :: i

        measure = 0.0_dp
        status = GC_CYL_INVALID_INPUT
        if (size(x) /= size(psi_star) .or. size(x) < 2) return
        if (.not. all(ieee_is_finite([x, psi_star]))) return
        do i = 1, size(x) - 1
            if (x(i + 1) <= x(i)) return
            measure = measure + abs(psi_star(i + 1) - psi_star(i))
        end do
        status = GC_CYL_SUCCESS
    end subroutine canonical_measure_from_curve

end module neort_gc_cylindrical_topology
