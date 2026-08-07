module neort_gc_full_resonance
    !! Status-aware resonance search for finite-width canonical orbits.
    !! Invalid samples split the scan into disjoint valid segments.  The scan
    !! is oversampled so that even roots and nearby sign changes do not depend
    !! on a coarse endpoint sign change.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

    implicit none
    private

    integer, parameter, public :: GC_RESONANCE_SUCCESS = 0
    integer, parameter, public :: GC_RESONANCE_INVALID_INPUT = 1
    integer, parameter, public :: GC_RESONANCE_NOT_CONVERGED = 2
    integer, parameter, public :: GC_RESONANCE_PARTIAL = 3
    integer, parameter, public :: GC_RESONANCE_BOUNDARY_INVALID = 4

    integer, parameter, public :: GC_RESONANCE_SAMPLE_VALID = 0
    integer, parameter, public :: GC_RESONANCE_SAMPLE_BOUNDARY = 1
    integer, parameter, public :: GC_RESONANCE_SAMPLE_UNCONFINED = 2
    !! Only this class is certified physical when the caller supplies a real
    !! wall-polygon crossing status.  The other non-valid classes remain
    !! fail-closed and do not contribute to physical coverage.
    integer, parameter, public :: GC_RESONANCE_SAMPLE_WALL = 3
    !! Computational-domain and timeout results are unresolved, not losses.
    integer, parameter, public :: GC_RESONANCE_SAMPLE_RADIAL_DOMAIN = 4
    integer, parameter, public :: GC_RESONANCE_SAMPLE_INVALID = 5

    integer, parameter :: scan_refinement_factor = 8
    integer, parameter :: maximum_scan_points = 4096
    integer, parameter :: maximum_refinement_iterations = 100

    type, public :: gc_resonance_diagnostics_t
        integer :: scan_samples = 0
        integer :: confined_samples = 0
        integer :: boundary_samples = 0
        integer :: lost_orbits = 0
        integer :: unconfined_samples = 0
        integer :: wall_orbits = 0
        integer :: radial_domain_orbits = 0
        integer :: numerical_samples = 0
        integer :: measure_failures = 0
        integer :: unknown_measure_cells = 0
        integer :: component_count = 0
        real(dp) :: canonical_scan_measure = 0.0_dp
        real(dp) :: canonical_confined_measure = 0.0_dp
        real(dp) :: canonical_physical_measure = 0.0_dp
        real(dp) :: canonical_boundary_measure = 0.0_dp
        real(dp) :: canonical_unresolved_measure = 0.0_dp
        real(dp) :: unknown_measure_coordinate_span = 0.0_dp
        real(dp) :: confined_coverage_fraction = 0.0_dp
        real(dp) :: physical_coverage_fraction = 0.0_dp
        real(dp) :: unresolved_coverage_fraction = 0.0_dp
        logical :: canonical_measure_certified = .false.
        logical :: component_identity_certified = .false.
    end type gc_resonance_diagnostics_t

    abstract interface
        subroutine gc_residual_i(eta, residual, status)
            import dp
            real(dp), intent(in) :: eta
            real(dp), intent(out) :: residual
            integer, intent(out) :: status
        end subroutine gc_residual_i

        integer function gc_residual_class_i(sample_status)
            integer, intent(in) :: sample_status
        end function gc_residual_class_i

        subroutine gc_residual_measure_i(eta, density, status)
            import dp
            real(dp), intent(in) :: eta
            real(dp), intent(out) :: density
            integer, intent(out) :: status
        end subroutine gc_residual_measure_i

        integer function gc_residual_component_i(eta, sample_status)
            import dp
            real(dp), intent(in) :: eta
            integer, intent(in) :: sample_status
        end function gc_residual_component_i
    end interface

    public :: find_gc_resonances, gc_residual_i, gc_residual_class_i, &
        gc_residual_measure_i, gc_residual_component_i

contains

    subroutine find_gc_resonances(evaluate, eta_min, eta_max, scan_intervals, &
            residual_tolerance, eta_tolerance, roots, derivatives, nroots, status, &
            diagnostics, classify_status, measure, component_id)
        procedure(gc_residual_i) :: evaluate
        real(dp), intent(in) :: eta_min, eta_max
        integer, intent(in) :: scan_intervals
        real(dp), intent(in) :: residual_tolerance, eta_tolerance
        real(dp), intent(out) :: roots(:), derivatives(:)
        integer, intent(out) :: nroots, status
        type(gc_resonance_diagnostics_t), intent(out), optional :: diagnostics
        procedure(gc_residual_class_i), optional :: classify_status
        procedure(gc_residual_measure_i), optional :: measure
        procedure(gc_residual_component_i), optional :: component_id

        real(dp), allocatable :: eta(:), residual(:)
        real(dp), allocatable :: sample_measure(:)
        integer, allocatable :: sample_status(:), sample_class(:), sample_component(:)
        integer :: nscan, k, component_start, active_component
        integer :: measure_status
        real(dp) :: measure_density
        logical :: partial_scan, capacity_exhausted, sample_valid, measure_known

        roots = 0.0_dp
        derivatives = 0.0_dp
        nroots = 0
        status = GC_RESONANCE_INVALID_INPUT
        capacity_exhausted = .false.
        if (present(diagnostics)) diagnostics = gc_resonance_diagnostics_t()
        if (present(diagnostics)) then
            diagnostics%canonical_measure_certified = present(measure)
            diagnostics%component_identity_certified = present(component_id)
        end if
        if (eta_max <= eta_min .or. scan_intervals < 1) return
        if (residual_tolerance <= 0.0_dp .or. eta_tolerance <= 0.0_dp) return
        if (size(derivatives) /= size(roots)) return

        nscan = min(maximum_scan_points, scan_intervals*scan_refinement_factor)
        allocate(eta(0:nscan), residual(0:nscan), sample_measure(0:nscan), &
            sample_status(0:nscan), sample_class(0:nscan), &
            sample_component(0:nscan))
        partial_scan = .false.
        do k = 0, nscan
            if (k == 0) then
                eta(k) = eta_min
            else if (k == nscan) then
                eta(k) = eta_max
            else
                eta(k) = eta_min + real(k, dp)*(eta_max - eta_min)/real(nscan, dp)
            end if
        end do
        do k = 0, nscan
            call evaluate(eta(k), residual(k), sample_status(k))
            sample_class(k) = resolve_sample_class(sample_status(k), classify_status)
            sample_component(k) = 0
            if (present(component_id)) then
                sample_component(k) = component_id(eta(k), sample_status(k))
            end if
            sample_valid = sample_class(k) == GC_RESONANCE_SAMPLE_VALID
            if (sample_valid) then
                if (.not. ieee_is_finite(residual(k))) then
                    sample_class(k) = GC_RESONANCE_SAMPLE_INVALID
                    sample_valid = .false.
                end if
            end if
            sample_measure(k) = sample_cell_width(eta, k)
            measure_known = .false.
            if (present(measure)) then
                call measure(eta(k), measure_density, measure_status)
                if (measure_status /= GC_RESONANCE_SUCCESS) then
                    partial_scan = .true.
                    if (present(diagnostics)) diagnostics%measure_failures = &
                        diagnostics%measure_failures + 1
                else if (.not. ieee_is_finite(measure_density) &
                        .or. measure_density < 0.0_dp) then
                    partial_scan = .true.
                    if (present(diagnostics)) diagnostics%measure_failures = &
                        diagnostics%measure_failures + 1
                else
                    sample_measure(k) = measure_density*sample_measure(k)
                    measure_known = .true.
                end if
            end if
            call record_sample_class(diagnostics, sample_class(k), &
                sample_measure(k), measure_known)
        end do

        component_start = -1
        active_component = 0
        do k = 0, nscan
            sample_valid = sample_class(k) == GC_RESONANCE_SAMPLE_VALID
            if (sample_valid) then
                if (component_start < 0) then
                    component_start = k
                    active_component = sample_component(k)
                    if (present(diagnostics)) diagnostics%component_count = &
                        diagnostics%component_count + 1
                else if (sample_component(k) /= active_component) then
                    call process_valid_component(evaluate, eta, residual, &
                        component_start, k - 1, residual_tolerance, eta_tolerance, &
                        roots, derivatives, nroots, capacity_exhausted, partial_scan)
                    component_start = k
                    active_component = sample_component(k)
                    if (present(diagnostics)) diagnostics%component_count = &
                        diagnostics%component_count + 1
                end if
                cycle
            end if
            if (component_start >= 0) then
                call process_valid_component(evaluate, eta, residual, component_start, &
                    k - 1, residual_tolerance, eta_tolerance, roots, derivatives, &
                    nroots, capacity_exhausted, partial_scan)
                component_start = -1
            end if
            if (sample_class(k) == GC_RESONANCE_SAMPLE_BOUNDARY) then
                if (k > 0) then
                    if (k < nscan) partial_scan = .true.
                end if
            else if (sample_class(k) /= GC_RESONANCE_SAMPLE_WALL) then
                partial_scan = .true.
            end if
        end do
        if (component_start >= 0) then
            call process_valid_component(evaluate, eta, residual, component_start, &
                nscan, residual_tolerance, eta_tolerance, roots, derivatives, &
                nroots, capacity_exhausted, partial_scan)
        end if

        if (capacity_exhausted) partial_scan = .true.
        call update_coverage(diagnostics)
        status = merge(GC_RESONANCE_PARTIAL, GC_RESONANCE_SUCCESS, partial_scan)
    end subroutine find_gc_resonances

    integer function resolve_sample_class(sample_status, classify_status)
        integer, intent(in) :: sample_status
        procedure(gc_residual_class_i), optional :: classify_status

        resolve_sample_class = GC_RESONANCE_SAMPLE_INVALID
        if (present(classify_status)) then
            resolve_sample_class = classify_status(sample_status)
        else if (sample_status == GC_RESONANCE_SUCCESS) then
            resolve_sample_class = GC_RESONANCE_SAMPLE_VALID
        else if (sample_status == GC_RESONANCE_BOUNDARY_INVALID) then
            resolve_sample_class = GC_RESONANCE_SAMPLE_BOUNDARY
        end if
        if (resolve_sample_class < GC_RESONANCE_SAMPLE_VALID .or. &
            resolve_sample_class > GC_RESONANCE_SAMPLE_INVALID) then
            resolve_sample_class = GC_RESONANCE_SAMPLE_INVALID
        end if
    end function resolve_sample_class

    real(dp) function sample_cell_width(eta, index)
        real(dp), intent(in) :: eta(0:)
        integer, intent(in) :: index

        integer :: last_index

        last_index = ubound(eta, 1)
        if (index == 0) then
            sample_cell_width = 0.5_dp*(eta(1) - eta(0))
        else if (index == last_index) then
            sample_cell_width = 0.5_dp*(eta(last_index) - eta(last_index - 1))
        else
            sample_cell_width = 0.5_dp*(eta(index + 1) - eta(index - 1))
        end if
    end function sample_cell_width

    subroutine record_sample_class(diagnostics, sample_class, sample_measure, &
            measure_known)
        type(gc_resonance_diagnostics_t), intent(inout), optional :: diagnostics
        integer, intent(in) :: sample_class
        real(dp), intent(in) :: sample_measure
        logical, intent(in) :: measure_known

        if (.not. present(diagnostics)) return
        diagnostics%scan_samples = diagnostics%scan_samples + 1
        select case (sample_class)
        case (GC_RESONANCE_SAMPLE_VALID)
            diagnostics%confined_samples = diagnostics%confined_samples + 1
        case (GC_RESONANCE_SAMPLE_BOUNDARY)
            diagnostics%boundary_samples = diagnostics%boundary_samples + 1
        case (GC_RESONANCE_SAMPLE_WALL)
            diagnostics%lost_orbits = diagnostics%lost_orbits + 1
            diagnostics%wall_orbits = diagnostics%wall_orbits + 1
        case (GC_RESONANCE_SAMPLE_UNCONFINED)
            diagnostics%unconfined_samples = diagnostics%unconfined_samples + 1
        case (GC_RESONANCE_SAMPLE_RADIAL_DOMAIN)
            diagnostics%radial_domain_orbits = diagnostics%radial_domain_orbits + 1
        case default
            diagnostics%numerical_samples = diagnostics%numerical_samples + 1
        end select
        if (.not. measure_known) then
            diagnostics%canonical_measure_certified = .false.
            diagnostics%unknown_measure_cells = diagnostics%unknown_measure_cells + 1
            diagnostics%unknown_measure_coordinate_span = &
                diagnostics%unknown_measure_coordinate_span + sample_measure
            return
        end if
        diagnostics%canonical_scan_measure = diagnostics%canonical_scan_measure &
            +sample_measure
        select case (sample_class)
        case (GC_RESONANCE_SAMPLE_VALID)
            diagnostics%canonical_confined_measure = &
                diagnostics%canonical_confined_measure + sample_measure
        case (GC_RESONANCE_SAMPLE_BOUNDARY)
            diagnostics%canonical_boundary_measure = &
                diagnostics%canonical_boundary_measure + sample_measure
        case (GC_RESONANCE_SAMPLE_WALL)
            diagnostics%canonical_physical_measure = &
                diagnostics%canonical_physical_measure + sample_measure
        case (GC_RESONANCE_SAMPLE_UNCONFINED, GC_RESONANCE_SAMPLE_RADIAL_DOMAIN)
            diagnostics%canonical_unresolved_measure = &
                diagnostics%canonical_unresolved_measure + sample_measure
        case default
            diagnostics%canonical_unresolved_measure = &
                diagnostics%canonical_unresolved_measure + sample_measure
        end select
    end subroutine record_sample_class

    subroutine update_coverage(diagnostics)
        type(gc_resonance_diagnostics_t), intent(inout), optional :: diagnostics

        real(dp) :: covered_measure

        if (.not. present(diagnostics)) return
        if (.not. diagnostics%canonical_measure_certified) then
            diagnostics%confined_coverage_fraction = 0.0_dp
            diagnostics%physical_coverage_fraction = 0.0_dp
            diagnostics%unresolved_coverage_fraction = 0.0_dp
            return
        end if
        covered_measure = diagnostics%canonical_scan_measure &
            -diagnostics%canonical_boundary_measure
        if (covered_measure <= 0.0_dp) then
            diagnostics%confined_coverage_fraction = 0.0_dp
            diagnostics%physical_coverage_fraction = 0.0_dp
            diagnostics%unresolved_coverage_fraction = 0.0_dp
            return
        end if
        diagnostics%confined_coverage_fraction = diagnostics%canonical_confined_measure &
            /covered_measure
        diagnostics%physical_coverage_fraction = ( &
            diagnostics%canonical_confined_measure &
            +diagnostics%canonical_physical_measure)/covered_measure
        diagnostics%canonical_unresolved_measure = max(0.0_dp, covered_measure &
            -diagnostics%canonical_confined_measure &
            -diagnostics%canonical_physical_measure)
        diagnostics%unresolved_coverage_fraction = &
            diagnostics%canonical_unresolved_measure/covered_measure
    end subroutine update_coverage

    subroutine process_valid_component(evaluate, eta, residual, first, last, &
            residual_tolerance, eta_tolerance, roots, derivatives, nroots, &
            capacity_exhausted, partial_scan)
        procedure(gc_residual_i) :: evaluate
        real(dp), intent(in) :: eta(0:), residual(0:)
        integer, intent(in) :: first, last
        real(dp), intent(in) :: residual_tolerance, eta_tolerance
        real(dp), intent(inout) :: roots(:), derivatives(:)
        integer, intent(inout) :: nroots
        logical, intent(inout) :: capacity_exhausted, partial_scan

        integer :: k, local_status
        logical :: failed

        if (last - first < 1) then
            partial_scan = .true.
            return
        end if
        do k = first, last
            if (abs(residual(k)) <= residual_tolerance) then
                call append_root(eta(k), evaluate, eta(first), eta(last), eta_tolerance, &
                    roots, derivatives, nroots, capacity_exhausted, local_status)
                if (local_status /= GC_RESONANCE_SUCCESS) partial_scan = .true.
            end if
            if (k > first) then
                if (opposite_signs(residual(k - 1), residual(k))) then
                    call bisect_valid_segment(evaluate, eta(k - 1), eta(k), &
                        residual(k - 1), residual(k), residual_tolerance, eta_tolerance, &
                        roots, derivatives, nroots, capacity_exhausted, local_status)
                    if (local_status /= GC_RESONANCE_SUCCESS) partial_scan = .true.
                end if
            end if
            if (k > first) then
                if (k < last) then
                    if (abs(residual(k)) <= abs(residual(k - 1))) then
                        if (abs(residual(k)) <= abs(residual(k + 1))) then
                            if (residual(k) /= residual(k - 1) .or. &
                                residual(k) /= residual(k + 1)) then
                                call refine_tangent_root(evaluate, eta(k - 1), eta(k + 1), &
                                    residual_tolerance, eta_tolerance, eta(k), residual(k), &
                                    roots, derivatives, nroots, capacity_exhausted, failed)
                                if (failed) partial_scan = .true.
                            end if
                        end if
                    end if
                end if
            end if
        end do
    end subroutine process_valid_component

    logical function opposite_signs(left_value, right_value)
        real(dp), intent(in) :: left_value, right_value

        opposite_signs = (left_value < 0.0_dp .and. right_value > 0.0_dp) &
            .or. (left_value > 0.0_dp .and. right_value < 0.0_dp)
    end function opposite_signs

    subroutine bisect_valid_segment(evaluate, eta_a, eta_b, residual_a, residual_b, &
            residual_tolerance, eta_tolerance, roots, derivatives, nroots, &
            capacity_exhausted, status)
        procedure(gc_residual_i) :: evaluate
        real(dp), intent(in) :: eta_a, eta_b, residual_a, residual_b
        real(dp), intent(in) :: residual_tolerance, eta_tolerance
        real(dp), intent(inout) :: roots(:), derivatives(:)
        integer, intent(inout) :: nroots
        logical, intent(inout) :: capacity_exhausted
        integer, intent(out) :: status

        real(dp) :: left, right, middle, fleft, fright, fmiddle, derivative
        integer :: k, middle_status, derivative_status
        logical :: converged

        status = GC_RESONANCE_NOT_CONVERGED
        left = eta_a
        right = eta_b
        fleft = residual_a
        fright = residual_b
        converged = .false.
        do k = 1, maximum_refinement_iterations
            if (abs(fleft) <= residual_tolerance) then
                middle = left
                converged = .true.
                exit
            end if
            if (abs(fright) <= residual_tolerance) then
                middle = right
                converged = .true.
                exit
            end if
            middle = 0.5_dp*(left + right)
            call evaluate(middle, fmiddle, middle_status)
            if (middle_status /= GC_RESONANCE_SUCCESS .or. &
                .not. ieee_is_finite(fmiddle)) return
            if (abs(fmiddle) <= residual_tolerance .or. &
                right - left <= eta_tolerance) then
                converged = .true.
                exit
            end if
            if (opposite_signs(fleft, fmiddle)) then
                right = middle
                fright = fmiddle
            else
                left = middle
                fleft = fmiddle
            end if
        end do
        if (.not. converged) return
        call estimate_derivative(evaluate, middle, eta_a, eta_b, eta_tolerance, &
            derivative, derivative_status)
        if (derivative_status /= GC_RESONANCE_SUCCESS) return
        call append_root_values(middle, derivative, eta_tolerance, roots, derivatives, &
            nroots, capacity_exhausted)
        status = GC_RESONANCE_SUCCESS
    end subroutine bisect_valid_segment

    subroutine refine_tangent_root(evaluate, eta_a, eta_b, residual_tolerance, &
            eta_tolerance, eta_initial, residual_initial, roots, derivatives, &
            nroots, capacity_exhausted, failed)
        procedure(gc_residual_i) :: evaluate
        real(dp), intent(in) :: eta_a, eta_b, residual_tolerance, eta_tolerance
        real(dp), intent(in) :: eta_initial, residual_initial
        real(dp), intent(inout) :: roots(:), derivatives(:)
        integer, intent(inout) :: nroots
        logical, intent(inout) :: capacity_exhausted
        logical, intent(out) :: failed

        real(dp), parameter :: golden_ratio = 0.6180339887498948482_dp
        real(dp) :: left, right, first, second, ffirst, fsecond
        real(dp) :: best_eta, best_residual, derivative
        integer :: k, first_status, second_status, derivative_status

        failed = .false.
        left = eta_a
        right = eta_b
        best_eta = eta_initial
        best_residual = abs(residual_initial)
        do k = 1, maximum_refinement_iterations
            if (right - left <= eta_tolerance) exit
            first = right - golden_ratio*(right - left)
            second = left + golden_ratio*(right - left)
            call evaluate(first, ffirst, first_status)
            call evaluate(second, fsecond, second_status)
            if (first_status /= GC_RESONANCE_SUCCESS .or. &
                second_status /= GC_RESONANCE_SUCCESS &
                .or. .not. ieee_is_finite(ffirst) &
                .or. .not. ieee_is_finite(fsecond)) then
                failed = .true.
                return
            end if
            if (abs(ffirst) < best_residual) then
                best_eta = first
                best_residual = abs(ffirst)
            end if
            if (abs(fsecond) < best_residual) then
                best_eta = second
                best_residual = abs(fsecond)
            end if
            if (abs(ffirst) <= abs(fsecond)) then
                right = second
            else
                left = first
            end if
        end do
        ! A resolved local minimum above tolerance is a valid no-root result,
        ! not an incomplete scan.  `failed` is reserved for invalid samples
        ! or refinement/derivative evaluations that prevent a decision.
        if (best_residual > residual_tolerance) return
        call estimate_derivative(evaluate, best_eta, eta_a, eta_b, eta_tolerance, &
            derivative, derivative_status)
        if (derivative_status /= GC_RESONANCE_SUCCESS) then
            failed = .true.
            return
        end if
        call append_root_values(best_eta, derivative, eta_tolerance, roots, derivatives, &
            nroots, capacity_exhausted)
    end subroutine refine_tangent_root

    subroutine estimate_derivative(evaluate, root, eta_min, eta_max, eta_tolerance, &
            derivative, status)
        procedure(gc_residual_i) :: evaluate
        real(dp), intent(in) :: root, eta_min, eta_max, eta_tolerance
        real(dp), intent(out) :: derivative
        integer, intent(out) :: status

        real(dp) :: h, eta_minus, eta_plus, fminus, fplus
        integer :: minus_status, plus_status

        derivative = 0.0_dp
        status = GC_RESONANCE_NOT_CONVERGED
        h = max(eta_tolerance, sqrt(epsilon(root))*max(1.0_dp, abs(root)))
        eta_minus = max(eta_min, root - h)
        eta_plus = min(eta_max, root + h)
        if (eta_plus <= eta_minus) return
        call evaluate(eta_minus, fminus, minus_status)
        call evaluate(eta_plus, fplus, plus_status)
        if (minus_status /= GC_RESONANCE_SUCCESS .or. &
            plus_status /= GC_RESONANCE_SUCCESS &
            .or. .not. ieee_is_finite(fminus) &
            .or. .not. ieee_is_finite(fplus)) return
        derivative = (fplus - fminus)/(eta_plus - eta_minus)
        status = GC_RESONANCE_SUCCESS
    end subroutine estimate_derivative

    subroutine append_root(root, evaluate, eta_min, eta_max, &
            eta_tolerance, roots, derivatives, nroots, capacity_exhausted, status)
        procedure(gc_residual_i) :: evaluate
        real(dp), intent(in) :: root, eta_min, eta_max, eta_tolerance
        real(dp), intent(inout) :: roots(:), derivatives(:)
        integer, intent(inout) :: nroots
        logical, intent(inout) :: capacity_exhausted
        integer, intent(out) :: status

        real(dp) :: derivative
        integer :: derivative_status

        call estimate_derivative(evaluate, root, eta_min, eta_max, eta_tolerance, &
            derivative, derivative_status)
        status = derivative_status
        if (status /= GC_RESONANCE_SUCCESS) return
        call append_root_values(root, derivative, eta_tolerance, roots, derivatives, &
            nroots, capacity_exhausted)
        status = GC_RESONANCE_SUCCESS
    end subroutine append_root

    subroutine append_root_values(root, derivative, eta_tolerance, roots, derivatives, &
            nroots, capacity_exhausted)
        real(dp), intent(in) :: root, derivative, eta_tolerance
        real(dp), intent(inout) :: roots(:), derivatives(:)
        integer, intent(inout) :: nroots
        logical, intent(inout) :: capacity_exhausted

        integer :: k, insertion_point

        if (nroots > 0) then
            if (any(abs(roots(1:nroots) - root) <= eta_tolerance)) return
        end if
        if (nroots == size(roots)) then
            capacity_exhausted = .true.
            return
        end if
        insertion_point = nroots + 1
        do k = 1, nroots
            if (root < roots(k)) then
                insertion_point = k
                exit
            end if
        end do
        do k = nroots, insertion_point, -1
            roots(k + 1) = roots(k)
            derivatives(k + 1) = derivatives(k)
        end do
        roots(insertion_point) = root
        derivatives(insertion_point) = derivative
        nroots = nroots + 1
    end subroutine append_root_values

end module neort_gc_full_resonance
