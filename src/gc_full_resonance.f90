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

    integer, parameter :: scan_refinement_factor = 8
    integer, parameter :: maximum_scan_points = 4096
    integer, parameter :: maximum_refinement_iterations = 100

    abstract interface
        subroutine gc_residual_i(eta, residual, status)
            import dp
            real(dp), intent(in) :: eta
            real(dp), intent(out) :: residual
            integer, intent(out) :: status
        end subroutine gc_residual_i
    end interface

    public :: find_gc_resonances, gc_residual_i

contains

    subroutine find_gc_resonances(evaluate, eta_min, eta_max, scan_intervals, &
            residual_tolerance, eta_tolerance, roots, derivatives, nroots, status)
        procedure(gc_residual_i) :: evaluate
        real(dp), intent(in) :: eta_min, eta_max
        integer, intent(in) :: scan_intervals
        real(dp), intent(in) :: residual_tolerance, eta_tolerance
        real(dp), intent(out) :: roots(:), derivatives(:)
        integer, intent(out) :: nroots, status

        real(dp), allocatable :: eta(:), residual(:)
        integer, allocatable :: sample_status(:)
        integer :: nscan, k, local_status, first_valid, last_valid
        logical :: partial_scan, failed, capacity_exhausted

        roots = 0.0_dp
        derivatives = 0.0_dp
        nroots = 0
        status = GC_RESONANCE_INVALID_INPUT
        capacity_exhausted = .false.
        if (eta_max <= eta_min .or. scan_intervals < 1) return
        if (residual_tolerance <= 0.0_dp .or. eta_tolerance <= 0.0_dp) return
        if (size(derivatives) /= size(roots)) return

        nscan = min(maximum_scan_points, scan_intervals*scan_refinement_factor)
        allocate(eta(0:nscan), residual(0:nscan), sample_status(0:nscan))
        partial_scan = .false.
        do k = 0, nscan
            if (k == 0) then
                eta(k) = eta_min
            else if (k == nscan) then
                eta(k) = eta_max
            else
                eta(k) = eta_min + real(k, dp)*(eta_max - eta_min)/real(nscan, dp)
            end if
            call evaluate(eta(k), residual(k), sample_status(k))
        end do

        ! An open orbit interval may be invalid exactly at its endpoints
        ! (most notably at the trapped/passing separatrix).  Those samples do
        ! not describe a missing interior interval: discard the invalid
        ! prefix/suffix and search only the valid open interval.  Any other
        ! invalid sample remains evidence that the scan could have hidden a
        ! root.
        first_valid = -1
        do k = 0, nscan
            if (sample_is_valid(sample_status(k), residual(k))) then
                first_valid = k
                exit
            end if
        end do
        last_valid = -1
        do k = nscan, 0, -1
            if (sample_is_valid(sample_status(k), residual(k))) then
                last_valid = k
                exit
            end if
        end do
        if (first_valid < 0 .or. last_valid < first_valid &
                .or. first_valid == last_valid) then
            status = GC_RESONANCE_PARTIAL
            return
        end if

        partial_scan = .false.
        do k = 0, first_valid - 1
            partial_scan = partial_scan .or. &
                sample_status(k) /= GC_RESONANCE_BOUNDARY_INVALID
        end do
        do k = first_valid, last_valid
            partial_scan = partial_scan .or. &
                .not. sample_is_valid(sample_status(k), residual(k))
        end do
        do k = last_valid + 1, nscan
            partial_scan = partial_scan .or. &
                sample_status(k) /= GC_RESONANCE_BOUNDARY_INVALID
        end do

        do k = first_valid, last_valid
            if (sample_is_valid(sample_status(k), residual(k)) &
                    .and. abs(residual(k)) <= residual_tolerance) then
                call append_root(eta(k), evaluate, eta(first_valid), eta(last_valid), &
                    eta_tolerance, roots, derivatives, nroots, &
                    capacity_exhausted, local_status)
                partial_scan = partial_scan .or. local_status /= GC_RESONANCE_SUCCESS
            end if
            if (k == first_valid) cycle
            if (sample_is_valid(sample_status(k - 1), residual(k - 1)) &
                    .and. sample_is_valid(sample_status(k), residual(k))) then
                if (opposite_signs(residual(k - 1), residual(k))) then
                    call bisect_valid_segment(evaluate, eta(k - 1), eta(k), &
                        residual(k - 1), residual(k), residual_tolerance, &
                        eta_tolerance, roots, derivatives, nroots, &
                        capacity_exhausted, local_status)
                    partial_scan = partial_scan .or. local_status /= GC_RESONANCE_SUCCESS
                end if
            end if
            if (k == last_valid) cycle
            if (sample_is_valid(sample_status(k - 1), residual(k - 1)) &
                    .and. sample_is_valid(sample_status(k), residual(k)) &
                    .and. sample_is_valid(sample_status(k + 1), residual(k + 1)) &
                    .and. abs(residual(k)) <= abs(residual(k - 1)) &
                    .and. abs(residual(k)) <= abs(residual(k + 1)) &
                    .and. (residual(k) /= residual(k - 1) &
                    .or. residual(k) /= residual(k + 1))) then
                call refine_tangent_root(evaluate, eta(k - 1), eta(k + 1), &
                    residual_tolerance, eta_tolerance, eta(k), residual(k), &
                    roots, derivatives, nroots, capacity_exhausted, failed)
                partial_scan = partial_scan .or. failed
            end if
        end do

        if (capacity_exhausted) partial_scan = .true.
        status = merge(GC_RESONANCE_PARTIAL, GC_RESONANCE_SUCCESS, partial_scan)
    end subroutine find_gc_resonances

    logical function sample_is_valid(sample_status, sample_residual)
        integer, intent(in) :: sample_status
        real(dp), intent(in) :: sample_residual

        sample_is_valid = sample_status == GC_RESONANCE_SUCCESS &
            .and. ieee_is_finite(sample_residual)
    end function sample_is_valid

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
            if (middle_status /= 0 .or. .not. ieee_is_finite(fmiddle)) return
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
            if (first_status /= 0 .or. second_status /= 0 &
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
        if (minus_status /= 0 .or. plus_status /= 0 &
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

        if (nroots > 0 .and. any(abs(roots(1:nroots) - root) <= eta_tolerance)) return
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
