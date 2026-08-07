module neort_gc_full_resonance
    !! Status-aware resonance search for finite-width canonical orbits.
    !! Invalid orbit samples split the scan into disjoint valid segments, so
    !! roots can never be manufactured across a lost/no-return interval.
    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none
    private

    integer, parameter, public :: GC_RESONANCE_SUCCESS = 0
    integer, parameter, public :: GC_RESONANCE_INVALID_INPUT = 1
    integer, parameter, public :: GC_RESONANCE_NOT_CONVERGED = 2

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

        real(dp) :: eta_left, eta_right, residual_left, residual_right
        real(dp) :: deta
        integer :: k, left_status, right_status, root_status

        roots = 0.0_dp
        derivatives = 0.0_dp
        nroots = 0
        status = GC_RESONANCE_INVALID_INPUT
        if (eta_max <= eta_min .or. scan_intervals < 1) return
        if (residual_tolerance <= 0.0_dp .or. eta_tolerance <= 0.0_dp) return
        if (size(derivatives) /= size(roots)) return

        deta = (eta_max - eta_min)/real(scan_intervals, dp)
        eta_left = eta_min
        call evaluate(eta_left, residual_left, left_status)
        do k = 1, scan_intervals
            eta_right = eta_min + real(k, dp)*deta
            call evaluate(eta_right, residual_right, right_status)
            if (left_status == 0 .and. right_status == 0) then
                if (residual_left == 0.0_dp .or. residual_right == 0.0_dp .or. &
                    sign(1.0_dp, residual_left) /= sign(1.0_dp, residual_right)) then
                    if (nroots == size(roots)) return
                    call bisect_valid_segment(evaluate, eta_left, eta_right, &
                        residual_left, residual_right, residual_tolerance, &
                        eta_tolerance, roots(nroots + 1), &
                        derivatives(nroots + 1), root_status)
                    if (root_status /= GC_RESONANCE_SUCCESS) then
                        status = root_status
                        return
                    end if
                    if (nroots == 0 .or. &
                        abs(roots(nroots + 1) - roots(nroots)) > eta_tolerance) then
                        nroots = nroots + 1
                    end if
                end if
            end if
            eta_left = eta_right
            residual_left = residual_right
            left_status = right_status
        end do
        status = GC_RESONANCE_SUCCESS
    end subroutine find_gc_resonances

    subroutine bisect_valid_segment(evaluate, eta_a, eta_b, residual_a, residual_b, &
            residual_tolerance, eta_tolerance, root, derivative, status)
        procedure(gc_residual_i) :: evaluate
        real(dp), intent(in) :: eta_a, eta_b, residual_a, residual_b
        real(dp), intent(in) :: residual_tolerance, eta_tolerance
        real(dp), intent(out) :: root, derivative
        integer, intent(out) :: status

        integer, parameter :: max_iterations = 100
        real(dp) :: left, right, middle, fleft, fright, fmiddle
        real(dp) :: eta_minus, eta_plus, fminus, fplus, h
        integer :: k, middle_status, minus_status, plus_status

        root = 0.0_dp
        derivative = 0.0_dp
        status = GC_RESONANCE_NOT_CONVERGED
        left = eta_a
        right = eta_b
        fleft = residual_a
        fright = residual_b
        if (fleft == 0.0_dp) then
            middle = left
            fmiddle = fleft
            k = 1
        else if (fright == 0.0_dp) then
            middle = right
            fmiddle = fright
            k = 1
        else
            do k = 1, max_iterations
                middle = 0.5_dp*(left + right)
                call evaluate(middle, fmiddle, middle_status)
                ! A topology loss inside an apparently valid coarse bracket makes
                ! the bracket unusable; never bridge it by selecting one side.
                if (middle_status /= 0) return
                if (abs(fmiddle) <= residual_tolerance .or. &
                    right - left <= eta_tolerance) exit
                if (sign(1.0_dp, fmiddle) == sign(1.0_dp, fleft)) then
                    left = middle
                    fleft = fmiddle
                else
                    right = middle
                    fright = fmiddle
                end if
            end do
        end if
        if (k > max_iterations) return

        root = middle
        h = max(eta_tolerance, sqrt(epsilon(root))*max(1.0_dp, abs(root)))
        eta_minus = max(eta_a, root - h)
        eta_plus = min(eta_b, root + h)
        if (eta_plus <= eta_minus) return
        call evaluate(eta_minus, fminus, minus_status)
        call evaluate(eta_plus, fplus, plus_status)
        if (minus_status /= 0 .or. plus_status /= 0) return
        derivative = (fplus - fminus)/(eta_plus - eta_minus)
        status = GC_RESONANCE_SUCCESS
    end subroutine bisect_valid_segment

end module neort_gc_full_resonance
