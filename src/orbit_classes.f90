module neort_orbit_classes
    ! Only a single periodic magnetic well is supported by the asymptotic model.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: require_single_well, find_periodic_extrema
contains
    subroutine require_single_well(bmod)
        real(dp), intent(in) :: bmod(:)
        logical :: minima(size(bmod)), maxima(size(bmod))
        integer :: nmin, nmax

        if (.not. all(ieee_is_finite(bmod))) error stop "nonfinite magnetic field"
        if (any(bmod <= 0.0_dp)) error stop "nonpositive magnetic field"
        call find_periodic_extrema(bmod, minima, maxima, nmin, nmax)
        if (nmin /= 1 .or. nmax /= 1) then
            error stop "single-well separatrix model: unsupported multi-well or degenerate geometry"
        end if
    end subroutine require_single_well

    pure subroutine find_periodic_extrema(bmod, minima, maxima, nmin, nmax)
        real(dp), intent(in) :: bmod(:)
        logical, intent(out) :: minima(:), maxima(:)
        integer, intent(out) :: nmin, nmax
        integer :: k, left, right, n

        n = size(bmod)
        minima = .false.
        maxima = .false.
        nmin = 0
        nmax = 0
        if (size(minima) /= n .or. size(maxima) /= n .or. n < 3) return
        do k = 1, n
            left = periodic_index(k - 1, n)
            right = periodic_index(k + 1, n)
            if (bmod(k) <= bmod(left) .and. bmod(k) < bmod(right)) then
                nmin = nmin + 1
                minima(k) = .true.
            end if
            if (bmod(k) >= bmod(left) .and. bmod(k) > bmod(right)) then
                nmax = nmax + 1
                maxima(k) = .true.
            end if
        end do
    end subroutine find_periodic_extrema

    pure integer function periodic_index(index, n)
        integer, intent(in) :: index, n
        periodic_index = modulo(index - 1, n) + 1
    end function periodic_index

end module neort_orbit_classes
