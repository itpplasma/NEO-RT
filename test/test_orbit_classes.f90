program test_orbit_classes
    ! Independent geometry oracle for periodic well detection.  The analytic
    ! field has two minima and two maxima on one period; the phase-shifted
    ! case checks that the periodic endpoint is not treated as a boundary.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_orbit_classes, only: find_periodic_extrema
    use util, only: pi
    implicit none

    integer, parameter :: n = 256
    real(dp) :: theta(n), bmod(n)
    logical :: minima(n), maxima(n)
    integer :: k, nmin, nmax

    do k = 1, n
        theta(k) = -pi + 2.0_dp*pi*real(k - 1, dp)/real(n, dp)
        bmod(k) = 3.0_dp + 0.7_dp*cos(2.0_dp*theta(k))
    end do
    call check_inventory(theta, bmod, minima, maxima, nmin, nmax)

    do k = 1, n
        theta(k) = -pi + (real(k - 1, dp) + 0.37_dp)*2.0_dp*pi/real(n, dp)
        bmod(k) = 3.0_dp + 0.7_dp*cos(2.0_dp*theta(k))
    end do
    call check_inventory(theta, bmod, minima, maxima, nmin, nmax)

    print *, "test_orbit_classes PASSED"

contains

    subroutine check_inventory(theta, bmod, minima, maxima, nmin, nmax)
        real(dp), intent(in) :: theta(:), bmod(:)
        logical, intent(out) :: minima(:), maxima(:)
        integer, intent(out) :: nmin, nmax
        integer :: k, left, right

        call find_periodic_extrema(bmod, minima, maxima, nmin, nmax)
        if (nmin /= 2 .or. nmax /= 2) error stop "periodic well count is wrong"
        do k = 1, size(bmod)
            if (minima(k)) then
                left = modulo(k - 2, size(bmod)) + 1
                right = modulo(k, size(bmod)) + 1
                if (bmod(k) > bmod(left) .or. bmod(k) > bmod(right)) then
                    error stop "classified minimum is not a local minimum"
                end if
            end if
            if (maxima(k)) then
                left = modulo(k - 2, size(bmod)) + 1
                right = modulo(k, size(bmod)) + 1
                if (bmod(k) < bmod(left) .or. bmod(k) < bmod(right)) then
                    error stop "classified maximum is not a local maximum"
                end if
            end if
        end do
    end subroutine check_inventory

end program test_orbit_classes
