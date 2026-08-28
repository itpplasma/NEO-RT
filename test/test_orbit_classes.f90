program test_orbit_classes
    ! Independent geometry oracle for periodic well detection.  The analytic
    ! field has two minima and two maxima on one period; the phase-shifted
    ! case checks that the periodic endpoint is not treated as a boundary.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_orbit_classes, only: find_periodic_well_barriers
    use util, only: pi
    implicit none

    integer, parameter :: n = 384
    real(dp) :: theta(n), bmod(n)
    logical :: minima(n), maxima(n)
    real(dp) :: barriers(n)
    integer :: k, nmin, nmax
    logical :: valid

    do k = 1, n
        theta(k) = -pi + 2.0_dp*pi*real(k - 1, dp)/real(n, dp)
        bmod(k) = 3.0_dp + 0.7_dp*cos(3.0_dp*theta(k)) + 0.18_dp*cos(theta(k))
    end do
    call check_inventory(theta, bmod, minima, maxima, barriers, nmin, nmax, valid)

    do k = 1, n
        theta(k) = -pi + (real(k - 1, dp) + 0.37_dp)*2.0_dp*pi/real(n, dp)
        bmod(k) = 3.0_dp + 0.7_dp*cos(3.0_dp*theta(k)) + 0.18_dp*cos(theta(k))
    end do
    call check_inventory(theta, bmod, minima, maxima, barriers, nmin, nmax, valid)

    print *, "test_orbit_classes PASSED"

contains

    subroutine check_inventory(theta, bmod, minima, maxima, barriers, nmin, nmax, valid)
        real(dp), intent(in) :: theta(:), bmod(:)
        logical, intent(out) :: minima(:), maxima(:)
        real(dp), intent(out) :: barriers(:)
        integer, intent(out) :: nmin, nmax
        logical, intent(out) :: valid
        integer :: k, left, right, left_max, right_max, step

        call find_periodic_well_barriers(bmod, minima, maxima, barriers, nmin, nmax, valid)
        if (.not. valid .or. nmin /= 3 .or. nmax /= 3) error stop "periodic well inventory is wrong"
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
            if (.not. minima(k)) cycle
            left_max = 0
            right_max = 0
            do step = 1, size(bmod)
                left = modulo(k - step - 1, size(bmod)) + 1
                if (maxima(left)) then
                    left_max = left
                    exit
                end if
            end do
            do step = 1, size(bmod)
                right = modulo(k + step - 1, size(bmod)) + 1
                if (maxima(right)) then
                    right_max = right
                    exit
                end if
            end do
            if (left_max == 0 .or. right_max == 0) error stop "barrier oracle found no neighbors"
            if (barriers(k) /= min(bmod(left_max), bmod(right_max))) then
                error stop "periodic well barrier is wrong"
            end if
        end do
    end subroutine check_inventory

end program test_orbit_classes
