program test_orbit_classes
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_orbit_classes, only: require_single_well
    use util, only: pi
    implicit none
    integer, parameter :: n = 384
    real(dp) :: theta(n), bmod(n)
    integer :: k, shifted
    character(32) :: mode

    call get_command_argument(1, mode)
    do shifted = 0, 1
        do k = 1, n
            theta(k) = -pi + (real(k - 1, dp) + 0.37_dp*shifted)*2*pi/n
            if (mode == 'composite') then
                ! Unequal maxima: energy between the barriers connects wells.
                bmod(k) = 3 + 0.7_dp*cos(3*theta(k)) + 0.18_dp*cos(theta(k))
            else
                bmod(k) = 3 + 0.7_dp*cos(theta(k))
            end if
        end do
        call require_single_well(bmod)
    end do
    print *, "test_orbit_classes PASSED"
end program test_orbit_classes
