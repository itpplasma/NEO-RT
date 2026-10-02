program test_bmod_pert_precision
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none

    integer, parameter :: nr = 13, nz = 11
    real(dp), parameter :: rmin = 100.0_dp, rmax = 200.0_dp
    real(dp), parameter :: zmin = -40.0_dp, zmax = 40.0_dp
    real(dp), parameter :: relative_tolerance = 1.0e-12_dp
    real(dp) :: r, z, relative_error, max_error
    complex(dp) :: actual, expected
    integer :: i, j

    external :: bmod_pert

    call write_grid()
    max_error = 0.0_dp
    do j = 0, 20
        z = zmin + (zmax - zmin)*real(j, dp)/20.0_dp
        if (j > 0 .and. j < 20) z = z + 0.23_dp
        do i = 0, 24
            r = rmin + (rmax - rmin)*real(i, dp)/24.0_dp
            if (i > 0 .and. i < 24) r = r + 0.37_dp
            call bmod_pert(r, z, actual)
            expected = analytic_field(r, z)
            relative_error = abs(actual - expected)/abs(expected)
            max_error = max(max_error, relative_error)
        end do
    end do
    print '(a, es12.4)', 'maximum complex relative error = ', max_error
    if (max_error > relative_tolerance) error stop 'complex precision lost'

contains

    pure function analytic_field(r, z) result(value)
        real(dp), intent(in) :: r, z
        complex(dp) :: value
        real(dp) :: re, im

        ! Low-degree polynomials have no spline truncation error, so this
        ! independently checks both components at nodes and between nodes.
        re = 0.123456789012345_dp + 0.003456789012345_dp*r &
             - 0.004567890123456_dp*z + 0.000012345678901_dp*r*z
        im = -1.414213562373095_dp + 0.000023456789012_dp*r**2 &
             + 0.000034567890123_dp*z**2
        value = cmplx(re, im, kind=dp)
    end function analytic_field

    subroutine write_grid()
        real(dp) :: rad(nr), zet(nz), bre(nr, nz), bim(nr, nz)
        complex(dp) :: value
        integer :: k, l, unit

        do k = 1, nr
            rad(k) = rmin + (rmax - rmin)*real(k - 1, dp)/real(nr - 1, dp)
        end do
        do l = 1, nz
            zet(l) = zmin + (zmax - zmin)*real(l - 1, dp)/real(nz - 1, dp)
        end do
        do l = 1, nz
            do k = 1, nr
                value = analytic_field(rad(k), zet(l))
                bre(k, l) = real(value, dp)
                bim(k, l) = aimag(value)
            end do
        end do
        open (newunit=unit, file='bmod_n.dat', form='unformatted', status='replace')
        write (unit) nr, nz
        write (unit) rad, zet
        write (unit) bre, bim
        close (unit)
    end subroutine write_grid

end program test_bmod_pert_precision
