program test_bmod_pert_bounds
! bmod_pert must reproduce a tabulated analytic field inside the bmod_n.dat
! R-Z grid and must stop with an error outside it instead of clamping.
  use, intrinsic :: iso_fortran_env, only : dp => real64
  implicit none

  integer, parameter :: nrad = 101, nzet = 81
  real(dp), parameter :: rmin = 100.0_dp, rmax = 200.0_dp
  real(dp), parameter :: zmin = -40.0_dp, zmax = 40.0_dp
  ! cmplx() in bmod_pert rounds to default real, so ~1e-7 relative is the floor.
  real(dp), parameter :: tol = 1.0e-5_dp
  character(len=32) :: mode
  complex(dp) :: bmod_n
  real(dp) :: r, z, err, errmax
  integer :: i, j

  external :: bmod_pert

  if (command_argument_count() /= 1) error stop 'usage: test_bmod_pert_bounds.x mode'
  call get_command_argument(1, mode)
  call write_grid()

  select case (trim(mode))
  case ('inside')
    errmax = 0.0_dp
    do j = 0, 40
      z = zmin + (zmax - zmin)*real(j, dp)/40.0_dp
      do i = 0, 50
        r = rmin + (rmax - rmin)*real(i, dp)/50.0_dp
        if (i > 0 .and. i < 50) r = r + 0.37_dp
        call bmod_pert(r, z, bmod_n)
        err = abs(bmod_n - cmplx(f_re(r, z), f_im(r, z), dp))
        errmax = max(errmax, err)
      end do
    end do
    print '(a, es10.3)', 'max |bmod_pert - analytic| = ', errmax
    if (errmax > tol) error stop 'interpolation error above tolerance'
  case ('outside_r')
    call probe(rmax + 0.5_dp, 0.0_dp)
  case ('outside_z')
    call probe(150.0_dp, zmin - 0.5_dp)
  case default
    error stop 'unknown mode'
  end select

contains

  subroutine probe(r0, z0)
    real(dp), intent(in) :: r0, z0
    complex(dp) :: b

    call bmod_pert(r0, z0, b)
    print '(a, 2es24.16)', 'bmod_pert returned outside the grid: ', b
    error stop 1
  end subroutine probe

  pure function f_re(r, z) result(f)
    real(dp), intent(in) :: r, z
    real(dp) :: f

    f = sin(0.05_dp*r)*cos(0.07_dp*z) + 0.01_dp*r
  end function f_re

  pure function f_im(r, z) result(f)
    real(dp), intent(in) :: r, z
    real(dp) :: f

    f = exp(-((r - 150.0_dp)**2 + z**2)/2000.0_dp)
  end function f_im

  subroutine write_grid()
    real(dp) :: rad(nrad), zet(nzet), bre(nrad, nzet), bim(nrad, nzet)
    integer :: k, l, u

    do k = 1, nrad
      rad(k) = rmin + (rmax - rmin)*real(k - 1, dp)/real(nrad - 1, dp)
    end do
    do l = 1, nzet
      zet(l) = zmin + (zmax - zmin)*real(l - 1, dp)/real(nzet - 1, dp)
    end do
    do l = 1, nzet
      do k = 1, nrad
        bre(k, l) = f_re(rad(k), zet(l))
        bim(k, l) = f_im(rad(k), zet(l))
      end do
    end do
    open (newunit=u, file='bmod_n.dat', form='unformatted', status='replace')
    write (u) nrad, nzet
    write (u) rad, zet
    write (u) bre, bim
    close (u)
  end subroutine write_grid

end program test_bmod_pert_bounds
