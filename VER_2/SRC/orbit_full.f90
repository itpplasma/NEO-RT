module orbit
  use common
  implicit none

  contains

  subroutine timestep(tau, z, vz)
  ! See velo for details
    real(8) :: tau   ! Time
    real(8) :: z(5)  ! Phase-position
    real(8) :: vz(5) ! Phase-velocity

    call velo(tau, z, vz)
  end subroutine timestep

  subroutine bounce_average(n, z, integrand, out)
  ! Returns bounce average of an quantity via orbit integration
    integer(4), intent(in) :: n  ! Number of components of the integrand
    real(8), intent(in)    :: z(neqm)    ! Starting position
    external               :: integrand  ! Subroutine
    real(8), intent(out)   :: out(n)     ! Bounce average output

    real(8) :: dtau, taub, delphi

    dtau = 2.0d0  ! TODO: time step is for strongly passing

    out = 0.0d0

    call find_bounce(n, timestep_ext, dtau, z, taub, delphi, out)
    out = out/taub
    print *, 'Out: '
    print *, out

    contains

    subroutine timestep_ext(t, y, dydt)
      real(8), intent(in) :: t
      real(8), intent(in) :: y(neqm+n)
      real(8), intent(out) :: dydt(neqm+n)

      call timestep(t, y(1:neqm), dydt(1:neqm))
      call integrand(y(1:neqm), dydt(neqm+1:(neqm+n)))

      write(999, *) y
      !print *, y

    end subroutine timestep_ext

  end subroutine bounce_average

  subroutine bounce_average_box
  end subroutine bounce_average_box
end module orbit
