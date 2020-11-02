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

  subroutine bounce_average(n, z, integrand, res)
  ! Returns bounce average of an quantity via orbit integration
    integer(4), intent(in) :: n  ! Number of components of the integrand
    real(8), intent(in)    :: z(neqm)    ! Starting position
    external               :: integrand  ! Subroutine
    real(8), intent(out)   :: res(n)     ! Bounce average output

    real(8) :: dtau, taub, delphi

    dtau = 2.0d0  ! TODO: time step is for strongly passing

    res = 0.0d0

    call find_bounce(n, timestep_ext, dtau, z, taub, delphi, res)
    res = res/taub
    print *, 'Out: '
    print *, res

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

  subroutine bounce_harmonic(next, z, mb, nph, res, fn)
    ! Computes bounce harmonic mb of fn
      integer(4), intent(in) :: next       ! Number of extra integrals
      real(8), intent(inout) :: z(5)       ! Position on orbit
      integer(4), intent(in) :: mb         ! Bounce harmonic number
      integer(4), intent(in) :: nph        ! Toroidal harmonic number
      real(8), intent(out)   :: res(next)  ! Complex harmonic of input fn
      external               :: fn         ! Subroutine fn(z, out) to treat

      ! TODO
      call bounce_average(2, z, fn, res)

    end subroutine bounce_harmonic
end module orbit
