module orbit
  use common
  implicit none

  integer(4), parameter :: neqm=5  ! Number of equations of motion

  interface
    module subroutine timestep(tau, z, vz)
      ! See velo for details
      real(8) :: tau   ! Time
      real(8) :: z(5)  ! Phase-position
      real(8) :: vz(5) ! Phase-velocity
    end subroutine timestep

    module subroutine bounce_average(n, z, integrand, taub, delphi, ret)
      ! Returns bounce average of an quantity via orbit integration
      integer(4), intent(in) :: n        ! Number of components of the integrand
      real(8), intent(in)    :: z(neqm)  ! Starting position
      external               :: integrand  ! Routine fcn(t, z, ret) to integrate
      real(8), intent(out)   :: taub       ! Bounce time output
      real(8), intent(out)   :: delphi     ! Change in varphi during bounce time
      real(8), intent(out)   :: ret(n)     ! Bounce average output
    end subroutine bounce_average

    module subroutine bounce_integral_box(n, z, integrand, sbox, tau, ret)
    ! Returns bounce average of an quantity via orbit integration
      integer(4), intent(in) :: n        ! Number of components of the integrand
      real(8), intent(in)    :: z(neqm)  ! Starting position
      real(8), intent(in)    :: sbox(:)  ! Box boundaries
      real(8), intent(out)   :: tau(size(sbox))    ! Time in each box
      real(8), intent(out)   :: ret(n,size(sbox))  ! Integral over each box

      external :: integrand   ! Subroutine f(z, ret) to integrate
    end subroutine bounce_integral_box

    module subroutine bounce_harmonic(next, z, fn, mb, nph, taub, delphi, ret)
      ! Computes bounce harmonic mb of fn
      integer(4), intent(in) :: next       ! Number of extra integrals
      real(8), intent(inout) :: z(neqm)    ! Position on orbit
      external               :: fn         ! Subroutine fn(z, out) to treat
      integer(4), intent(in) :: mb         ! Bounce harmonic number
      integer(4), intent(in) :: nph        ! Toroidal harmonic number
      real(8), intent(out)   :: ret(next)  ! Complex harmonic of input fn
      real(8), intent(out)   :: taub     ! Bounce time
      real(8), intent(out)   :: delphi   ! Change in varphi during bounce time
    end subroutine bounce_harmonic
  end interface
end module orbit
