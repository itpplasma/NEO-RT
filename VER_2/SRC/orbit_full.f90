module orbit
use common
implicit none

integer(4), parameter :: neqm=5  ! Number of equations of motion

contains

subroutine timestep(tau, z, vz)
  ! See velo for details
  real(8), intent(in)    :: tau   ! Time
  real(8), intent(in)    :: z(neqm)  ! Phase-position
  real(8), intent(in)    :: vz(neqm) ! Phase-velocity
  call velo(tau, z, vz)
end subroutine timestep


subroutine bounce_average(n, z, integrand, taub, delphi, ret)
  ! Returns bounce average of an quantity via orbit integration
  integer(4), intent(in) :: n        ! Number of components of the integrand
  real(8), intent(in)    :: z(neqm)  ! Starting position
  external               :: integrand  ! Routine fcn(t, z, ret) to integrate
  real(8), intent(out)   :: taub       ! Bounce time output
  real(8), intent(out)   :: delphi     ! Change in varphi during bounce time
  real(8), intent(out)   :: ret(n)     ! Bounce average output
  real(8) :: dtau

  dtau = 2.0d0  ! TODO: time step is for strongly passing

  ret = 0.0d0

  call find_bounce(n, timestep_ext, dtau, z, taub, delphi, ret)

  contains

  subroutine timestep_ext(t, y, dydt)
    real(8), intent(in)    :: t
    real(8), intent(in)    :: y(neqm+n)
    real(8), intent(out)   :: dydt(neqm+n)

    call timestep(t, y(1:neqm), dydt(1:neqm))
    call integrand(t, y(1:neqm), dydt(neqm+1:(neqm+n)))
  end subroutine timestep_ext

end subroutine bounce_average


subroutine bounce_integral_box(n, z, integrand, sbox, tau, ret)
! Returns bounce average of an quantity via orbit integration
  use dvode_f90_m, only: vode_opts, set_normal_opts, dvode_f90, get_stats
  use field_eq_mod, only: psif, psi_axis, psi_sep

  integer(4), intent(in) :: n        ! Number of components of the integrand
  real(8), intent(in)    :: z(neqm)  ! Starting position
  real(8), intent(in)    :: sbox(:)  ! Box boundaries
  real(8), intent(out)   :: tau(size(sbox))    ! Time in each box
  real(8), intent(out)   :: ret(n,size(sbox))  ! Integral over each box
  external :: integrand   ! Subroutine f(z, ret) to integrate

  real(8) :: dtau, taub, delphi
  real(8) :: sprev, snext ! Previous and next flux radius box
  real(8) :: y(neqm+n)

  integer(4) :: nmax  ! Maximum loop iterations

  integer(4) :: k
  real(8) :: ti

  real(8) :: atol(neqm+n), rtol, tout, rstats(22)
  integer(4) :: neq, itask, istate, istats(31), numevents
  type (vode_opts) :: options
  real(8) :: bmod, phi_elec, s
  real(8) :: sold, told, yold(neqm+n)
  integer(4) :: sind, sind0 ! s index

  integer(4) :: jroots(2)

  ! TODO make this more efficient
  dtau = 2.0d0  ! TODO: time step is for strongly passing
  ret = 0.0d0
  y = 0.0d0
  y(1:neqm) = z
  call find_bounce(n, timestep_ext_forbounce, dtau, z, taub, delphi, ret)

  neq = neqm + n
  rtol = 1d-12
  atol = 1d-13
  itask = 1
  istate = 1
  numevents = 2
  options = set_normal_opts(abserr_vector=atol, relerr=rtol, nevents=numevents)

  nmax = 3*size(sbox)
  tau = 0d0
  ti = 0d0

  call get_bmod_and_Phi(z(1:3), bmod, phi_elec)

  s = abs((psif-psi_axis)/(psi_sep-psi_axis))
  sind = size(sbox)+1
  sprev = -1e5
  snext = 1e5
  do k = 1,size(sbox)
      if(sbox(k)>s) then
        sind = k
        snext = sbox(k)
        if(k>1) sprev = sbox(k-1)
        exit
      end if
  enddo
  sind0 = sind

  told = 0d0
  do k = 1,nmax
      yold = y
      sold = s
      tout = taub
      call dvode_f90( &
        timestep_ext, neq, y, ti, tout, itask, istate, options, g_fcn = sroots &
      )
      if (istate == 2) exit
      if (istate == 3) then
        tau(sind) = tau(sind) + ti-told
        told = ti
        call get_stats(rstats, istats, numevents, jroots)
        if (jroots(2).ne.0) then
            sind = sind + 1 ! moving outwards
            sprev = snext
            if (sind == size(sbox)+1) then
              snext = 1e5
            else
              snext = sbox(sind)
            end if
        end if
        if (jroots(1).ne.0) then
            sind = sind - 1 ! moving inwards
            snext = sprev
            if (sind == 1) then
              sprev = -1e5
            else
              sprev = sbox(sind-1)
            end if
        end if
      end if
  end do

  tau(sind) = tau(sind) + taub-told

  contains

  subroutine timestep_ext_forbounce(t, yext, dydt)
    real(8), intent(in)    :: t
    real(8), intent(in)    :: yext(neqm+n)
    real(8), intent(out)   :: dydt(neqm+n)

    call timestep(t, yext(1:neqm), dydt(1:neqm))
    call integrand(yext(1:neqm), dydt(neqm+1:(neqm+n)))
  end subroutine timestep_ext_forbounce

  subroutine timestep_ext(neqext, t, yext, dydt)
    integer(4), intent(in) :: neqext
    real(8), intent(in)    :: t
    real(8), intent(in)    :: yext(neqext)
    real(8), intent(out)   :: dydt(neqext)

    call timestep(t, yext(1:neqm), dydt(1:neqm))
    call integrand(yext(1:neqm), dydt(neqm+1:neqext))
  end subroutine timestep_ext

  subroutine sroots(neqext, t, yext, ng, gout)
    ! For finding roots between boxes

    integer, intent(in) :: neqext, ng
    real(8), intent(in) :: t, yext(neqext)
    real(8), intent(out) :: gout(ng)

    call get_bmod_and_Phi(yext(1:3), bmod, phi_elec)

    s = abs((psif-psi_axis)/(psi_sep-psi_axis))

    gout(1) = s - sprev
    gout(2) = s - snext
  end subroutine sroots
end subroutine bounce_integral_box


subroutine bounce_harmonic(next, z, fn, mb, nph, taub, delphi, ret)
  ! Computes bounce harmonic mb of fn
  integer(4), intent(in) :: next       ! Number of extra integrals
  real(8), intent(inout) :: z(neqm)    ! Position on orbit
  external               :: fn         ! Subroutine fn(z, out) to treat
  integer(4), intent(in) :: mb         ! Bounce harmonic number
  integer(4), intent(in) :: nph        ! Toroidal harmonic number
  real(8), intent(out)   :: ret(next)  ! Complex harmonic of input fn
  real(8), intent(out)   :: taub     ! Bounce time
  real(8), intent(out)   :: delphi   ! Change in varphi during bounce time
  real(8) :: omb, omphi
  ! TODO
  call bounce_average(next, z, fn, taub, delphi, ret)

  omb = 2d0*pi/taub
  omphi = delphi/taub

  call bounce_average(next, z, integrand, taub, delphi, ret)

  contains

  subroutine integrand(t, y, res)
  ! Integrand for Fourier integral: real and imaginary part
    real(8), intent(in)  :: t           ! Orbit time parameter
    real(8), intent(in)  :: y(5)        ! Orbit phase-space variables
    real(8), intent(out) :: res(2)      ! Output

    real(8) :: bmod, sqrtg
    real(8), dimension(3) :: bder,hcovar,hctrvr,hcurl

    real(8) :: fnres(2)
    complex(8) :: fnval, expfac, resval

    call fn(t, y, fnres)
    fnval = cmplx(fnres(1), fnres(2))
    expfac = exp(-imun*(mb*omb*t + nph*(y(2) - omphi*t)))
    resval = fnval*expfac

    res(1) = real(resval)
    res(2) = aimag(resval)
  end subroutine integrand
end subroutine bounce_harmonic

end module orbit
