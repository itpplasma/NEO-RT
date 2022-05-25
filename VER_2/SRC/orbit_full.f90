module orbit
use util
implicit none

integer(4), parameter :: neqm=5  ! Number of equations of motion
real(8) :: dtau=0.0

contains

subroutine orbit_init(a_dtau)
  real(8) :: a_dtau
  dtau = a_dtau
end subroutine orbit_init

subroutine timestep(tau, z, vz)
  ! See velo for details
  real(8), intent(in)    :: tau   ! Time
  real(8), intent(in)    :: z(neqm)  ! Phase-position
  real(8), intent(out)   :: vz(neqm) ! Phase-velocity
  call velo(tau, z, vz)
end subroutine timestep


subroutine timestep_vode(n, tau, z, vz)
  ! See velo for details
  integer(4), intent(in) :: n     ! number of equations
  real(8), intent(in)    :: tau   ! Time
  real(8), intent(in)    :: z(n)  ! Phase-position
  real(8), intent(out)   :: vz(n) ! Phase-velocity
  call timestep(tau, z, vz)
end subroutine timestep_vode

subroutine bounce_average(n, z, integrand, taub, delphi, ret)
  ! Returns bounce average of an quantity via orbit integration
  integer(4), intent(in) :: n        ! Number of components of the integrand
  real(8), intent(in)    :: z(neqm)  ! Starting position
  external               :: integrand  ! Routine fcn(t, z, ret) to integrate
  real(8), intent(out)   :: taub       ! Bounce time output
  real(8), intent(out)   :: delphi     ! Change in varphi during bounce time
  real(8), intent(out)   :: ret(n)     ! Bounce average output

  ret = 0.0d0

  call find_bounce(n, timestep_ext, dtau, z, taub, delphi, ret)
  ret = ret/taub

  contains

  subroutine timestep_ext(t, y, dydt)
    real(8), intent(in)    :: t
    real(8), intent(in)    :: y(neqm+n)
    real(8), intent(out)   :: dydt(neqm+n)

    call timestep(t, y(1:neqm), dydt(1:neqm))
    call integrand(t, y(1:neqm), dydt((neqm+1):(neqm+n)))
  end subroutine timestep_ext

end subroutine bounce_average


subroutine time_in_box(z, sbox, taub, tau)
! Returns time spent in boxes
  use dvode_f90_m, only: vode_opts, set_normal_opts, dvode_f90, get_stats
  use field_eq_mod, only: psif, psi_axis, psi_sep

  real(8), intent(in)    :: z(neqm)  ! Starting position
  real(8), intent(in)    :: sbox(:)  ! Box boundaries
  real(8), intent(in)    :: taub     ! Bounce time
  real(8), intent(out)   :: tau(size(sbox))    ! Time in each box

  real(8) :: delphi
  real(8) :: sprev, snext ! Previous and next flux radius box
  real(8) :: y(neqm)

  integer(4) :: nmax  ! Maximum loop iterations

  integer(4) :: k
  real(8) :: ti

  real(8) :: atol(neqm), rtol, tout, rstats(22)
  integer(4) :: itask, istate, istats(31), numevents
  type (vode_opts) :: options
  real(8) :: bmod, phi_elec, s
  real(8) :: sold, told, yold(neqm)
  integer(4) :: sind, sind0 ! s index

  integer(4) :: jroots(2)

  ! TODO make this more efficient
  y = z

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
        timestep_vode, neqm, y, ti, tout, itask, istate, options, &
        g_fcn = sroots )
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
end subroutine time_in_box


subroutine bounce_harmonic(z, fn, mb, nph, taub, deltaphi, ret)
  ! Computes bounce harmonic mb of complex fn
  integer(4), parameter  :: next=3       ! Number of extra integrals
  real(8), intent(inout) :: z(neqm)      ! Position on orbit
  external               :: fn           ! Subroutine fn(z, out) to treat
  integer(4), intent(in) :: mb           ! Bounce harmonic number
  integer(4), intent(in) :: nph          ! Toroidal harmonic number
  real(8), intent(out)   :: ret(next-1)  ! Complex harmonic of input fn
  real(8), intent(out)   :: taub     ! Bounce time
  real(8), intent(out)   :: deltaphi   ! Change in varphi during bounce time
  real(8) :: omb, omphi, ret_full(next)

  ! This includes also variable to track orbit time
  ret_full = 0.0d0
  call find_bounce(0, timestep, dtau, z, taub, deltaphi, ret_full)
  omb = 2d0*pi/taub
  omphi = deltaphi/taub

  ret_full = 0.0d0
  call find_bounce(next, timestep_ext, dtau, z, taub, deltaphi, ret_full)
  ret = ret_full(2:next)/taub

  contains

  subroutine timestep_ext(t, y, res)
  ! Integrand for Fourier integral: real and imaginary part
    real(8), intent(in)  :: t           ! Orbit time parameter
    real(8), intent(in)  :: y(neqm+next)      ! Orbit phase-space variables
    real(8), intent(out) :: res(neqm+next)    ! Output

    real(8) :: fnres(2)
    complex(8) :: fnval, expfac, resval

    call timestep(t, y(1:neqm), res(1:neqm))

    call fn(t, y(1:neqm), fnres)
    fnval = cmplx(fnres(1), fnres(2), 8)
    !expfac = exp(-imun*mb*omb*y(6))
    expfac = exp(-imun*((mb*omb+nph*omphi)*y(6) - nph*y(2)))
    resval = fnval*expfac

    ! print *, y(6)*omb/(2d0*pi), y(6)*omphi/(2d0*pi), y(2)

    res(6) = 1d0
    res(7) = real(resval)
    res(8) = aimag(resval)
  end subroutine timestep_ext
end subroutine bounce_harmonic

end module orbit
