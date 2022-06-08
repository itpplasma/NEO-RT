subroutine linspace(a, b, cnt, out)
  implicit none

  integer :: cnt
  real(8) :: out(cnt)
  integer :: i
  real(8) :: a, b, delta

  delta = (b-a)/(cnt-1)
  out = a + delta*[(i,i=0,cnt-1)]
end subroutine linspace

subroutine timestep_vode(n, tau, z, vz)
  implicit none
  ! See velo for details
  integer(4), intent(in) :: n     ! number of equations
  real(8), intent(in)    :: tau   ! Time
  real(8), intent(in)    :: z(n)  ! Phase-position
  real(8), intent(out)   :: vz(n) ! Phase-velocity
  call velo(tau, z, vz)
end subroutine timestep_vode

subroutine time_in_box(z, sbox, taub, tau)
  ! Returns time spent in boxes
  use dvode_f90_m, only: vode_opts, set_normal_opts, dvode_f90, get_stats
  use field_eq_mod, only: psif, psi_axis, psi_sep
  use orbit_dim_mod, only: neqm

  implicit none

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

  external timestep_vode

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
