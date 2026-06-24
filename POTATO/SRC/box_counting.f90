subroutine linspace(a, b, cnt, out)
  implicit none

  real(8), intent(in) :: a, b
  integer :: i, cnt
  real(8), intent(inout) :: out(cnt)
  real(8) :: delta

  delta = (b-a)/(cnt-1)
  do i=1,cnt
    out(i) = a + delta*(i-1)
  end do
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

subroutine time_in_box(z, cnt, sbox, taub, taub_max, tau)
  ! Returns time spent in boxes
  use dvode_f90_m, only: vode_opts, set_opts, dvode_f90, get_stats
  use field_sub, only: psif
  use field_eq_mod, only: psi_axis, psi_sep
  use orbit_dim_mod, only: neqm

  implicit none

  real(8), intent(in)    :: z(neqm)    ! Starting position
  integer, intent(in)    :: cnt        ! Box boundary count
  real(8), intent(in)    :: sbox(cnt)  ! Box boundaries
  real(8), intent(in)    :: taub       ! Bounce time
  real(8), intent(in)    :: taub_max   ! Largest integrable bounce time (find_bounce cap)
  real(8), intent(out)   :: tau(cnt)   ! Time in each box

  real(8) :: smid
  real(8) :: y(neqm)

  integer(4) :: k, nsample
  real(8) :: ti

  real(8) :: atol(neqm), rtol, tout
  integer(4) :: itask, istate, method_flag
  type (vode_opts) :: options
  real(8) :: bmod, phi_elec, s
  real(8) :: sold, told
  integer(4) :: sind ! s index

  external timestep_vode

  tau = 0d0

  ! A near-separatrix resonance root carries an interpolated Omega_b ~ 0, so its
  ! recovered bounce time |taub| runs orders of magnitude past taub_max (the
  ! find_bounce acceptance cap, 200*dtau). VODE cannot integrate such an orbit --
  ! it grinds to mxstep on every subinterval and returns a zero residence
  ! distribution anyway -- so skip it up front, the same way find_bounce skips
  ! the matching Omega_b=0 grazer. taub<=0 (interpolant overshoot through zero)
  ! is likewise unintegrable. Same zero contribution, no grind, no log flood.
  if (taub <= 0d0 .or. taub > taub_max) return

  y = z

  rtol = 1d-12
  atol = 1d-13
  itask = 1
  istate = 1
  method_flag = 10
  options = set_opts(method_flag=method_flag, abserr_vector=atol, &
    relerr=rtol, mxstep=2000)

  nsample = min(256, max(64, size(sbox)))
  ti = 0d0

  call get_bmod_and_Phi(z(1:3), bmod, phi_elec)

  s = abs((psif-psi_axis)/(psi_sep-psi_axis))
  sold = s
  told = 0d0
  do k = 1,nsample
    tout = taub*dble(k)/dble(nsample)
    call dvode_f90( &
        timestep_vode, neqm, y, ti, tout, itask, istate, options)
    if (istate == -1) then
      print *, 'time_in_box: VODE exceeded mxstep at k =', k, &
        ' ti =', ti, ' taub =', taub
      tau = 0d0
      return
    endif
    if (istate /= 2) then
      print *, 'time_in_box: VODE failed with istate =', istate, &
        ' k =', k, ' ti =', ti, ' taub =', taub
      tau = 0d0
      return
    end if

    call get_bmod_and_Phi(y(1:3), bmod, phi_elec)

    s = abs((psif-psi_axis)/(psi_sep-psi_axis))
    smid = 0.5d0*(sold+s)
    sind = radial_box_index(smid)
    tau(sind) = tau(sind) + ti-told
    sold = s
    told = ti
  end do

contains

  integer function radial_box_index(sval)
    real(8), intent(in) :: sval
    integer :: ibox

    radial_box_index = size(sbox)
    do ibox = 1,size(sbox)
      if (sval <= sbox(ibox)) then
        radial_box_index = ibox
        exit
      endif
    enddo
  end function radial_box_index
end subroutine time_in_box
