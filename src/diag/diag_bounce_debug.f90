module diag_bounce_debug
  use iso_fortran_env, only: real64
  use neort, only: read_and_set_control, init, check_magfie, runname => runname, &
                   set_to_passing_region, set_to_trapped_region, vsteps
  use neort_profiles, only: read_and_init_profile_input, read_and_init_plasma_input, init_profiles, vth, Om_tE
  use neort_freq, only: Om_th
  use neort_transport, only: timestep_transport
  use neort_orbit, only: nvar
  use neort_resonance, only: driftorbit_coarse, driftorbit_root
  use driftorbit, only: mth, mph, mi, pertfile, etamin, etamax, sign_vpar, &
                        etatp, etadt, efac
  use do_magfie_mod, only: R0, s, q, bfac, do_magfie_init
  use do_magfie_pert_mod, only: do_magfie_pert_init
  use dvode_f90_m, only: dvode_f90, set_normal_opts, vode_opts, get_stats
  implicit none

contains

  subroutine run_bounce_debug(arg_runname)
    character(*), intent(in) :: arg_runname
    logical :: file_exists
    integer :: i, j, mmin, mmax, nm, u
    real(real64) :: vminp, vmaxp, vmint, vmaxt
    real(real64) :: du, ux, v
    real(real64) :: roots(100,3), eta_res(2), eta
    integer :: nroots, kr
    real(real64) :: Omth, dOmthdv, dOmthdeta, taub
    integer :: istate
    integer :: istats(50)
    real(real64) :: rstats(50)

    ! Initialize environment just like main
    runname = trim(arg_runname)
    call read_and_set_control
    call do_magfie_init()
    if (pertfile) call do_magfie_pert_init()
    call init_profiles(R0)

    inquire(file="plasma.in", exist=file_exists)
    if (file_exists) call read_and_init_plasma_input("plasma.in", s)

    inquire(file="profile.in", exist=file_exists)
    if (file_exists) call read_and_init_profile_input("profile.in", s, R0, efac, bfac)

    call init
    call check_magfie

    vminp = 1.0d-6*vth
    vmaxp = 3.0d0*vth
    vmint = vminp
    vmaxt = vmaxp

    mmin = -ceiling(2.0_real64*abs(mph*q))
    mmax =  ceiling(2.0_real64*abs(mph*q))
    nm = mmax - mmin + 1

    open(newunit=u, file=trim(arg_runname)//'_bounce_debug.txt', status='replace', action='write')
    write(u,'(A)') '# bounce debug: logs where integration stalls (MXSTEP)'
    write(u,'(A,F10.6)') '# s = ', s
    write(u,'(A)') '# cols: branch mth ux eta Omth taub istate nst accrej last_h'

    do j = 1, nm
      mth = mmin + (j - 1)

      ! Passing co
      sign_vpar = 1.0_real64
      call set_to_passing_region(etamin, etamax)
      call scan_branch('cop', vminp, vmaxp, vsteps, u)

      ! Passing ctr
      sign_vpar = -1.0_real64
      call set_to_passing_region(etamin, etamax)
      call scan_branch('ctr', vminp, vmaxp, vsteps, u)

      ! Trapped
      sign_vpar = 1.0_real64
      call set_to_trapped_region(etamin, etamax)
      call scan_branch('trap', vmint, vmaxt, vsteps, u)
    end do
    close(u)
  end subroutine run_bounce_debug

  subroutine scan_branch(label, vmin, vmax, vsteps_loc, u)
    character(*), intent(in) :: label
    real(real64), intent(in) :: vmin, vmax
    integer, intent(in) :: vsteps_loc
    integer, intent(in) :: u
    real(real64) :: du, ux, v
    real(real64) :: roots(100,3), eta_res(2), eta
    integer :: nroots, kr
    real(real64) :: Omth, dOmthdv, dOmthdeta, taub
    integer :: istate
    integer :: istats(50)
    real(real64) :: rstats(50)

    du = (vmax - vmin)/(real(vsteps_loc, real64)*vth)
    ux = vmin/vth + du/2.0_real64
    do while (ux*vth <= vmax + 1d-12)
      v = ux*vth
      call driftorbit_coarse(v, etamin, etamax, roots, nroots)
      if (nroots > 0) then
        do kr = 1, nroots
          eta_res = driftorbit_root(v, 1d-8*abs(Om_tE), roots(kr,1), roots(kr,2))
          eta = eta_res(1)
          call Om_th(v, eta, Omth, dOmthdv, dOmthdeta)
          taub = 2.0_real64*acos(-1.0_real64)/abs(Omth)
          call probe_bounce(v, eta, taub, istate, istats, rstats)
          write(u,'(A,1X,I4,1X,F8.4,1X,ES12.5,1X,ES12.5,1X,ES12.5,1X,I4,1X,I9,1X,ES12.5)') &
               trim(label), mth, ux, eta, Omth, taub, istate, istats(1), rstats(6)
        end do
      end if
      ux = ux + du
    end do
  end subroutine scan_branch

  subroutine probe_bounce(v, eta, taub, istate, istats, rstats)
    real(real64), intent(in) :: v, eta, taub
    integer, intent(out) :: istate
    integer, intent(out) :: istats(:)
    real(real64), intent(out) :: rstats(:)

    real(real64) :: t1, t2
    real(real64) :: y(nvar)
    real(real64) :: atol(nvar), rtol
    integer :: neq, itask
    type(vode_opts) :: options

    t1 = 0.0_real64
    t2 = taub
    y = 1.0e-15_real64
    y(1) = 0.0_real64
    y(2) = 0.0_real64
    y(3:6) = 0.0_real64
    neq = nvar
    rtol = 1.0e-8_real64
    atol = 1.0e-9_real64
    itask = 1
    istate = 1
    options = set_normal_opts(abserr_vector=atol, relerr=rtol)

    call dvode_f90(timestep_wrapper, neq, y, t1, t2, itask, istate, options)
    call get_stats(rstats, istats)

  contains
    subroutine timestep_wrapper(neq_, t_, y_, ydot_)
      integer, intent(in) :: neq_
      real(real64), intent(in) :: t_
      real(real64), intent(in) :: y_(neq_)
      real(real64), intent(out) :: ydot_(neq_)
      call timestep_transport(v, eta, neq_, t_, y_, ydot_)
    end subroutine timestep_wrapper
  end subroutine probe_bounce

end module diag_bounce_debug

