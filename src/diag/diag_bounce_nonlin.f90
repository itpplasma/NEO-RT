module diag_bounce_nonlin
  use iso_fortran_env, only: real64
  use fortplot, only: figure, plot, title, xlabel, ylabel, legend, savefig
  use neort, only: read_control, init, check_magfie, runname => runname
  use neort_profiles, only: init_profile_input, init_plasma_input, init_profiles, vth
  use neort_nonlin, only: nonlinear_attenuation
  use neort_freq, only: Om_th
  use neort_transport, only: timestep_transport
  use neort_orbit, only: bounce_fast, nvar
  use driftorbit, only: nonlin, mth, mph, etatp, etadt, epsst_spl
  use do_magfie_mod, only: R0, s
  use do_magfie_mod, only: do_magfie_init
  use do_magfie_pert_mod, only: do_magfie_pert_init
  implicit none

contains

  subroutine run_bounce_nonlin_diag(arg_runname)
    character(*), intent(in) :: arg_runname
    logical :: file_exists
    real(real64) :: v, eta_min, eta_max
    integer :: i, npts
    real(real64), allocatable :: eta(:), att_lin(:), att_nonlin(:)
    real(real64) :: taub, bounceavg(nvar), Omth, dOmthdv, dOmthdeta
    real(real64) :: ux, Hmn2

    ! Initialize NEO-RT environment similar to neort:main
    runname = trim(arg_runname)
    call read_control
    call do_magfie_init()
    call do_magfie_pert_init()
    call init_profiles(R0)

    inquire(file="plasma.in", exist=file_exists)
    if (file_exists) call init_plasma_input(s)

    inquire(file="profile.in", exist=file_exists)
    if (file_exists) call init_profile_input(s, R0, 1.0_real64, 1.0_real64)

    call init
    call check_magfie

    ! Sample attenuation over trapped pitch range at v=vth (ux=1)
    v = vth
    ux = 1.0_real64

    eta_min = etatp*(1.0_real64 + 1.0e-4_real64)
    eta_max = etatp + (etadt - etatp)*(1.0_real64 - epsst_spl)
    npts = 100
    allocate(eta(npts), att_lin(npts), att_nonlin(npts))
    do i = 1, npts
      eta(i) = eta_min + (eta_max-eta_min)*(real(i-1,real64)/real(npts-1,real64))
    end do

    do i = 1, npts
      call Om_th(v, eta(i), Omth, dOmthdv, dOmthdeta)
      taub = 2.0_real64*acos(-1.0_real64)/abs(Omth)
      call bounce_fast(v, eta(i), taub, bounceavg, timestep_transport)
      Hmn2 = (bounceavg(3)**2 + bounceavg(4)**2)*(0.5_real64*v**2)**2

      nonlin = .false.
      att_lin(i) = nonlinear_attenuation(ux, eta(i), bounceavg, Omth, dOmthdv, dOmthdeta, Hmn2)

      nonlin = .true.
      att_nonlin(i) = nonlinear_attenuation(ux, eta(i), bounceavg, Omth, dOmthdv, dOmthdeta, Hmn2)
    end do

    call figure()
    call plot(eta, att_lin, label="linear (nonlin=F)")
    call plot(eta, att_nonlin, label="nonlinear (nonlin=T)")
    call title("Nonlinear attenuation vs eta")
    call xlabel("eta")
    call ylabel("attenuation factor")
    call legend()
    call savefig(trim(arg_runname)//"_bounce_nonlin.pdf")

  end subroutine run_bounce_nonlin_diag

end module diag_bounce_nonlin

