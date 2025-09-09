module diag_dnorm_breakdown
  use iso_fortran_env, only: real64
  use fortplot, only: figure, plot, title, xlabel, ylabel, legend, savefig
  use util, only: qe, mu, c
  use collis_alp, only: coleff
  use neort, only: read_control, init, check_magfie, runname => runname
  use neort_profiles, only: init_profile_input, init_plasma_input, init_profiles, vth
  use neort_nonlin, only: nonlinear_attenuation, omega_prime
  use neort_freq, only: Om_th, Om_ph, d_Om_ds
  use neort_transport, only: timestep_transport
  use neort_orbit, only: bounce_fast, nvar
  use driftorbit, only: nonlin, mth, mph, etatp, etadt, epsst_spl, mi, qi, iota, psi_pr, pertfile
  use do_magfie_mod, only: R0, s, sign_theta, q
  use do_magfie_mod, only: do_magfie_init
  use do_magfie_pert_mod, only: do_magfie_pert_init
  implicit none

contains

  subroutine run_dnorm_breakdown(arg_runname)
    character(*), intent(in) :: arg_runname
    logical :: file_exists
    real(real64) :: eta, eta_min, eta_max
    integer :: j, i, nm, u
    integer :: mth_min, mth_max
    real(real64), parameter :: ux_list(3) = [1.0_real64, 2.0_real64, 3.0_real64]
    real(real64) :: ux, v
    real(real64) :: Omth, dOmthdv, dOmthdeta
    real(real64) :: Omph, dOmphdv, dOmphdeta
    real(real64) :: dOmthds, dOmphds
    real(real64) :: dOmdv, dOmdeta, dOmdpph, Ompr
    real(real64) :: taub, bounceavg(nvar)
    real(real64) :: Hmn2, dpp, dhh, fpeff, dres, dnorm, theta

    ! Initialize like main
    runname = trim(arg_runname)
    call read_control
    call do_magfie_init()
    if (pertfile) call do_magfie_pert_init()
    call init_profiles(R0)
    inquire(file="plasma.in", exist=file_exists)
    if (file_exists) call init_plasma_input(s)
    inquire(file="profile.in", exist=file_exists)
    if (file_exists) call init_profile_input(s, R0, 1.0_real64, 1.0_real64)
    call init
    call check_magfie

    ! Use the same eta_max definition as bounce_nonlin
    eta_min = etatp*(1.0_real64 + 1.0e-4_real64)
    eta_max = etatp + (etadt - etatp)*(1.0_real64 - epsst_spl)
    eta = eta_max

    mth_min = -ceiling(2.0_real64*abs(mph*q))
    mth_max =  ceiling(2.0_real64*abs(mph*q))
    nm = mth_max - mth_min + 1

    open(newunit=u, file=trim(arg_runname)//'_dnorm_breakdown.txt', status='replace', action='write')
    write(u,'(A)') '# dnorm breakdown at eta_max (near trapped boundary)'
    write(u,'(A,F10.6)') '# s = ', s
    write(u,'(A,ES12.5)') '# eta_min = ', eta_min
    write(u,'(A,ES12.5)') '# eta_max = ', eta_max
    write(u,'(A,ES12.5)') '# eta_used = ', eta
    write(u,'(A)') '# columns: mth ux Ompr dOmdv dOmdeta dOmdpph Omth dOmthdv dOmthdeta Omph dOmphdv dOmphdeta dOmthds dOmphds dpp dhh Hmn2 dres dnorm theta'

    ! enable nonlin just to compute bounceavg(5:6) = <1/B>, <B>
    nonlin = .true.
    do j = mth_min, mth_max
      mth = j
      do i = 1, size(ux_list)
        ux = ux_list(i)
        v = ux*vth

        call Om_th(v, eta, Omth, dOmthdv, dOmthdeta)
        call Om_ph(v, eta, Omph, dOmphdv, dOmphdeta)
        call d_Om_ds(v, eta, dOmthds, dOmphds)

        dOmdv   = mth*dOmthdv   + mph*dOmphdv
        dOmdeta = mth*dOmthdeta + mph*dOmphdeta
        dOmdpph = -(qi/c*iota*sign_theta*psi_pr)**(-1) * (mth*dOmthds + mph*dOmphds)

        taub = 2.0_real64*acos(-1.0_real64)/abs(Omth)
        call bounce_fast(v, eta, taub, bounceavg, timestep_transport)
        Hmn2 = (bounceavg(3)**2 + bounceavg(4)**2)*(mi*(v*v/2.0_real64))**2

        call coleff(ux, dpp, dhh, fpeff)
        ! Match scaling used in nonlinear_attenuation
        dhh = vth*dhh
        dpp = vth**3*dpp

        Ompr = omega_prime(ux, eta, bounceavg, Omth, dOmdv, dOmdeta, dOmdpph)
        dres = dpp*(dOmdv/Ompr)**2 + dhh*eta*(bounceavg(5) - eta)*(dOmdeta/Ompr)**2
        dnorm = dres*sqrt(abs(Ompr)) / (abs(Hmn2)**(3.0_real64/2.0_real64))
        call attenuation_factor(dnorm, theta)

        write(u,'(I4,1X,F5.2,1X,18(1X,ES14.6))') mth, ux, Ompr, dOmdv, dOmdeta, dOmdpph, &
              Omth, dOmthdv, dOmthdeta, Omph, dOmphdv, dOmphdeta, dOmthds, dOmphds, &
              dpp, dhh, Hmn2, dres, dnorm, theta
      end do
    end do
    close(u)

  end subroutine run_dnorm_breakdown

end module diag_dnorm_breakdown
