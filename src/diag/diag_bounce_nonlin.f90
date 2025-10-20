module diag_bounce_nonlin
  use iso_fortran_env, only: real64
  use fortplot, only: figure, plot, title, xlabel, ylabel, legend, savefig, xlim, ylim
  use neort, only: read_and_set_control, init, check_magfie, runname => runname
  use neort_profiles, only: read_and_init_profile_input, read_and_init_plasma_input, init_profiles, vth
  use neort_nonlin, only: nonlinear_attenuation
  use neort_freq, only: Om_th
  use neort_transport, only: timestep_transport
  use neort_orbit, only: bounce_fast, nvar
  use driftorbit, only: nonlin, mth, mph, etatp, etadt, epsst_spl, mi, pertfile
  use do_magfie_mod, only: R0, s
  use do_magfie_mod, only: do_magfie_init
  use do_magfie_pert_mod, only: do_magfie_pert_init
  implicit none

contains

  subroutine run_bounce_nonlin_diag(arg_runname)
    character(*), intent(in) :: arg_runname
    logical :: file_exists
    real(real64) :: v, eta_min, eta_max
    integer :: i, j, npts
    real(real64), allocatable :: eta(:)
    real(real64), parameter :: ux_list(3) = [1.0_real64, 2.0_real64, 3.0_real64]
    real(real64), allocatable :: att_nonlin(:, :)
    real(real64) :: taub, bounceavg(nvar), Omth, dOmthdv, dOmthdeta
    real(real64) :: ux, Hmn2

    ! Initialize NEO-RT environment similar to neort:main
    runname = trim(arg_runname)
    call read_and_set_control(runname)
    call do_magfie_init("in_file")
    if (pertfile) call do_magfie_pert_init("in_file_pert")
    call init_profiles(R0)

    inquire(file="plasma.in", exist=file_exists)
    if (file_exists) call read_and_init_plasma_input("plasma.in", s)

    inquire(file="profile.in", exist=file_exists)
    if (file_exists) call read_and_init_profile_input("profile.in", s, R0, 1.0_real64, 1.0_real64)

    call init
    call check_magfie(runname)

    ! Sample attenuation over trapped pitch range for several ux values
    eta_min = etatp*(1.0_real64 + 1.0e-4_real64)
    eta_max = etatp + (etadt - etatp)*(1.0_real64 - epsst_spl)
    npts = 100
    allocate(eta(npts))
    allocate(att_nonlin(npts, size(ux_list)))
    do i = 1, npts
      eta(i) = eta_min + (eta_max-eta_min)*(real(i-1,real64)/real(npts-1,real64))
    end do

    do j = 1, size(ux_list)
      ux = ux_list(j)
      v = ux*vth
      do i = 1, npts
        call Om_th(v, eta(i), Omth, dOmthdv, dOmthdeta)
        taub = 2.0_real64*acos(-1.0_real64)/abs(Omth)
        call bounce_fast(v, eta(i), taub, bounceavg, timestep_transport)
        Hmn2 = (bounceavg(3)**2 + bounceavg(4)**2)*(mi*(v*v/2.0_real64))**2

        nonlin = .true.
        att_nonlin(i, j) = nonlinear_attenuation(ux, eta(i), bounceavg, Omth, dOmthdv, dOmthdeta, Hmn2)
      end do
    end do

    call figure()
    do j = 1, size(ux_list)
      select case (j)
      case (1)
        call plot(eta, att_nonlin(:, j), label="ux="//trim(adjustl(to_str(ux_list(j)))), linestyle='-')
      case (2)
        call plot(eta, att_nonlin(:, j), label="ux="//trim(adjustl(to_str(ux_list(j)))), linestyle='--')
      case (3)
        call plot(eta, att_nonlin(:, j), label="ux="//trim(adjustl(to_str(ux_list(j)))), linestyle='-.')
      end select
    end do
    call title("Nonlinear attenuation vs eta")
    call xlabel("eta")
    call ylabel("attenuation factor")
    call xlim(minval(eta), maxval(eta))
    call ylim(0.0_real64, 1.2_real64)
    call legend()
    call savefig(trim(arg_runname)//"_bounce_nonlin.png")

    ! Also write data to a text file with header
    call write_data_file(trim(arg_runname)//"_bounce_nonlin.txt", eta, att_nonlin, ux_list)

    ! Print simple stats
    do j = 1, size(ux_list)
      write(*,'(A,F5.2,A,1x,ES12.4,1x,ES12.4)') 'ux=', ux_list(j), ' min/max:', minval(att_nonlin(:,j)), maxval(att_nonlin(:,j))
    end do

  end subroutine run_bounce_nonlin_diag

  pure function to_str(x) result(s)
    real(real64), intent(in) :: x
    character(len=16) :: s
    write(s,'(F4.1)') x
    s = adjustl(s)
  end function to_str

  subroutine write_data_file(fname, eta, att, uxvals)
    character(*), intent(in) :: fname
    real(real64), intent(in) :: eta(:)
    real(real64), intent(in) :: att(:, :)
    real(real64), intent(in) :: uxvals(:)
    integer :: i, u

    open(newunit=u, file=fname, status='replace', action='write')
    write(u,'(A)') '# NEO-RT diagnostic: nonlinear attenuation vs eta'
    write(u,'(A,F10.6)') '# s = ', s
    write(u,'(A,I0)')    '# mth = ', mth
    write(u,'(A,ES12.5)') '# mph = ', mph
    write(u,'(A,*(F5.2,1X))') '# ux values = ', uxvals
    write(u,'(A)') '# columns: eta  att(ux='//trim(adjustl(to_str(uxvals(1))))//')  '&
                   //'att(ux='//trim(adjustl(to_str(uxvals(2))))//')  '&
                   //'att(ux='//trim(adjustl(to_str(uxvals(3))))//')'
    do i = 1, size(eta)
      write(u,'(ES14.6,3(1X,ES14.6))') eta(i), att(i,1), att(i,2), att(i,3)
    end do
    close(u)
  end subroutine write_data_file

end module diag_bounce_nonlin
