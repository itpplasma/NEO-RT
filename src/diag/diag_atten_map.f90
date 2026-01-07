module diag_atten_map
  use iso_fortran_env, only: real64
  use fortplot, only: figure, pcolormesh, title, xlabel, ylabel, savefig
  use neort, only: read_and_set_control, init, check_magfie, write_magfie_data_to_files, runname
  use neort_datatypes, only: magfie_data_t
  use neort_profiles, only: read_and_init_profile_input, read_and_init_plasma_input, init_profiles, vth
  use neort_nonlin, only: nonlinear_attenuation
  use neort_freq, only: Om_th
  use neort_transport, only: timestep_transport
  use neort_orbit, only: bounce_fast, nvar
  use driftorbit, only: nonlin, mth, etatp, etadt, epsst_spl, mi, pertfile
  use do_magfie_mod, only: R0, s, q
  use do_magfie_mod, only: do_magfie_init
  use do_magfie_pert_mod, only: do_magfie_pert_init, mph
  implicit none

contains

  subroutine run_atten_map_diag(arg_runname)
    character(*), intent(in) :: arg_runname
    logical :: file_exists
    integer :: i, j, nu, nm
    real(real64) :: eta_max
    real(real64), allocatable :: ux(:), ux_edges(:), mth_vals(:), mth_edges(:)
    real(real64), allocatable :: z(:, :)
    real(real64) :: v, Omth, dOmthdv, dOmthdeta, taub, bounceavg(nvar)
    real(real64) :: Hmn2
    integer :: mth_min, mth_max
    integer :: u
    type(magfie_data_t) :: magfie_data

    ! Initialize like main
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
    call check_magfie(magfie_data)
    call write_magfie_data_to_files(magfie_data, runname)

    ! Max trapped eta value consistent with earlier diagnostic
    eta_max = etatp + (etadt - etatp)*(1.0_real64 - epsst_spl)

    ! Grids
    nu = 80
    nm = 0
    mth_min = -ceiling(2.0_real64*abs(mph*q))
    mth_max =  ceiling(2.0_real64*abs(mph*q))
    nm = mth_max - mth_min + 1
    allocate(ux(nu), ux_edges(nu+1))
    allocate(mth_vals(nm), mth_edges(nm+1))
    allocate(z(nm, nu))

    do i = 1, nu
      ux(i) = 0.2_real64 + (3.0_real64 - 0.2_real64) * real(i-1,real64) / real(nu-1,real64)
    end do
    do j = 1, nm
      mth_vals(j) = real(mth_min + (j-1), real64)
    end do

    ! Compute cell edges for pcolormesh
    call edges_from_centers(ux, ux_edges)
    call edges_from_centers(mth_vals, mth_edges)

    ! Evaluate attenuation at eta_max for each (mth, ux)
    nonlin = .true.
    do j = 1, nm
      mth = nint(mth_vals(j))
      do i = 1, nu
        v = ux(i)*vth
        call Om_th(v, eta_max, Omth, dOmthdv, dOmthdeta)
        taub = 2.0_real64*acos(-1.0_real64)/abs(Omth)
        call bounce_fast(v, eta_max, taub, bounceavg, timestep_transport)
        Hmn2 = (bounceavg(3)**2 + bounceavg(4)**2)*(mi*(v*v/2.0_real64))**2
        z(j, i) = nonlinear_attenuation(ux(i), eta_max, bounceavg, Omth, dOmthdv, dOmthdeta, Hmn2)
      end do
    end do

    call figure()
    call pcolormesh(ux_edges, mth_edges, z, colormap='viridis')
    call title('Attenuation factor at eta_max')
    call xlabel('ux')
    call ylabel('mth')
    call savefig(trim(arg_runname)//'_atten_map.png')

    ! Write a text dump for further analysis
    open(newunit=u, file=trim(arg_runname)//'_atten_map.txt', status='replace', action='write')
    write(u,'(A)') '# NEO-RT diagnostic: attenuation map at eta_max'
    write(u,'(A,F10.6)') '# s = ', s
    write(u,'(A,ES12.5)') '# eta_max = ', eta_max
    write(u,'(A,ES12.5)') '# mph = ', mph
    write(u,'(A,ES12.5)') '# q = ', q
    write(u,'(A)') '# rows: mth from mth_min to mth_max; columns: ux grid values'
    write(u,'(A,I0,A,I0)') '# mth_min=', mth_min, ' mth_max=', mth_max
    write(u,'(A,*(1X,ES12.5))') '# ux =', ux
    do j = 1, nm
      write(u,'(*(ES12.5,1X))') z(j, :)
    end do
    close(u)

  end subroutine run_atten_map_diag

  subroutine edges_from_centers(c, e)
    real(real64), intent(in) :: c(:)
    real(real64), intent(out) :: e(:)
    integer :: n, i
    n = size(c)
    if (size(e) /= n+1) stop 'edges_from_centers: size mismatch'
    e(1) = c(1) - 0.5_real64*(c(2)-c(1))
    do i = 2, n
      e(i) = 0.5_real64*(c(i-1)+c(i))
    end do
    e(n+1) = c(n) + 0.5_real64*(c(n)-c(n-1))
  end subroutine edges_from_centers

end module diag_atten_map
