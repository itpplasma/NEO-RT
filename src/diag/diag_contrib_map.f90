module diag_contrib_map
  use iso_fortran_env, only: real64
  use fortplot, only: figure, plot, pcolormesh, title, xlabel, ylabel, legend, savefig
  use neort, only: read_control, init, check_magfie, runname => runname, set_to_passing_region, set_to_trapped_region
  use neort_profiles, only: init_profile_input, init_plasma_input, init_profiles, vth, Om_tE
  use neort_nonlin, only: nonlinear_attenuation
  use neort_freq, only: Om_th
  use neort_transport, only: timestep_transport, Tphi_int
  use neort_orbit, only: bounce_fast, nvar
  use neort_resonance, only: driftorbit_coarse, driftorbit_root
  use driftorbit, only: nonlin, mth, mph, etatp, etadt, epsst_spl, mi, pertfile, etamin, etamax, sign_vpar, sign_vpar_htheta
  use do_magfie_mod, only: R0, s, q
  use do_magfie_mod, only: do_magfie_init
  use do_magfie_pert_mod, only: do_magfie_pert_init
  implicit none

contains

  subroutine run_contrib_diag(arg_runname)
    character(*), intent(in) :: arg_runname
    logical :: file_exists
    integer :: i, j, nu, nm, u
    real(real64) :: v
    real(real64), allocatable :: ux(:), ux_edges(:), mth_vals(:), mth_edges(:)
    real(real64), allocatable :: z(:, :), y(:)
    integer :: mth_min, mth_max

    ! Temp variables for resonance calculation
    real(real64) :: Omth, dOmthdv, dOmthdeta, taub, bounceavg(nvar)
    real(real64) :: Hmn2, eta, du
    real(real64) :: roots(100, 3)
    integer :: nroots, kr
    real(real64) :: eta_res(2)
    real(real64) :: contrib

    ! Initialize environment
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

    ! Grids
    nu = 400
    mth_min = -ceiling(2.0_real64*abs(mph*q))
    mth_max =  ceiling(2.0_real64*abs(mph*q))
    nm = mth_max - mth_min + 1
    allocate(ux(nu), ux_edges(nu+1))
    allocate(mth_vals(nm), mth_edges(nm+1))
    allocate(z(nm, nu))
    allocate(y(nu))
    z = 0.0_real64
    y = 0.0_real64

    do i = 1, nu
      ux(i) = 0.3_real64 + (2.5_real64 - 0.3_real64) * real(i-1,real64) / real(nu-1,real64)
    end do
    do j = 1, nm
      mth_vals(j) = real(mth_min + (j-1), real64)
    end do
    call edges_from_centers(ux, ux_edges)
    call edges_from_centers(mth_vals, mth_edges)

    ! Compute contribution density vs ux and mth (torque integrand)
    ! Force linear contributions: temporarily set nonlin = .false.
    nonlin = .false.
    do j = 1, nm
      mth = nint(mth_vals(j))
      do i = 1, nu
        v = ux(i)*vth
        contrib = 0.0_real64

        ! Passing (co- and counter-)
        call accumulate_class(v, ux(i), +1.0_real64, .true., contrib)
        call accumulate_class(v, ux(i), -1.0_real64, .true., contrib)

        ! Trapped
        call accumulate_class(v, ux(i), +1.0_real64, .false., contrib)

        z(j, i) = contrib
        y(i) = y(i) + contrib
      end do
    end do
    ! Note: leave nonlin as originally set for subsequent diagnostics

    ! Plot sum over mth vs ux
    call figure()
    call plot(ux, y, label='sum over mth')
    call title('Torque integrand vs ux (nonlinear)')
    call xlabel('ux')
    call ylabel('dT/dx (arbitrary units)')
    call legend()
    call savefig(trim(arg_runname)//'_contrib_ux.png')

    ! Heatmap over mth vs ux
    call figure()
    call pcolormesh(ux_edges, mth_edges, z, colormap='plasma')
    call title('Torque integrand map (nonlinear)')
    call xlabel('ux')
    call ylabel('mth')
    call savefig(trim(arg_runname)//'_contrib_map.png')

    ! Dumps
    open(newunit=u, file=trim(arg_runname)//'_contrib_ux.txt', status='replace', action='write')
    write(u,'(A)') '# torque integrand density vs ux, summed over mth'
    write(u,'(A,F10.6)') '# s = ', s
    write(u,'(A,*(1X,ES12.5))') '# ux =', ux
    write(u,'(*(ES12.5,1X))') y
    close(u)

    open(newunit=u, file=trim(arg_runname)//'_contrib_map.txt', status='replace', action='write')
    write(u,'(A)') '# torque integrand map (rows=mth, cols=ux)'
    write(u,'(A,F10.6)') '# s = ', s
    write(u,'(A,I0,A,I0)') '# mth_min=', mth_min, ' mth_max=', mth_max
    write(u,'(A,*(1X,ES12.5))') '# ux =', ux
    do j = 1, nm
      write(u,'(*(ES12.5,1X))') z(j, :)
    end do
    close(u)

  contains

    subroutine accumulate_contrib(v, ux, contrib_out)
      real(real64), intent(in) :: v, ux
      real(real64), intent(inout) :: contrib_out
      real(real64) :: Omth, dOmthdv, dOmthdeta, taub, bounceavg(nvar), Hmn2
      real(real64) :: eta
      real(real64) :: eta_res(2)
      real(real64) :: roots(100, 3)
      integer :: nroots, kr
      real(real64) :: att

      call driftorbit_coarse(v, etamin, etamax, roots, nroots)
      if (nroots == 0) return
      do kr = 1, nroots
        eta_res = driftorbit_root(v, 1d-8*abs(Om_tE), roots(kr, 1), roots(kr, 2))
        eta = eta_res(1)

        call Om_th(v, eta, Omth, dOmthdv, dOmthdeta)
        taub = 2.0_real64*acos(-1.0_real64)/abs(Omth)
        call bounce_fast(v, eta, taub, bounceavg, timestep_transport)
        Hmn2 = (bounceavg(3)**2 + bounceavg(4)**2)*(mi*(v*v/2.0_real64))**2
        att = nonlinear_attenuation(ux, eta, bounceavg, Omth, dOmthdv, dOmthdeta, Hmn2)
        contrib_out = contrib_out + Tphi_int(ux, taub, Hmn2)/abs(eta_res(2)) * att
      end do
    end subroutine accumulate_contrib

    subroutine accumulate_class(v, ux, sign_vpar_in, is_passing, contrib_out)
      real(real64), intent(in) :: v, ux, sign_vpar_in
      logical, intent(in) :: is_passing
      real(real64), intent(inout) :: contrib_out
      sign_vpar = sign_vpar_in
      if (is_passing) then
        call set_to_passing_region(etamin, etamax)
      else
        call set_to_trapped_region(etamin, etamax)
      end if
      call accumulate_contrib(v, ux, contrib_out)
    end subroutine accumulate_class

  end subroutine run_contrib_diag

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

end module diag_contrib_map
