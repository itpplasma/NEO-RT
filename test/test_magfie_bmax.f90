program test_magfie_bmax
    ! etatp = 1/Bmax must bound the passing orbits. The golden-record field
    ! has its maximum of B between the theta grid samples; the grid maximum
    ! put etatp 5e-8 to 3e-7 above 1/max(B), so orbits classified as passing
    ! reflected. Compare with a dense theta scan on several surfaces.
    use iso_fortran_env, only: dp => real64
    use neort_lib, only: neort_init, neort_prepare_splines, neort_setup_at_s
    use do_magfie_mod, only: do_magfie, s
    use driftorbit, only: etatp
    implicit none
    integer, parameter :: samples = 200000
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: surfaces(3) = [0.1_dp, 0.5_dp, 0.9_dp]
    character(len=4096) :: config
    real(dp) :: x(3), bmod, sqrtg, hder(3), hcovar(3), hctrvr(3), hcurl(3)
    real(dp) :: scan_maximum, excess
    integer :: j, k, unit

    call get_command_argument(1, config)
    open(newunit=unit, file=trim(config), status='replace', action='write')
    write(unit, '(a)') '&params', ' s = 0.5', ' m_t = 0.036', ' qs = 1.0', &
        ' ms = 2.014', ' epsmn = 1.0', ' m0 = -1', ' mph = 1', ' magdrift = .true.', &
        ' nopassing = .false.', ' noshear = .false.', ' pertfile = .false.', &
        ' nonlin = .false.', ' bfac = 1.0', ' efac = 1.0', ' inp_swi = 9', &
        ' vsteps = 16', ' log_level = -1', '/'
    close(unit)
    call neort_init(trim(config), 'in_file')
    call neort_prepare_splines('plasma.in', 'profile.in')
    do j = 1, size(surfaces)
        call neort_setup_at_s(surfaces(j))
        scan_maximum = 0.0_dp
        do k = 0, samples - 1
            x = [s, 0.0_dp, 2.0_dp*pi*real(k, dp)/real(samples, dp)]
            call do_magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
            scan_maximum = max(scan_maximum, bmod)
        end do
        excess = etatp*scan_maximum - 1.0_dp
        print *, 's, etatp*max(B) - 1:', surfaces(j), excess
        if (excess > 1.0e-13_dp) error stop 'etatp exceeds 1/max(B): passing orbits reflect'
        if (excess < -1.0e-9_dp) error stop 'etatp below 1/max(B) beyond scan resolution'
    end do
    print *, 'test_magfie_bmax PASSED'
end program test_magfie_bmax
