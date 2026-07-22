program probe_neort_eqdsk_drift
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use do_magfie_mod, only: inp_swi, bfac, read_boozer_file, set_s, &
        init_magfie_at_s, do_magfie, q
    use neort_orbit, only: timestep
    use neort_profiles, only: Om_tE
    use util, only: mi, qi

    implicit none

    integer, parameter :: nsample = 8
    character(len=1024) :: eqdsk_file
    real(dp), parameter :: sample_s(nsample) = &
        [0.10_dp, 0.10_dp, 0.35_dp, 0.35_dp, 0.65_dp, 0.65_dp, 0.90_dp, 0.90_dp]
    real(dp), parameter :: sample_phi(nsample) = &
        [0.00_dp, 0.71_dp, 1.43_dp, 2.14_dp, 2.86_dp, 3.57_dp, 4.29_dp, 5.00_dp]
    real(dp), parameter :: sample_theta(nsample) = &
        [0.00_dp, 0.79_dp, 1.57_dp, 2.36_dp, 3.14_dp, 3.93_dp, 4.71_dp, 5.50_dp]
    real(dp), parameter :: pitch_b = 0.30_dp
    real(dp), parameter :: v_probe = 1.0_dp

    real(dp) :: x(3), y(7), ydot_plus(7), ydot_minus(7), ydot_eplus(7), ydot_eminus(7)
    real(dp) :: bmod, sqrtg, bder(3), hcov(3), hcon(3), hcurl(3)
    real(dp) :: eta, mi_saved, qi_saved, omte_saved
    integer :: k, status

    call get_command_argument(1, eqdsk_file, status=status)
    if (status /= 0 .or. len_trim(eqdsk_file) == 0) then
        error stop 'usage: probe_neort_eqdsk_drift.x EQDSK_FILE'
    end if

    inp_swi = 11
    bfac = 1.0_dp
    call read_boozer_file(trim(eqdsk_file))

    mi_saved = mi
    qi_saved = qi
    omte_saved = Om_tE
    mi = 1.0_dp
    do k = 1, nsample
        x = [sample_s(k), sample_phi(k), sample_theta(k)]
        call set_s(x(1))
        call init_magfie_at_s()
        call do_magfie(x, bmod, sqrtg, bder, hcov, hcon, hcurl)

        eta = pitch_b/bmod
        y = 0.0_dp
        y(1) = x(3)
        y(2) = v_probe*sqrt(1.0_dp - pitch_b)

        qi = 1.0_dp
        Om_tE = 0.0_dp
        call timestep(v_probe, eta, size(y), 0.0_dp, y, ydot_plus)
        qi = -1.0_dp
        call timestep(v_probe, eta, size(y), 0.0_dp, y, ydot_minus)
        qi = 1.0_dp
        Om_tE = 1.0e9_dp
        call timestep(v_probe, eta, size(y), 0.0_dp, y, ydot_eplus)
        Om_tE = -1.0e9_dp
        call timestep(v_probe, eta, size(y), 0.0_dp, y, ydot_eminus)

        write(*, '(a,1x,i0,1x,*(es24.16e3,1x))') 'NEORT_DRIFT_SAMPLE', k, &
            x, bmod, sqrtg, eta, q, bder, hcov, hcon, hcurl, &
            ydot_plus(3), ydot_minus(3), ydot_eplus(3), ydot_eminus(3)
    end do
    mi = mi_saved
    qi = qi_saved
    Om_tE = omte_saved

end program probe_neort_eqdsk_drift
