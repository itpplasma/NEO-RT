program test_fow_radial
    ! First-order radial orbit width (fow_radial).
    ! 1. Size: the displacement at theta = 0 matches the large-aspect-ratio
    !    poloidal-gyroradius estimate Delta s = 2 q R0 v_par / (a**2 omega_c)
    !    (Delta r = q R0 v_par/(r omega_c), s = r**2/a**2) to O(eps).
    ! 2. Labelling: the displacement averages to zero over trapped and passing
    !    orbits, so s stays the orbit-averaged radius.
    use iso_fortran_env, only: dp => real64
    use do_magfie_mod, only: do_magfie_init, do_magfie, s, a, R0, q
    use do_magfie_pert_mod, only: do_magfie_pert_init
    use neort, only: init_profiles, init, set_to_trapped_region
    use neort_config, only: read_and_set_config
    use neort_profiles, only: read_and_init_plasma_input, read_and_init_profile_input
    use neort_orbit, only: nvar, bounce_fast, bounce_time
    use neort_freq, only: Om_th
    use neort_transport, only: timestep_transport, fow_displacement, Omth
    use util, only: pi, c, qi
    use driftorbit

    implicit none

    real(dp) :: v, d1, d2
    integer :: failures

    failures = 0
    call read_and_set_config("driftorbit.in")
    call do_magfie_init("in_file")
    call init_profiles(R0)
    call read_and_init_plasma_input("plasma.in", s)
    call read_and_init_profile_input("profile.in", s, R0, efac, bfac)
    call init
    sign_vpar = 1
    mth = 0
    call set_to_trapped_region(etamin, etamax)
    v = 0.9_dp*vth

    call check_size(0.5_dp*(etamin + etamax))
    call check_mean("passing", 0.5_dp*etatp)
    call check_mean("trapped", 0.5_dp*(etamin + etamax))
    sign_vpar = -1
    call check_mean("passing counter", 0.5_dp*etatp)
    sign_vpar = 1

    if (failures > 0) error stop "test_fow_radial failed"
    print *, "test_fow_radial passed"

contains

    subroutine check_size(eta)
        real(dp), intent(in) :: eta
        real(dp) :: x(3), bmod, sqrtg, hder(3), hcovar(3), hctrvr(3), hcurl(3)
        real(dp) :: v_par, ds, ds_ref, omc

        x = [s, 0.0_dp, 0.0_dp]
        call do_magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
        v_par = v*sqrt(1.0_dp - eta*bmod)
        omc = qi*bmod/(mi*c)
        ds = fow_displacement(v_par, hcovar(2), eta)
        ds_ref = 2.0_dp*q*R0*v_par/(a**2*omc)
        print '(a,es12.4,a,es12.4,a,f8.4)', "trapped theta=0: Delta s =", ds, &
            "  large-aspect =", ds_ref, "  |ratio| =", abs(ds/ds_ref)
        if (abs(abs(ds/ds_ref) - 1.0_dp) > 0.25_dp) failures = failures + 1
    end subroutine check_size

    subroutine check_mean(label, eta)
        character(*), intent(in) :: label
        real(dp), intent(in) :: eta
        real(dp) :: taub, ba(nvar), ba_rms(nvar)

        ! Exact orbit period (bounce_time locates co-directed orbits only), so
        ! that the orbit closes and the mean formula is checked without the
        ! frequency-spline error
        call Om_th(v, eta, Omth, d1, d2)
        taub = 2.0_dp*pi/abs(Omth)
        if (sign_vpar > 0) then
            taub = bounce_time(v, eta, taub)
            Omth = sign(2.0_dp*pi/taub, Omth)
        end if
        call bounce_fast(v, eta, taub, ba, ts_mean)
        call bounce_fast(v, eta, taub, ba_rms, ts_rms)
        print '(a,a,es12.4,a,es12.4)', label, ": <Delta s> =", ba(3), &
            "  rms Delta s =", sqrt(ba_rms(3))
        if (.not. (sqrt(ba_rms(3)) > 0.0_dp)) failures = failures + 1
        if (abs(ba(3)) > 1.0e-3_dp*sqrt(ba_rms(3))) failures = failures + 1
    end subroutine check_mean

    subroutine ts_mean(v_, eta_, neq, t, y, ydot)
        real(dp), intent(in) :: v_, eta_
        integer, intent(in) :: neq
        real(dp), intent(in) :: t
        real(dp), intent(in) :: y(neq)
        real(dp), intent(out) :: ydot(neq)

        call timestep_transport(v_, eta_, neq, t, y, ydot)
        ydot(3) = fow_displacement(y(2), hphi_at(y(1)), eta_)
        ydot(4:) = 0.0_dp
    end subroutine ts_mean

    subroutine ts_rms(v_, eta_, neq, t, y, ydot)
        real(dp), intent(in) :: v_, eta_
        integer, intent(in) :: neq
        real(dp), intent(in) :: t
        real(dp), intent(in) :: y(neq)
        real(dp), intent(out) :: ydot(neq)

        call timestep_transport(v_, eta_, neq, t, y, ydot)
        ydot(3) = fow_displacement(y(2), hphi_at(y(1)), eta_)**2
        ydot(4:) = 0.0_dp
    end subroutine ts_rms

    function hphi_at(theta) result(hphi)
        real(dp), intent(in) :: theta
        real(dp) :: hphi
        real(dp) :: x(3), bmod, sqrtg, hder(3), hcovar(3), hctrvr(3), hcurl(3)

        x = [s, 0.0_dp, theta]
        call do_magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
        hphi = hcovar(2)
    end function hphi_at

end program test_fow_radial
