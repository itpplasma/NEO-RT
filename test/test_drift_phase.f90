program test_drift_phase
    ! The orbit Fourier coefficient of the perturbation with the oscillating
    ! drift phase mph*int Om_osc dt, Om_osc = A cos(w t), must equal the
    ! Jacobi-Anger sum  H(mth) = sum_k J_k(beta) H0(mth - k)  of the
    ! unmodulated coefficients H0, beta = mph*A/w.  Passing and trapped orbits.
    ! Also checks that the real local drift integrates to ~zero over one orbit.
    use iso_fortran_env, only: dp => real64
    use do_magfie_mod, only: do_magfie_init, s
    use do_magfie_pert_mod, only: do_magfie_pert_init
    use neort, only: init_profiles, init, set_to_trapped_region
    use neort_config, only: read_and_set_config
    use neort_profiles, only: read_and_init_plasma_input, read_and_init_profile_input
    use neort_orbit, only: nvar, bounce_fast
    use neort_freq, only: Om_th
    use neort_transport, only: timestep_transport, drift_osc_override, Omth
    use util, only: pi
    use driftorbit

    implicit none

    integer, parameter :: kmax = 20
    real(dp), parameter :: beta = 1.5_dp
    real(dp) :: amp, w, d1, d2, v
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
    drift_phase = .true.
    drift_osc_override => osc

    call run_case("passing", 0.5_dp*etatp)
    call run_case("trapped", 0.5_dp*(etamin + etamax))
    call run_real_field("passing", 0.5_dp*etatp)
    call run_real_field("trapped", 0.5_dp*(etamin + etamax))

    if (failures > 0) error stop "test_drift_phase failed"
    print *, "test_drift_phase passed"

contains

    function osc(t) result(Om_osc)
        real(dp), intent(in) :: t
        real(dp) :: Om_osc
        Om_osc = amp*cos(w*t)
    end function osc

    function jk(k, x) result(j)
        ! Bessel function of integer order k, J_-k = (-1)^k J_k
        integer, intent(in) :: k
        real(dp), intent(in) :: x
        real(dp) :: j
        j = bessel_jn(abs(k), x)
        if (k < 0 .and. mod(k, 2) /= 0) j = -j
    end function jk

    function coeff(eta, taub, mth_) result(H)
        real(dp), intent(in) :: eta, taub
        integer, intent(in) :: mth_
        complex(dp) :: H
        real(dp) :: ba(nvar)
        mth = mth_
        call bounce_fast(v, eta, taub, ba, timestep_transport)
        H = cmplx(ba(3), ba(4), dp)
    end function coeff

    subroutine run_case(label, eta)
        character(*), intent(in) :: label
        real(dp), intent(in) :: eta
        integer :: mths(4), i, k
        complex(dp) :: H, Href
        real(dp) :: taub, scale, err, diff0

        call Om_th(v, eta, Omth, d1, d2)
        taub = 2.0_dp*pi/abs(Omth)
        w = Omth
        mths = [0, 1, -1, 2]
        do i = 1, size(mths)
            amp = beta*w/mph
            H = coeff(eta, taub, mths(i))
            Href = (0.0_dp, 0.0_dp)
            amp = 0.0_dp
            do k = -kmax, kmax
                Href = Href + jk(k, beta)*coeff(eta, taub, mths(i) - k)
            end do
            scale = abs(coeff(eta, taub, 0))
            err = abs(H - Href)/scale
            diff0 = abs(H - coeff(eta, taub, mths(i)))/scale
            print '(a,a,i3,a,es10.3,a,es10.3)', label, " mth=", mths(i), &
                " err=", err, " modulation effect=", diff0
            if (err > 1.0e-4_dp) failures = failures + 1
            if (diff0 < 1.0e-2_dp .and. mths(i) == 0) failures = failures + 1
        end do
    end subroutine run_case

    subroutine run_real_field(label, eta)
        ! The local drift minus its bounce average must integrate to zero
        ! over a closed orbit: y(7)/taub small against the bounce frequency.
        character(*), intent(in) :: label
        real(dp), intent(in) :: eta
        real(dp) :: taub, ba(nvar)

        drift_osc_override => null()
        call Om_th(v, eta, Omth, d1, d2)
        taub = 2.0_dp*pi/abs(Omth)
        mth = 0
        call bounce_fast(v, eta, taub, ba, timestep_transport)
        print '(a,a,es10.3,a,es10.3)', label, " real-field mean osc drift =", &
            ba(7), " |Omth|=", abs(Omth)
        if (abs(ba(7)) > 1.0e-3_dp*abs(Omth)/mph) failures = failures + 1
        drift_osc_override => osc
    end subroutine run_real_field

end program test_drift_phase
