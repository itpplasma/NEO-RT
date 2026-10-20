program test_passing_drift_derivative
    use iso_fortran_env, only: dp => real64
    use do_magfie_mod, only: do_magfie_init, R0, s
    use driftorbit, only: efac, bfac, magdrift, magdrift_passing
    use neort, only: init, set_to_passing_region
    use neort_config, only: read_and_set_config
    use neort_freq, only: Om_th, d_Om_ds
    use neort_magfie, only: init_flux_surface_average
    use neort_profiles, only: init_profile_at_s, init_profiles, read_and_init_plasma_input, &
        read_and_init_profile_input, vth
    use util, only: pi

    implicit none

    real(dp), parameter :: RTOL = 1.0e-10_dp
    real(dp) :: eta, eta_min, eta_max, v, taub
    real(dp) :: d_no_drift, d_passing_off, scale, s0

    call read_and_set_config("driftorbit.in")
    magdrift = .true.
    magdrift_passing = 0

    call do_magfie_init("in_file")
    call init_profiles(R0)
    call read_and_init_plasma_input("plasma.in", s)
    call read_and_init_profile_input("profile.in", s, R0, efac, bfac)
    call init

    s0 = s
    call set_to_passing_region(eta_min, eta_max)
    eta = 0.5_dp * (eta_min + eta_max)
    v = 0.9_dp * vth

    magdrift = .false.
    d_no_drift = evaluate_d_omphds(s0, eta, v)

    magdrift = .true.
    magdrift_passing = 0
    d_passing_off = evaluate_d_omphds(s0, eta, v)
    scale = max(1.0_dp, abs(d_no_drift), abs(d_passing_off))

    if (abs(d_no_drift - d_passing_off) > RTOL * scale) then
        write (*, '(A,2ES16.8)') "passing drift-off differs from no-drift control: ", &
            d_passing_off, d_no_drift
        error stop "passing magnetic-drift derivative regression failed"
    end if

    write (*, '(A,2ES16.8)') "passing drift-off and no-drift control: ", &
        d_passing_off, d_no_drift

contains

    function evaluate_d_omphds(s_eval, eta_eval, v_eval) result(d_omphds)
        real(dp), intent(in) :: s_eval, eta_eval, v_eval
        real(dp) :: d_omphds
        real(dp) :: omth, domthdv, domthdeta, domthds

        s = s_eval
        call init_flux_surface_average(s)
        call init_profile_at_s(R0, efac, bfac)
        call Om_th(v_eval, eta_eval, omth, domthdv, domthdeta)
        taub = 2.0_dp * pi / abs(omth)
        call d_Om_ds(v_eval, eta_eval, taub, domthds, d_omphds)
    end function evaluate_d_omphds

end program test_passing_drift_derivative
