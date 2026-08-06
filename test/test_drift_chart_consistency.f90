program test_drift_chart_consistency
    !! Oracle for the toroidal magnetic drift.
    !!
    !! NEO-RT carries two expressions for Om_tB/v**2, selected by the chart:
    !! an analytic one written with Boozer covariant components, and a
    !! geometric one built from curvature and grad-B and projected across the
    !! field-line label.  Om_tB is the entire toroidal precession of a trapped
    !! particle, so with mth=0 the resonance condition is Om_tE + Om_tB = 0 and
    !! any error in it silently removes the trapped contribution while leaving
    !! passing particles, which resonate through Omth/iota, almost untouched.
    !!
    !! The geometric expression is chart independent by construction, so on a
    !! Boozer equilibrium it has to reproduce the analytic one.  Comparing them
    !! there needs no new equilibrium and no reference data: each is a complete
    !! statement of the same physical quantity, so they check each other.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use do_magfie_mod, only: inp_swi, read_boozer_file, set_s, &
        init_magfie_at_s, do_magfie
    use neort_orbit, only: drift_toroidal_geometric, drift_toroidal_boozer
    use util, only: pi

    implicit none

    !> The two forms differ only by discretisation of the equilibrium splines.
    real(dp), parameter :: tol_relative = 1.0e-3_dp
    integer, parameter :: ntheta = 64
    real(dp), parameter :: s_probe(3) = [0.06_dp, 0.25_dp, 0.49_dp]
    real(dp), parameter :: eta_fraction(3) = [0.2_dp, 0.5_dp, 0.9_dp]

    character(len=1024) :: boozer_file
    real(dp) :: x(3), bmod, sqrtg, hder(3), hcovar(3), hctrvr(3), hcurl(3)
    real(dp) :: eta, geometric, analytic, scale, worst, worst_s, worst_theta
    integer :: is, k, ie

    call get_environment_variable('BOOZER_FILE', boozer_file)
    if (len_trim(boozer_file) == 0) then
        write (*, '(a)') 'test_drift_chart_consistency skipped: BOOZER_FILE unset'
        stop
    end if

    inp_swi = 9
    call read_boozer_file(trim(boozer_file))

    ! The comparison is only meaningful where both expressions have their
    ! inputs.  do_magfie leaves hcurl at zero for the Boozer and chartmap
    ! charts (marked TODO in do_magfie_standalone), and without curl(h) the
    ! geometric form loses its curvature term entirely, so it would fail here
    ! for want of an input rather than because it is wrong.
    call set_s(s_probe(1))
    call init_magfie_at_s()
    x = [s_probe(1), 0.0_dp, 0.3_dp]
    call do_magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
    if (maxval(abs(hcurl)) <= 0.0_dp) then
        write (*, '(a)') 'test_drift_chart_consistency skipped: this chart '// &
            'does not supply curl(h), so the geometric drift cannot be formed'
        stop
    end if

    worst = 0.0_dp
    worst_s = 0.0_dp
    worst_theta = 0.0_dp
    scale = 0.0_dp

    do is = 1, size(s_probe)
        call set_s(s_probe(is))
        call init_magfie_at_s()
        do k = 1, ntheta
            x(1) = s_probe(is)
            x(2) = 0.0_dp
            x(3) = 2.0_dp*pi*real(k - 1, dp)/real(ntheta, dp)
            call do_magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
            do ie = 1, size(eta_fraction)
                eta = eta_fraction(ie)/bmod
                geometric = drift_toroidal_geometric(eta, bmod, sqrtg, hder, &
                                                     hcovar, hctrvr, hcurl)
                analytic = drift_toroidal_boozer(eta, bmod, hder, hctrvr)
                scale = max(scale, abs(analytic))
                if (abs(geometric - analytic) > worst) then
                    worst = abs(geometric - analytic)
                    worst_s = s_probe(is)
                    worst_theta = x(3)
                end if
            end do
        end do
    end do

    write (*, '(a,es12.4,a,es12.4)') &
        'drift chart consistency: worst absolute ', worst, ' against scale ', scale
    write (*, '(a,f8.4,a,f8.4)') '   worst at s = ', worst_s, ', theta = ', worst_theta
    if (worst > tol_relative*scale) then
        write (*, '(a,es12.4,a,es12.4)') &
            'FAIL: geometric and Boozer magnetic drift disagree: ', &
            worst/max(scale, tiny(scale)), ' >', tol_relative
        error stop
    end if

end program test_drift_chart_consistency
