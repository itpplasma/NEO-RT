program test_separatrix_model
    ! Independent manufactured oracle for the singular-fit algebra.  The
    ! production orbit integrator supplies the sampled values; this test pins
    ! the fit and its physical-distance derivatives against the known formula.
    !
    ! The oracle is written out here from the endpoint form directly, without
    ! reusing the generated kernel: both orbit integrals are
    ! A log(1/x) + B + x (C log(1/x) + D), the bounce average is their ratio,
    ! and the derivatives are the quotient rule applied to that ratio.  The
    ! manufactured coefficients are chosen so that the separatrix limit
    ! num_a/tau_a is a nonzero number of order one.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
    use neort_separatrix, only: separatrix_model_t, fit_separatrix_model, &
        evaluate_separatrix_model
    implicit none

    real(dp), parameter :: tau_coeff(4) = [2.3_dp, 0.7_dp, -0.4_dp, 1.1_dp]
    real(dp), parameter :: num_coeff(4) = [-0.8_dp, 0.2_dp, 0.6_dp, -0.3_dp]
    real(dp), parameter :: separatrix_limit = num_coeff(1)/tau_coeff(1)
    real(dp), parameter :: sample_x(12) = [ &
        1.0e-6_dp, 2.0e-6_dp, 4.0e-6_dp, 8.0e-6_dp, 1.6e-5_dp, 3.2e-5_dp, &
        6.4e-5_dp, 1.28e-4_dp, 2.56e-4_dp, 5.12e-4_dp, 1.024e-3_dp, &
        2.048e-3_dp]
    real(dp) :: tau_values(size(sample_x)), drift_values(size(sample_x))
    real(dp) :: filtered_x(size(sample_x))
    real(dp) :: narrow_x(8), narrow_tau(8), narrow_drift(8)
    real(dp) :: adaptive_x(20), adaptive_tau(20), adaptive_drift(20)
    real(dp) :: tau, dtau_dx, omega, domega_dx, drift, ddrift_dx
    type(separatrix_model_t) :: model
    logical :: ok
    integer :: k

    do k = 1, size(sample_x)
        tau_values(k) = series_value(tau_coeff, sample_x(k))
        drift_values(k) = series_value(num_coeff, sample_x(k))/tau_values(k)
    end do

    ! A failed orbit sample must not poison the boundary fit when enough
    ! resolved samples remain.  The expected result is the same asymptotic
    ! model, with the invalid point omitted from the selected boundary window.
    filtered_x = sample_x
    filtered_x(1) = ieee_value(filtered_x(1), ieee_quiet_nan)
    call fit_separatrix_model(filtered_x, tau_values, drift_values, model, ok)
    if (.not. ok) error stop "finite separatrix samples were not selected"
    if (model%tau(1) <= 0.0_dp) error stop "filtered fit lost positive log coefficient"

    call fit_separatrix_model(sample_x, tau_values, drift_values, model, ok)
    if (.not. ok) error stop "manufactured separatrix fit was rejected"
    if (model%tau_relative_rms > 1.0e-10_dp .or. &
            model%numerator_relative_rms > 1.0e-10_dp) then
        error stop "manufactured separatrix fit residual is too large"
    end if

    ! The bounce average must reach its known finite separatrix limit, and it
    ! must do so from a distance far below any sampled one.  A model that
    ! diverges logarithmically or that truncates the 1/log series fails here.
    if (abs(model%numerator(1)/model%tau(1) - separatrix_limit) > &
            1.0e-9_dp*abs(separatrix_limit)) then
        error stop "fitted separatrix limit disagrees with num_a/tau_a"
    end if
    call evaluate_separatrix_model(model, 1.0e-300_dp, tau, dtau_dx, omega, &
        domega_dx, drift, ddrift_dx, ok)
    if (.not. ok) error stop "machine-adjacent separatrix evaluation was rejected"
    if (abs(drift - separatrix_limit) > 1.0e-2_dp*abs(separatrix_limit)) then
        error stop "bounce average does not approach its separatrix limit"
    end if

    ! A narrow resolved window cannot distinguish the x*log(x) correction
    ! columns from the leading logarithm, but it still resolves that leading
    ! pair.  The production fit must reduce its model order instead of
    ! rejecting a usable passing surface.
    do k = 1, size(narrow_x)
        narrow_x(k) = 1.0e-3_dp * (1.0_dp - real(size(narrow_x) - k, dp)*1.0e-5_dp)
        narrow_tau(k) = series_value([2.0_dp, 0.5_dp, 0.0_dp, 0.0_dp], narrow_x(k))
        narrow_drift(k) = 0.0_dp
    end do
    call fit_separatrix_model(narrow_x, narrow_tau, narrow_drift, model, ok)
    if (.not. ok) error stop "narrow separatrix window was rejected"
    if (abs(model%tau(3)) > 0.0_dp .or. abs(model%tau(4)) > 0.0_dp) then
        error stop "narrow separatrix window did not reduce model order"
    end if
    if (model%tau(1) <= 0.0_dp) error stop "narrow fit lost positive log coefficient"

    ! A vanishing drift integrand must stay identically zero rather than
    ! producing a spurious finite limit from a zero-over-zero ratio.
    call evaluate_separatrix_model(model, 1.0e-12_dp, tau, dtau_dx, omega, &
        domega_dx, drift, ddrift_dx, ok)
    if (.not. ok) error stop "narrow-window evaluation was rejected"
    if (drift /= 0.0_dp .or. ddrift_dx /= 0.0_dp) then
        error stop "zero drift samples did not give a zero bounce average"
    end if

    ! The closest samples can carry a reversed numerical trend.  The adaptive
    ! boundary-window search must move outward until the positive leading
    ! coefficient is resolved by the independent manufactured data.
    do k = 1, size(adaptive_x)
        adaptive_x(k) = 1.0e-6_dp * 2.0_dp**real(k - 1, dp)
        if (k <= 8) then
            adaptive_tau(k) = series_value([-1.0_dp, 0.5_dp, 0.0_dp, 0.0_dp], &
                adaptive_x(k))
        else
            adaptive_tau(k) = series_value([2.0_dp, 0.5_dp, 0.0_dp, 0.0_dp], &
                adaptive_x(k))
        end if
        adaptive_drift(k) = 0.0_dp
    end do
    call fit_separatrix_model(adaptive_x, adaptive_tau, adaptive_drift, model, ok)
    if (.not. ok .or. model%tau(1) <= 0.0_dp) then
        error stop "adaptive separatrix window did not recover positive trend"
    end if

    ! A rejected fit must not leave a window behind: callers gate on
    ! x <= xscale and would otherwise enter an invalid model.
    call fit_separatrix_model(sample_x(1:3), tau_values(1:3), drift_values(1:3), &
        model, ok)
    if (ok) error stop "under-determined separatrix fit was accepted"
    if (model%valid .or. model%xscale /= 0.0_dp) then
        error stop "rejected separatrix fit left a usable window behind"
    end if

    call fit_separatrix_model(sample_x, tau_values, drift_values, model, ok)
    if (.not. ok) error stop "manufactured separatrix fit was not reproducible"

    ! Probe the whole validity range, from the outer edge of the fitted
    ! window down to eight decades below the closest sampled orbit.
    do k = 0, 8
        call check_point(model, model%xscale*10.0_dp**real(-k, dp), &
            tau_coeff, num_coeff)
    end do
    write (*, '(a)') "PASS separatrix asymptotic fit"

contains

    subroutine check_point(model, x, tau_coeff, num_coeff)
        type(separatrix_model_t), intent(in) :: model
        real(dp), intent(in) :: x, tau_coeff(4), num_coeff(4)
        real(dp) :: tau, dtau_dx, omega, domega_dx, drift, ddrift_dx
        real(dp) :: want_tau, want_dtau, want_omega, want_domega
        real(dp) :: want_num, want_dnum, want_drift, want_ddrift
        logical :: good

        call evaluate_separatrix_model(model, x, tau, dtau_dx, omega, &
            domega_dx, drift, ddrift_dx, good)
        if (.not. good) error stop "separatrix evaluation was rejected"
        want_tau = series_value(tau_coeff, x)
        want_dtau = series_derivative(tau_coeff, x)
        want_omega = 2.0_dp*acos(-1.0_dp)/want_tau
        want_domega = -2.0_dp*acos(-1.0_dp)*want_dtau/want_tau**2
        want_num = series_value(num_coeff, x)
        want_dnum = series_derivative(num_coeff, x)
        want_drift = want_num/want_tau
        want_ddrift = want_dnum/want_tau - want_num*want_dtau/want_tau**2
        if (abs(tau - want_tau) > 2.0e-9_dp*max(1.0_dp, abs(want_tau))) &
            error stop "tau asymptotic value disagrees with oracle"
        if (abs(dtau_dx - want_dtau) > 2.0e-6_dp*max(1.0_dp, abs(want_dtau))) &
            error stop "tau asymptotic derivative disagrees with oracle"
        if (abs(omega - want_omega) > 2.0e-9_dp*max(1.0_dp, abs(want_omega))) &
            error stop "frequency asymptotic value disagrees with oracle"
        if (abs(domega_dx - want_domega) > 2.0e-6_dp*max(1.0_dp, abs(want_domega))) &
            error stop "frequency asymptotic derivative disagrees with oracle"
        if (abs(drift - want_drift) > 2.0e-9_dp*max(1.0_dp, abs(want_drift))) &
            error stop "drift asymptotic value disagrees with oracle"
        if (abs(ddrift_dx - want_ddrift) > 2.0e-6_dp*max(1.0_dp, abs(want_ddrift))) &
            error stop "drift asymptotic derivative disagrees with oracle"
    end subroutine check_point

    pure function series_value(coeff, x) result(value)
        real(dp), intent(in) :: coeff(4), x
        real(dp) :: value, log_distance

        log_distance = log(1.0_dp/x)
        value = coeff(1)*log_distance + coeff(2) + x * &
            (coeff(3)*log_distance + coeff(4))
    end function series_value

    pure function series_derivative(coeff, x) result(value)
        real(dp), intent(in) :: coeff(4), x
        real(dp) :: value, log_distance

        log_distance = log(1.0_dp/x)
        value = -coeff(1)/x + coeff(3)*log_distance - coeff(3) + coeff(4)
    end function series_derivative

end program test_separatrix_model
