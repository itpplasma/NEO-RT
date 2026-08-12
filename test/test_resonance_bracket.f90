program test_resonance_bracket
    ! Regression test for the driftorbit_root bracket-failure sentinel.
    !
    ! When driftorbit_root is handed an eta bracket that contains no
    ! resonance (driftorbit_nroot == 0), it must return a *defined* result
    ! rather than leaving driftorbit_root(1:2) unassigned. The contract is a
    ! negative sentinel in eta_res(1) (eta is a non-negative pitch parameter)
    ! and a defined eta_res(2). The production caller only ever passes
    ! sign-change brackets, so this path is inactive for production inputs;
    ! this test pins the defined failure behaviour directly.
    use iso_fortran_env, only: dp => real64
    use logger, only: set_log_level
    use neort_lib, only: neort_init, neort_prepare_splines, neort_setup_at_s
    use driftorbit, only: mth, mph, nlev, sign_vpar, vth, etatp, etadt, &
        etamin, etamax
    use neort_freq, only: Om_ph, Om_th
    use neort_profiles, only: Om_tE
    use neort_resonance, only: driftorbit_coarse, driftorbit_root, solve_bracketed_root
    use neort, only: set_to_passing_region, set_to_trapped_region

    implicit none

    real(dp) :: eta_res(2), eta_mid, v, eta_min, eta_max
    real(dp) :: Omph, dOmphdv, dOmphdeta
    real(dp) :: Omth, dOmthdv, dOmthdeta
    real(dp) :: residual, residual_scale, residual_tol, dresdeta
    real(dp) :: roots(nlev, 3)
    integer :: nroots, trial_mth, trial_v, status
    real(dp) :: analytic_eta, analytic_derivative
    logical :: found_root

    call set_log_level(-1)  ! silence the expected bracket-failure warning

    call neort_init("driftorbit.in", "in_file")
    call neort_prepare_splines("plasma.in", "profile.in")
    call neort_setup_at_s(0.5_dp)

    call set_to_trapped_region(eta_min, eta_max)
    if (.not. (eta_min > etatp .and. eta_min < eta_max)) then
        error stop "trapped region does not reach the open boundary"
    end if
    if ((eta_min - etatp)/etatp > 8.0_dp*epsilon(1.0_dp)) then
        error stop "trapped boundary retains a physical exclusion layer"
    end if
    call set_to_passing_region(eta_min, eta_max)
    if (.not. (eta_min < eta_max .and. eta_max < etatp)) then
        error stop "passing region does not reach the open boundary"
    end if
    if ((etatp - eta_max)/etatp > 8.0_dp*epsilon(1.0_dp)) then
        error stop "passing boundary retains a physical exclusion layer"
    end if

    call solve_bracketed_root(0.0_dp, 1.0e-14_dp, 0.0_dp, 1.0_dp, &
        analytic_residual, analytic_eta, analytic_derivative, status)
    if (status /= 1 .or. abs(analytic_eta - 0.375_dp) > 1.0e-14_dp .or. &
        abs(analytic_derivative - 1.0_dp) > 1.0e-14_dp) then
        error stop "analytic safeguarded root oracle failed"
    end if

    mth = 1
    v = vth

    ! Degenerate (zero-width) bracket inside the trapped domain: no sign
    ! change, so driftorbit_nroot returns 0 and the failure path is taken.
    eta_mid = 0.5_dp*(etatp + etadt)
    eta_res = driftorbit_root(v, 1.0e-8_dp, eta_mid, eta_mid)

    if (.not. (eta_res(1) < 0.0_dp)) then
        error stop "bracket failure must return a negative eta sentinel"
    end if
    if (eta_res(2) /= 0.0_dp) then
        error stop "bracket failure must return a defined eta_res(2)"
    end if

    ! A zero electric-rotation lane previously passed tol=0 to the bisection
    ! solver. The strict abs(res)<tol predicate could then never succeed, and
    ! the last iterate was returned as if it were a usable resonance. Locate
    ! an actual magnetic-drift resonance, solve it with the production zero
    ! tolerance, and verify the independently recomputed native residual.
    Om_tE = 0.0_dp
    sign_vpar = 1
    found_root = .false.
    do trial_v = 1, 24
        v = (0.25_dp * trial_v) * vth
        do trial_mth = -12, 12
            if (trial_mth == 0) cycle
            mth = trial_mth
            call set_to_trapped_region(etamin, etamax)
            call driftorbit_coarse(v, etamin, etamax, roots, nroots)
            if (nroots > 0) then
                found_root = .true.
                exit
            end if
        end do
        if (found_root) exit
    end do
    if (.not. found_root) error stop "test fixture has no zero-rotation resonance"

    eta_res = driftorbit_root(v, 0.0_dp, roots(1, 1), roots(1, 2))
    if (eta_res(1) < roots(1, 1) .or. eta_res(1) > roots(1, 2)) then
        error stop "zero-rotation resonance escaped its bracket"
    end if

    call Om_ph(v, eta_res(1), Omph, dOmphdv, dOmphdeta)
    call Om_th(v, eta_res(1), Omth, dOmthdv, dOmthdeta)
    residual = mph * Omph + mth * Omth
    dresdeta = mph*dOmphdeta + mth*dOmthdeta
    residual_scale = max(1.0_dp, abs(mph * Omph), abs(mth * Omth))
    residual_tol = max(1.0e-9_dp * residual_scale, &
        64.0_dp*epsilon(1.0_dp)*abs(dresdeta)*max(abs(eta_res(1)), tiny(1.0_dp)))
    if (abs(residual) > residual_tol) then
        write(*,*) "zero-rotation residual and tolerance:", residual, residual_tol
        error stop "zero-rotation resonance residual is not converged"
    end if

    print *, "test_resonance_bracket PASSED"

contains

    subroutine analytic_residual(v_in, eta_in, residual_out, derivative_out)
        real(dp), intent(in) :: v_in, eta_in
        real(dp), intent(out) :: residual_out, derivative_out

        residual_out = eta_in - 0.375_dp
        derivative_out = 1.0_dp
    end subroutine analytic_residual

end program test_resonance_bracket
