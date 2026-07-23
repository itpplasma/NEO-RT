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
        epst, epsp, epst_spl, epsp_spl
    use neort_freq, only: Om_ph, Om_th
    use neort_profiles, only: Om_tE
    use neort_resonance, only: driftorbit_coarse, driftorbit_root

    implicit none

    real(dp) :: eta_res(2), eta_mid, v
    real(dp) :: Omph, dOmphdv, dOmphdeta
    real(dp) :: Omth, dOmthdv, dOmthdeta
    real(dp) :: residual, residual_scale
    real(dp) :: roots(nlev, 3)
    integer :: nroots, trial_mth, trial_v
    logical :: converged, found_root

    call set_log_level(-1) ! silence the expected bracket-failure warning

    call neort_init("driftorbit.in", "in_file")
    call neort_prepare_splines("plasma.in", "profile.in")
    call neort_setup_at_s(0.5_dp)

    if (epst <= epst_spl .or. epsp <= epsp_spl) then
        error stop "transport roots need an explicit-orbit margin beyond the spline edge"
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
            call driftorbit_coarse(v, etatp * (1.0_dp + epst), &
                etadt * (1.0_dp - epst), roots, nroots)
            if (nroots > 0) then
                found_root = .true.
                exit
            end if
        end do
        if (found_root) exit
    end do
    if (.not. found_root) error stop "test fixture has no zero-rotation resonance"

    eta_res = driftorbit_root(v, 0.0_dp, roots(1, 1), roots(1, 2), converged)
    if (.not. converged) error stop "zero-rotation resonance must converge"
    if (eta_res(1) < roots(1, 1) .or. eta_res(1) > roots(1, 2)) then
        error stop "zero-rotation resonance escaped its bracket"
    end if

    call Om_ph(v, eta_res(1), Omph, dOmphdv, dOmphdeta)
    call Om_th(v, eta_res(1), Omth, dOmthdv, dOmthdeta)
    residual = mph * Omph + mth * Omth
    residual_scale = max(1.0_dp, abs(mph * Omph), abs(mth * Omth))
    if (abs(residual) > 1.0e-12_dp * residual_scale) then
        error stop "zero-rotation resonance residual is not converged"
    end if

    print *, "test_resonance_bracket PASSED"
end program test_resonance_bracket
