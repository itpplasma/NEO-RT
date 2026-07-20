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
    use driftorbit, only: mth, vth, etatp, etadt
    use neort_resonance, only: driftorbit_root

    implicit none

    real(dp) :: eta_res(2), eta_mid, v

    call set_log_level(-1)  ! silence the expected bracket-failure warning

    call neort_init("driftorbit.in", "in_file")
    call neort_prepare_splines("plasma.in", "profile.in")
    call neort_setup_at_s(0.5_dp)

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

    print *, "test_resonance_bracket PASSED"
end program test_resonance_bracket
