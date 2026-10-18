program test_omph_shaing_units
    ! Shaing's analytic bounce-averaged toroidal drift frequency (omph_shaing,
    ! Shaing 2009 PPCF 51 035009 Eq. 8) against NEO-RT's numerically bounce-
    ! averaged magnetic precession (Om_tB, the supban = .false. path) on
    ! examples/circ.bc.  Shaing's expression is the large-aspect-ratio,
    ! circular limit of the numerical one, so the two must agree to O(eps).
    ! The radial derivative must be taken with respect to s and divided by the
    ! poloidal-flux derivative iota*psi_pr = psi_pr/q (chi' of Eq. 8); using
    ! deps/dr [1/cm] with psi_pr = dpsi_tor/ds is too small by q*dr/ds [cm]
    ! (~2e2 at s = 0.09).
    use iso_fortran_env, only: dp => real64
    use do_magfie_mod, only: eps, s
    use neort_lib, only: neort_init, neort_prepare_splines, neort_setup_at_s
    use neort_freq, only: Om_tB
    use neort_profiles, only: Om_tE
    use shaing, only: omph_shaing
    use driftorbit, only: vth, etadt, etatp, sign_vpar

    implicit none

    real(dp), parameter :: s_test(2) = [0.05_dp, 0.12_dp]
    real(dp), parameter :: tol_eps = 5.0_dp   ! allowed |ratio-1| in units of eps
    real(dp), parameter :: frac(3) = [0.5_dp, 0.7_dp, 0.9_dp]
    real(dp) :: v, eta, OmtB, dOmtBdv, dOmtBdeta, shaing, ratio, worst
    integer :: i, j

    call neort_init("driftorbit.in", "circ.bc", "in_file_pert")
    call neort_prepare_splines("plasma.in", "profile.in")
    sign_vpar = 1
    worst = 0.0_dp
    do i = 1, size(s_test)
        call neort_setup_at_s(s_test(i))
        v = vth
        do j = 1, size(frac)
            eta = etatp + frac(j)*(etadt - etatp)
            call Om_tB(v, eta, OmtB, dOmtBdv, dOmtBdeta)
            shaing = omph_shaing(v, eta) - Om_tE
            ratio = shaing/OmtB
            print '(A,F6.3,A,F4.2,A,ES11.3,A,ES11.3,A,F8.4,A,F7.4)', ' s=', s, &
                ' f=', frac(j), ' Om_tB=', OmtB, ' shaing=', shaing, &
                ' ratio=', ratio, ' eps=', eps
            worst = max(worst, abs(ratio - 1.0_dp)/eps)
        end do
    end do
    print *, "max |ratio-1|/eps =", worst
    ! Observed 2-4 eps (the circular large-aspect-ratio model neglects the
    ! O(eps) shaping and current-gradient terms of the numerical precession);
    ! the unit defect gave |ratio-1| ~ 1, i.e. ~20 eps, and the wrong sign.
    if (worst > tol_eps) then
        print *, "FAIL: omph_shaing deviates from the numerical precession"
        error stop 1
    end if
    print *, "PASS"
end program test_omph_shaing_units
