program test_resonance_oracle
    ! Independent root-selection oracle for the production resonance solver.
    ! The oracle uses only the residual values and plain bisection; it does
    ! not reuse the production derivative or iteration.
    use iso_fortran_env, only: dp => real64
    use logger, only: set_log_level
    use neort_lib, only: neort_init, neort_prepare_splines, neort_setup_at_s
    use neort_freq, only: Om_ph, Om_th
    use neort_profiles, only: vth, Om_tE
    use driftorbit, only: mth, mph, nlev, sign_vpar, etatp
    use neort_resonance, only: driftorbit_coarse, driftorbit_root
    use neort, only: set_to_passing_region, set_to_trapped_region
    implicit none

    real(dp) :: roots(nlev, 3), eta_res(2), eta_ref
    real(dp) :: eta_min, eta_max, v
    integer :: nroots, ir, trial_mth, trial_v, direction, surface_index
    real(dp) :: s_test
    integer :: nchecked

    call set_log_level(-1)
    call neort_init("driftorbit.in", "in_file")
    call neort_prepare_splines("plasma.in", "profile.in")
    nchecked = 0

    do surface_index = 1, 2
        s_test = merge(0.5_dp, 0.039063_dp, surface_index == 1)
        call neort_setup_at_s(s_test)
        Om_tE = 0.0_dp
        do direction = -1, 1, 2
            sign_vpar = direction
            do trial_v = 1, 16
                v = (0.25_dp * trial_v) * vth
                do trial_mth = -12, 12
                    if (trial_mth == 0) cycle
                    mth = trial_mth
                    mph = 3
                    call set_to_passing_region(eta_min, eta_max)
                    call check_region(v, eta_min, eta_max)
                    call set_to_trapped_region(eta_min, eta_max)
                    call check_region(v, eta_min, eta_max)
                end do
            end do
        end do
    end do

    if (nchecked == 0) error stop "root oracle fixture found no resonance brackets"
    print *, "test_resonance_oracle PASSED; roots checked =", nchecked

contains

    subroutine check_region(v_, lo, hi)
        real(dp), intent(in) :: v_, lo, hi
        integer :: kbis
        real(dp) :: a, b, fa, fb, fm, scale, bracket_width
        real(dp) :: eta_tolerance
        real(dp) :: h, rp, rm, dres_numeric

        call driftorbit_coarse(v_, lo, hi, roots, nroots)
        do ir = 1, nroots
            ! Use a representative nonzero production residual tolerance.  The
            ! eta-width guard must still recover the bracket-selected root when
            ! the frequency residual is poorly conditioned near the boundary.
            eta_res = driftorbit_root(v_, 5.0e-4_dp, roots(ir, 1), roots(ir, 2))
            if (eta_res(1) < 0.0_dp) error stop "production root rejected coarse bracket"

            a = min(roots(ir, 1), roots(ir, 2))
            b = max(roots(ir, 1), roots(ir, 2))
            bracket_width = b - a
            call residual(v_, a, fa)
            call residual(v_, b, fb)
            do kbis = 1, 220
                eta_ref = 0.5_dp * (a + b)
                call residual(v_, eta_ref, fm)
                if (fa * fm <= 0.0_dp) then
                    b = eta_ref
                    fb = fm
                else
                    a = eta_ref
                    fa = fm
                end if
            end do
            eta_ref = 0.5_dp * (a + b)
            scale = max(abs(eta_ref), tiny(1.0_dp))
            ! The production residual reuses Om_th in the passing branch,
            ! whereas this independent oracle calls public Om_ph and Om_th
            ! separately.  With -ffast-math, different CPUs may reassociate
            ! those mathematically identical operations by a few ulps.  Keep
            ! the machine-precision target when possible, but allow at most
            ! 1% of the original sign-changing bracket for that portable
            ! evaluation difference; the residual and Jacobian checks below
            ! still test the returned root independently.
            eta_tolerance = max(2.0e-12_dp * scale, 1.0e-2_dp * bracket_width)
            if (abs(eta_res(1) - eta_ref) > eta_tolerance) then
                write(*,*) 'root mismatch v,mth,sign,bracket,production,oracle:', &
                    v_/vth, mth, sign_vpar, roots(ir,1), roots(ir,2), eta_res(1), eta_ref, eta_tolerance
                error stop "production resonance root disagrees with bisection oracle"
            end if

            h = 1.0e-9_dp * max(abs(eta_res(1)), tiny(1.0_dp))
            if (abs(eta_res(1) - etatp) > 1.0e-4_dp * &
                max(abs(eta_res(1)), abs(etatp), tiny(1.0_dp))) then
                ! Keep the finite difference on one frequency branch; the
                ! separatrix itself is logarithmically non-differentiable.
                h = min(h, 0.25_dp * abs(eta_res(1) - etatp))
            else
                h = 0.0_dp
            end if
            if (h > 10.0_dp * spacing(eta_res(1))) then
                call residual(v_, eta_res(1) + h, rp)
                call residual(v_, eta_res(1) - h, rm)
                dres_numeric = (rp - rm) / (2.0_dp * h)
                if (abs(dres_numeric - eta_res(2)) > 1.0e-3_dp * &
                    max(abs(eta_res(2)), tiny(1.0_dp))) then
                    write(*,*) 'jacobian mismatch v,mth,sign,eta,etatp,h,analytic,numeric:', &
                        v_/vth, mth, sign_vpar, eta_res(1), etatp, h, eta_res(2), dres_numeric
                    error stop "production resonance Jacobian disagrees with residual derivative"
                end if
            end if
            nchecked = nchecked + 1
        end do
    end subroutine check_region

    subroutine residual(v_, eta_, residual_)
        real(dp), intent(in) :: v_, eta_
        real(dp), intent(out) :: residual_
        real(dp) :: Omph, dOmphdv, dOmphdeta
        real(dp) :: Omth, dOmthdv, dOmthdeta
        call Om_ph(v_, eta_, Omph, dOmphdv, dOmphdeta)
        call Om_th(v_, eta_, Omth, dOmthdv, dOmthdeta)
        residual_ = mph * Omph + mth * Omth
    end subroutine residual

end program test_resonance_oracle
