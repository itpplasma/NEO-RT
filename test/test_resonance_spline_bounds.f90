program test_resonance_spline_bounds
    use iso_fortran_env, only: dp => real64
    use logger, only: set_log_level
    use neort_lib, only: neort_init, neort_prepare_splines, neort_setup_at_s
    use neort, only: set_to_passing_region, set_to_trapped_region
    use driftorbit, only: etatp, etadt, epsp_spl, epsst_spl, epst_spl

    implicit none

    real(dp) :: eta_min, eta_max

    call set_log_level(-1)
    call neort_init("driftorbit.in", "in_file")
    call neort_prepare_splines("plasma.in", "profile.in")
    call neort_setup_at_s(0.5_dp)

    call set_to_trapped_region(eta_min, eta_max)
    if (eta_min < (1.0_dp + epst_spl) * etatp) then
        error stop "trapped lower bound reaches the spline separatrix"
    end if

    call set_to_passing_region(eta_min, eta_max)
    if (eta_max > (1.0_dp - epsp_spl) * etatp) then
        error stop "passing upper bound reaches the spline separatrix"
    end if

    print *, "test_resonance_spline_bounds PASSED"
end program test_resonance_spline_bounds
