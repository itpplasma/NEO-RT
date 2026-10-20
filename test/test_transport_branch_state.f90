program test_transport_branch_state
    ! Identical physical inputs must not produce transport that depends on the
    ! cached parallel direction left by an earlier orbit class.
    use iso_fortran_env, only: dp => real64
    use logger, only: set_log_level
    use neort, only: set_to_passing_region
    use neort_lib, only: neort_init, neort_prepare_splines, neort_setup_at_s
    use neort_profiles, only: vth, Om_tE
    use neort_transport, only: compute_transport_integral
    use driftorbit, only: mth, mph, etamin, etamax, sign_vpar, &
        sign_vpar_htheta
    implicit none

    real(dp) :: D_from_negative_cache(2), D_from_positive_cache(2)
    real(dp) :: T_from_negative_cache, T_from_positive_cache
    integer :: direction
    real(dp) :: scale

    call set_log_level(-1)
    call neort_init('driftorbit.in', 'in_file')
    call neort_prepare_splines('plasma.in', 'profile.in')
    call neort_setup_at_s(0.5_dp)
    Om_tE = 0.0_dp

    mph = 3
    call set_to_passing_region(etamin, etamax)

    ! Different harmonics provide nonzero resonances in the two directions
    ! for the base equilibrium at s=0.5 and zero electric precession.
    do direction = 1, -1, -2
        mth = merge(-4, -3, direction == 1)
        sign_vpar = real(direction, dp)

        sign_vpar_htheta = -1.0_dp
        call compute_transport_integral(1.0e-6_dp*vth, 4.0_dp*vth, 32, -1.0_dp, &
            D_from_negative_cache, T_from_negative_cache)

        sign_vpar_htheta = 1.0_dp
        call compute_transport_integral(1.0e-6_dp*vth, 4.0_dp*vth, 32, -1.0_dp, &
            D_from_positive_cache, T_from_positive_cache)

        if (any(D_from_negative_cache <= 1.0e-14_dp) .or. &
            abs(T_from_negative_cache) <= 1.0e-2_dp) then
            error stop 'branch-state fixture no longer has a transport signal'
        end if

        scale = max(maxval(abs(D_from_negative_cache)), &
            maxval(abs(D_from_positive_cache)), tiny(1.0_dp))
        if (maxval(abs(D_from_negative_cache-D_from_positive_cache)) > &
            1.0e-12_dp*scale) then
            error stop 'transport depends on stale parallel-direction cache'
        end if

        scale = max(abs(T_from_negative_cache), abs(T_from_positive_cache), &
            tiny(1.0_dp))
        if (abs(T_from_negative_cache-T_from_positive_cache) > &
            1.0e-12_dp*scale) then
            error stop 'torque depends on stale parallel-direction cache'
        end if

        print *, 'direction, mth, D(-cache), D(+cache), T(-cache), T(+cache):', &
            direction, mth, D_from_negative_cache, D_from_positive_cache, &
            T_from_negative_cache, T_from_positive_cache
    end do
    print *, 'test_transport_branch_state PASSED'
end program test_transport_branch_state
