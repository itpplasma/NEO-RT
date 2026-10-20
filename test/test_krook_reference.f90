program test_krook_reference
    ! Retain the sensitive T co-passing m=4, velocity node 3 of the MARS L5/8
    ! reproducer. The fixture keeps 64 original equilibrium surfaces around
    ! s=.25 and a three-surface linear snapshot of all 641 perturbation modes.
    ! Reference: origin/main b501a00, native Fast on the same input, using the
    ! forward equilibrium recurrence. Equilibrium Horner changes this node by
    ! 1.99e-6; the cached forward evaluation agrees to 2.2e-12.
    use iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use neort_lib, only: neort_init, neort_setup_at_s
    use neort_species, only: selected_species
    use neort, only: set_to_passing_region
    use driftorbit, only: mth, sign_vpar, etamin, etamax
    use neort_transport, only: compute_krook_velocity
    implicit none
    real(dp), parameter :: expected(3) = [ &
        -2.55679990099455903e-9_dp, -8.98875624368494869e-9_dp, &
        3.59361902247261969e4_dp]
    real(dp) :: values(3), ux

    call neort_init('r0p50000.in', 'in_file', 'in_file_pert')
    selected_species = 2
    sign_vpar = 1
    call neort_setup_at_s(0.25_dp)
    mth = 4
    call set_to_passing_region(etamin, etamax)
    ux = 1.0e-6_dp + 2.5_dp*(6.0_dp - 1.0e-6_dp)/8
    call compute_krook_velocity(ux, values, -1.0_dp)
    if (.not. all(ieee_is_finite(values))) error stop 'nonfinite Krook reference'
    print *, 'Krook reference relative deviation:', abs(values/expected - 1.0_dp)
    if (any(abs(values - expected) > 1.0e-6_dp*abs(expected))) &
        error stop 'Krook transport exceeds the main reference bound'
end program test_krook_reference
