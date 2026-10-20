program test_krook_c12_curvature_cancellation
    ! Exact opposite/constant/0.01 C12 deck, D counter-passing m=+11,
    ! first thermal speed node. Its resonant Fourier harmonic is below the
    ! orbit accuracy. A zero absolute tolerance used to exhaust 16384 pieces.
    use iso_fortran_env, only: dp => real64
    use ieee_arithmetic, only: ieee_is_finite
    use neort_lib, only: neort_init, neort_setup_at_s
    use neort_species, only: selected_species
    use neort, only: prepare_transport, select_component, component_ctr
    use neort_transport, only: compute_krook_velocity
    use driftorbit, only: pertfile_scale
    use do_magfie_mod, only: s
    implicit none
    real(dp), parameter :: ux = 0.01562599609375_dp
    ! Native zero-atol adaptive estimate, rounded to 0.1%.
    ! This tiny channel is tested to 1%, not advertised as relatively resolved.
    real(dp), parameter :: reference = 3.48e-34_dp
    real(dp), parameter :: scales(3) = [1.0_dp, 1.0e-3_dp, 1.0e3_dp]
    real(dp) :: vmin, vmax, values(3), normalized(3), baseline(3)
    integer :: lower, upper, k
    call neort_init('driftorbit.in', 'in_file', 'in_file_pert', 'in_file_pert_kappa')
    selected_species = 1
    call neort_setup_at_s(s)
    call prepare_transport(lower, upper)
    call select_component(11, component_ctr, vmin, vmax)
    do k = 1, size(scales)
        pertfile_scale = scales(k)
        call compute_krook_velocity(ux, values, -1.0_dp)
        if (.not. all(ieee_is_finite(values))) error stop 'nonfinite C12 sample'
        normalized = values/scales(k)**2
        print '(a,es12.4,3es26.17)', 'C12 scale and values:', scales(k), normalized
        if (abs(normalized(1)/reference - 1.0_dp) > 0.01_dp) &
            error stop 'C12 cancelling channel is outside its bounded reference'
        if (abs(normalized(2)/(ux**2*normalized(1)) - 1.0_dp) > 1.0e-12_dp) &
            error stop 'C12 energy moment changed'
        if (k == 1) baseline = normalized
        if (any(abs(normalized/baseline - 1.0_dp) > 1.0e-3_dp)) &
            error stop 'C12 absolute error budget does not scale quadratically'
    end do
end program test_krook_c12_curvature_cancellation
