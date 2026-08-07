program test_gc_full_resonance
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_gc_full_resonance, only: GC_RESONANCE_SUCCESS, GC_RESONANCE_PARTIAL, &
        find_gc_resonances

    implicit none

    real(dp) :: roots(4), derivatives(4)
    integer :: nroots, status

    call find_gc_resonances(two_roots, 0.0_dp, 1.0_dp, 40, 1.0e-12_dp, &
        1.0e-12_dp, roots, derivatives, nroots, status)
    if (status /= GC_RESONANCE_SUCCESS .or. nroots /= 2) &
        error stop "manufactured roots not recovered"
    if (maxval(abs(roots(1:2) - [0.25_dp, 0.75_dp])) > 1.0e-10_dp) &
        error stop "manufactured root positions mismatch"
    if (abs(derivatives(1) + 0.5_dp) > 1.0e-6_dp .or. &
        abs(derivatives(2) - 0.5_dp) > 1.0e-6_dp) &
        error stop "manufactured root derivatives mismatch"

    call find_gc_resonances(root_hidden_by_loss, 0.0_dp, 1.0_dp, 20, &
        1.0e-12_dp, 1.0e-12_dp, roots, derivatives, nroots, status)
    if (status /= GC_RESONANCE_SUCCESS .or. nroots /= 0) &
        error stop "search bridged a no-return interval"

    call find_gc_resonances(one_bad_one_good, 0.0_dp, 1.0_dp, 4, &
        1.0e-12_dp, 1.0e-12_dp, roots, derivatives, nroots, status)
    if (status /= GC_RESONANCE_PARTIAL .or. nroots /= 1) &
        error stop "valid root was discarded after a failed bracket"
    if (abs(roots(1) - 0.8_dp) > 1.0e-10_dp) &
        error stop "retained root position mismatch"

    write (*, '(A)') "test_gc_full_resonance OK"

contains

    subroutine two_roots(eta, residual, local_status)
        real(dp), intent(in) :: eta
        real(dp), intent(out) :: residual
        integer, intent(out) :: local_status

        residual = (eta - 0.25_dp)*(eta - 0.75_dp)
        local_status = 0
    end subroutine two_roots

    subroutine root_hidden_by_loss(eta, residual, local_status)
        real(dp), intent(in) :: eta
        real(dp), intent(out) :: residual
        integer, intent(out) :: local_status

        residual = eta - 0.5_dp
        local_status = merge(1, 0, eta >= 0.45_dp .and. eta <= 0.55_dp)
    end subroutine root_hidden_by_loss

    subroutine one_bad_one_good(eta, residual, local_status)
        real(dp), intent(in) :: eta
        real(dp), intent(out) :: residual
        integer, intent(out) :: local_status

        residual = (eta - 0.2_dp)*(eta - 0.8_dp)
        local_status = merge(1, 0, eta > 0.1_dp .and. eta < 0.15_dp)
    end subroutine one_bad_one_good

end program test_gc_full_resonance
