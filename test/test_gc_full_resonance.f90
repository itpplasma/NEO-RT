program test_gc_full_resonance
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_gc_full_resonance, only: GC_RESONANCE_SUCCESS, GC_RESONANCE_PARTIAL, &
        GC_RESONANCE_BOUNDARY_INVALID, find_gc_resonances

    implicit none

    real(dp) :: roots(8), derivatives(8)
    integer :: nroots, status

    call find_gc_resonances(invalid_open_endpoints, 0.0_dp, 1.0_dp, 8, &
        1.0e-12_dp, 1.0e-12_dp, roots, derivatives, nroots, status)
    if (status /= GC_RESONANCE_SUCCESS .or. nroots /= 1) &
        error stop "invalid open endpoints were treated as a partial scan"
    if (abs(roots(1) - 0.35_dp) > 1.0e-9_dp) &
        error stop "root in the shrunken open interval was not recovered"

    call find_gc_resonances(positive_local_minimum, 0.0_dp, 1.0_dp, 3, &
        1.0e-12_dp, 1.0e-10_dp, roots, derivatives, nroots, status)
    if (status /= GC_RESONANCE_SUCCESS .or. nroots /= 0) &
        error stop "positive local minimum was reported as a failed scan/root"

    call find_gc_resonances(two_roots, 0.0_dp, 1.0_dp, 40, 1.0e-12_dp, &
        1.0e-12_dp, roots, derivatives, nroots, status)
    if (status /= GC_RESONANCE_SUCCESS .or. nroots /= 2) &
        error stop "manufactured roots not recovered"
    if (maxval(abs(roots(1:2) - [0.25_dp, 0.75_dp])) > 1.0e-10_dp) &
        error stop "manufactured root positions mismatch"
    if (abs(derivatives(1) + 0.5_dp) > 1.0e-6_dp .or. &
        abs(derivatives(2) - 0.5_dp) > 1.0e-6_dp) &
        error stop "manufactured root derivatives mismatch"

    call find_gc_resonances(tangent_root, 0.0_dp, 1.0_dp, 3, &
        1.0e-12_dp, 1.0e-10_dp, roots, derivatives, nroots, status)
    if (status /= GC_RESONANCE_SUCCESS .or. nroots /= 1) &
        error stop "tangent root was not recovered"
    if (abs(roots(1) - 0.37_dp) > 1.0e-9_dp) &
        error stop "tangent root position mismatch"
    if (abs(derivatives(1)) > 1.0e-7_dp) &
        error stop "tangent root derivative mismatch"

    call find_gc_resonances(close_sign_roots, 0.0_dp, 1.0_dp, 4, &
        1.0e-12_dp, 1.0e-10_dp, roots, derivatives, nroots, status)
    if (status /= GC_RESONANCE_SUCCESS .or. nroots /= 2) &
        error stop "closely spaced roots were not both recovered"
    if (maxval(abs(roots(1:2) - [0.401_dp, 0.409_dp])) > 1.0e-9_dp) &
        error stop "closely spaced root positions mismatch"

    call find_gc_resonances(root_hidden_by_loss, 0.0_dp, 1.0_dp, 20, &
        1.0e-12_dp, 1.0e-12_dp, roots, derivatives, nroots, status)
    if (status /= GC_RESONANCE_PARTIAL .or. nroots /= 0) &
        error stop "no-return interval was not exposed"

    call find_gc_resonances(one_bad_one_good, 0.0_dp, 1.0_dp, 4, &
        1.0e-12_dp, 1.0e-12_dp, roots, derivatives, nroots, status)
    if (status /= GC_RESONANCE_PARTIAL .or. nroots /= 2) &
        error stop "valid root was discarded after a failed bracket"
    if (maxval(abs(roots(1:2) - [0.2_dp, 0.8_dp])) > 1.0e-10_dp) &
        error stop "retained root position mismatch"

    write (*, '(A)') "test_gc_full_resonance OK"

contains

    subroutine invalid_open_endpoints(eta, residual, local_status)
        real(dp), intent(in) :: eta
        real(dp), intent(out) :: residual
        integer, intent(out) :: local_status

        residual = eta - 0.35_dp
        local_status = merge(GC_RESONANCE_BOUNDARY_INVALID, 0, &
            eta == 0.0_dp .or. eta == 1.0_dp)
    end subroutine invalid_open_endpoints

    subroutine positive_local_minimum(eta, residual, local_status)
        real(dp), intent(in) :: eta
        real(dp), intent(out) :: residual
        integer, intent(out) :: local_status

        residual = (eta - 0.37_dp)**2 + 1.0e-4_dp
        local_status = 0
    end subroutine positive_local_minimum

    subroutine two_roots(eta, residual, local_status)
        real(dp), intent(in) :: eta
        real(dp), intent(out) :: residual
        integer, intent(out) :: local_status

        residual = (eta - 0.25_dp)*(eta - 0.75_dp)
        local_status = 0
    end subroutine two_roots

    subroutine tangent_root(eta, residual, local_status)
        real(dp), intent(in) :: eta
        real(dp), intent(out) :: residual
        integer, intent(out) :: local_status

        residual = (eta - 0.37_dp)**2
        local_status = 0
    end subroutine tangent_root

    subroutine close_sign_roots(eta, residual, local_status)
        real(dp), intent(in) :: eta
        real(dp), intent(out) :: residual
        integer, intent(out) :: local_status

        residual = (eta - 0.401_dp)*(eta - 0.409_dp)
        local_status = 0
    end subroutine close_sign_roots

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
