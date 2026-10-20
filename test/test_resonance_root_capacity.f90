module root_capacity_fixture
    use iso_fortran_env, only: dp => real64
    implicit none
contains
    subroutine cubic(v, eta, res, dresdeta)
        real(dp), intent(in) :: v, eta
        real(dp), intent(out) :: res, dresdeta

        if (v < 0.0_dp) error stop 'unused velocity'
        res = (eta - 0.1_dp) * (eta - 0.2_dp) * (eta - 0.8_dp)
        dresdeta = (eta - 0.2_dp) * (eta - 0.8_dp) &
            + (eta - 0.1_dp) * (eta - 0.8_dp) &
            + (eta - 0.1_dp) * (eta - 0.2_dp)
    end subroutine cubic
end module root_capacity_fixture

program test_resonance_root_capacity
    use iso_fortran_env, only: dp => real64
    use driftorbit, only: vth
    use neort_resonance, only: driftorbit_coarse, driftorbit_root, &
        resonance_override
    use root_capacity_fixture, only: cubic
    implicit none
    real(dp) :: roots(2, 3), eta_res(2)
    real(dp), parameter :: expected(2) = [0.1_dp, 0.2_dp]
    integer :: nroots, k

    vth = 1.0_dp
    resonance_override => cubic
    ! The first cell [0, 0.5] fills both rows with a fold pair. The ordinary
    ! sign change in [0.5, 1] must not write a third row or inflate nroots.
    call driftorbit_coarse(1.0_dp, 0.0_dp, 1.0_dp, roots, nroots)
    if (nroots /= 2) error stop 'root count exceeds available storage'
    do k = 1, nroots
        eta_res = driftorbit_root(1.0_dp, 0.0_dp, roots(k, 1), roots(k, 2))
        if (abs(eta_res(1) - expected(k)) > 1.0e-10_dp) &
            error stop 'capacity exhaustion corrupted the retained fold pair'
    end do
    nullify(resonance_override)
    print *, 'test_resonance_root_capacity PASSED'
end program test_resonance_root_capacity
