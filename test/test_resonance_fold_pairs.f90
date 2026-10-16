module fold_pair_fixture
    use iso_fortran_env, only: dp => real64
    implicit none
    ! Manufactured resonance R = (eta - c)**2 - w**2 (h: coarse interval).
    real(dp) :: c = 0.0_dp, w = 0.0_dp, h = 0.0_dp

contains


    subroutine parabola(v, eta, res, dresdeta)
        real(dp), intent(in) :: v, eta
        real(dp), intent(out) :: res, dresdeta
        if (v < 0.0_dp) error stop 'unused velocity'
        res = (eta - c)**2 - w**2
        dresdeta = 2.0_dp * (eta - c)
    end subroutine parabola

    subroutine shifted_parabola(v, eta, res, dresdeta)
        real(dp), intent(in) :: v, eta
        real(dp), intent(out) :: res, dresdeta
        if (v < 0.0_dp) error stop 'unused velocity'
        res = (eta - c)**2 + 1.0e-3_dp * h**2
        dresdeta = 2.0_dp * (eta - c)
    end subroutine shifted_parabola

    subroutine line(v, eta, res, dresdeta)
        real(dp), intent(in) :: v, eta
        real(dp), intent(out) :: res, dresdeta
        if (v < 0.0_dp) error stop 'unused velocity'
        res = eta - 1.234567_dp
        dresdeta = 1.0_dp
    end subroutine line
end module fold_pair_fixture

program test_resonance_fold_pairs
    ! Two resonances closer than one coarse eta interval must both be found.
    ! Near a fold of R(v, eta) (an extremum in eta), a root pair is born at
    ! one velocity and separates like sqrt(v - v_fold).  The pair first lies
    ! inside one coarse interval, where R has the same sign at both ends; the
    ! transport sum then drops the pair, although its weight 1/|dR/deta| is
    ! largest there.  Oracle: R = (eta - c)**2 - w**2 has the exact roots c - w
    ! and c + w, each simple with |dR/deta| = 2 w.
    use iso_fortran_env, only: dp => real64
    use driftorbit, only: nlev
    use neort_resonance, only: driftorbit_coarse, driftorbit_root, &
        resonance_override
    use fold_pair_fixture, only: c, w, h, parabola, shifted_parabola, line
    implicit none

    real(dp), parameter :: eta_min = 1.0_dp, eta_max = 2.0_dp
    real(dp) :: roots(nlev, 3), eta_res(2), expected(2)
    integer :: nroots, k, kcase
    real(dp), parameter :: offsets(3) = [0.5_dp, 0.3_dp, 0.77_dp]
    real(dp), parameter :: widths(3) = [0.2_dp, 1.0e-3_dp, 0.45_dp]

    resonance_override => parabola
    h = (eta_max - eta_min) / nlev
    do kcase = 1, size(offsets)
        ! Centre inside interval 37, half-separation a fraction of h, so both
        ! roots lie strictly inside the same coarse interval.
        c = eta_min + (36.0_dp + offsets(kcase)) * h
        w = widths(kcase) * min(offsets(kcase), 1.0_dp - offsets(kcase)) * h
        expected = [c - w, c + w]
        call driftorbit_coarse(1.0_dp, eta_min, eta_max, roots, nroots)
        if (nroots /= 2) then
            write(*,*) 'case, nroots, c, w:', kcase, nroots, c, w
            error stop 'root pair inside one coarse interval was not bracketed'
        end if
        do k = 1, 2
            eta_res = driftorbit_root(1.0_dp, 0.0_dp, roots(k, 1), roots(k, 2))
            if (abs(eta_res(1) - expected(k)) > 4.0e-12_dp) then
                write(*,*) 'case, root, found, exact:', kcase, k, eta_res(1), expected(k)
                error stop 'bracketed root differs from the exact root'
            end if
            if (abs(abs(eta_res(2)) / (2.0_dp * w) - 1.0_dp) > 1.0e-5_dp) then
                write(*,*) 'case, root, slope, exact:', kcase, k, eta_res(2), 2.0_dp * w
                error stop 'root slope differs from the exact slope'
            end if
        end do
    end do

    ! Same-sign extremum without a root, and a simple root: unchanged counts.
    c = eta_min + 36.5_dp * h
    w = 0.0_dp
    call driftorbit_coarse(1.0_dp, eta_min, eta_max, roots, nroots)
    if (nroots /= 0) error stop 'extremum touching zero must not produce roots'
    resonance_override => shifted_parabola
    call driftorbit_coarse(1.0_dp, eta_min, eta_max, roots, nroots)
    if (nroots /= 0) error stop 'positive extremum must not produce roots'
    resonance_override => line
    call driftorbit_coarse(1.0_dp, eta_min, eta_max, roots, nroots)
    if (nroots /= 1) error stop 'simple root must give one bracket'
    nullify(resonance_override)
    print *, 'test_resonance_fold_pairs PASSED'

end program test_resonance_fold_pairs
