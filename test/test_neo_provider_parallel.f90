program test_neo_provider_parallel
    ! Shared fixture is produced by the independent handedness oracle. Each
    ! thread changes the prepared surface, so a stale global cache cannot pass.
    use iso_fortran_env, only: dp => real64
    use do_magfie_mod, only: read_boozer_file, set_s, init_magfie_at_s, &
        do_magfie, Bthcov, Bphcov, q, inp_swi, bfac
    implicit none
    integer :: j
    inp_swi = 9
    bfac = 1.0_dp
    call read_boozer_file('manufactured_left.bc')
    !$omp parallel do num_threads(2) schedule(static, 1)
    do j = 1, 8
        call check_surface(0.25_dp + 0.05_dp*j)
    end do
    !$omp end parallel do
contains
    subroutine check_surface(surface)
        real(dp), intent(in) :: surface
        real(dp) :: x(3), bmod, sqrtg, bder(3), hc(3), ht(3), curl(3)
        real(dp) :: expected_theta
        call set_s(surface)
        call init_magfie_at_s()
        x = [surface, 0.0_dp, 0.7_dp]
        call do_magfie(x, bmod, sqrtg, bder, hc, ht, curl)
        expected_theta = -2.0e4_dp*(50.0_dp**2*surface)*0.2_dp/500.0_dp
        if (abs(Bthcov/expected_theta - 1.0_dp) > 1.0e-10_dp) &
            error stop 'thread-local surface or stale NEO-2 cache'
        if (abs(Bphcov/1.0e7_dp - 1.0_dp) > 1.0e-10_dp) &
            error stop 'thread-local covariant field'
        if (abs(q + 5.0_dp) > 1.0e-10_dp) error stop 'thread-local safety factor'
    end subroutine check_surface
end program test_neo_provider_parallel
