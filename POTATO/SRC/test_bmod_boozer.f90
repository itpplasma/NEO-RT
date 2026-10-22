program test_bmod_boozer
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use potato_input_mod, only: read_potato_input, n_tor
    use do_magfie_pert_mod, only: do_magfie_pert, do_magfie_pert_amp, &
        init_magfie_pert_at_s
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp), phi = 0.37_dp
    real(dp) :: s, theta, nu, r, z, x(3), field, err, worst
    complex(dp) :: got, amp, expected
    integer :: j

    call read_potato_input('potato.in')
    ! First call loads shared coefficients, then workers allocate their own caches.
    call bmod_pert(320.0_dp, 0.0_dp, got)
    worst = 0.0_dp
    !$omp parallel do private(s,theta,nu,r,z,x,got,amp,expected,field,err) &
    !$omp reduction(max:worst)
    do j = 0, 191
        s = 0.15_dp + 0.7_dp*real(mod(j, 23), dp)/22.0_dp
        theta = 2*pi*real(mod(j, 64), dp)/64.0_dp
        nu = 0.13_dp + 0.02_dp*sin(theta)
        r = 300.0_dp + 100.0_dp*s*cos(theta)
        z = 100.0_dp*s*sin(theta)
        call bmod_pert(r, z, got)
        call init_magfie_pert_at_s()
        x = [s, phi + nu, theta]
        call do_magfie_pert_amp(x, amp)
        expected = amp*exp(cmplx(0.0_dp, n_tor*nu, dp))
        call do_magfie_pert(x, field)
        err = max(abs(got - expected), &
            abs(real(got*exp(cmplx(0.0_dp, n_tor*phi, dp)), dp) - field))
        worst = max(worst, err/abs(expected))
    end do
    !$omp end parallel do
    print *, 'maximum relative complex/real field error:', worst
    if (worst > 1.0e-10_dp) error stop 'Boozer perturbation mismatch'
end program test_bmod_boozer
