program test_phi1_reader
    use iso_fortran_env, only: dp => real64
    use neort_phi1, only: phi1, phi1file, phi1_scale, read_phi1_file, &
                          do_elec_pert_amp, phi1_amplitude_bound
    use neort_config, only: config_t, set_config, read_and_set_config
    use do_magfie_pert_mod, only: set_mph
    use do_magfie_mod, only: bfac
    use util, only: imun
    implicit none
    character(32) :: mode
    integer :: unit, ns, k, j, n, rows
    real(dp) :: surface, x(3), bound
    complex(dp) :: actual, expected
    type(config_t) :: config
    call get_command_argument(1, mode)
    ns = 3
    if (mode == 'single') ns = 1
    if (mode == 'two') ns = 2
    rows = 3
    if (mode == 'one-mode') rows = 1
    call set_mph(-4)
    open (newunit=unit, file='phi1_fixture.bc', status='replace')
    write (unit, '(a)') 'CC synthetic affine Phi coefficients, test_phi1_reader.f90'
    if (mode /= 'quantity') &
        write (unit, '(a)') 'CC quantity: Phi_1 [statvolt] in bmnc/bmns columns'
    write (unit, '(a)') ' m0b n0b nsurf nper flux a R'
    write (unit, *) rows - 1, 0, ns, 2, 0.0_dp, 0.0_dp, 0.0_dp
    do k = 1, ns
        surface = 0.2_dp*k
        write (unit, '(a)') ' s iota Jpol Itor pprime sqrtg'
        write (unit, '(a)') ' units'
        write (unit, *) surface, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp
        write (unit, '(a)') ' m n rmnc rmns zmnc zmns vmnc vmns bmnc bmns'
        do j = 1, rows
            n = -2
            if (mode == 'n' .or. (mode == 'mixed-n' .and. j == 2)) n = 2
            write (unit, *) j - 2, n, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                0.0_dp, 0.0_dp, j*(1 + surface), -0.3_dp*j*(1 + surface)
        end do
    end do
    close (unit)
    phi1 = .true.
    phi1file = 'phi1_fixture.bc'
    call read_phi1_file()
    phi1_scale = 2.0_dp
    bfac = 7.0_dp
    x = [0.35_dp, 0.7_dp, -0.42_dp]
    surface = x(1)
    if (ns == 1) surface = 0.2_dp
    expected = (0.0_dp, 0.0_dp)
    do j = 1, rows
        expected = expected + 2*j*(1 + surface)*cmplx(1.0_dp, 0.3_dp, dp)* &
                   exp(imun*(j - 2)*x(3))
    end do
    call do_elec_pert_amp(x, actual)
  if (abs(actual - expected) > 1.0e-12_dp) error stop 'Phi spline/Horner units or phase'
    if (abs(real(actual*exp(-4*imun*x(2)), dp) &
            - real(expected*exp(-4*imun*x(2)), dp)) > 1.0e-12_dp) &
        error stop 'real Phi angle convention'
    call phi1_amplitude_bound(x, bound)
    if (bound < abs(actual)) error stop 'Phi Fourier triangle bound'
    !$omp parallel private(actual)
    call do_elec_pert_amp(x, actual)
    if (abs(actual - expected) > 1.0e-12_dp) error stop 'Phi thread state'
    !$omp end parallel
    ! Both configuration interfaces must reset the default-off switch.
    call set_config(config)
    call do_elec_pert_amp(x, actual)
    if (actual /= (0.0_dp, 0.0_dp)) error stop 'config default does not disable Phi'
    open (newunit=unit, file='phi1_config.in', status='replace')
    write (unit, '(a)') '&params phi1=T, phi_mn=(0.2,-0.1), m0_phi=2, mph=-4 /'
    close (unit)
    call read_and_set_config('phi1_config.in')
    call do_elec_pert_amp(x, actual)
    expected = cmplx(0.2_dp, -0.1_dp, dp)*exp(2*imun*x(3))
    if (abs(actual - expected) > 1.0e-14_dp) error stop 'analytic Phi namelist'
    print *, 'test_phi1_reader PASSED ', trim(mode)
end program test_phi1_reader
