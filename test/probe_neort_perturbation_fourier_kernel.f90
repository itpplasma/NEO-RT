program probe_neort_perturbation_fourier_kernel
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_config, only: config_t, read_and_set_config, set_config
    use neort, only: check_magfie
    use neort_datatypes, only: magfie_data_t
    use do_magfie_mod, only: bfac, do_magfie, init_magfie_at_s, inp_swi, &
        magfie_thread_init, read_boozer_file, set_s
    use do_magfie_pert_mod, only: do_magfie_pert, do_magfie_pert_amp, &
        init_magfie_pert_at_s, read_boozer_pert_file
    use driftorbit, only: epsmn, m0, pertfile, pertfile_scale
    use neort_magfie, only: B0

    implicit none

    character(len=*), parameter :: fixture = 'manufactured_perturbation.bc'
    character(len=*), parameter :: equilibrium_fixture = 'manufactured_equilibrium.bc'
    real(dp), parameter :: theta = 0.6_dp, phi = -0.4_dp
    real(dp), parameter :: cosine_amplitude_t = 2.0e-4_dp
    real(dp), parameter :: sine_amplitude_t = 0.7e-4_dp
    complex(dp) :: amplitude_g, expected_consumer
    real(dp) :: bmod, sqrtg, bder(3), hcovar(3), hctrvr(3), hcurl(3)
    real(dp) :: full_field_g, expected_g, wrong_sign_g, x(3), theta_sample
    integer :: sample
    type(config_t) :: config
    type(magfie_data_t) :: magfie_data

    call write_fixture(fixture)
    call write_equilibrium_fixture(equilibrium_fixture)

    call write_config("mixed_readers.in")
    call read_and_set_config("mixed_readers.in")
    call check_manufactured_field("namelist mixed readers")

    config%s = 0.5_dp
    config%inp_swi = 9
    config%pertfile = .true.
    call set_config(config)
    call check_manufactured_field("legacy inherited reader")

    config = config_t()
    config%s = 0.5_dp
    config%qs = 1.0_dp
    config%ms = 1.0_dp
    config%inp_swi = 9
    config%inp_swi_pert = 9
    config%pertfile = .true.
    config%pertfile_scale = 2.0_dp
    call set_config(config)
    bfac = 1.0_dp
    B0 = 1.0_dp
    call magfie_thread_init()
    call read_boozer_file(equilibrium_fixture)
    call set_s(0.5_dp)
    call init_magfie_at_s()
    call read_boozer_pert_file(fixture)
    call init_magfie_pert_at_s()
    call check_magfie(magfie_data)

    sample = 13
    theta_sample = -acos(-1.0_dp) + (sample - 1) * 2.0_dp * acos(-1.0_dp) / 49.0_dp
    x = [0.5_dp, 0.0_dp, theta_sample]
    call do_magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
    call do_magfie_pert_amp(x, amplitude_g)
    expected_consumer = pertfile_scale * amplitude_g / bmod
    if (abs(magfie_data%tensors(sample)%bn - expected_consumer) > 1.0e-12_dp) then
        error stop 'pertfile_scale does not scale file perturbations'
    end if

    config%pertfile = .false.
    config%epsmn = 0.25_dp
    config%m0 = 2
    call set_config(config)
    call check_magfie(magfie_data)
    expected_consumer = epsmn * exp((0.0_dp, 1.0_dp) * real(m0, dp) * theta_sample)
    if (abs(magfie_data%tensors(sample)%bn - expected_consumer) > 1.0e-12_dp) then
        error stop 'epsmn does not control analytic perturbations independently'
    end if

    call delete_file(fixture)
    call delete_file(equilibrium_fixture)
    call delete_file("mixed_readers.in")

contains

    subroutine check_manufactured_field(label)
        character(len=*), intent(in) :: label

        call read_boozer_pert_file(fixture)
        call init_magfie_pert_at_s()

        x = [0.5_dp, phi, theta]
        call do_magfie_pert_amp(x, amplitude_g)
        call do_magfie_pert(x, full_field_g)
        expected_g = expected_value(theta + 3.0_dp*phi)
        wrong_sign_g = expected_value(theta - 3.0_dp*phi)

        write (*, '(*(g0,1x))') 'NEORT_PERT_KERNEL', 'case', trim(label), &
            'theta', theta, 'phi_ccw', phi, 'amplitude_real_g', real(amplitude_g), &
            'amplitude_imag_g', aimag(amplitude_g), 'full_field_g', full_field_g, &
            'expected_plus_n_g', expected_g, 'wrong_minus_n_g', wrong_sign_g
        if (abs(full_field_g - expected_g) > 1.0e-12_dp) then
            error stop 'perturbation reader does not reproduce the manufactured field'
        end if
    end subroutine check_manufactured_field

    real(dp) function expected_value(phase)
        real(dp), intent(in) :: phase

        expected_value = 1.0e4_dp*(cosine_amplitude_t*cos(phase) &
            + sine_amplitude_t*sin(phase))
    end function expected_value

    subroutine write_fixture(path)
        character(len=*), intent(in) :: path

        integer :: unit, surface
        real(dp), parameter :: surfaces(3) = [0.2_dp, 0.5_dp, 0.8_dp]

        open (newunit=unit, file=path, status='replace', action='write')
        write (unit, '(a)') 'CC manufactured perturbation kernel'
        write (unit, '(a)') 'CC physical phi is CCW'
        write (unit, '(a)') 'CC kernel m*theta+n*phi'
        write (unit, '(a)') 'CC complex coefficient bmnc-i*bmns'
        write (unit, '(a)') 'm0b n0b nsurf nper flux a R'
        write (unit, '(*(g0,1x))') 2, 0, 3, 1, 1.0_dp, 1.0_dp, 1.0_dp
        do surface = 1, size(surfaces)
            call write_surface(unit, surfaces(surface))
        end do
        close (unit)
    end subroutine write_fixture

    subroutine write_equilibrium_fixture(path)
        character(len=*), intent(in) :: path

        integer :: unit, surface
        real(dp), parameter :: surfaces(3) = [0.2_dp, 0.5_dp, 0.8_dp]

        open (newunit=unit, file=path, status='replace', action='write')
        write (unit, '(a)') 'CC manufactured equilibrium field'
        write (unit, '(a)') 'CC independent synthetic equilibrium oracle'
        write (unit, '(a)') 'CC m0b n0b nsurf nper flux a R'
        write (unit, '(a)') 'CC'
        write (unit, '(a)') 'CC'
        write (unit, '(*(g0,1x))') 1, 0, 3, 1, 1.0_dp, 0.5_dp, 5.0_dp
        do surface = 1, size(surfaces)
            write (unit, '(a)') 'CC s iota Jpol Itor pprime sqrtg00'
            write (unit, '(a)') 'CC units A A Pa m3'
            write (unit, '(*(g0,1x))') surfaces(surface), 0.2_dp, 100.0_dp, &
                10.0_dp, 0.0_dp, 1.0_dp
            write (unit, '(a)') 'CC m n rmnc rmns zmnc zmns vmnc vmns bmnc bmns'
            write (unit, '(*(g0,1x))') 0, 0, 5.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                0.0_dp, 0.0_dp, 2.0_dp, 0.0_dp
            write (unit, '(*(g0,1x))') 1, 0, 0.5_dp, 0.0_dp, 0.0_dp, 0.5_dp, &
                0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp
        end do
        close (unit)
    end subroutine write_equilibrium_fixture

    subroutine write_config(path)
        character(len=*), intent(in) :: path

        integer :: unit

        open (newunit=unit, file=path, status='replace', action='write')
        write (unit, '(a)') '&params'
        write (unit, '(a)') '  s = 0.5'
        write (unit, '(a)') '  qs = 1.0'
        write (unit, '(a)') '  ms = 1.0'
        write (unit, '(a)') '  mph = 3'
        write (unit, '(a)') '  pertfile = .true.'
        write (unit, '(a)') '  inp_swi = 10'
        write (unit, '(a)') '  inp_swi_pert = 9'
        write (unit, '(a)') '/'
        close (unit)
    end subroutine write_config

    subroutine delete_file(path)
        character(len=*), intent(in) :: path

        integer :: unit

        open (newunit=unit, file=path, status='old')
        close (unit, status='delete')
    end subroutine delete_file

    subroutine write_surface(unit, radial_coordinate)
        integer, intent(in) :: unit
        real(dp), intent(in) :: radial_coordinate

        write (unit, '(a)') 's iota Jpol Itor pprime sqrtg'
        write (unit, '(a)') 'units A A Pa m3'
        write (unit, '(*(es24.16e3,1x))') radial_coordinate, 1.0_dp, &
            0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp
        write (unit, '(a)') 'm n rmnc rmns zmnc zmns vmnc vmns bmnc bmns'
        call write_mode(unit, -1, 0.0_dp, 0.0_dp)
        call write_mode(unit, 0, 0.0_dp, 0.0_dp)
        call write_mode(unit, 1, cosine_amplitude_t, sine_amplitude_t)
    end subroutine write_surface

    subroutine write_mode(unit, poloidal_mode, bmnc, bmns)
        integer, intent(in) :: unit, poloidal_mode
        real(dp), intent(in) :: bmnc, bmns

        write (unit, '(*(g0,1x))') poloidal_mode, 3, 0.0_dp, 0.0_dp, &
            0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, bmnc, bmns
    end subroutine write_mode

end program probe_neort_perturbation_fourier_kernel
