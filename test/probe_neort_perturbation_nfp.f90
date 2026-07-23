program probe_neort_perturbation_nfp
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use do_magfie_mod, only: bfac, set_s
    use do_magfie_pert_mod, only: do_magfie_pert, do_magfie_pert_amp, &
        init_magfie_pert_at_s, inp_swi_pert, mph, read_boozer_pert_file

    implicit none

    character(len=*), parameter :: fixture = 'manufactured_perturbation_nfp.bc'
    integer, parameter :: field_periods = 3, reduced_n = 2
    real(dp), parameter :: theta = 0.6_dp, phi = -0.4_dp
    real(dp), parameter :: cosine_amplitude_t = 2.0e-4_dp
    real(dp), parameter :: sine_amplitude_t = 0.7e-4_dp
    real(dp), parameter :: tolerance_g = 1.0e-12_dp
    complex(dp) :: amplitude_g
    real(dp) :: field_g, period_field_g, expected_g, wrong_reduced_g, x(3)

    call write_fixture(fixture)
    inp_swi_pert = 9
    bfac = 1.0_dp
    call set_s(0.5_dp)
    call read_boozer_pert_file(fixture)
    call init_magfie_pert_at_s()
    x = [0.5_dp, phi, theta]
    call do_magfie_pert_amp(x, amplitude_g)
    call do_magfie_pert(x, field_g)
    x(2) = phi + 2.0_dp*acos(-1.0_dp)/real(field_periods, dp)
    call do_magfie_pert(x, period_field_g)
    expected_g = expected_value(theta + field_periods*reduced_n*phi)
    wrong_reduced_g = expected_value(theta + reduced_n*phi)
    write (*, '(*(g0,1x))') 'NEORT_PERT_NFP', 'nfp', field_periods, &
        'reduced_n', reduced_n, 'physical_n', mph, 'field_g', field_g, &
        'period_field_g', period_field_g, 'expected_g', expected_g, &
        'wrong_reduced_g', wrong_reduced_g, 'amplitude_real_g', real(amplitude_g), &
        'amplitude_imag_g', aimag(amplitude_g)
    if (abs(field_g - expected_g) > tolerance_g) error stop 'wrong physical n'
    if (abs(period_field_g - field_g) > tolerance_g) error stop 'wrong NFP period'

contains

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
        write (unit, '(a)') 'CC manufactured NFP perturbation kernel'
        write (unit, '(a)') 'CC physical n equals nper times stored n'
        write (unit, '(a)') 'CC kernel m*theta+n*phi_CCW'
        write (unit, '(a)') 'CC complex coefficient bmnc-i*bmns'
        write (unit, '(a)') 'm0b n0b nsurf nper flux a R'
        write (unit, '(*(g0,1x))') 2, 0, 3, field_periods, 1.0_dp, 1.0_dp, 1.0_dp
        do surface = 1, size(surfaces)
            call write_surface(unit, surfaces(surface))
        end do
        close (unit)
    end subroutine write_fixture

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

        write (unit, '(*(g0,1x))') poloidal_mode, reduced_n, 0.0_dp, 0.0_dp, &
            0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, bmnc, bmns
    end subroutine write_mode

end program probe_neort_perturbation_nfp
