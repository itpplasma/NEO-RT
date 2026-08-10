program test_boozer_handedness
    ! Independent large-aspect-ratio circular oracle:
    ! R = R0 + a*sqrt(s)*cos(theta), Z = sigma_theta*a*sqrt(s)*sin(theta).
    ! The metric gives B_phi = sigma_phi*B0*R0,
    ! B_theta = r**2*iota*sigma_phi*B0/R0, and
    ! sqrt(g) = handedness*R0*a**2/2.
    use iso_fortran_env, only: dp => real64
    use do_magfie_mod, only: Bphcov, Bthcov, bfac, do_magfie, &
        init_magfie_at_s, inp_swi, magfie_thread_init, &
        read_boozer_file, set_s
    use logger, only: set_log_level
    use util, only: pi

    implicit none

    real(dp), parameter :: field_t = 2.0_dp
    real(dp), parameter :: major_radius_m = 5.0_dp
    real(dp), parameter :: minor_radius_m = 0.5_dp
    real(dp), parameter :: iota_abs = 0.2_dp
    real(dp), parameter :: evaluation_s = 0.5_dp
    real(dp), parameter :: current_to_covar = 0.2_dp
    real(dp), parameter :: meter_to_cm = 100.0_dp
    real(dp), parameter :: tesla_to_gauss = 1.0e4_dp
    real(dp), parameter :: tolerance = 1.0e-12_dp

    call set_log_level(-1)
    call write_chart("manufactured_left.bc", 1, -1, -1)
    call write_chart("manufactured_left_positive_iota.bc", 1, 1, -1)
    call write_chart("manufactured_right.bc", -1, 1, 1)
    call check_chart("manufactured_left.bc", 1, -1, -1)
    call check_chart("manufactured_left_positive_iota.bc", 1, 1, -1)
    call check_chart("manufactured_right.bc", -1, 1, 1)

contains

    subroutine check_chart(path, orientation, iota_sign, handedness)
        character(len=*), intent(in) :: path
        integer, intent(in) :: orientation, iota_sign, handedness

        real(dp) :: bder(3), bmod, hcovar(3), hctrvr(3), hcurl(3)
        real(dp) :: expected_Bphcov, expected_Bthcov, expected_sqrtg
        real(dp) :: iota_chart, radius_cm, sigma_phi, sqrtg, x(3)

        sigma_phi = -real(handedness, dp)/real(orientation, dp)
        iota_chart = real(iota_sign, dp)*iota_abs
        radius_cm = minor_radius_m*meter_to_cm*sqrt(evaluation_s)
        expected_Bphcov = sigma_phi*field_t*tesla_to_gauss*major_radius_m*meter_to_cm
        expected_Bthcov = field_t*tesla_to_gauss*radius_cm**2*iota_chart*sigma_phi &
            /(major_radius_m*meter_to_cm)
        expected_sqrtg = real(handedness, dp)*major_radius_m*meter_to_cm &
            *(minor_radius_m*meter_to_cm)**2/2.0_dp

        inp_swi = 9
        bfac = 1.0_dp
        call magfie_thread_init()
        call read_boozer_file(path)
        call set_s(evaluation_s)
        call init_magfie_at_s()
        x = [evaluation_s, 0.0_dp, 0.0_dp]
        call do_magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)

        call assert_close(path//" B_theta", Bthcov, expected_Bthcov)
        call assert_close(path//" B_phi", Bphcov, expected_Bphcov)
        call assert_close(path//" h_theta", hcovar(3)*bmod, expected_Bthcov)
        call assert_close(path//" h_phi", hcovar(2)*bmod, expected_Bphcov)
        call assert_close(path//" sqrt(g)", sqrtg, expected_sqrtg)
    end subroutine check_chart

    subroutine write_chart(path, orientation, iota_sign, handedness)
        character(len=*), intent(in) :: path
        integer, intent(in) :: orientation, iota_sign, handedness

        integer :: file_unit, surface_index
        real(dp) :: Bph_si, Bth_si, flux, iota_chart, radius_m
        real(dp) :: raw_Jpol, raw_Itor, s_surface, sigma_phi
        real(dp), parameter :: surfaces(3) = [0.25_dp, 0.5_dp, 0.75_dp]

        sigma_phi = -real(handedness, dp)/real(orientation, dp)
        iota_chart = real(iota_sign, dp)*iota_abs
        flux = -real(handedness, dp)*sigma_phi*pi*minor_radius_m**2*field_t
        Bph_si = sigma_phi*field_t*major_radius_m

        open (newunit=file_unit, file=path, status="replace", action="write")
        call write_header(file_unit, flux)
        do surface_index = 1, size(surfaces)
            s_surface = surfaces(surface_index)
            radius_m = minor_radius_m*sqrt(s_surface)
            Bth_si = field_t*radius_m**2*iota_chart*sigma_phi/major_radius_m
            raw_Jpol = Bph_si*tesla_to_gauss*meter_to_cm &
                /(current_to_covar*real(handedness, dp))
            raw_Itor = Bth_si*tesla_to_gauss*meter_to_cm &
                /(current_to_covar*real(handedness, dp))
            call write_surface(file_unit, s_surface, iota_chart, raw_Jpol, &
                raw_Itor, radius_m, orientation)
        end do
        close (file_unit)
    end subroutine write_chart

    subroutine write_header(file_unit, flux)
        integer, intent(in) :: file_unit
        real(dp), intent(in) :: flux

        write (file_unit, "(a)") "CC manufactured large-aspect-ratio circular chart"
        write (file_unit, "(a)") "CC independent analytic field and geometry"
        write (file_unit, "(a)") "CC m0b n0b nsurf nper flux a R"
        write (file_unit, "(a)") "CC"
        write (file_unit, "(a)") "CC"
        write (file_unit, *) 1, 0, 3, 1, flux, minor_radius_m, major_radius_m
    end subroutine write_header

    subroutine write_surface(file_unit, s_surface, iota_chart, Jpol, Itor, &
            radius_m, orientation)
        integer, intent(in) :: file_unit
        integer, intent(in) :: orientation
        real(dp), intent(in) :: s_surface, iota_chart, Jpol, Itor, radius_m

        real(dp) :: field_magnitude

        field_magnitude = field_t*sqrt(1.0_dp + (radius_m*iota_abs/major_radius_m)**2)
        write (file_unit, "(a)") "CC s iota Jpol Itor pprime sqrtg00"
        write (file_unit, "(a)") "CC units A A Pa m3"
        write (file_unit, *) s_surface, iota_chart, Jpol, Itor, 0.0_dp, 0.0_dp
        write (file_unit, "(a)") "CC m n rmnc rmns zmnc zmns vmnc vmns bmnc bmns"
        write (file_unit, *) 0, 0, major_radius_m, 0.0_dp, 0.0_dp, 0.0_dp, &
            0.0_dp, 0.0_dp, field_magnitude, 0.0_dp
        write (file_unit, *) 1, 0, radius_m, 0.0_dp, 0.0_dp, &
            real(orientation, dp)*radius_m, &
            0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp
    end subroutine write_surface

    subroutine assert_close(label, actual, expected)
        character(len=*), intent(in) :: label
        real(dp), intent(in) :: actual, expected

        real(dp) :: relative_error

        relative_error = abs(actual - expected)/max(abs(expected), tiny(expected))
        if (relative_error > tolerance) then
            print *, trim(label), " expected", expected, "got", actual
            error stop "manufactured Boozer oracle failed"
        end if
    end subroutine assert_close

end program test_boozer_handedness
