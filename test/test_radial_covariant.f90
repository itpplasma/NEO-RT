program test_radial_covariant
    ! Independent oracle for the source-bound .bc radial covariant hook.
    ! The manufactured chart is circular with v_shift=s, so the geometric
    ! dot-product term vanishes analytically and
    !     B_s = Bphcov * 2*pi * d(v_shift)/ds = Bphcov*2*pi.
    ! Jpol and Itor are intentionally different: only Jpol/nper maps to the
    ! toroidal field function Bphcov.  The native do_magfie hcovar(1) must
    ! remain zero on this diagnostic branch.
    use iso_fortran_env, only: dp => real64
    use do_magfie_mod, only: Bphcov, bfac, do_magfie, &
        finite_orbit_toroidal_phase, init_magfie_at_s, inp_swi, iota, &
        magfie_thread_init, psi_pr, radial_covariant_component, &
        read_boozer_file, set_s
    use logger, only: set_log_level
    use util, only: c, mi, pi, qi

    implicit none

    real(dp), parameter :: field_t = 2.0_dp
    real(dp), parameter :: major_radius_m = 5.0_dp
    real(dp), parameter :: minor_radius_m = 0.5_dp
    real(dp), parameter :: iota_chart = -0.2_dp
    real(dp), parameter :: current_to_covar = 0.2_dp
    real(dp), parameter :: meter_to_cm = 100.0_dp
    real(dp), parameter :: tesla_to_gauss = 1.0e4_dp
    real(dp), parameter :: evaluation_s = 0.5_dp
    real(dp), parameter :: tolerance = 1.0e-12_dp
    real(dp), parameter :: theta_samples(3) = [0.0_dp, 0.7_dp, 2.3_dp]

    real(dp) :: bder(3), bmod, hcovar(3), hctrvr(3), hcurl(3)
    real(dp) :: a_phi_prime, delta_phi_h, expected, expected_a_phi_prime
    real(dp) :: expected_delta_phi_h, h_s, v_parallel, x(3)
    integer :: index

    call set_log_level(-1)
    call write_chart("manufactured_radial_covariant.bc")
    inp_swi = 9
    bfac = 1.0_dp
    call magfie_thread_init()
    call read_boozer_file("manufactured_radial_covariant.bc")
    call set_s(evaluation_s)
    call init_magfie_at_s()

    do index = 1, size(theta_samples)
        x = [evaluation_s, 0.0_dp, theta_samples(index)]
        call do_magfie(x, bmod, expected, bder, hcovar, hctrvr, hcurl)
        call radial_covariant_component(x, bmod, h_s)
        ! Bphcov is in G*cm and bmod in G, so h_s is cm.
        expected = 2.0_dp*pi*Bphcov/bmod
        call assert_close("h_s", theta_samples(index), h_s, expected)
        v_parallel = (-1.0_dp)**index*1.25e8_dp
        call finite_orbit_toroidal_phase(v_parallel, h_s, a_phi_prime, delta_phi_h)
        expected_a_phi_prime = psi_pr*iota
        expected_delta_phi_h = -c*mi*v_parallel*h_s/(qi*expected_a_phi_prime)
        call assert_close("A_phi'", theta_samples(index), a_phi_prime, &
            expected_a_phi_prime)
        call assert_close("delta_phi_H", theta_samples(index), delta_phi_h, &
            expected_delta_phi_h)
        if (abs(hcovar(1)) > 0.0_dp) then
            error stop "native hcovar(1) changed on diagnostic hook branch"
        end if
    end do

    print *, "test_radial_covariant: PASS"

contains

    subroutine write_chart(path)
        character(len=*), intent(in) :: path

        integer :: file_unit, surface_index
        real(dp) :: flux, radius_m, s_surface
        real(dp) :: Bph_si, raw_Jpol
        real(dp), parameter :: surfaces(4) = [0.2_dp, 0.4_dp, 0.6_dp, 0.8_dp]
        integer, parameter :: handedness = -1

        flux = -pi*minor_radius_m**2*field_t
        Bph_si = field_t*major_radius_m
        raw_Jpol = Bph_si*tesla_to_gauss*meter_to_cm / &
            (current_to_covar*real(handedness, dp))

        open (newunit=file_unit, file=path, status="replace", action="write")
        write (file_unit, "(a)") "CC manufactured circular radial-covariant chart"
        write (file_unit, "(a)") "CC independent source-formula oracle"
        write (file_unit, "(a)") "CC m0b n0b nsurf nper flux a R"
        write (file_unit, "(a)") "CC"
        write (file_unit, "(a)") "CC"
        write (file_unit, *) 1, 0, size(surfaces), 1, flux, minor_radius_m, &
            major_radius_m

        do surface_index = 1, size(surfaces)
            s_surface = surfaces(surface_index)
            radius_m = minor_radius_m*sqrt(s_surface)
            write (file_unit, "(a)") "CC s iota Jpol Itor pprime sqrtg00"
            write (file_unit, "(a)") "CC units A A Pa m3"
            write (file_unit, *) s_surface, iota_chart, raw_Jpol, 0.0_dp, &
                0.0_dp, 0.0_dp
            write (file_unit, "(a)") &
                "CC m n rmnc rmns zmnc zmns vmnc vmns bmnc bmns"
            write (file_unit, *) 0, 0, major_radius_m, 0.0_dp, 0.0_dp, &
                0.0_dp, s_surface, 0.0_dp, field_t, 0.0_dp
            write (file_unit, *) 1, 0, radius_m, 0.0_dp, 0.0_dp, radius_m, &
                0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp
        end do
        close (file_unit)
    end subroutine write_chart

    subroutine assert_close(name, theta, found, reference)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: theta, found, reference
        real(dp) :: deviation, scale

        scale = max(abs(found), abs(reference), tiny(1.0_dp))
        deviation = abs(found - reference)/scale
        if (deviation > tolerance) then
            print *, "FAIL ", name, " at theta=", theta
            print *, "  found=", found, " reference=", reference
            print *, "  relative deviation=", deviation
            error stop "radial covariant oracle failed"
        end if
    end subroutine assert_close

end program test_radial_covariant
