program test_hcurl_axisymmetric
    ! Independent oracle for curl of the unit field direction in the standalone
    ! Boozer path.
    !
    ! The routine computes the radial and poloidal components analytically from
    ! bder and dBphcovds. This test never uses those formulas. It differentiates
    ! hcovar(2) numerically, with central differences over separate do_magfie
    ! evaluations, and requires
    !
    !     hcurl(1) = -(1/sqrtg) d_theta h_phi
    !     hcurl(3) = +(1/sqrtg) d_s     h_phi
    !
    ! The chart is manufactured with a theta-varying field magnitude, because a
    ! chart with a flux-function |B| would make the radial component vanish and
    ! the test vacuous. That non-vanishing is asserted: the previous
    ! implementation returned zero for all three components, and an independent
    ! axisymmetric identity in the ITER TC24 evidence packet
    ! results/neort_x3_drift_vector_20260905 rejected that zero as a physical
    ! curl.
    !
    ! The toroidal component stays zero, because it needs the covariant radial
    ! component of the field and a Boozer file does not store it. The test pins
    ! that it is advertised as unavailable rather than silently returned.
    use iso_fortran_env, only: dp => real64
    use do_magfie_mod, only: bfac, do_magfie, hcurl_toroidal_is_available, &
        init_magfie_at_s, inp_swi, magfie_thread_init, read_boozer_file, set_s
    use logger, only: set_log_level
    use util, only: pi

    implicit none

    real(dp), parameter :: field_t = 2.0_dp
    real(dp), parameter :: major_radius_m = 5.0_dp
    real(dp), parameter :: minor_radius_m = 0.5_dp
    real(dp), parameter :: iota_abs = 0.2_dp
    real(dp), parameter :: mirror_fraction = 0.15_dp
    real(dp), parameter :: current_to_covar = 0.2_dp
    real(dp), parameter :: meter_to_cm = 100.0_dp
    real(dp), parameter :: tesla_to_gauss = 1.0e4_dp
    ! Central differences are second order in the step. Comparing against a
    ! single step would mean choosing a tolerance that hides its truncation
    ! error, so both derivatives are Richardson-extrapolated from two steps,
    ! which cancels the leading term and leaves O(h^4). The tolerance then
    ! measures the formula rather than the differencing.
    real(dp), parameter :: theta_step = 1.0e-3_dp
    real(dp), parameter :: s_step = 1.0e-3_dp
    real(dp), parameter :: tolerance = 1.0e-8_dp

    character(len=*), parameter :: chart = "manufactured_hcurl.bc"
    real(dp) :: theta_samples(5)
    integer :: index

    call set_log_level(-1)
    theta_samples = [0.0_dp, 0.7_dp, 1.9_dp, 3.4_dp, 5.1_dp]

    call write_chart(chart)
    call load_chart(chart)

    do index = 1, size(theta_samples)
        call check_radial_component(theta_samples(index))
        call check_poloidal_component(theta_samples(index))
    end do
    call check_radial_component_is_not_trivially_zero()
    call check_toroidal_component_is_declared_unavailable()

    print *, "test_hcurl_axisymmetric: PASS"

contains

    subroutine load_chart(path)
        character(len=*), intent(in) :: path

        inp_swi = 9
        bfac = 1.0_dp
        call magfie_thread_init()
        call read_boozer_file(path)
    end subroutine load_chart

    subroutine evaluate(s_value, theta, bmod, sqrtg, hcovar, hcurl)
        real(dp), intent(in) :: s_value, theta
        real(dp), intent(out) :: bmod, sqrtg
        real(dp), intent(out) :: hcovar(3), hcurl(3)

        real(dp) :: bder(3), hctrvr(3), x(3)

        call set_s(s_value)
        call init_magfie_at_s()
        x = [s_value, 0.0_dp, theta]
        call do_magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
    end subroutine evaluate

    subroutine check_radial_component(theta)
        ! hcurl(1) against a central difference of h_phi in theta.
        real(dp), intent(in) :: theta

        real(dp) :: bmod, sqrtg, hcovar(3), hcurl(3)
        real(dp) :: derivative, expected

        call evaluate(evaluation_s(), theta, bmod, sqrtg, hcovar, hcurl)
        derivative = richardson_derivative(evaluation_s(), theta, 3, theta_step)
        expected = -derivative/sqrtg
        call assert_close("hcurl(1)", theta, hcurl(1), expected, tolerance)
    end subroutine check_radial_component

    subroutine check_poloidal_component(theta)
        ! hcurl(3) against a central difference of h_phi in s.
        real(dp), intent(in) :: theta

        real(dp) :: bmod, sqrtg, hcovar(3), hcurl(3)
        real(dp) :: derivative, expected

        call evaluate(evaluation_s(), theta, bmod, sqrtg, hcovar, hcurl)
        derivative = richardson_derivative(evaluation_s(), theta, 1, s_step)
        expected = derivative/sqrtg
        call assert_close("hcurl(3)", theta, hcurl(3), expected, tolerance)
    end subroutine check_poloidal_component

    subroutine check_radial_component_is_not_trivially_zero()
        ! A chart with a flux-function |B| would satisfy the checks above with
        ! every component zero, which is what the routine used to return.
        real(dp) :: bmod, sqrtg, hcovar(3), hcurl(3), largest
        integer :: sample

        largest = 0.0_dp
        do sample = 1, size(theta_samples)
            call evaluate(evaluation_s(), theta_samples(sample), bmod, sqrtg, &
                hcovar, hcurl)
            largest = max(largest, abs(hcurl(1)))
        end do
        if (largest <= 0.0_dp) then
            print *, "FAIL: hcurl(1) vanishes on every sample; the test is vacuous"
            error stop 1
        end if
    end subroutine check_radial_component_is_not_trivially_zero

    subroutine check_toroidal_component_is_declared_unavailable()
        real(dp) :: bmod, sqrtg, hcovar(3), hcurl(3)

        if (hcurl_toroidal_is_available) then
            print *, "FAIL: hcurl(2) is advertised as available but is not computed"
            error stop 1
        end if
        call evaluate(evaluation_s(), 1.1_dp, bmod, sqrtg, hcovar, hcurl)
        if (hcurl(2) /= 0.0_dp) then
            print *, "FAIL: hcurl(2) is nonzero while advertised unavailable"
            error stop 1
        end if
    end subroutine check_toroidal_component_is_declared_unavailable

    function richardson_derivative(s_value, theta, direction, step) result(derivative)
        ! d h_phi / d(coordinate `direction`), Richardson-extrapolated from
        ! central differences at `step` and `step/2`, so the leading O(h^2)
        ! truncation term cancels.
        real(dp), intent(in) :: s_value, theta, step
        integer, intent(in) :: direction
        real(dp) :: derivative

        real(dp) :: coarse, fine

        coarse = central_difference(s_value, theta, direction, step)
        fine = central_difference(s_value, theta, direction, 0.5_dp*step)
        derivative = (4.0_dp*fine - coarse)/3.0_dp
    end function richardson_derivative

    function central_difference(s_value, theta, direction, step) result(derivative)
        real(dp), intent(in) :: s_value, theta, step
        integer, intent(in) :: direction
        real(dp) :: derivative

        real(dp) :: bmod, sqrtg, h_plus(3), h_minus(3), scratch(3)
        real(dp) :: s_plus, s_minus, theta_plus, theta_minus

        s_plus = s_value
        s_minus = s_value
        theta_plus = theta
        theta_minus = theta
        if (direction == 1) then
            s_plus = s_value + step
            s_minus = s_value - step
        else
            theta_plus = theta + step
            theta_minus = theta - step
        end if
        call evaluate(s_plus, theta_plus, bmod, sqrtg, h_plus, scratch)
        call evaluate(s_minus, theta_minus, bmod, sqrtg, h_minus, scratch)
        derivative = (h_plus(2) - h_minus(2))/(2.0_dp*step)
    end function central_difference

    pure function evaluation_s() result(value)
        real(dp) :: value
        value = 0.5_dp
    end function evaluation_s

    subroutine assert_close(name, theta, found, expected, tolerance)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: theta, found, expected, tolerance

        real(dp) :: scale, deviation

        scale = max(abs(expected), abs(found), tiny(1.0_dp))
        deviation = abs(found - expected)/scale
        if (deviation > tolerance) then
            print *, "FAIL ", name, " at theta =", theta
            print *, "  computed =", found, " finite difference =", expected
            print *, "  relative deviation =", deviation
            error stop 1
        end if
    end subroutine assert_close

    subroutine write_chart(path)
        ! Large-aspect-ratio circular chart with |B| = B0*(1 - eps*cos(theta)),
        ! so the field magnitude varies poloidally and the radial curl component
        ! is nonzero.
        character(len=*), intent(in) :: path

        integer :: file_unit, surface_index
        real(dp) :: Bph_si, Bth_si, flux, radius_m, s_surface
        real(dp) :: raw_Jpol, raw_Itor, field_magnitude
        real(dp), parameter :: surfaces(5) = &
            [0.2_dp, 0.35_dp, 0.5_dp, 0.65_dp, 0.8_dp]

        flux = pi*minor_radius_m**2*field_t
        Bph_si = field_t*major_radius_m

        open (newunit=file_unit, file=path, status="replace", action="write")
        write (file_unit, "(a)") "CC manufactured circular chart with a mirror term"
        write (file_unit, "(a)") "CC independent oracle for curl of the field direction"
        write (file_unit, "(a)") "CC m0b n0b nsurf nper flux a R"
        write (file_unit, "(a)") "CC"
        write (file_unit, "(a)") "CC"
        write (file_unit, *) 1, 0, size(surfaces), 1, flux, minor_radius_m, &
            major_radius_m

        do surface_index = 1, size(surfaces)
            s_surface = surfaces(surface_index)
            radius_m = minor_radius_m*sqrt(s_surface)
            Bth_si = field_t*radius_m**2*iota_abs/major_radius_m
            raw_Jpol = Bph_si*tesla_to_gauss*meter_to_cm/current_to_covar
            raw_Itor = Bth_si*tesla_to_gauss*meter_to_cm/current_to_covar
            field_magnitude = field_t &
                *sqrt(1.0_dp + (radius_m*iota_abs/major_radius_m)**2)
            write (file_unit, "(a)") "CC s iota Jpol Itor pprime sqrtg00"
            write (file_unit, "(a)") "CC units A A Pa m3"
            write (file_unit, *) s_surface, iota_abs, raw_Jpol, raw_Itor, &
                0.0_dp, 0.0_dp
            write (file_unit, "(a)") &
                "CC m n rmnc rmns zmnc zmns vmnc vmns bmnc bmns"
            write (file_unit, *) 0, 0, major_radius_m, 0.0_dp, 0.0_dp, 0.0_dp, &
                0.0_dp, 0.0_dp, field_magnitude, 0.0_dp
            ! The m=1 bmnc term is what makes |B| vary with theta.
            write (file_unit, *) 1, 0, radius_m, 0.0_dp, 0.0_dp, radius_m, &
                0.0_dp, 0.0_dp, -mirror_fraction*field_magnitude*sqrt(s_surface), &
                0.0_dp
        end do
        close (file_unit)
    end subroutine write_chart

end program test_hcurl_axisymmetric
