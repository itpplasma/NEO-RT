program probe_neort_eqdsk_coordinate
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use do_magfie_mod, only: inp_swi, bfac, read_boozer_file, set_s, &
        init_magfie_at_s, do_magfie, q, iota, psi_pr, R0, a
    use geoflux_coordinates, only: geoflux_to_cyl
    use field_sub, only: field_eq

    implicit none

    integer, parameter :: nsample = 8
    character(len=1024) :: eqdsk_file
    real(dp), parameter :: sample_s(nsample) = &
        [0.10_dp, 0.10_dp, 0.35_dp, 0.35_dp, 0.65_dp, 0.65_dp, 0.90_dp, 0.90_dp]
    real(dp), parameter :: sample_phi(nsample) = &
        [0.00_dp, 0.71_dp, 1.43_dp, 2.14_dp, 2.86_dp, 3.57_dp, 4.29_dp, 5.00_dp]
    real(dp), parameter :: sample_theta(nsample) = &
        [0.00_dp, 0.79_dp, 1.57_dp, 2.36_dp, 3.14_dp, 3.93_dp, 4.71_dp, 5.50_dp]

    real(dp) :: x(3), xgeo(3), xcyl_cm(3), jac_cm(3, 3)
    real(dp) :: bmod, sqrtg, bder(3), hcov(3), hcon(3), hcurl(3)
    real(dp) :: b_native_g(3), b_from_hcon_g(3)
    real(dp) :: derivs(12), det_rz_cm2, sqrtg_oriented_cm3
    real(dp) :: phi, cphi, sphi
    integer :: k, status

    call get_command_argument(1, eqdsk_file, status=status)
    if (status /= 0 .or. len_trim(eqdsk_file) == 0) then
        error stop 'usage: probe_neort_eqdsk_coordinate.x EQDSK_FILE'
    end if

    inp_swi = 11
    bfac = 1.0_dp
    call read_boozer_file(trim(eqdsk_file))

    write(*, '(a,1x,3(es24.16e3,1x))') 'NEORT_EQDSK_META', R0, a, bfac
    do k = 1, nsample
        x = [sample_s(k), sample_phi(k), sample_theta(k)]
        xgeo = [x(1), x(3), x(2)]
        call geoflux_to_cyl(xgeo, xcyl_cm, jac_cm)

        call set_s(x(1))
        call init_magfie_at_s()
        call do_magfie(x, bmod, sqrtg, bder, hcov, hcon, hcurl)

        call evaluate_field(xcyl_cm(1), xcyl_cm(2), xcyl_cm(3), b_native_g, derivs)

        phi = xcyl_cm(2)
        cphi = cos(phi)
        sphi = sin(phi)
        b_native_g = [b_native_g(1)*cphi - b_native_g(2)*sphi, &
            b_native_g(1)*sphi + b_native_g(2)*cphi, b_native_g(3)]
        b_from_hcon_g = bmod*( &
            hcon(1)*[jac_cm(1, 1)*cphi, jac_cm(1, 1)*sphi, jac_cm(3, 1)] &
            + hcon(2)*[-xcyl_cm(1)*sphi, xcyl_cm(1)*cphi, 0.0_dp] &
            + hcon(3)*[jac_cm(1, 2)*cphi, jac_cm(1, 2)*sphi, jac_cm(3, 2)])

        det_rz_cm2 = jac_cm(1, 1)*jac_cm(3, 2) - jac_cm(1, 2)*jac_cm(3, 1)
        sqrtg_oriented_cm3 = xcyl_cm(1)*det_rz_cm2

        write(*, '(a,1x,i0,1x,*(es24.16e3,1x))') 'NEORT_EQDSK_SAMPLE', k, &
            x, xcyl_cm, b_native_g, b_from_hcon_g, bmod, sqrtg, &
            sqrtg_oriented_cm3, q, iota, psi_pr, hcov, hcon
    end do

contains

    subroutine evaluate_field(r, phi_in, z, b_cyl, field_derivs)
        real(dp), intent(in) :: r, phi_in, z
        real(dp), intent(out) :: b_cyl(3), field_derivs(12)

        call field_eq(r, phi_in, z, b_cyl(1), b_cyl(2), b_cyl(3), &
            field_derivs(1), field_derivs(2), field_derivs(3), &
            field_derivs(4), field_derivs(5), field_derivs(6), &
            field_derivs(7), field_derivs(8), field_derivs(9))
        field_derivs(10:) = 0.0_dp
    end subroutine evaluate_field

end program probe_neort_eqdsk_coordinate
