module neort_eqdsk_field
    !! Axisymmetric field read directly from a GEQDSK (inp_swi = 11).
    !!
    !! The chart is libneo's geoflux chart: normalized toroidal flux s, the
    !! geometric poloidal angle theta about the magnetic axis and the
    !! cylindrical toroidal angle phi.  NEO-RT orders it (s, phi, theta), which
    !! at the outboard midplane has the orientation of (R, phi, Z), so the
    !! Jacobian is positive.  The field itself is libneo's cylindrical field_eq
    !! of the same file; everything here is in CGS and without bfac.
    !!
    !! Shared read-only state is set once by init_eqdsk_field on the main thread.
    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none
    private

    public :: init_eqdsk_field, eqdsk_field, eqdsk_local_pitch
    public :: eqdsk_dpitch_ds, eqdsk_flux_profiles, eqdsk_axis, eqdsk_minor_radius
    public :: eqdsk_flux_sign

    ! Sign of the local pitch B^phi/B^theta and of B^phi, fixed by the file.
    real(dp), save :: pitch_sign = 1.0_dp
    real(dp), save :: toroidal_field_sign = 1.0_dp
    ! Radial step of the centred difference for d(pitch)/ds.
    real(dp), parameter :: pitch_ds = 1.0e-4_dp

contains

    subroutine init_eqdsk_field(path)
        use field_eq_mod, only: reset_field_eq_state, use_fpol, nwindow_r, nwindow_z
        use geoflux_coordinates, only: init_geoflux_coordinates
        use input_files, only: gfile, ieqfile

        character(len=*), intent(in) :: path
        real(dp) :: bmod, sqrtg, bder(3), hcovar(3), hctrvr(3), hcurl(3)

        call reset_field_eq_state()
        gfile = trim(path)
        ieqfile = 1
        use_fpol = .true.
        nwindow_r = 0
        nwindow_z = 0
        call init_geoflux_coordinates(path)

        ! The first field_eq call reads and splines the file; do it here so that
        ! worker threads only ever read the shared tables.
        call eqdsk_field([0.5_dp, 0.0_dp, 0.0_dp], bmod, sqrtg, bder, hcovar, &
                         hctrvr, hcurl)
        if (abs(hctrvr(3)) <= tiny(1.0_dp)) then
            error stop "GEQDSK field has no poloidal component at s=0.5"
        end if
        pitch_sign = sign(1.0_dp, hctrvr(2)/hctrvr(3))
        toroidal_field_sign = sign(1.0_dp, hctrvr(2))
    end subroutine init_eqdsk_field

    subroutine eqdsk_field(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
        !! Field in the chart at x = (s, phi, theta): |B|, Jacobian, covariant
        !! d(ln B), covariant and contravariant h = B/|B|, contravariant curl h.
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: bmod, sqrtg
        real(dp), dimension(3), intent(out) :: bder, hcovar, hctrvr, hcurl

        real(dp) :: R, jac(3, 3), bder_cyl(3), hcov_cyl(3), hcon_cyl(3)
        real(dp) :: hcurl_cyl(3)

        call cylindrical_field(x, R, jac, bmod, bder_cyl, hcov_cyl, hcon_cyl, &
                               hcurl_cyl)
        call covariant_to_chart(jac, bder_cyl, bder)
        call covariant_to_chart(jac, hcov_cyl, hcovar)
        call contravariant_to_chart(jac, hcon_cyl, hctrvr)
        call contravariant_to_chart(jac, hcurl_cyl, hcurl)
        sqrtg = R*meridional_det(jac)
    end subroutine eqdsk_field

    subroutine cylindrical_field(x, R, jac, bmod, bder_cyl, hcov_cyl, hcon_cyl, &
                                 hcurl_cyl)
        !! Cylindrical (R, phi, Z) components at the chart point x.  With
        !! sqrt(g) = R and h_phi = R*h^phi, (curl h)^i = eps^ijk d_j h_k / R and
        !! d_j h_k = d_j B_k / B - h_k d_j(ln B).
        use field_sub, only: field_eq
        use geoflux_coordinates, only: geoflux_to_cyl

        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: R, jac(3, 3), bmod
        real(dp), dimension(3), intent(out) :: bder_cyl, hcov_cyl, hcon_cyl, hcurl_cyl

        real(dp) :: xcyl(3), br, bf, bz, brr, brf, brz, bfr, bff, bfz, bzr, bzf, bzz
        real(dp) :: h(3)

        ! libneo geoflux order is (s, theta, phi).
        call geoflux_to_cyl([x(1), x(3), x(2)], xcyl, jac)
        R = xcyl(1)
        call field_eq(xcyl(1), xcyl(2), xcyl(3), br, bf, bz, &
                      brr, brf, brz, bfr, bff, bfz, bzr, bzf, bzz)
        bmod = sqrt(br**2 + bf**2 + bz**2)
        h = [br, bf, bz]/bmod
        bder_cyl(1) = (brr*h(1) + bfr*h(2) + bzr*h(3))/bmod
        bder_cyl(2) = (brf*h(1) + bff*h(2) + bzf*h(3))/bmod
        bder_cyl(3) = (brz*h(1) + bfz*h(2) + bzz*h(3))/bmod
        hcov_cyl = [h(1), R*h(2), h(3)]
        hcon_cyl = [h(1), h(2)/R, h(3)]
        ! d_phi h_Z - d_Z h_phi, with d_Z(R B_phi) = R dB_phi/dZ
        hcurl_cyl(1) = ((bzf - R*bfz)/bmod - hcov_cyl(3)*bder_cyl(2) &
                        + hcov_cyl(2)*bder_cyl(3))/R
        ! d_Z h_R - d_R h_Z
        hcurl_cyl(2) = ((brz - bzr)/bmod - hcov_cyl(1)*bder_cyl(3) &
                        + hcov_cyl(3)*bder_cyl(1))/R
        ! d_R h_phi - d_phi h_R, with d_R(R B_phi) = B_phi + R dB_phi/dR
        hcurl_cyl(3) = ((bf + R*bfr - brf)/bmod - hcov_cyl(2)*bder_cyl(1) &
                        + hcov_cyl(1)*bder_cyl(2))/R
    end subroutine cylindrical_field

    pure real(dp) function meridional_det(jac)
        !! det d(R, Z)/d(s, theta); jac(i, j) = d xcyl(i) / d xgeo(j).
        real(dp), intent(in) :: jac(3, 3)

        meridional_det = jac(1, 1)*jac(3, 2) - jac(1, 2)*jac(3, 1)
    end function meridional_det

    pure subroutine covariant_to_chart(jac, v_cyl, v)
        !! v_j = sum_i v_cyl,i d xcyl(i)/d x(j), reordered to (s, phi, theta).
        real(dp), intent(in) :: jac(3, 3), v_cyl(3)
        real(dp), intent(out) :: v(3)

        v(1) = jac(1, 1)*v_cyl(1) + jac(3, 1)*v_cyl(3)
        v(2) = v_cyl(2)
        v(3) = jac(1, 2)*v_cyl(1) + jac(3, 2)*v_cyl(3)
    end subroutine covariant_to_chart

    pure subroutine contravariant_to_chart(jac, v_cyl, v)
        !! v^j = sum_i d x(j)/d xcyl(i) v_cyl^i via the inverse meridional map.
        real(dp), intent(in) :: jac(3, 3), v_cyl(3)
        real(dp), intent(out) :: v(3)
        real(dp) :: det

        det = meridional_det(jac)
        v(1) = (jac(3, 2)*v_cyl(1) - jac(1, 2)*v_cyl(3))/det
        v(2) = v_cyl(2)
        v(3) = (-jac(3, 1)*v_cyl(1) + jac(1, 1)*v_cyl(3))/det
    end subroutine contravariant_to_chart

    real(dp) function eqdsk_local_pitch(s, theta) result(pitch)
        !! Field-line pitch d(phi)/d(theta) = B^phi/B^theta at fixed s.
        real(dp), intent(in) :: s, theta
        real(dp) :: bmod, sqrtg, bder(3), hcovar(3), hctrvr(3), hcurl(3)

        call eqdsk_field([s, 0.0_dp, theta], bmod, sqrtg, bder, hcovar, hctrvr, &
                         hcurl)
        if (abs(hctrvr(3)) <= tiny(1.0_dp)) then
            error stop "direct GEQDSK: B^theta vanishes, field-line pitch undefined"
        end if
        pitch = hctrvr(2)/hctrvr(3)
    end function eqdsk_local_pitch

    real(dp) function eqdsk_dpitch_ds(s, theta) result(dpitch)
        !! Radial derivative of the local pitch at fixed geometric theta,
        !! centred difference (one-sided at the ends of the s range).
        real(dp), intent(in) :: s, theta
        real(dp) :: s_lo, s_hi

        s_lo = max(s - pitch_ds, pitch_ds)
        s_hi = min(s + pitch_ds, 1.0_dp)
        dpitch = (eqdsk_local_pitch(s_hi, theta) - eqdsk_local_pitch(s_lo, theta)) &
                 /(s_hi - s_lo)
    end function eqdsk_dpitch_ds

    subroutine eqdsk_flux_profiles(s, q, dqds, psi_tor)
        !! Safety factor signed like the local pitch, so that a passing orbit
        !! advances phi by 2*pi*q per poloidal turn, and the toroidal flux per
        !! radian signed like B^phi, i.e. the theta average of sqrt(g)*B^phi
        !! with the positive (s, phi, theta) chart Jacobian.
        use geoflux_coordinates, only: geoflux_get_flux_profiles

        real(dp), intent(in) :: s
        real(dp), intent(out) :: q, dqds, psi_tor
        real(dp) :: q_file, dq_file, psi_pol, dpsi_pol, psi_tor_edge

        call geoflux_get_flux_profiles(s, q_file, dq_file, psi_pol, dpsi_pol, &
                                       psi_tor_edge)
        q = pitch_sign*abs(q_file)
        dqds = pitch_sign*sign(1.0_dp, q_file)*dq_file
        psi_tor = toroidal_field_sign*abs(psi_tor_edge)
    end subroutine eqdsk_flux_profiles

    real(dp) function eqdsk_flux_sign()
        !! Sign of B^phi, and of sqrt(g)*B^phi in the (s, phi, theta) chart.
        eqdsk_flux_sign = toroidal_field_sign
    end function eqdsk_flux_sign

    subroutine eqdsk_axis(R_axis, Z_axis)
        use geoflux_coordinates, only: geoflux_get_axis

        real(dp), intent(out) :: R_axis, Z_axis

        call geoflux_get_axis(R_axis, Z_axis)
    end subroutine eqdsk_axis

    real(dp) function eqdsk_minor_radius() result(a)
        !! Outboard midplane distance of the boundary surface from the axis.
        use geoflux_coordinates, only: geoflux_to_cyl

        real(dp) :: xcyl(3), R_axis, Z_axis

        call eqdsk_axis(R_axis, Z_axis)
        call geoflux_to_cyl([1.0_dp, 0.0_dp, 0.0_dp], xcyl)
        a = abs(xcyl(1) - R_axis)
        if (a <= 0.0_dp) error stop "GEQDSK has no usable minor radius"
    end function eqdsk_minor_radius

end module neort_eqdsk_field
