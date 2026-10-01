program test_eqdsk_field
    ! Direct GEQDSK field (inp_swi=11) against the closed forms of the analytic
    ! circular equilibrium written by test/gen_analytic_eqdsk.py (circ):
    ! chart map, |B|, local pitch and its radial derivative, q, toroidal flux.
    ! curl(h) is checked against a centred difference of the returned covariant
    ! components in the chart, which does not use the cylindrical curl formula.
    use iso_fortran_env, only: dp => real64
    use do_magfie_mod, only: inp_swi, bfac, read_boozer_file, do_magfie, q, &
        psi_pr, sign_theta
    use neort_eqdsk_field, only: eqdsk_local_pitch, eqdsk_dpitch_ds
    use geoflux_coordinates, only: geoflux_to_cyl
    use util, only: pi

    implicit none

    real(dp), parameter :: R0 = 165.0_dp, a = 50.0_dp, B0 = 2.0e4_dp
    real(dp), parameter :: Q0 = 1.2_dp, QA = 3.5_dp
    real(dp), parameter :: C = (QA - Q0)/a**2
    real(dp), parameter :: s_test(3) = [0.2_dp, 0.5_dp, 0.8_dp]
    integer, parameter :: ntheta = 8
    character(len=1024) :: path
    integer :: is, it, nfail
    real(dp) :: s, theta, rho, R, err_max(9)

    call get_environment_variable("EQDSK_CIRC_FILE", path)
    if (len_trim(path) == 0) error stop "EQDSK_CIRC_FILE not set"
    inp_swi = 11
    bfac = 1.0_dp
    call read_boozer_file(trim(path))

    err_max = 0.0_dp
    do is = 1, size(s_test)
        s = s_test(is)
        rho = rho_of_s(s)
        do it = 0, ntheta - 1
            theta = 2.0_dp*pi*it/ntheta + 0.1_dp
            R = R0 + rho*cos(theta)
            call check_point(s, theta, rho, R, err_max)
        end do
        call check_surface(s, rho, err_max)
    end do

    nfail = 0
    call report("chart map |dx|/a", err_max(1), 2.0e-4_dp, nfail)
    call report("|B| rel", err_max(2), 2.0e-4_dp, nfail)
    call report("local pitch rel", err_max(3), 2.0e-4_dp, nfail)
    ! The radial derivatives inherit the error of libneo's s(psi) table, which
    ! is ~1e-3 in d/ds although positions agree to 3e-5.
    call report("d pitch/ds rel", err_max(4), 3.0e-3_dp, nfail)
    call report("q rel", err_max(5), 2.0e-4_dp, nfail)
    call report("toroidal flux rel", err_max(6), 1.0e-4_dp, nfail)
    call report("<sqrt(g) B^phi> / sign_theta psi_pr - 1", err_max(9), 3.0e-3_dp, nfail)
    call report("curl h vs chart difference", err_max(7), 2.0e-3_dp, nfail)
    call report("h_i h^i - 1, h^s", err_max(8), 1.0e-6_dp, nfail)
    if (nfail > 0) error stop "test_eqdsk_field failed"
    print *, "PASS test_eqdsk_field"

contains

    pure real(dp) function rho_of_s(s_val)
        ! Toroidal flux per radian F*(R0 - sqrt(R0**2 - rho**2)), normalized.
        real(dp), intent(in) :: s_val
        real(dp) :: d

        d = R0 - sqrt(R0**2 - a**2)
        rho_of_s = sqrt(R0**2 - (R0 - s_val*d)**2)
    end function rho_of_s

    pure real(dp) function qm(r)
        real(dp), intent(in) :: r
        qm = Q0 + C*r**2
    end function qm

    subroutine check_point(s_val, th, r, Rmaj, err)
        real(dp), intent(in) :: s_val, th, r, Rmaj
        real(dp), intent(inout) :: err(:)
        real(dp) :: x(3), xcyl(3), bmod, sqrtg, bder(3), hcov(3), hcon(3), hcurl(3)
        real(dp) :: b_exact, p_exact, dp_exact, d, xgeo(3)

        x = [s_val, 0.3_dp, th]
        xgeo = [s_val, th, 0.0_dp]
        call geoflux_to_cyl(xgeo, xcyl)
        err(1) = max(err(1), hypot(xcyl(1) - Rmaj, xcyl(3) - r*sin(th))/a)
        call do_magfie(x, bmod, sqrtg, bder, hcov, hcon, hcurl)
        b_exact = sqrt((R0*B0)**2 + (B0*r/qm(r))**2)/Rmaj
        err(2) = max(err(2), abs(bmod/b_exact - 1.0_dp))
        p_exact = qm(r)*R0/Rmaj
        err(3) = max(err(3), abs(abs(hcon(2)/hcon(3))/p_exact - 1.0_dp))
        d = R0 - sqrt(R0**2 - a**2)
        dp_exact = R0*(2.0_dp*C*r*Rmaj - qm(r)*cos(th))/Rmaj**2 &
            *(R0 - s_val*d)*d/r
        err(4) = max(err(4), abs(abs(eqdsk_dpitch_ds(s_val, th)) &
            /abs(dp_exact) - 1.0_dp))
        err(7) = max(err(7), curl_mismatch(x, hcurl))
        err(8) = max(err(8), abs(sum(hcov*hcon) - 1.0_dp), abs(hcon(1))*a)
        if (sign(1.0_dp, q) /= sign(1.0_dp, eqdsk_local_pitch(s_val, th))) then
            err(5) = huge(1.0_dp)
        end if
    end subroutine check_point

    subroutine check_surface(s_val, r, err)
        ! q and the toroidal flux per radian, sign_theta*psi_pr being the theta
        ! average of sqrtg*B^phi with sqrtg the (s, theta, phi) Jacobian.
        real(dp), intent(in) :: s_val, r
        real(dp), intent(inout) :: err(:)
        integer, parameter :: nth = 256
        real(dp) :: bmod, sqrtg, bder(3), hcov(3), hcon(3), hcurl(3), avg, flux
        real(dp) :: x(3)
        integer :: k

        avg = 0.0_dp
        do k = 0, nth - 1
            x = [s_val, 0.0_dp, 2.0_dp*pi*k/nth]
            call do_magfie(x, bmod, sqrtg, bder, hcov, hcon, hcurl)
            avg = avg + sqrtg*hcon(2)*bmod/nth
        end do
        err(5) = max(err(5), abs(abs(q)*sqrt(1.0_dp - (r/R0)**2)/qm(r) - 1.0_dp))
        flux = R0*B0*(R0 - sqrt(R0**2 - a**2))
        err(6) = max(err(6), abs(abs(sign_theta*psi_pr)/flux - 1.0_dp))
        err(9) = max(err(9), abs(avg/(sign_theta*psi_pr) - 1.0_dp))
    end subroutine check_surface

    real(dp) function curl_mismatch(x, hcurl) result(mismatch)
        ! (curl h)^i = eps^ijk d_j h_k / sqrt(g) in (s, phi, theta); axisymmetry
        ! removes the phi derivatives.
        real(dp), intent(in) :: x(3), hcurl(3)
        real(dp), parameter :: ds = 1.0e-4_dp, dth = 1.0e-4_dp
        real(dp) :: hs_p(3), hs_m(3), ht_p(3), ht_m(3), sqrtg, curl_fd(3)

        call covariant_h(x + [ds, 0.0_dp, 0.0_dp], hs_p, sqrtg)
        call covariant_h(x - [ds, 0.0_dp, 0.0_dp], hs_m, sqrtg)
        call covariant_h(x + [0.0_dp, 0.0_dp, dth], ht_p, sqrtg)
        call covariant_h(x - [0.0_dp, 0.0_dp, dth], ht_m, sqrtg)
        call covariant_h(x, curl_fd, sqrtg)
        curl_fd(1) = -(ht_p(2) - ht_m(2))/(2.0_dp*dth)/sqrtg
        curl_fd(2) = ((ht_p(1) - ht_m(1))/(2.0_dp*dth) &
            - (hs_p(3) - hs_m(3))/(2.0_dp*ds))/sqrtg
        curl_fd(3) = (hs_p(2) - hs_m(2))/(2.0_dp*ds)/sqrtg
        ! Compare physical magnitudes: |v|^2 = g_ij v^i v^j ~ (a v^s)^2 + ...
        mismatch = maxval(abs(curl_fd - hcurl)*[a, R0, a]) &
            /maxval(abs(hcurl)*[a, R0, a])
    end function curl_mismatch

    subroutine covariant_h(x_in, hcov, sqrtg)
        ! sqrtg is returned as the (s, phi, theta) Jacobian.
        real(dp), intent(in) :: x_in(3)
        real(dp), intent(out) :: hcov(3), sqrtg
        real(dp) :: bmod, bder(3), hcon(3), hcurl(3), x(3)

        x = x_in
        call do_magfie(x, bmod, sqrtg, bder, hcov, hcon, hcurl)
        sqrtg = -sqrtg
    end subroutine covariant_h

    subroutine report(label, value, tol, nfail)
        character(len=*), intent(in) :: label
        real(dp), intent(in) :: value, tol
        integer, intent(inout) :: nfail

        if (value <= tol) then
            print '(a,": ",es10.3," <= ",es9.2)', label, value, tol
        else
            print '("FAIL ",a,": ",es10.3," > ",es9.2)', label, value, tol
            nfail = nfail + 1
        end if
    end subroutine report

end program test_eqdsk_field
