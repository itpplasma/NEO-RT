program test_eqdsk_local_pitch
    !> The direct-EQDSK field-line label must carry the LOCAL pitch.
    !>
    !> In geoflux coordinates theta is the geometric poloidal angle, so
    !> alpha = phi - q*theta with the GLOBAL q is not a field-line label and the
    !> precession built from it is not the canonical one (NEO-RT #85).  Along a
    !> field line dphi/dtheta = B^phi/B^theta, so the label gradient carries
    !> hctrvr(2)/hctrvr(3), which reduces to q only when theta is already a
    !> straight-field-line angle.
    !>
    !> ORACLE.  This test does NOT compare against the tabulated q.  It compares
    !> against a field line integrated in CYLINDRICAL coordinates from field_eq
    !> alone:
    !>
    !>     dR/dphi = R B_R/B_phi ,   dZ/dphi = R B_Z/B_phi ,
    !>
    !> advanced until the geometric angle about the magnetic axis has completed
    !> one circuit, giving the rotational transform directly as Delta_phi/(2 pi).
    !> That shares no code with do_magfie_eqdsk's coordinate transformation, so it
    !> is an independent behavioral oracle rather than a restatement of the patch.
    !>
    !> Deliberately NOT the tabulated q: on the circ.eqdsk fixture the tabulated
    !> q disagrees with the rotational transform of the fixture's own field by a
    !> factor of about 1.55 (see NEO-RT #93).  A test written against q would
    !> therefore fail for a reason that has nothing to do with the field-line
    !> label, which is exactly why test_eqdsk_axis had to use a permissive
    !> factor-of-two bound.
    !>
    !> SCOPE.  This validates the substituted quantity, not the full canonical
    !> precession.  The straight-field-line angle map also depends on s, so its
    !> gradient has a radial term that the inp_swi=11 branch still omits; see
    !> NEO-RT #85.  A complete test of Om_tB has to wait for that term.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use do_magfie_mod, only: inp_swi, read_boozer_file, set_s, &
        init_magfie_at_s, do_magfie
    use geoflux_coordinates, only: geoflux_to_cyl
    use field_sub, only: field_eq
    use util, only: pi
    use util_for_test, only: pass_test, fail_test

    implicit none

    integer, parameter :: ntheta = 2048
    integer, parameter :: nstep = 400000
    !> The two constructions are independent discretizations of the same
    !> continuous quantity, so a few tenths of a percent is the honest bar; the
    !> observed agreement is 6e-4.
    real(dp), parameter :: pitch_tol = 5.0e-3_dp
    !> Below this the local pitch is indistinguishable from a constant and the
    !> test would be vacuous: an implementation returning the global q would pass.
    real(dp), parameter :: min_variation = 1.0e-2_dp

    real(dp) :: s_values(3) = [0.15_dp, 0.35_dp, 0.60_dp]
    character(len=1024) :: eqdsk_file
    real(dp) :: x(3), bmod, sqrtg, bder(3), hcov(3), hcon(3), hcurl(3)
    real(dp) :: xgeo(3), xcyl(3), jac(3, 3), raxis, zaxis
    real(dp) :: rr, zz, phi, dphi, th, thprev, unwrapped
    real(dp) :: theta, dtheta, pitch, mean_pitch, pmin, pmax
    real(dp) :: transform, relative_error, variation, s
    integer :: k, j

    call get_environment_variable('EQDSK_FILE', eqdsk_file)
    if (len_trim(eqdsk_file) == 0) then
        write (*, *) 'test_eqdsk_local_pitch: EQDSK_FILE is not set'
        call fail_test
    end if

    inp_swi = 11
    call read_boozer_file(trim(eqdsk_file))

    ! Magnetic axis from the chart, so the poloidal angle used by the oracle is
    ! measured about the same point the chart uses.
    xgeo = [1.0e-8_dp, 0.0_dp, 0.0_dp]
    call geoflux_to_cyl(xgeo, xcyl, jac)
    raxis = xcyl(1)
    zaxis = xcyl(3)

    dtheta = 2.0_dp*pi/real(ntheta, dp)

    do j = 1, size(s_values)
        s = s_values(j)
        call set_s(s)
        call init_magfie_at_s()

        ! ---- the metric's field-line pitch, theta-averaged ------------------
        mean_pitch = 0.0_dp
        pmin = huge(pmin)
        pmax = -huge(pmax)
        do k = 0, ntheta - 1
            theta = -pi + (real(k, dp) + 0.5_dp)*dtheta
            x = [s, 0.0_dp, theta]
            call do_magfie(x, bmod, sqrtg, bder, hcov, hcon, hcurl)
            if (.not. all(ieee_is_finite([bmod, sqrtg, hcon]))) then
                write (*, *) 'test_eqdsk_local_pitch: non-finite metric, s =', s
                call fail_test
            end if
            ! .and. does not short-circuit in Fortran, so guard separately.
            if (abs(hcon(3)) <= tiny(hcon(3))) then
                write (*, *) 'test_eqdsk_local_pitch: B^theta vanishes, s =', s
                call fail_test
            end if
            pitch = hcon(2)/hcon(3)
            mean_pitch = mean_pitch + pitch*dtheta/(2.0_dp*pi)
            pmin = min(pmin, pitch)
            pmax = max(pmax, pitch)
        end do

        ! ---- the independent oracle: integrate the actual field line ---------
        xgeo = [s, 0.0_dp, 0.0_dp]
        call geoflux_to_cyl(xgeo, xcyl, jac)
        rr = xcyl(1)
        zz = xcyl(3)
        phi = 0.0_dp
        dphi = 16.0_dp*pi/real(nstep, dp)
        thprev = atan2(zz - zaxis, rr - raxis)
        unwrapped = 0.0_dp
        do k = 1, nstep
            call rk4_step(rr, zz, dphi)
            phi = phi + dphi
            th = atan2(zz - zaxis, rr - raxis)
            unwrapped = unwrapped + wrapped_delta(th - thprev)
            thprev = th
            if (abs(unwrapped) >= 2.0_dp*pi) exit
        end do
        if (abs(unwrapped) < 2.0_dp*pi) then
            write (*, *) 'test_eqdsk_local_pitch: field line did not close a ', &
                'poloidal circuit, s =', s
            call fail_test
        end if
        ! Scale out the small overshoot of the final step.
        transform = phi/(2.0_dp*pi)*(2.0_dp*pi/abs(unwrapped))

        if (abs(transform) <= tiny(transform)) then
            write (*, *) 'test_eqdsk_local_pitch: zero rotational transform, s =', s
            call fail_test
        end if

        relative_error = abs(abs(mean_pitch) - abs(transform))/abs(transform)
        variation = (pmax - pmin)/abs(mean_pitch)

        write (*, '(a,f6.3,a,es12.5,a,es12.5,a,es10.3)') &
            ' s =', s, '  <B^phi/B^theta> =', mean_pitch, &
            '  field-line transform =', transform, '  rel err =', relative_error

        if (relative_error > pitch_tol) then
            write (*, *) 'test_eqdsk_local_pitch: the metric field-line pitch does'
            write (*, *) '    not reproduce the integrated rotational transform'
            call fail_test
        end if
        if (variation < min_variation) then
            write (*, *) 'test_eqdsk_local_pitch: local pitch is effectively'
            write (*, *) '    constant, so this fixture cannot distinguish it from'
            write (*, *) '    the global q; variation/mean =', variation
            call fail_test
        end if
    end do

    call pass_test

contains

    pure function wrapped_delta(d) result(w)
        real(dp), intent(in) :: d
        real(dp) :: w
        w = d
        do while (w > pi)
            w = w - 2.0_dp*pi
        end do
        do while (w < -pi)
            w = w + 2.0_dp*pi
        end do
    end function wrapped_delta

    subroutine derivs(rr_, zz_, dr, dz)
        real(dp), intent(in) :: rr_, zz_
        real(dp), intent(out) :: dr, dz
        real(dp) :: br, bf, bz
        real(dp) :: brr, brf, brz, bfr, bff, bfz, bzr, bzf, bzz
        call field_eq(rr_, 0.0_dp, zz_, br, bf, bz, &
            brr, brf, brz, bfr, bff, bfz, bzr, bzf, bzz)
        dr = rr_*br/bf
        dz = rr_*bz/bf
    end subroutine derivs

    subroutine rk4_step(rr_, zz_, h)
        real(dp), intent(inout) :: rr_, zz_
        real(dp), intent(in) :: h
        real(dp) :: k1r, k1z, k2r, k2z, k3r, k3z, k4r, k4z
        call derivs(rr_, zz_, k1r, k1z)
        call derivs(rr_ + 0.5_dp*h*k1r, zz_ + 0.5_dp*h*k1z, k2r, k2z)
        call derivs(rr_ + 0.5_dp*h*k2r, zz_ + 0.5_dp*h*k2z, k3r, k3z)
        call derivs(rr_ + h*k3r, zz_ + h*k3z, k4r, k4z)
        rr_ = rr_ + h/6.0_dp*(k1r + 2.0_dp*k2r + 2.0_dp*k3r + k4r)
        zz_ = zz_ + h/6.0_dp*(k1z + 2.0_dp*k2z + 2.0_dp*k3z + k4z)
    end subroutine rk4_step

end program test_eqdsk_local_pitch
