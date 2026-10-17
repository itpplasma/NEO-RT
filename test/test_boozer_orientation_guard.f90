program test_boozer_orientation_guard
    ! Usage: test_boozer_orientation_guard.x accept|flip|scale|drift file inp_swi
    !   accept: read the file; it must load (prints "accepted").
    !   flip:   write a copy with the header toroidal flux negated and read it;
    !           the loader must stop with "Inconsistent Boozer orientation".
    ! The oracle is the file's own (dV/ds)/nper column, computed by the
    ! equilibrium code independently of NEO-RT's orientation convention.
    use iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use do_magfie_mod, only: Bphcov, bfac, do_magfie, init_magfie_at_s, inp_swi, &
        magfie_thread_init, modes0, nmode, q, read_boozer_file, set_s, spl_coeff2
    use logger, only: set_log_level
    use neort_orbit, only: timestep
    use spline, only: spline_val_0

    implicit none

    character(len=1024) :: mode, path, swi_arg
    character(len=64) :: flipped_path

    call get_command_argument(1, mode)
    call get_command_argument(2, path)
    call get_command_argument(3, swi_arg)
    read (swi_arg, *) inp_swi
    ! One copy per input switch, so concurrent ctest cases do not share a file.
    write (flipped_path, "(a,i0,a)") "orientation_guard_flipped_swi", inp_swi, ".bc"
    call set_log_level(-1)
    bfac = 1.0_dp
    call magfie_thread_init()

    select case (trim(mode))
    case ("accept")
        call read_boozer_file(trim(path))
        print *, "accepted ", trim(path)
    case ("flip")
        call write_flux_changed_copy(trim(path), trim(flipped_path), -1.0_dp)
        call read_boozer_file(trim(flipped_path))
        error stop "flux-flipped Boozer file was accepted"
    case ("scale")
        call write_flux_changed_copy(trim(path), trim(flipped_path), 1.2_dp)
        call read_boozer_file(trim(flipped_path))
        error stop "flux-scaled Boozer file was accepted"
    case ("drift")
        call read_boozer_file(trim(path))
        call check_deeply_trapped_drift()
    case default
        error stop "usage: test_boozer_orientation_guard.x accept|flip|scale|drift file inp_swi"
    end select

contains

    subroutine write_flux_changed_copy(source, target, flux_factor)
        character(len=*), intent(in) :: source, target
        real(dp), intent(in) :: flux_factor

        character(len=4096) :: line
        integer :: in_unit, out_unit, line_number, ios
        integer :: m0b, n0b, nflux, nfp
        real(dp) :: flux, minor_radius, major_radius

        open (newunit=in_unit, file=source, status="old", action="read")
        open (newunit=out_unit, file=target, status="replace", action="write")
        line_number = 0
        do
            read (in_unit, "(a)", iostat=ios) line
            if (ios /= 0) exit
            line_number = line_number + 1
            if (line_number == 6) then
                read (line, *) m0b, n0b, nflux, nfp, flux, minor_radius, &
                    major_radius
                write (out_unit, "(4(1x,i5),3(1x,es24.16))") m0b, n0b, nflux, &
                    nfp, flux_factor*flux, minor_radius, major_radius
            else
                write (out_unit, "(a)") trim(line)
            end if
        end do
        close (in_unit)
        close (out_unit)
    end subroutine write_flux_changed_copy

    subroutine check_deeply_trapped_drift()
        ! At the outboard midplane, grad-B drift has v_Z proportional to
        ! -B_phi*dB/dR. The field-line label alpha=phi-q*theta therefore has
        ! alpha_dot with sign q*B_phi*(dB/dR)/(dZ/dtheta). This independent
        ! cylindrical oracle uses the file's R,Z geometry, never its flux sign.
        real(dp), parameter :: surface = 0.3_dp
        real(dp) :: bmod, sqrtg, bder(3), hcovar(3), hctrvr(3), hcurl(3)
        real(dp) :: x(3), y(3), ydot(3), drds, dzdtheta, expected_sign

        if (inp_swi /= 8) error stop "circular drift oracle requires inp_swi=8"
        call set_s(surface)
        call init_magfie_at_s()
        x = [surface, 0.0_dp, 0.0_dp]
        call do_magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
        call circular_geometry_derivatives(surface, drds, dzdtheta)
        if (drds == 0.0_dp .or. dzdtheta == 0.0_dp) error stop "degenerate drift oracle"
        expected_sign = q*Bphcov*bder(1)/(drds*dzdtheta)
        y = 0.0_dp
        call timestep(1.0e8_dp, 1.0_dp/bmod, 3, 0.0_dp, y, ydot)
        print *, "deeply trapped grad-B sign oracle:", expected_sign, ydot(3)
        if (.not. ieee_is_finite(expected_sign)) error stop "nonfinite drift oracle"
        if (.not. ieee_is_finite(ydot(3))) error stop "nonfinite trapped precession"
        if (expected_sign*ydot(3) <= 0.0_dp) error stop "wrong trapped precession sign"
    end subroutine check_deeply_trapped_drift

    subroutine circular_geometry_derivatives(surface, drds, dzdtheta)
        real(dp), intent(in) :: surface
        real(dp), intent(out) :: drds, dzdtheta

        real(dp), parameter :: ds = 1.0e-5_dp
        real(dp) :: val(3), rminus, rplus
        integer :: j

        rminus = 0.0_dp
        rplus = 0.0_dp
        dzdtheta = 0.0_dp
        do j = 1, nmode
            val = spline_val_0(spl_coeff2(:, :, 1, j), surface - ds)
            rminus = rminus + val(1)
            val = spline_val_0(spl_coeff2(:, :, 1, j), surface + ds)
            rplus = rplus + val(1)
            val = spline_val_0(spl_coeff2(:, :, 2, j), surface)
            dzdtheta = dzdtheta + modes0(1, j, 1)*val(1)
        end do
        drds = (rplus - rminus)/(2.0_dp*ds)
    end subroutine circular_geometry_derivatives

end program test_boozer_orientation_guard
