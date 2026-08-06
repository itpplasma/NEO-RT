program test_boozer_angle_map
    !! Oracle for the geometric-to-Boozer angle map.
    !!
    !! The map's defining property is that it inverts the Boozer file's own
    !! forward parametrisation: taking theta_B to a physical point, measuring
    !! the geometric angle of that point about the chart axis, and asking the
    !! map for the Boozer angle must return the theta_B one started from.  That
    !! is checked directly against the file's Fourier series, independently of
    !! how the map is built, so a wrong inversion cannot pass by agreeing with
    !! itself.
    use, intrinsic :: iso_fortran_env, only: dp => real64, output_unit
    use math_constants, only: pi
    use boozer_angle_map, only: build_boozer_angle_map, boozer_angle_at, &
        boozer_angle_map_ready

    implicit none

    !> One part in 1e4 of a full turn.  The map is tabulated on 256 nodes and
    !> read back with a quintic spline, so the inversion error is set by the
    !> tabulation, not by the arithmetic.
    real(dp), parameter :: tol_theta = 2.0_dp*pi/1.0e4_dp
    integer, parameter :: n_probe = 37

    character(len=1024) :: path
    integer :: status, length
    real(dp) :: R_axis, Z_axis
    real(dp) :: worst_theta, worst_nu_jump

    call get_command_argument(1, value=path, length=length, status=status)
    if (status /= 0 .or. length == 0) then
        write (output_unit, '(a)') &
            'test_boozer_angle_map skipped: no axisymmetric Boozer file given'
        stop
    end if

    call axis_from_file(trim(path), R_axis, Z_axis)
    call build_boozer_angle_map(trim(path), R_axis, Z_axis)
    if (.not. boozer_angle_map_ready()) then
        write (output_unit, '(a)') 'FAIL: map reports not ready after build'
        error stop
    end if

    call check_round_trip(trim(path), R_axis, Z_axis, worst_theta)
    call check_nu_continuity(worst_nu_jump)

    write (output_unit, '(a,es12.4,a,es12.4)') &
        'boozer_angle_map: worst inversion error ', worst_theta, &
        ' rad, worst nu jump ', worst_nu_jump
    if (worst_theta > tol_theta) then
        write (output_unit, '(a,es12.4,a,es12.4)') &
            'FAIL: angle map does not invert the Boozer parametrisation: ', &
            worst_theta, ' >', tol_theta
        error stop
    end if

contains

    subroutine check_round_trip(file_path, R_axis_in, Z_axis_in, worst)
        character(len=*), intent(in) :: file_path
        real(dp), intent(in) :: R_axis_in, Z_axis_in
        real(dp), intent(out) :: worst

        real(dp), allocatable :: s_grid(:), mpol(:, :)
        real(dp), allocatable :: rmnc(:, :), rmns(:, :), zmnc(:, :), zmns(:, :)
        real(dp), allocatable :: vmnc(:, :), vmns(:, :)
        real(dp) :: theta_b_in, theta_geo, theta_b_out, nu_out, err
        real(dp) :: R_pt, Z_pt, nu_pt
        integer :: nflux, nmode, nper, k, i, ksel(5), kk

        call read_axisymmetric_boozer_public(file_path, s_grid, mpol, rmnc, &
                                             rmns, zmnc, zmns, vmnc, vmns, &
                                             nflux, nmode, nper)

        ksel = [nflux/8, nflux/4, nflux/2, (3*nflux)/4, (9*nflux)/10]
        worst = 0.0_dp
        do kk = 1, size(ksel)
            k = max(1, min(nflux, ksel(kk)))
            do i = 1, n_probe
                theta_b_in = 2.0_dp*pi*real(i - 1, dp)/real(n_probe, dp)
                call evaluate_surface(theta_b_in, mpol(k, :), rmnc(k, :), &
                                      rmns(k, :), zmnc(k, :), zmns(k, :), &
                                      vmnc(k, :), vmns(k, :), nmode, &
                                      R_pt, Z_pt, nu_pt)
                theta_geo = atan2(Z_pt - Z_axis_in, R_pt - R_axis_in)
                call boozer_angle_at(s_grid(k), theta_geo, theta_b_out, nu_out)
                err = abs(wrap_pi(theta_b_out - theta_b_in))
                worst = max(worst, err)
                ! The stream function is carried alongside the angle and must
                ! agree with the file at the same point, or the toroidal phase
                ! applied to the perturbation is wrong even when the poloidal
                ! angle is right.
                if (abs(nu_out - nu_pt) > 1.0e-3_dp*max(1.0_dp, abs(nu_pt))) then
                    write (output_unit, '(a,i0,a,es12.4,a,es12.4)') &
                        'FAIL: stream function mismatch on surface ', k, &
                        ': map ', nu_out, ' file ', nu_pt
                    error stop
                end if
            end do
        end do
    end subroutine check_round_trip

    subroutine check_nu_continuity(worst)
        !! nu must come back periodic: a jump across the seam would put a
        !! discontinuity into the toroidal phase of the perturbation.
        real(dp), intent(out) :: worst

        real(dp) :: theta, tb1, nu1, tb2, nu2
        integer :: i

        worst = 0.0_dp
        do i = 1, 64
            theta = 2.0_dp*pi*real(i - 1, dp)/64.0_dp
            call boozer_angle_at(0.25_dp, theta, tb1, nu1)
            call boozer_angle_at(0.25_dp, theta + 2.0_dp*pi, tb2, nu2)
            worst = max(worst, abs(nu1 - nu2))
        end do
        if (worst > 1.0e-8_dp) then
            write (output_unit, '(a,es12.4)') &
                'FAIL: stream function is not periodic in theta_geo: ', worst
            error stop
        end if
    end subroutine check_nu_continuity

    pure real(dp) function wrap_pi(angle)
        real(dp), intent(in) :: angle

        wrap_pi = modulo(angle + pi, 2.0_dp*pi) - pi
    end function wrap_pi

    subroutine evaluate_surface(theta, m, rmnc, rmns, zmnc, zmns, vmnc, vmns, &
                                nmode, R_pt, Z_pt, nu_pt)
        real(dp), intent(in) :: theta, m(:), rmnc(:), rmns(:)
        real(dp), intent(in) :: zmnc(:), zmns(:), vmnc(:), vmns(:)
        integer, intent(in) :: nmode
        real(dp), intent(out) :: R_pt, Z_pt, nu_pt

        real(dp) :: angle, cs, sn
        integer :: j

        R_pt = 0.0_dp
        Z_pt = 0.0_dp
        nu_pt = 0.0_dp
        do j = 1, nmode
            angle = m(j)*theta
            cs = cos(angle)
            sn = sin(angle)
            R_pt = R_pt + rmnc(j)*cs + rmns(j)*sn
            Z_pt = Z_pt + zmnc(j)*cs + zmns(j)*sn
            nu_pt = nu_pt + vmnc(j)*cs + vmns(j)*sn
        end do
    end subroutine evaluate_surface

    subroutine axis_from_file(file_path, R_axis_out, Z_axis_out)
        character(len=*), intent(in) :: file_path
        real(dp), intent(out) :: R_axis_out, Z_axis_out

        real(dp), allocatable :: s_grid(:), mpol(:, :)
        real(dp), allocatable :: rmnc(:, :), rmns(:, :), zmnc(:, :), zmns(:, :)
        real(dp), allocatable :: vmnc(:, :), vmns(:, :)
        integer :: nflux, nmode, nper, j

        call read_axisymmetric_boozer_public(file_path, s_grid, mpol, rmnc, &
                                             rmns, zmnc, zmns, vmnc, vmns, &
                                             nflux, nmode, nper)
        R_axis_out = 0.0_dp
        Z_axis_out = 0.0_dp
        do j = 1, nmode
            if (nint(mpol(1, j)) == 0) then
                R_axis_out = R_axis_out + rmnc(1, j)
                Z_axis_out = Z_axis_out + zmnc(1, j)
            end if
        end do
    end subroutine axis_from_file

    subroutine read_axisymmetric_boozer_public(file_path, s_grid, mpol, rmnc, &
                                               rmns, zmnc, zmns, vmnc, vmns, &
                                               nflux, nmode, nper)
        character(len=*), intent(in) :: file_path
        real(dp), allocatable, intent(out) :: s_grid(:), mpol(:, :)
        real(dp), allocatable, intent(out) :: rmnc(:, :), rmns(:, :)
        real(dp), allocatable, intent(out) :: zmnc(:, :), zmns(:, :)
        real(dp), allocatable, intent(out) :: vmnc(:, :), vmns(:, :)
        integer, intent(out) :: nflux, nmode, nper

        integer, parameter :: ncol1 = 5, ncol2 = 8
        integer :: unit, m0b, n0b, ksurf, kmode
        real(dp) :: flux, a_minor, R0
        real(dp) :: params(ncol1 + 1), row(ncol2 + 2)

        open (newunit=unit, file=file_path, action='read', status='old')
        read (unit, '(////)')
        read (unit, *) m0b, n0b, nflux, nper, flux, a_minor, R0
        nmode = (m0b + 1)*(n0b + 1)
        allocate (s_grid(nflux), mpol(nflux, nmode))
        allocate (rmnc(nflux, nmode), rmns(nflux, nmode))
        allocate (zmnc(nflux, nmode), zmns(nflux, nmode))
        allocate (vmnc(nflux, nmode), vmns(nflux, nmode))
        do ksurf = 1, nflux
            read (unit, '(/)')
            read (unit, *) params
            read (unit, *)
            do kmode = 1, nmode
                read (unit, *) row
                mpol(ksurf, kmode) = row(1)
                rmnc(ksurf, kmode) = row(3)
                rmns(ksurf, kmode) = row(4)
                zmnc(ksurf, kmode) = row(5)
                zmns(ksurf, kmode) = row(6)
                vmnc(ksurf, kmode) = row(7)
                vmns(ksurf, kmode) = row(8)
            end do
            s_grid(ksurf) = params(1)
        end do
        close (unit)
    end subroutine read_axisymmetric_boozer_public

end program test_boozer_angle_map
