module boozer_angle_map
    !! Boozer angles as functions of the geometric angle of the geoflux chart.
    !!
    !! A perturbation stored as Boozer-Fourier coefficients is a function of
    !! (s, theta_B, phi_B).  An orbit integrated in the direct-GEQDSK chart
    !! reports (s, theta_geo, phi_geom) instead, and the two poloidal angles
    !! differ by the straight-field-line shift while the two toroidal angles
    !! differ by the Boozer stream function.  Evaluating the Fourier sum
    !! against theta_geo is not an approximation, it is the wrong series, so
    !! the mapping has to be supplied before the two can be combined.
    !!
    !! The axisymmetric Boozer file parametrises each surface *by* theta_B, so
    !! the forward map theta_B -> (R, Z, nu) is a Fourier sum and the inverse
    !! needs no root finding: sampling theta_B on a uniform grid yields
    !! theta_geo at each node, and theta_geo is monotonic in theta_B on a
    !! star-shaped surface.  What is tabulated is the difference
    !! theta_B - theta_geo, which is 2*pi-periodic and smooth, rather than the
    !! angle itself, which is not.
    !!
    !! The geometric angle is measured about the axis reported by the geoflux
    !! chart, not about the Boozer file's own axis, so that the abscissa here
    !! is the same quantity the orbit hands over.

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use math_constants, only: pi
    use interpolate, only: BatchSplineData2D, construct_batch_splines_2d, &
        destroy_batch_splines_2d, evaluate_batch_splines_2d

    implicit none
    private

    public :: build_boozer_angle_map, boozer_angle_at, boozer_angle_map_ready
    public :: destroy_boozer_angle_map, boozer_angle_map_nper

    integer, parameter :: n_theta_map = 256

    logical, save :: map_ready = .false.
    integer, save :: map_nper = 1
    real(dp), save :: map_s_min = 0.0_dp, map_s_max = 1.0_dp
    type(BatchSplineData2D), save :: map_spline

contains

    pure logical function boozer_angle_map_ready()
        boozer_angle_map_ready = map_ready
    end function boozer_angle_map_ready

    pure integer function boozer_angle_map_nper()
        boozer_angle_map_nper = map_nper
    end function boozer_angle_map_nper

    subroutine destroy_boozer_angle_map()
        if (map_ready) then
            call destroy_batch_splines_2d(map_spline)
            map_ready = .false.
        end if
    end subroutine destroy_boozer_angle_map

    subroutine build_boozer_angle_map(path, R_axis, Z_axis)
        !! Read an axisymmetric Boozer file and tabulate theta_B - theta_geo
        !! and nu against (s, theta_geo).  R_axis and Z_axis define the
        !! geometric angle and must come from the chart the orbit uses.
        character(len=*), intent(in) :: path
        real(dp), intent(in) :: R_axis, Z_axis

        real(dp), allocatable :: s_grid(:), mpol(:, :)
        real(dp), allocatable :: rmnc(:, :), rmns(:, :), zmnc(:, :), zmns(:, :)
        real(dp), allocatable :: vmnc(:, :), vmns(:, :)
        real(dp), allocatable :: data(:, :, :)
        real(dp) :: theta_b(n_theta_map), theta_geo(n_theta_map)
        real(dp) :: shift(n_theta_map), nu(n_theta_map)
        real(dp) :: theta_uniform, delta, nu_at
        integer :: nflux, nmode, k, j, i

        call read_axisymmetric_boozer(path, s_grid, mpol, rmnc, rmns, &
                                      zmnc, zmns, vmnc, vmns, nflux, nmode, &
                                      map_nper)

        allocate (data(nflux, n_theta_map + 1, 2))

        do k = 1, nflux
            do i = 1, n_theta_map
                theta_b(i) = 2.0_dp*pi*real(i - 1, dp)/real(n_theta_map, dp)
            end do
            call surface_samples(theta_b, mpol(k, :), rmnc(k, :), rmns(k, :), &
                                 zmnc(k, :), zmns(k, :), vmnc(k, :), vmns(k, :), &
                                 nmode, R_axis, Z_axis, theta_geo, nu)

            ! theta_geo is monotonic in theta_B on a star-shaped surface, so the
            ! inverse is an interpolation of the periodic difference rather than
            ! a root find.  Fail loudly if the surface is not star-shaped about
            ! this axis: silently interpolating a folded abscissa would place
            ! the perturbation on the wrong part of the surface.
            call require_monotonic(theta_geo, k)
            shift = theta_b - theta_geo

            do i = 1, n_theta_map
                theta_uniform = 2.0_dp*pi*real(i - 1, dp)/real(n_theta_map, dp)
                call interp_periodic(theta_geo, shift, nu, theta_uniform, &
                                     delta, nu_at)
                data(k, i, 1) = delta
                data(k, i, 2) = nu_at
            end do
            data(k, n_theta_map + 1, :) = data(k, 1, :)
        end do

        map_s_min = s_grid(1)
        map_s_max = s_grid(nflux)
        call destroy_boozer_angle_map()
        call construct_batch_splines_2d([map_s_min, 0.0_dp], &
                                        [map_s_max, 2.0_dp*pi], data, [5, 5], &
                                        [.false., .true.], map_spline)
        map_ready = .true.
        deallocate (data)
    end subroutine build_boozer_angle_map

    subroutine boozer_angle_at(s, theta_geo, theta_b, nu)
        !! theta_b is the Boozer poloidal angle at the orbit point; nu is the
        !! normalized stream function, so phi_B - phi_geom = 2*pi*nu/nper.
        real(dp), intent(in) :: s, theta_geo
        real(dp), intent(out) :: theta_b, nu

        real(dp) :: x(2), y(2), s_use, theta_use

        if (.not. map_ready) error stop 'boozer_angle_map: not built'

        s_use = min(max(s, map_s_min), map_s_max)
        theta_use = modulo(theta_geo, 2.0_dp*pi)
        x = [s_use, theta_use]
        call evaluate_batch_splines_2d(map_spline, x, y)
        theta_b = theta_geo + y(1)
        nu = y(2)
    end subroutine boozer_angle_at

    subroutine surface_samples(theta_b, m, rmnc, rmns, zmnc, zmns, vmnc, vmns, &
                               nmode, R_axis, Z_axis, theta_geo, nu)
        real(dp), intent(in) :: theta_b(:), m(:)
        real(dp), intent(in) :: rmnc(:), rmns(:), zmnc(:), zmns(:)
        real(dp), intent(in) :: vmnc(:), vmns(:)
        integer, intent(in) :: nmode
        real(dp), intent(in) :: R_axis, Z_axis
        real(dp), intent(out) :: theta_geo(:), nu(:)

        real(dp) :: angle, cs, sn, R_pt, Z_pt, nu_pt
        integer :: i, j

        do i = 1, size(theta_b)
            R_pt = 0.0_dp
            Z_pt = 0.0_dp
            nu_pt = 0.0_dp
            do j = 1, nmode
                angle = m(j)*theta_b(i)
                cs = cos(angle)
                sn = sin(angle)
                R_pt = R_pt + rmnc(j)*cs + rmns(j)*sn
                Z_pt = Z_pt + zmnc(j)*cs + zmns(j)*sn
                nu_pt = nu_pt + vmnc(j)*cs + vmns(j)*sn
            end do
            theta_geo(i) = modulo(atan2(Z_pt - Z_axis, R_pt - R_axis), &
                                  2.0_dp*pi)
            nu(i) = nu_pt
        end do
    end subroutine surface_samples

    subroutine require_monotonic(theta_geo, ksurf)
        real(dp), intent(in) :: theta_geo(:)
        integer, intent(in) :: ksurf

        integer :: i, n
        real(dp) :: step, total

        n = size(theta_geo)
        total = 0.0_dp
        do i = 1, n
            step = modulo(theta_geo(modulo(i, n) + 1) - theta_geo(i), &
                          2.0_dp*pi)
            if (step <= 0.0_dp .or. step >= pi) then
                write (*, *) 'boozer_angle_map: surface ', ksurf, &
                    ' is not star-shaped about the chart axis at node ', i
                error stop
            end if
            total = total + step
        end do
        if (abs(total - 2.0_dp*pi) > 1.0e-8_dp) then
            write (*, *) 'boozer_angle_map: surface ', ksurf, &
                ' does not wind once about the chart axis: ', total
            error stop
        end if
    end subroutine require_monotonic

    subroutine interp_periodic(abscissa, first, second, at, first_out, second_out)
        !! Linear interpolation on a strictly increasing periodic abscissa.
        real(dp), intent(in) :: abscissa(:), first(:), second(:), at
        real(dp), intent(out) :: first_out, second_out

        integer :: n, lo, hi, mid
        real(dp) :: span, weight, target_value

        n = size(abscissa)
        target_value = modulo(at, 2.0_dp*pi)

        if (target_value < abscissa(1) .or. target_value >= abscissa(n)) then
            ! Wrap interval between the last and the first node.
            span = modulo(abscissa(1) - abscissa(n), 2.0_dp*pi)
            weight = modulo(target_value - abscissa(n), 2.0_dp*pi)/span
            ! theta_B and theta_geo both advance by 2*pi over one circuit, so
            ! their difference is periodic and needs no branch correction.
            first_out = first(n) + weight*(first(1) - first(n))
            second_out = second(n) + weight*(second(1) - second(n))
            return
        end if

        lo = 1
        hi = n
        do while (hi - lo > 1)
            mid = (lo + hi)/2
            if (abscissa(mid) <= target_value) then
                lo = mid
            else
                hi = mid
            end if
        end do
        weight = (target_value - abscissa(lo))/(abscissa(hi) - abscissa(lo))
        first_out = first(lo) + weight*(first(hi) - first(lo))
        second_out = second(lo) + weight*(second(hi) - second(lo))
    end subroutine interp_periodic

    subroutine read_axisymmetric_boozer(path, s_grid, mpol, rmnc, rmns, zmnc, &
                                        zmns, vmnc, vmns, nflux, nmode, nper)
        !! Private reader.  The shared boozer_read state belongs to the
        !! axisymmetric field, which under inp_swi=11 is the GEQDSK; reusing it
        !! here would overwrite the equilibrium the orbit is integrated in.
        character(len=*), intent(in) :: path
        real(dp), allocatable, intent(out) :: s_grid(:), mpol(:, :)
        real(dp), allocatable, intent(out) :: rmnc(:, :), rmns(:, :)
        real(dp), allocatable, intent(out) :: zmnc(:, :), zmns(:, :)
        real(dp), allocatable, intent(out) :: vmnc(:, :), vmns(:, :)
        integer, intent(out) :: nflux, nmode, nper

        integer, parameter :: ncol1 = 5, ncol2 = 8
        integer :: unit, m0b, n0b, ksurf, kmode
        real(dp) :: flux, a_minor, R0
        real(dp) :: params(ncol1 + 1), row(ncol2 + 2)

        open (newunit=unit, file=path, action='read', status='old')
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
    end subroutine read_axisymmetric_boozer

end module boozer_angle_map
