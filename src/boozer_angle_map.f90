module neort_boozer_angle_map
    !! Boozer angles at a point of the direct GEQDSK chart (inp_swi = 11).
    !!
    !! A Boozer-Fourier perturbation is a function of (s_B, theta_B, phi_B);
    !! the direct orbit reports (s, theta_geo, phi).  An axisymmetric Boozer
    !! file of the same equilibrium parametrizes each surface by theta_B, so
    !! theta_B -> (R, Z, nu) is a Fourier sum and theta_geo(theta_B) follows
    !! about the GEQDSK axis without root finding.  The tabulated quantities
    !! are the periodic difference theta_B - theta_geo and nu, against
    !! (s, theta_geo).  Surfaces are matched by toroidal flux,
    !! s_B = s*flux_scale, since the two files need not share the edge flux.
    !!
    !! Shared read-only state, built once on the main thread.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use interpolate, only: BatchSplineData2D, construct_batch_splines_2d, &
        destroy_batch_splines_2d, evaluate_batch_splines_2d

    implicit none
    private

    public :: build_boozer_angle_map, boozer_angles, boozer_angle_map_ready

    real(dp), parameter :: pi = acos(-1.0_dp)
    integer, parameter :: n_theta_map = 512
    ! Relative tolerance for a uniform s grid, which the batch spline assumes.
    real(dp), parameter :: uniform_tol = 1.0e-6_dp

    logical, save :: map_ready = .false.
    integer, save :: map_nper = 1
    real(dp), save :: map_s_min = 0.0_dp, map_s_max = 1.0_dp
    real(dp), save :: flux_scale = 1.0_dp, map_ds = 0.0_dp
    type(BatchSplineData2D), save :: map_spline

contains

    logical function boozer_angle_map_ready()
        boozer_angle_map_ready = map_ready
    end function boozer_angle_map_ready

    subroutine build_boozer_angle_map(path, R_axis, Z_axis, psi_tor)
        !! path: axisymmetric Boozer file in the 8-column (inp_swi = 9) format.
        !! R_axis, Z_axis [m]: axis of the chart the orbit uses.
        !! psi_tor [G cm^2]: toroidal flux per radian of the GEQDSK edge.
        character(len=*), intent(in) :: path
        real(dp), intent(in) :: R_axis, Z_axis, psi_tor
        real(dp), allocatable :: s_grid(:), mpol(:), rz(:, :, :), data(:, :, :)
        real(dp) :: flux
        integer :: nflux, k

        call read_axisymmetric_boozer(path, s_grid, mpol, rz, flux, map_nper)
        nflux = size(s_grid)
        call require_uniform(s_grid)
        flux_scale = abs(psi_tor)/abs(1.0e8_dp*flux/(2.0_dp*pi))
        allocate (data(nflux, n_theta_map + 1, 2))
        do k = 1, nflux
            call tabulate_surface(mpol, rz(:, :, k), R_axis, Z_axis, k, &
                                  data(k, :, :))
        end do
        map_s_min = s_grid(1)
        map_s_max = s_grid(nflux)
        map_ds = (map_s_max - map_s_min)/real(nflux - 1, dp)
        if (map_ready) call destroy_batch_splines_2d(map_spline)
        call construct_batch_splines_2d([map_s_min, 0.0_dp], [map_s_max, 2.0_dp*pi], &
                                        data, [5, 5], [.false., .true.], map_spline)
        map_ready = .true.
    end subroutine build_boozer_angle_map

    subroutine boozer_angles(s, theta_geo, s_b, theta_b, dphi)
        !! Boozer flux label, poloidal angle and phi_B - phi at (s, theta_geo).
        real(dp), intent(in) :: s, theta_geo
        real(dp), intent(out) :: s_b, theta_b, dphi
        real(dp) :: x(2), y(2)

        if (.not. map_ready) error stop "boozer_angle_map: not built"
        s_b = s*flux_scale
        ! Fail closed beyond the outermost tabulated surface rather than
        ! extrapolating; inside the first surface the series is clamped, as in
        ! the Boozer path.
        if (s_b > map_s_max + map_ds) then
            error stop "boozer_angle_map: orbit outside the pert_angle_map surfaces"
        end if
        x = [min(max(s_b, map_s_min), map_s_max), modulo(theta_geo, 2.0_dp*pi)]
        call evaluate_batch_splines_2d(map_spline, x, y)
        theta_b = theta_geo + y(1)
        ! libneo's efit_to_boozer and vmec_to_boozer write the stream-function
        ! column as (phi - phi_B)*nper/(2*pi); test_boozer_angle_map checks this
        ! sign against the direct field line.  Some legacy file headers label
        ! the column with the opposite sign.
        dphi = -2.0_dp*pi*y(2)/real(map_nper, dp)
    end subroutine boozer_angles

    subroutine tabulate_surface(mpol, coeff, R_axis, Z_axis, ksurf, table)
        !! table(:, 1) = theta_B - theta_geo, table(:, 2) = nu on a uniform
        !! theta_geo grid including the periodic end point.
        real(dp), intent(in) :: mpol(:), coeff(:, :), R_axis, Z_axis
        integer, intent(in) :: ksurf
        real(dp), intent(out) :: table(:, :)
        real(dp) :: theta_b(n_theta_map), theta_geo(n_theta_map), nu(n_theta_map)
        real(dp) :: R, Z, angle(size(mpol))
        integer :: i

        do i = 1, n_theta_map
            theta_b(i) = 2.0_dp*pi*real(i - 1, dp)/real(n_theta_map, dp)
            angle = mpol*theta_b(i)
            R = sum(coeff(:, 1)*cos(angle) + coeff(:, 2)*sin(angle))
            Z = sum(coeff(:, 3)*cos(angle) + coeff(:, 4)*sin(angle))
            nu(i) = sum(coeff(:, 5)*cos(angle) + coeff(:, 6)*sin(angle))
            theta_geo(i) = modulo(atan2(Z - Z_axis, R - R_axis), 2.0_dp*pi)
        end do
        call require_monotonic(theta_geo, ksurf)
        do i = 1, n_theta_map
            call interp_periodic(theta_geo, theta_b - theta_geo, nu, &
                                 2.0_dp*pi*real(i - 1, dp)/real(n_theta_map, dp), &
                                 table(i, 1), table(i, 2))
        end do
        table(n_theta_map + 1, :) = table(1, :)
    end subroutine tabulate_surface

    subroutine require_uniform(s_grid)
        real(dp), intent(in) :: s_grid(:)
        real(dp) :: ds

        ds = (s_grid(size(s_grid)) - s_grid(1))/real(size(s_grid) - 1, dp)
        if (any(abs(s_grid(2:) - s_grid(:size(s_grid) - 1) - ds) > uniform_tol*ds)) then
            error stop "boozer_angle_map: s grid of pert_angle_map is not uniform"
        end if
    end subroutine require_uniform

    subroutine require_monotonic(theta_geo, ksurf)
        !! theta_geo must wind once, increasing, as theta_B goes round: the
        !! surface is star-shaped about the axis and both angles have the same
        !! orientation.  Otherwise the interpolation abscissa would fold.
        real(dp), intent(in) :: theta_geo(:)
        integer, intent(in) :: ksurf
        real(dp) :: step, total
        integer :: i, n

        n = size(theta_geo)
        total = 0.0_dp
        do i = 1, n
            step = modulo(theta_geo(modulo(i, n) + 1) - theta_geo(i), 2.0_dp*pi)
            if (step <= 0.0_dp .or. step >= pi) then
                print *, "boozer_angle_map: surface", ksurf, "node", i
                error stop "theta_geo is not monotonic in theta_B"
            end if
            total = total + step
        end do
        if (abs(total - 2.0_dp*pi) > 1.0e-8_dp) then
            error stop "boozer_angle_map: surface does not wind once about the axis"
        end if
    end subroutine require_monotonic

    subroutine interp_periodic(abscissa, first, second, at, first_out, second_out)
        !! Linear interpolation on an increasing periodic abscissa.  theta_B and
        !! theta_geo both advance by 2*pi per turn, so first is periodic.
        real(dp), intent(in) :: abscissa(:), first(:), second(:), at
        real(dp), intent(out) :: first_out, second_out
        real(dp) :: weight, target_value
        integer :: n, lo, hi, mid

        n = size(abscissa)
        target_value = modulo(at, 2.0_dp*pi)
        if (target_value < abscissa(1) .or. target_value >= abscissa(n)) then
            weight = modulo(target_value - abscissa(n), 2.0_dp*pi) &
                     /modulo(abscissa(1) - abscissa(n), 2.0_dp*pi)
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

    subroutine read_axisymmetric_boozer(path, s_grid, mpol, coeff, flux, nper)
        !! Private reader: the shared boozer_read state of do_magfie_mod belongs
        !! to the axisymmetric field, which here is the GEQDSK.
        !! coeff(:, :, k) = rmnc rmns zmnc zmns vmnc vmns on surface k [m, 1].
        character(len=*), intent(in) :: path
        real(dp), allocatable, intent(out) :: s_grid(:), mpol(:), coeff(:, :, :)
        real(dp), intent(out) :: flux
        integer, intent(out) :: nper
        real(dp) :: params(6), row(10), a_minor, R0
        integer :: unit, m0b, n0b, nflux, nmode, ksurf, kmode

        open (newunit=unit, file=path, action='read', status='old')
        read (unit, '(////)')
        read (unit, *) m0b, n0b, nflux, nper, flux, a_minor, R0
        if (n0b /= 0) error stop "boozer_angle_map: pert_angle_map must be axisymmetric"
        nmode = m0b + 1
        allocate (s_grid(nflux), mpol(nmode), coeff(nmode, 6, nflux))
        do ksurf = 1, nflux
            read (unit, '(/)')
            read (unit, *) params
            read (unit, *)
            do kmode = 1, nmode
                read (unit, *) row
                mpol(kmode) = row(1)
                coeff(kmode, :, ksurf) = row(3:8)
            end do
            s_grid(ksurf) = params(1)
        end do
        close (unit)
    end subroutine read_axisymmetric_boozer

end module neort_boozer_angle_map
