module neort_krook
    use iso_fortran_env, only: dp => real64
    use ieee_arithmetic, only: ieee_is_finite
    use fortnum_quadrature, only: gauss_legendre
    implicit none
    private
    public :: krook_rate, lorentzian, integrate_krook
    public :: krook_collisions, krook_nu, krook_model, krook_harmonic
    public :: krook_geometry_trapped, krook_geometry_passing
    public :: krook_eta_steps, krook_rtol, krook_atol, krook_max_depth
    public :: validate_krook, reset_krook

    real(dp), parameter :: pi = acos(-1.0_dp)
    logical :: krook_collisions = .false.
    real(dp) :: krook_nu = 0.0_dp
    integer :: krook_model = 0, krook_harmonic = 0
    real(dp) :: krook_geometry_trapped = 1.0_dp, krook_geometry_passing = 1.0_dp
    integer :: krook_eta_steps = 32, krook_max_depth = 20
    real(dp) :: krook_rtol = 1.0e-5_dp, krook_atol = 0.0_dp

    type :: resonance_interval_t
        real(dp) :: a, b, ra, rb, bracket_ra, bracket_rb, nu, angle
        logical :: direct
    end type resonance_interval_t

    abstract interface
        subroutine residual_i(eta, res, dres)
            import :: dp
            real(dp), intent(in) :: eta
            real(dp), intent(out) :: res, dres
        end subroutine residual_i
        subroutine weight_i(eta, values)
            import :: dp
            real(dp), intent(in) :: eta
            real(dp), intent(out) :: values(:)
        end subroutine weight_i
    end interface

contains

    function krook_rate(base, speed, harmonic, model, geometry) result(rate)
        real(dp), intent(in) :: base, speed, harmonic, geometry
        integer, intent(in) :: model
        real(dp) :: rate, kernel, energy_rate, harmonic_rate
        if (base == 0.0_dp) then
            rate = 0.0_dp
            return
        end if
        call krook_leaf(0.0_dp, 1.0_dp, 1.0_dp, pi, base, speed, harmonic, &
                        geometry, kernel, energy_rate, harmonic_rate)
        select case (model)
        case (0)
            rate = base*geometry
        case (1)
            rate = harmonic_rate
        case (2)
            rate = energy_rate
        case default
            error stop 'invalid Krook model'
        end select
    end function krook_rate

    function lorentzian(res, nu) result(value)
        real(dp), intent(in) :: res, nu
        real(dp) :: value, scale, energy_rate, harmonic_rate
        scale = max(abs(res), nu)
        call krook_leaf(res/scale, nu/scale, scale, pi, 0.0_dp, &
                        1.0_dp, 0.0_dp, 1.0_dp, value, energy_rate, harmonic_rate)
    end function lorentzian

    subroutine validate_krook()
        if (.not. ieee_is_finite(krook_nu) .or. krook_nu < 0.0_dp) &
            error stop 'krook_nu must be finite and nonnegative [s^-1]'
        if (krook_model < 0 .or. krook_model > 2) &
            error stop 'krook_model must be 0, 1 or 2'
        if (krook_harmonic < 0 .or. krook_harmonic > 1) &
            error stop 'krook_harmonic must be 0 or 1'
        if (.not. ieee_is_finite(krook_geometry_trapped) .or. &
            .not. ieee_is_finite(krook_geometry_passing)) &
            error stop 'Krook geometry factors must be finite'
        if (min(krook_geometry_trapped, krook_geometry_passing) <= 0.0_dp) &
            error stop 'Krook geometry factors must be positive'
        if (krook_eta_steps < 2 .or. krook_max_depth < 1) &
            error stop 'Krook eta steps >= 2 and max depth >= 1 required'
        if (.not. ieee_is_finite(krook_rtol) .or. krook_rtol <= 0.0_dp) &
            error stop 'krook_rtol must be finite and positive'
        if (.not. ieee_is_finite(krook_atol) .or. krook_atol < 0.0_dp) &
            error stop 'krook_atol must be finite and nonnegative'
    end subroutine validate_krook

    subroutine reset_krook()
        krook_collisions = .false.
        krook_nu = 0.0_dp
        krook_model = 0
        krook_harmonic = 0
        krook_geometry_trapped = 1.0_dp
        krook_geometry_passing = 1.0_dp
        krook_eta_steps = 32
        krook_max_depth = 20
        krook_rtol = 1.0e-5_dp
        krook_atol = 0.0_dp
    end subroutine reset_krook

    subroutine integrate_krook(residual, weight, a, b, nu, nseed, rtol, atol, &
                               value, status, max_depth)
        procedure(residual_i) :: residual
        procedure(weight_i) :: weight
        real(dp), intent(in) :: a, b, nu, rtol, atol
        integer, intent(in) :: nseed
        real(dp), intent(out) :: value(:)
        integer, intent(out) :: status
        integer, intent(in), optional :: max_depth
        real(dp) :: knots(4*nseed + 2), x(16), w(16), total_error(size(value))
        integer :: count, depth_limit
        value = 0.0_dp
        total_error = 0.0_dp
        status = 0
        if (nu <= 0.0_dp .or. b <= a .or. nseed < 2) then
            status = 1
            return
        end if
        depth_limit = 20
        if (present(max_depth)) depth_limit = max_depth
        call gauss_legendre(16, x, w)
        call monotonic_knots(residual, a, b, nseed, knots, count)
        call integrate_monotonic_branches(residual, weight, knots(:count), nu, &
                              rtol, atol, x, w, depth_limit, value, total_error, status)
        if (status == 0) then
            if (any(total_error > atol + rtol*abs(value))) status = 2
        end if
    end subroutine integrate_krook

    subroutine integrate_monotonic_branches(residual, weight, knots, nu, rtol, &
                                            atol, x, w, depth, value, errors, status)
        procedure(residual_i) :: residual
        procedure(weight_i) :: weight
        real(dp), intent(in) :: knots(:), nu, rtol, atol, x(:), w(:)
        integer, intent(in) :: depth
        real(dp), intent(inout) :: value(:), errors(:)
        integer, intent(inout) :: status
        type(resonance_interval_t) :: interval
        real(dp) :: piece(size(value)), piece_error(size(value))
        integer :: k
        do k = 1, size(knots) - 1
            call make_interval(residual, knots(k), knots(k + 1), nu, interval)
            call integrate_branch(residual, weight, interval, rtol, &
                atol/real(size(knots) - 1, dp), x, w, depth, piece, piece_error, status)
            value = value + piece
            errors = errors + piece_error
        end do
    end subroutine integrate_monotonic_branches

    subroutine integrate_branch(residual, weight, interval, rtol, atol, x, w, &
                                depth, value, error_estimate, status)
        procedure(residual_i) :: residual
        procedure(weight_i) :: weight
        type(resonance_interval_t), intent(in) :: interval
        real(dp), intent(in) :: rtol, atol, x(:), w(:)
        integer, intent(in) :: depth
        real(dp), intent(out) :: value(:), error_estimate(:)
        integer, intent(inout) :: status
        type(resonance_interval_t) :: band
        real(dp) :: limits(1100), piece(size(value)), piece_error(size(value))
        integer :: count, k
        if (interval%direct) then
            call adaptive_interval(residual, weight, interval, 0.0_dp, 1.0_dp, &
                                 rtol, atol, x, w, depth, value, error_estimate, status)
            return
        end if
        call residual_limits(interval, limits, count)
        value = 0.0_dp
        error_estimate = 0.0_dp
        do k = 1, count - 1
            call make_band(interval, limits(k:k + 1), band)
            call adaptive_interval(residual, weight, band, 0.0_dp, 1.0_dp, &
                rtol, atol/real(count - 1, dp), x, w, depth, piece, piece_error, status)
            value = value + piece
            error_estimate = error_estimate + piece_error
        end do
    end subroutine integrate_branch

    subroutine make_band(interval, limits, band)
        type(resonance_interval_t), intent(in) :: interval
        real(dp), intent(in) :: limits(2)
        type(resonance_interval_t), intent(out) :: band
        band = interval
        band%ra = limits(1)
        band%rb = limits(2)
        call set_angle(band)
    end subroutine make_band

    subroutine residual_limits(interval, limits, count)
        type(resonance_interval_t), intent(in) :: interval
        real(dp), intent(out) :: limits(:)
        integer, intent(out) :: count
        real(dp), parameter :: ratio = 16.0_dp
        real(dp) :: levels(550), lower, upper, largest
        integer :: k, n
        lower = min(interval%ra, interval%rb)
        upper = max(interval%ra, interval%rb)
        largest = max(abs(lower), abs(upper))
        levels(1) = interval%nu
        n = 1
        do while (levels(n) < largest/ratio)
            n = n + 1
            levels(n) = levels(n - 1)*ratio
        end do
        count = 1
        limits(count) = lower
        do k = n, 1, -1
            call append_limit(-levels(k), lower, upper, limits, count)
        end do
        call append_limit(0.0_dp, lower, upper, limits, count)
        do k = 1, n
            call append_limit(levels(k), lower, upper, limits, count)
        end do
        count = count + 1
        limits(count) = upper
    end subroutine residual_limits

    subroutine append_limit(value, lower, upper, limits, count)
        real(dp), intent(in) :: value, lower, upper
        real(dp), intent(inout) :: limits(:)
        integer, intent(inout) :: count
        if (value <= lower .or. value >= upper) return
        count = count + 1
        limits(count) = value
    end subroutine append_limit

    subroutine monotonic_knots(residual, a, b, nseed, knots, count)
        procedure(residual_i) :: residual
        real(dp), intent(in) :: a, b
        integer, intent(in) :: nseed
        real(dp), intent(out) :: knots(:)
        integer, intent(out) :: count
        real(dp) :: left, right, rl, rr, dl, dr, stationary(2)
        integer :: k, j, found
        count = 1
        knots(count) = a
        left = a
        call residual(left, rl, dl)
        do k = 1, 2*nseed
            right = a + (b - a)*real(k, dp)/real(2*nseed, dp)
            call residual(right, rr, dr)
        call stationary_points(residual, left, right, rl, rr, dl, dr, stationary, found)
            do j = 1, found
                if (stationary(j) <= knots(count) .or. stationary(j) >= b) cycle
                count = count + 1
                knots(count) = stationary(j)
            end do
            left = right
            rl = rr
            dl = dr
        end do
        count = count + 1
        knots(count) = b
    end subroutine monotonic_knots

    subroutine stationary_points(residual, a, b, ra, rb, da, db, points, count)
        procedure(residual_i) :: residual
        real(dp), intent(in) :: a, b, ra, rb, da, db
        real(dp), intent(out) :: points(2)
        integer, intent(out) :: count
        real(dp) :: quadratic, linear, fraction, probe, rp, derivative
        count = 0
        if (different_sign(da, db)) then
            count = 1
            call find_extremum(residual, a, b, da, points(1))
            return
        end if
        call krook_hermite_leaf(da, db, (rb - ra)/(b - a), quadratic, linear)
        if (quadratic /= 0.0_dp) then
            fraction = -linear/(2.0_dp*quadratic)
            if (fraction > 0.0_dp .and. fraction < 1.0_dp) then
                probe = a + (b - a)*fraction
                call residual(probe, rp, derivative)
                call append_extremum(residual, a, probe, da, derivative, points, count)
                call append_extremum(residual, probe, b, derivative, db, points, count)
            end if
        end if
        if (count == 0 .and. db == 0.0_dp .and. da /= 0.0_dp) then
            count = 1
            points(1) = b
        end if
    end subroutine stationary_points

    subroutine append_extremum(residual, a, b, da, db, points, count)
        procedure(residual_i) :: residual
        real(dp), intent(in) :: a, b, da, db
        real(dp), intent(inout) :: points(:)
        integer, intent(inout) :: count
        if (.not. different_sign(da, db)) return
        count = count + 1
        call find_extremum(residual, a, b, da, points(count))
    end subroutine append_extremum

    logical function different_sign(a, b)
        real(dp), intent(in) :: a, b
        different_sign = (a < 0.0_dp .and. b > 0.0_dp) .or. &
                         (a > 0.0_dp .and. b < 0.0_dp)
    end function different_sign

    subroutine find_extremum(residual, left, right, dleft, eta)
        procedure(residual_i) :: residual
        real(dp), intent(in) :: left, right, dleft
        real(dp), intent(out) :: eta
        real(dp) :: a, b, res, derivative
        integer :: k
        a = left
        b = right
        do k = 1, 80
            eta = a + (b - a)/2.0_dp
            if (eta == a .or. eta == b) exit
            call residual(eta, res, derivative)
            if (derivative == 0.0_dp) exit
            if ((derivative > 0.0_dp) .eqv. (dleft > 0.0_dp)) then
                a = eta
            else
                b = eta
            end if
        end do
    end subroutine find_extremum

    subroutine make_interval(residual, a, b, nu, interval)
        procedure(residual_i) :: residual
        real(dp), intent(in) :: a, b, nu
        type(resonance_interval_t), intent(out) :: interval
        real(dp) :: derivative
        interval%a = a
        interval%b = b
        interval%nu = nu
        call residual(a, interval%ra, derivative)
        call residual(b, interval%rb, derivative)
        interval%bracket_ra = interval%ra
        interval%bracket_rb = interval%rb
        interval%direct = abs(interval%rb - interval%ra) <= nu
        call set_angle(interval)
    end subroutine make_interval

    subroutine set_angle(interval)
        type(resonance_interval_t), intent(inout) :: interval
        real(dp) :: scale, ra, rb, width
        scale = max(interval%nu, abs(interval%ra), abs(interval%rb))
        ra = interval%ra/scale
        rb = interval%rb/scale
        width = interval%nu/scale
        interval%angle = atan2(width*(rb - ra), width**2 + ra*rb)
    end subroutine set_angle

    recursive subroutine adaptive_interval(residual, weight, interval, a, b, &
                                 rtol, atol, x, w, depth, value, error_estimate, status)
        procedure(residual_i) :: residual
        procedure(weight_i) :: weight
        type(resonance_interval_t), intent(in) :: interval
        real(dp), intent(in) :: a, b, rtol, atol, x(:), w(:)
        integer, intent(in) :: depth
        real(dp), intent(out) :: value(:)
        real(dp), intent(out) :: error_estimate(:)
        integer, intent(inout) :: status
        real(dp) :: coarse(size(value)), fine(size(value)), other(size(value)), mid
        real(dp) :: left_error(size(value)), right_error(size(value))
        mid = (a + b)/2.0_dp
        call fixed_interval(residual, weight, interval, a, b, x, w, coarse, status)
        call fixed_interval(residual, weight, interval, a, mid, x, w, fine, status)
        call fixed_interval(residual, weight, interval, mid, b, x, w, other, status)
        value = fine + other
        error_estimate = abs(value - coarse)
        if (status /= 0) return
        if (all(abs(value - coarse) <= atol + rtol*abs(value))) return
        if (depth <= 0) return
        call adaptive_interval(residual, weight, interval, a, mid, rtol, atol/2, &
                               x, w, depth - 1, fine, left_error, status)
        call adaptive_interval(residual, weight, interval, mid, b, rtol, atol/2, &
                               x, w, depth - 1, other, right_error, status)
        value = fine + other
        error_estimate = left_error + right_error
    end subroutine adaptive_interval

    subroutine fixed_interval(residual, weight, interval, a, b, x, w, value, status)
        procedure(residual_i) :: residual
        procedure(weight_i) :: weight
        type(resonance_interval_t), intent(in) :: interval
        real(dp), intent(in) :: a, b, x(:), w(:)
        real(dp), intent(out) :: value(:)
        integer, intent(inout) :: status
        real(dp) :: node, eta, jacobian, values(size(value))
        integer :: k
        value = 0.0_dp
        do k = 1, size(x)
            node = (a + b)/2.0_dp + (b - a)*x(k)/2.0_dp
            call mapped_node(residual, interval, node, eta, jacobian, status)
            if (status /= 0) return
            call weight(eta, values)
            value = value + w(k)*jacobian*values
        end do
        value = value*(b - a)/2.0_dp
    end subroutine fixed_interval

    subroutine mapped_node(residual, interval, node, eta, jacobian, status)
        procedure(residual_i) :: residual
        type(resonance_interval_t), intent(in) :: interval
        real(dp), intent(in) :: node
        real(dp), intent(out) :: eta, jacobian
        integer, intent(inout) :: status
        real(dp) :: scale, ra, width, target, res, derivative
        if (interval%direct) then
            eta = interval%a + node*(interval%b - interval%a)
            call residual(eta, res, derivative)
            jacobian = (interval%b - interval%a)*lorentzian(res, interval%nu)
            return
        end if
        scale = max(abs(interval%ra), interval%nu)
        ra = interval%ra/scale
        width = interval%nu/scale
        call krook_angle_leaf(ra, width, interval%nu, pi, node, interval%angle, &
                              abs(interval%angle), target, jacobian)
        target = max(min(interval%ra, interval%rb), &
                     min(max(interval%ra, interval%rb), target))
        call inverse_residual(residual, interval, target, eta, derivative, status)
        if (derivative == 0.0_dp) status = 3
        if (status /= 0) return
        jacobian = jacobian/abs(derivative)
    end subroutine mapped_node

    subroutine inverse_residual(residual, interval, target, eta, derivative, status)
        procedure(residual_i) :: residual
        type(resonance_interval_t), intent(in) :: interval
        real(dp), intent(in) :: target
        real(dp), intent(out) :: eta, derivative
        integer, intent(inout) :: status
        real(dp) :: a, b, res, next
        logical :: increasing
        integer :: k
        a = interval%a
        b = interval%b
        increasing = interval%bracket_rb > interval%bracket_ra
        eta = a + (b - a)*(target - interval%bracket_ra)/ &
              (interval%bracket_rb - interval%bracket_ra)
        do k = 1, 100
            call residual(eta, res, derivative)
            call check_inverse_sample(res, derivative, interval, status)
            if (status /= 0) return
            if (res == target) return
            if ((res < target) .eqv. increasing) then
                a = eta
            else
                b = eta
            end if
            next = a + (b - a)/2.0_dp
            if (next == a .or. next == b) return
            call next_inverse_step(eta, res - target, derivative, a, b, next)
            eta = next
        end do
    end subroutine inverse_residual

    subroutine check_inverse_sample(residual, derivative, interval, status)
        real(dp), intent(in) :: residual, derivative
        type(resonance_interval_t), intent(in) :: interval
        integer, intent(inout) :: status
        if (.not. ieee_is_finite(residual) .or. &
            .not. ieee_is_finite(derivative)) status = 3
        if (different_sign(derivative, interval%bracket_rb-interval%bracket_ra)) &
            status = 3
    end subroutine check_inverse_sample

    subroutine next_inverse_step(eta, residual, derivative, a, b, next)
        real(dp), intent(in) :: eta, residual, derivative, a, b
        real(dp), intent(inout) :: next
        if (derivative == 0.0_dp) return
        next = eta - residual/derivative
        if (next <= a .or. next >= b) next = a + (b - a)/2.0_dp
    end subroutine next_inverse_step

    include 'generated/krook_leaf.inc'
end module neort_krook
