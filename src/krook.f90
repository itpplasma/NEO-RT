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
    integer :: krook_eta_steps = 32, krook_max_depth = 40
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
        krook_max_depth = 40
        krook_rtol = 1.0e-5_dp
        krook_atol = 0.0_dp
    end subroutine reset_krook

    subroutine integrate_krook(residual, weight, a, b, nu, nseed, rtol, atol, &
                               value, status, max_depth)
        ! Globally adaptive quadrature of weight(eta)*L(R(eta), nu) on [a, b].
        ! The domain is split into monotonic branches of R and, away from the
        ! direct regime, into Lorentzian residual bands mapped by the angle
        ! substitution.  Pieces are bisected in order of their error relative
        ! to the requested tolerance of the whole integral, until the summed
        ! error estimate satisfies atol + rtol*|value| in every component.
        ! status: 0 converged, 1 invalid input, 2 not converged within
        ! max_depth bisections (or the piece budget), 3 invalid residual.
        ! The seed scan can miss a stationary point of R between two samples
        ! (e.g. a spline wiggle next to the trapped-passing boundary). Such a
        ! branch is detected when the inverse map meets a derivative of the
        ! wrong sign; its extrema are then located around that pitch, added
        ! as knots and the integration restarts (at most max_repairs times).
        integer, parameter :: max_repairs = 32
        procedure(residual_i) :: residual
        procedure(weight_i) :: weight
        real(dp), intent(in) :: a, b, nu, rtol, atol
        integer, intent(in) :: nseed
        real(dp), intent(out) :: value(:)
        integer, intent(out) :: status
        integer, intent(in), optional :: max_depth
        real(dp) :: x(16), w(16), turn
        real(dp), allocatable :: knots(:)
        type(resonance_interval_t), allocatable :: regions(:)
        integer :: count, depth_limit, nregion, repair
        value = 0.0_dp
        status = 0
        if (nu <= 0.0_dp .or. b <= a .or. nseed < 2) then
            status = 1
            return
        end if
        depth_limit = 40
        if (present(max_depth)) depth_limit = max_depth
        call gauss_legendre(16, x, w)
        allocate(knots(4*nseed + 2))
        call monotonic_knots(residual, a, b, nseed, knots, count)
        do repair = 0, max_repairs
            status = 0
            call collect_regions(residual, knots(:count), nu, regions, nregion)
            call integrate_regions(residual, weight, regions(:nregion), rtol, atol, &
                                   x, w, depth_limit, value, status, turn)
            if (status /= 4) return
            ! Non-monotonic branch: split it at the extrema around turn.
            value = 0.0_dp
            status = 3
            if (repair == max_repairs) return
            call insert_extrema(residual, turn, knots, count)
        end do
    end subroutine integrate_krook

    subroutine insert_extrema(residual, turn, knots, count)
        ! R'(turn) has the opposite sign of its branch between two knots.
        ! Samples graded geometrically towards turn from both knots bracket
        ! the nearby sign changes of R', whose extrema become new knots. If
        ! none is bracketed, turn itself splits the branch.
        procedure(residual_i) :: residual
        real(dp), intent(in) :: turn
        real(dp), allocatable, intent(inout) :: knots(:)
        integer, intent(inout) :: count
        integer, parameter :: levels = 60
        real(dp) :: samples(2*levels + 3), found(2*(2*levels + 2) + 1)
        real(dp) :: res(2*levels + 3), der(2*levels + 3), points(2), left, right
        real(dp), allocatable :: merged(:)
        integer :: k, j, cell, n, nfound, npoints
        cell = 0
        do k = 1, count - 1
            if (knots(k) < turn .and. turn < knots(k + 1)) cell = k
        end do
        if (cell == 0) return
        left = knots(cell)
        right = knots(cell + 1)
        n = 1
        samples(n) = left
        do j = 1, levels
            n = n + 1
            samples(n) = turn - (turn - left)*0.5_dp**j
        end do
        n = n + 1
        samples(n) = turn
        do j = levels, 1, -1
            n = n + 1
            samples(n) = turn + (right - turn)*0.5_dp**j
        end do
        n = n + 1
        samples(n) = right
        do k = 1, n
            call residual(samples(k), res(k), der(k))
        end do
        nfound = 0
        do k = 1, n - 1
            if (samples(k + 1) <= samples(k)) cycle
            call stationary_points(residual, samples(k), samples(k + 1), res(k), &
                                   res(k + 1), der(k), der(k + 1), points, npoints)
            do j = 1, npoints
                if (points(j) <= left .or. points(j) >= right) cycle
                if (nfound > 0) then
                    if (points(j) <= found(nfound)) cycle
                end if
                nfound = nfound + 1
                found(nfound) = points(j)
            end do
        end do
        if (nfound == 0) then
            nfound = 1
            found(1) = turn
        end if
        allocate(merged(count + nfound))
        merged(:cell) = knots(:cell)
        merged(cell + 1:cell + nfound) = found(:nfound)
        merged(cell + nfound + 1:) = knots(cell + 1:count)
        count = count + nfound
        call move_alloc(merged, knots)
    end subroutine insert_extrema

    subroutine collect_regions(residual, knots, nu, regions, nregion)
        procedure(residual_i) :: residual
        real(dp), intent(in) :: knots(:), nu
        type(resonance_interval_t), allocatable, intent(out) :: regions(:)
        integer, intent(out) :: nregion
        type(resonance_interval_t) :: interval
        type(resonance_interval_t), allocatable :: copy(:)
        real(dp) :: limits(1100)
        integer :: k, j, count
        allocate(regions(2*size(knots)))
        nregion = 0
        do k = 1, size(knots) - 1
            call make_interval(residual, knots(k), knots(k + 1), nu, interval)
            count = 2
            if (.not. interval%direct) call residual_limits(interval, limits, count)
            if (nregion + count > size(regions)) then
                allocate(copy(2*(nregion + count)))
                copy(:nregion) = regions(:nregion)
                call move_alloc(copy, regions)
            end if
            if (interval%direct) then
                nregion = nregion + 1
                regions(nregion) = interval
                cycle
            end if
            do j = 1, count - 1
                nregion = nregion + 1
                call make_band(interval, limits(j:j + 1), regions(nregion))
            end do
        end do
    end subroutine collect_regions

    subroutine integrate_regions(residual, weight, regions, rtol, atol, x, w, &
                                 depth_limit, value, status, turn)
        procedure(residual_i) :: residual
        procedure(weight_i) :: weight
        type(resonance_interval_t), intent(in) :: regions(:)
        real(dp), intent(in) :: rtol, atol, x(:), w(:)
        integer, intent(in) :: depth_limit
        real(dp), intent(out) :: value(:)
        integer, intent(inout) :: status
        real(dp), intent(out) :: turn
        ! A budget of pieces bounds the work when the integrand is not
        ! resolvable (e.g. noisy) to the requested tolerance.
        integer, parameter :: max_pieces = 2**14
        integer, allocatable :: owner(:), depth(:), heap(:)
        real(dp), allocatable :: lower(:), upper(:), left(:, :), right(:, :), &
                                 error(:, :), priority(:)
        real(dp) :: scale(size(value)), total_error(size(value)), &
                    frozen_error(size(value)), coarse(size(value)), mid
        integer :: n, k, piece, child, nheap, capacity
        n = size(value)
        capacity = max(64, 2*size(regions))
        allocate(owner(capacity), depth(capacity), lower(capacity), &
                 upper(capacity), left(n, capacity), right(n, capacity), &
                 error(n, capacity), priority(capacity), heap(capacity))
        do k = 1, size(regions)
            owner(k) = k
            depth(k) = 0
            lower(k) = 0.0_dp
            upper(k) = 1.0_dp
            call fixed_interval(residual, weight, regions(k), 0.0_dp, 1.0_dp, &
                                x, w, coarse, status, turn)
            if (status /= 0) return
            call refine_piece(k, coarse)
            if (status /= 0) return
        end do
        piece = size(regions)
        value = sum(left(:, :piece) + right(:, :piece), dim=2)
        total_error = sum(error(:, :piece), dim=2)
        scale = atol + rtol*abs(value)
        nheap = 0
        do k = 1, piece
            call heap_push(k)
        end do
        frozen_error = 0.0_dp
        do
            if (all(total_error <= atol + rtol*abs(value))) then
                ! Confirm the running sums with an exact summation.
                value = sum(left(:, :piece) + right(:, :piece), dim=2)
                total_error = sum(error(:, :piece), dim=2)
                if (all(total_error <= atol + rtol*abs(value))) return
            end if
            if (any(frozen_error > atol + rtol*abs(value)) .or. nheap == 0 &
                .or. piece + 2 > max_pieces) then
                status = 2
                return
            end if
            k = heap_pop()
            if (depth(k) >= depth_limit) then
                frozen_error = frozen_error + error(:, k)
                cycle
            end if
            if (piece + 2 > capacity) call grow()
            mid = (lower(k) + upper(k))/2.0_dp
            do child = 1, 2
                piece = piece + 1
                owner(piece) = owner(k)
                depth(piece) = depth(k) + 1
                if (child == 1) then
                    lower(piece) = lower(k)
                    upper(piece) = mid
                    coarse = left(:, k)
                else
                    lower(piece) = mid
                    upper(piece) = upper(k)
                    coarse = right(:, k)
                end if
                call refine_piece(piece, coarse)
                if (status /= 0) return
                value = value + left(:, piece) + right(:, piece)
                total_error = total_error + error(:, piece)
                call heap_push(piece)
            end do
            ! The parent is replaced by its children.
            value = value - left(:, k) - right(:, k)
            total_error = total_error - error(:, k)
            left(:, k) = 0.0_dp
            right(:, k) = 0.0_dp
            error(:, k) = 0.0_dp
        end do
    contains
        subroutine refine_piece(index, coarse_value)
            integer, intent(in) :: index
            real(dp), intent(in) :: coarse_value(:)
            real(dp) :: half
            half = (lower(index) + upper(index))/2.0_dp
            call fixed_interval(residual, weight, regions(owner(index)), &
                lower(index), half, x, w, left(:, index), status, turn)
            if (status /= 0) return
            call fixed_interval(residual, weight, regions(owner(index)), &
                half, upper(index), x, w, right(:, index), status, turn)
            if (status /= 0) return
            error(:, index) = abs(left(:, index) + right(:, index) - coarse_value)
        end subroutine refine_piece

        real(dp) function piece_priority(index) result(p)
            integer, intent(in) :: index
            integer :: c
            p = 0.0_dp
            do c = 1, n
                if (error(c, index) == 0.0_dp) cycle
                if (scale(c) > 0.0_dp) then
                    p = max(p, error(c, index)/scale(c))
                else
                    p = huge(1.0_dp)
                end if
            end do
        end function piece_priority

        subroutine heap_push(index)
            integer, intent(in) :: index
            integer :: i, parent
            priority(index) = piece_priority(index)
            nheap = nheap + 1
            i = nheap
            heap(i) = index
            do while (i > 1)
                parent = i/2
                if (priority(heap(parent)) >= priority(heap(i))) exit
                call swap(i, parent)
                i = parent
            end do
        end subroutine heap_push

        integer function heap_pop() result(index)
            integer :: i, largest, c
            index = heap(1)
            heap(1) = heap(nheap)
            nheap = nheap - 1
            i = 1
            do
                largest = i
                do c = 2*i, min(2*i + 1, nheap)
                    if (priority(heap(c)) > priority(heap(largest))) largest = c
                end do
                if (largest == i) exit
                call swap(i, largest)
                i = largest
            end do
        end function heap_pop

        subroutine swap(i, j)
            integer, intent(in) :: i, j
            integer :: t
            t = heap(i)
            heap(i) = heap(j)
            heap(j) = t
        end subroutine swap

        subroutine grow()
            integer :: new_capacity
            new_capacity = 2*capacity
            call grow_int(owner, new_capacity)
            call grow_int(depth, new_capacity)
            call grow_int(heap, new_capacity)
            call grow_real(lower, new_capacity)
            call grow_real(upper, new_capacity)
            call grow_real(priority, new_capacity)
            call grow_matrix(left, new_capacity)
            call grow_matrix(right, new_capacity)
            call grow_matrix(error, new_capacity)
            capacity = new_capacity
        end subroutine grow
    end subroutine integrate_regions

    subroutine grow_int(array, capacity)
        integer, allocatable, intent(inout) :: array(:)
        integer, intent(in) :: capacity
        integer, allocatable :: copy(:)
        allocate(copy(capacity))
        copy(:size(array)) = array
        call move_alloc(copy, array)
    end subroutine grow_int

    subroutine grow_real(array, capacity)
        real(dp), allocatable, intent(inout) :: array(:)
        integer, intent(in) :: capacity
        real(dp), allocatable :: copy(:)
        allocate(copy(capacity))
        copy(:size(array)) = array
        call move_alloc(copy, array)
    end subroutine grow_real

    subroutine grow_matrix(array, capacity)
        real(dp), allocatable, intent(inout) :: array(:, :)
        integer, intent(in) :: capacity
        real(dp), allocatable :: copy(:, :)
        allocate(copy(size(array, 1), capacity))
        copy(:, :size(array, 2)) = array
        call move_alloc(copy, array)
    end subroutine grow_matrix

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

    subroutine fixed_interval(residual, weight, interval, a, b, x, w, value, status, &
                              turn)
        procedure(residual_i) :: residual
        procedure(weight_i) :: weight
        type(resonance_interval_t), intent(in) :: interval
        real(dp), intent(in) :: a, b, x(:), w(:)
        real(dp), intent(out) :: value(:)
        integer, intent(inout) :: status
        real(dp), intent(inout) :: turn
        real(dp) :: node, eta, jacobian, values(size(value))
        integer :: k
        value = 0.0_dp
        do k = 1, size(x)
            node = (a + b)/2.0_dp + (b - a)*x(k)/2.0_dp
            call mapped_node(residual, interval, node, eta, jacobian, status)
            if (status == 4) turn = eta
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
        if (derivative == 0.0_dp .and. status == 0) status = 3
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
            if (status == 4) return
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
        ! status 4: R' opposes the branch, which is then not monotonic.
        if (different_sign(derivative, interval%bracket_rb-interval%bracket_ra)) &
            status = 4
        if (.not. ieee_is_finite(residual) .or. &
            .not. ieee_is_finite(derivative)) status = 3
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
