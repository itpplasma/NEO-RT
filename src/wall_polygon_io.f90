module neort_wall_io
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64

    implicit none

    private

    integer, parameter, public :: WALL_IO_OK = 0
    integer, parameter, public :: WALL_IO_MISSING_PATH = 1
    integer, parameter, public :: WALL_IO_MISSING_FILE = 2
    integer, parameter, public :: WALL_IO_OPEN_ERROR = 3
    integer, parameter, public :: WALL_IO_READ_ERROR = 4
    integer, parameter, public :: WALL_IO_MALFORMED_ROW = 5
    integer, parameter, public :: WALL_IO_NONFINITE = 6
    integer, parameter, public :: WALL_IO_TOO_FEW_VERTICES = 7
    integer, parameter, public :: WALL_IO_INVALID_UNITS = 8
    integer, parameter, public :: WALL_IO_NONPOSITIVE_RADIUS = 9
    integer, parameter, public :: WALL_IO_ZERO_AREA = 10
    integer, parameter, public :: WALL_IO_SELF_INTERSECTION = 11
    integer, parameter, public :: WALL_IO_TOO_LONG_RECORD = 12
    integer, parameter, public :: WALL_IO_CLOSE_ERROR = 13

    integer, parameter :: MAX_WALL_RECORD = 4096

    type, public :: wall_polygon_t
        ! Vertices are stored in backend CGS centimetres as (R/Z, vertex).
        ! Input records may be metres or centimetres, but the caller must make
        ! that choice explicitly at the load boundary.
        real(dp), allocatable :: vertices(:, :)
        character(len=16) :: input_units = ''
        character(len=64) :: hash = ''
    contains
        procedure, public :: clear => clear_wall_polygon
    end type wall_polygon_t

    public :: load_wall_polygon
    public :: require_readable_wall_polygon
    public :: wall_io_status_message

contains

    subroutine clear_wall_polygon(wall)
        class(wall_polygon_t), intent(inout) :: wall

        if (allocated(wall%vertices)) deallocate (wall%vertices)
        wall%input_units = ''
        wall%hash = ''
    end subroutine clear_wall_polygon

    subroutine load_wall_polygon(path, wall, status, message, input_units)
        character(len=*), intent(in) :: path
        type(wall_polygon_t), intent(out) :: wall
        integer, intent(out) :: status
        character(len=*), intent(out), optional :: message
        character(len=*), intent(in), optional :: input_units

        character(len=4096) :: line
        character(len=256) :: parse_message
        integer :: capacity, count, ios, close_status, unit
        integer :: n_tokens, parse_status
        logical :: exists
        real(dp) :: r, z, unit_scale
        real(dp), allocatable :: data(:, :), grown(:, :)
        character(len=16) :: units

        call wall%clear()
        wall%input_units = ''
        status = WALL_IO_OK
        call set_message(message, '')

        units = 'm'
        if (present(input_units)) units = normalize_units(input_units)
        call wall_unit_scale(units, unit_scale, status)
        if (status /= WALL_IO_OK) then
            call set_message(message, wall_io_status_message(status))
            return
        end if

        if (len_trim(path) == 0) then
            status = WALL_IO_MISSING_PATH
            call set_message(message, wall_io_status_message(status))
            return
        end if

        inquire (file=trim(path), exist=exists, iostat=ios)
        if (ios /= 0) then
            status = WALL_IO_OPEN_ERROR
            call set_message(message, wall_io_status_message(status))
            return
        end if
        if (.not. exists) then
            status = WALL_IO_MISSING_FILE
            call set_message(message, wall_io_status_message(status))
            return
        end if

        open (newunit=unit, file=trim(path), status='old', action='read', &
            form='formatted', iostat=ios)
        if (ios /= 0) then
            status = WALL_IO_OPEN_ERROR
            call set_message(message, wall_io_status_message(status))
            return
        end if

        capacity = 0
        count = 0
        do
            read (unit, '(A)', iostat=ios) line
            if (ios < 0) exit
            if (ios > 0) then
                status = WALL_IO_READ_ERROR
                exit
            end if
            if (len_trim(line) >= len(line)) then
                status = WALL_IO_TOO_LONG_RECORD
                exit
            end if

            call parse_wall_line(line, r, z, n_tokens, parse_status, parse_message)
            if (parse_status /= WALL_IO_OK) then
                status = parse_status
                exit
            end if
            if (n_tokens == 0) cycle
            if (n_tokens /= 2) then
                status = WALL_IO_MALFORMED_ROW
                exit
            end if

            if (count == capacity) then
                if (capacity == 0) then
                    capacity = 16
                else
                    capacity = 2 * capacity
                end if
                allocate (grown(2, capacity))
                if (count > 0) grown(:, 1:count) = data(:, 1:count)
                call move_alloc(grown, data)
            end if
            count = count + 1
            data(1, count) = r
            data(2, count) = z
        end do
        close (unit, iostat=close_status)
        if (status == WALL_IO_OK .and. close_status /= 0) then
            status = WALL_IO_CLOSE_ERROR
        end if

        if (status /= WALL_IO_OK) then
            call wall%clear()
            call set_message(message, wall_io_status_message(status))
            return
        end if
        call finalize_wall_polygon(data, count, unit_scale, units, wall, status)
        call set_message(message, wall_io_status_message(status))
    end subroutine load_wall_polygon

    subroutine require_readable_wall_polygon(path, input_units)
        character(len=*), intent(in) :: path
        character(len=*), intent(in), optional :: input_units

        type(wall_polygon_t) :: wall
        integer :: status

        call load_wall_polygon(path, wall, status, input_units=input_units)
        if (status /= WALL_IO_OK) then
            error stop 'frequency_model=2 requires a readable wall_file polygon'
        end if
    end subroutine require_readable_wall_polygon

    subroutine finalize_wall_polygon(data, count, unit_scale, units, wall, status)
        real(dp), intent(in) :: data(:, :), unit_scale
        integer, intent(in) :: count
        character(len=*), intent(in) :: units
        type(wall_polygon_t), intent(inout) :: wall
        integer, intent(out) :: status

        real(dp), allocatable :: normalized(:, :)
        real(dp) :: tolerance, area, coordinate_scale
        integer :: i, j, normalized_count, unique_count
        logical :: duplicate

        status = WALL_IO_OK
        call wall%clear()
        wall%input_units = normalize_units(units)
        wall%hash = ''
        if (count < 1) then
            status = WALL_IO_TOO_FEW_VERTICES
            return
        end if

        coordinate_scale = max(1.0_dp, maxval(abs(data(:, 1:count)))*unit_scale)
        tolerance = 1.0e-11_dp*coordinate_scale
        allocate(normalized(2, count))
        normalized_count = 0
        do i = 1, count
            if (normalized_count == 0) then
                normalized_count = 1
                normalized(:, normalized_count) = unit_scale*data(:, i)
            else if (distance_squared(normalized(:, normalized_count), &
                    unit_scale*data(:, i)) > tolerance**2) then
                normalized_count = normalized_count + 1
                normalized(:, normalized_count) = unit_scale*data(:, i)
            end if
        end do
        if (normalized_count >= 2) then
            if (distance_squared(normalized(:, 1), &
                    normalized(:, normalized_count)) <= tolerance**2) then
                normalized_count = normalized_count - 1
            end if
        end if
        if (normalized_count < 3) then
            status = WALL_IO_TOO_FEW_VERTICES
            return
        end if

        unique_count = 0
        do i = 1, normalized_count
            duplicate = .false.
            do j = 1, i - 1
                if (distance_squared(normalized(:, i), normalized(:, j)) &
                        <= tolerance**2) then
                    duplicate = .true.
                    exit
                end if
            end do
            if (.not. duplicate) unique_count = unique_count + 1
        end do
        if (unique_count < 3) then
            status = WALL_IO_TOO_FEW_VERTICES
            return
        end if
        if (any(normalized(1, 1:normalized_count) <= 0.0_dp)) then
            status = WALL_IO_NONPOSITIVE_RADIUS
            return
        end if

        area = polygon_signed_area(normalized(:, 1:normalized_count), normalized_count)
        if (abs(area) <= tolerance**2) then
            status = WALL_IO_ZERO_AREA
            return
        end if
        if (has_self_intersection(normalized(:, 1:normalized_count), &
                normalized_count, tolerance)) then
            status = WALL_IO_SELF_INTERSECTION
            return
        end if

        allocate(wall%vertices(2, normalized_count))
        wall%vertices = normalized(:, 1:normalized_count)
        call compute_wall_hash(wall)
        status = WALL_IO_OK
    end subroutine finalize_wall_polygon

    subroutine compute_wall_hash(wall)
        type(wall_polygon_t), intent(inout) :: wall

        character(len=:), allocatable :: payload
        character(len=64) :: record
        character(len=64) :: digest
        integer :: i

        wall%hash = ''
        if (.not. allocated(wall%vertices)) return
        payload = ''
        do i = 1, size(wall%vertices, 2)
            write (record, '(ES24.16E3,1X,ES24.16E3)') &
                wall%vertices(1, i), wall%vertices(2, i)
            payload = payload//trim(record)//achar(10)
        end do
        call sha256_digest(payload, digest)
        wall%hash = digest
    end subroutine compute_wall_hash

    subroutine sha256_digest(message, digest)
        character(len=*), intent(in) :: message
        character(len=64), intent(out) :: digest

        integer(int64), parameter :: modulus = 4294967296_int64
        integer(int64), parameter :: mask = 4294967295_int64
        integer(int64), parameter :: initial(8) = [ &
            1779033703_int64, 3144134277_int64, 1013904242_int64, &
            2773480762_int64, 1359893119_int64, 2600822924_int64, &
            528734635_int64, 1541459225_int64]
        integer(int64), parameter :: round_constant(64) = [ &
            1116352408_int64, 1899447441_int64, 3049323471_int64, &
            3921009573_int64, 961987163_int64, 1508970993_int64, &
            2453635748_int64, 2870763221_int64, 3624381080_int64, &
            310598401_int64, 607225278_int64, 1426881987_int64, &
            1925078388_int64, 2162078206_int64, 2614888103_int64, &
            3248222580_int64, 3835390401_int64, 4022224774_int64, &
            264347078_int64, 604807628_int64, 770255983_int64, &
            1249150122_int64, 1555081692_int64, 1996064986_int64, &
            2554220882_int64, 2821834349_int64, 2952996808_int64, &
            3210313671_int64, 3336571891_int64, 3584528711_int64, &
            113926993_int64, 338241895_int64, 666307205_int64, &
            773529912_int64, 1294757372_int64, 1396182291_int64, &
            1695183700_int64, 1986661051_int64, 2177026350_int64, &
            2456956037_int64, 2730485921_int64, 2820302411_int64, &
            3259730800_int64, 3345764771_int64, 3516065817_int64, &
            3600352804_int64, 4094571909_int64, 275423344_int64, &
            430227734_int64, 506948616_int64, 659060556_int64, &
            883997877_int64, 958139571_int64, 1322822218_int64, &
            1537002063_int64, 1747873779_int64, 1955562222_int64, &
            2024104815_int64, 2227730452_int64, 2361852424_int64, &
            2428436474_int64, 2756734187_int64, 3204031479_int64, &
            3329325298_int64]
        integer(int64), allocatable :: padded(:), words(:)
        integer(int64) :: state(8), working(8), t1, t2, bit_length
        integer :: message_length, block_count, base, i, t

        message_length = len(message)
        block_count = (message_length + 9 + 63)/64
        allocate (padded(0:block_count*64 - 1), words(0:63))
        padded = 0_int64
        do i = 1, message_length
            padded(i - 1) = int(iachar(message(i:i)), int64)
        end do
        padded(message_length) = 128_int64
        bit_length = int(message_length, int64)*8_int64
        do i = 0, 7
            padded(block_count*64 - 1 - i) = &
                iand(shiftr(bit_length, 8*i), 255_int64)
        end do

        state = initial
        do base = 0, block_count*64 - 64, 64
            do t = 0, 15
                words(t) = iand(padded(base + 4*t), 255_int64)*16777216_int64 &
                    +iand(padded(base + 4*t + 1), 255_int64)*65536_int64 &
                    +iand(padded(base + 4*t + 2), 255_int64)*256_int64 &
                    +iand(padded(base + 4*t + 3), 255_int64)
            end do
            do t = 16, 63
                words(t) = modulo32(words(t - 16) + small_sigma0(words(t - 15)) &
                    +words(t - 7) + small_sigma1(words(t - 2)), modulus)
            end do
            working = state
            do t = 0, 63
                t1 = modulo32(working(8) + big_sigma1(working(5)) &
                    +choose_word(working(5), working(6), working(7)) &
                    +round_constant(t + 1) + words(t), modulus)
                t2 = modulo32(big_sigma0(working(1)) &
                    +majority_word(working(1), working(2), working(3)), modulus)
                working(8) = working(7)
                working(7) = working(6)
                working(6) = working(5)
                working(5) = modulo32(working(4) + t1, modulus)
                working(4) = working(3)
                working(3) = working(2)
                working(2) = working(1)
                working(1) = modulo32(t1 + t2, modulus)
            end do
            do i = 1, 8
                state(i) = modulo32(state(i) + working(i), modulus)
            end do
        end do
        write (digest, '(8(Z8.8))') state
        deallocate (padded, words)
    end subroutine sha256_digest

    integer(int64) function modulo32(value, modulus)
        integer(int64), intent(in) :: value, modulus

        modulo32 = mod(value, modulus)
    end function modulo32

    integer(int64) function rotate_right(value, amount)
        integer(int64), intent(in) :: value
        integer, intent(in) :: amount
        integer(int64), parameter :: mask = 4294967295_int64

        rotate_right = ior(shiftr(iand(value, mask), amount), &
            iand(ishft(iand(value, mask), 32 - amount), mask))
    end function rotate_right

    integer(int64) function choose_word(x, y, z)
        integer(int64), intent(in) :: x, y, z
        integer(int64), parameter :: mask = 4294967295_int64

        choose_word = iand(ieor(iand(x, y), iand(not(x), z)), mask)
    end function choose_word

    integer(int64) function majority_word(x, y, z)
        integer(int64), intent(in) :: x, y, z
        integer(int64), parameter :: mask = 4294967295_int64

        majority_word = iand(ieor(ieor(iand(x, y), iand(x, z)), iand(y, z)), mask)
    end function majority_word

    integer(int64) function big_sigma0(value)
        integer(int64), intent(in) :: value
        integer(int64), parameter :: mask = 4294967295_int64

        big_sigma0 = iand(ieor(ieor(rotate_right(value, 2), &
            rotate_right(value, 13)), rotate_right(value, 22)), mask)
    end function big_sigma0

    integer(int64) function big_sigma1(value)
        integer(int64), intent(in) :: value
        integer(int64), parameter :: mask = 4294967295_int64

        big_sigma1 = iand(ieor(ieor(rotate_right(value, 6), &
            rotate_right(value, 11)), rotate_right(value, 25)), mask)
    end function big_sigma1

    integer(int64) function small_sigma0(value)
        integer(int64), intent(in) :: value
        integer(int64), parameter :: mask = 4294967295_int64

        small_sigma0 = iand(ieor(ieor(rotate_right(value, 7), &
            rotate_right(value, 18)), shiftr(value, 3)), mask)
    end function small_sigma0

    integer(int64) function small_sigma1(value)
        integer(int64), intent(in) :: value
        integer(int64), parameter :: mask = 4294967295_int64

        small_sigma1 = iand(ieor(ieor(rotate_right(value, 17), &
            rotate_right(value, 19)), shiftr(value, 10)), mask)
    end function small_sigma1

    subroutine wall_unit_scale(units, scale, status)
        character(len=*), intent(in) :: units
        real(dp), intent(out) :: scale
        integer, intent(out) :: status

        character(len=16) :: normalized

        normalized = normalize_units(units)
        scale = 0.0_dp
        status = WALL_IO_INVALID_UNITS
        select case (trim(normalized))
        case ('m')
            scale = 100.0_dp
        case ('cm')
            scale = 1.0_dp
        case default
            return
        end select
        status = WALL_IO_OK
    end subroutine wall_unit_scale

    function normalize_units(units) result(normalized)
        character(len=*), intent(in) :: units
        character(len=16) :: normalized
        integer :: i, code

        normalized = adjustl(units)
        do i = 1, len(normalized)
            code = iachar(normalized(i:i))
            if (code >= iachar('A') .and. code <= iachar('Z')) then
                normalized(i:i) = achar(code + iachar('a') - iachar('A'))
            end if
        end do
        normalized = trim(normalized)
        select case (trim(normalized))
        case ('meter', 'meters', 'metre', 'metres')
            normalized = 'm'
        case ('centimeter', 'centimeters', 'centimetre', 'centimetres', 'cgs')
            normalized = 'cm'
        end select
    end function normalize_units

    pure real(dp) function distance_squared(first, second)
        real(dp), intent(in) :: first(2), second(2)

        distance_squared = sum((first - second)**2)
    end function distance_squared

    pure real(dp) function polygon_signed_area(vertices, count)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: count
        integer :: i, j

        polygon_signed_area = 0.0_dp
        j = count
        do i = 1, count
            polygon_signed_area = polygon_signed_area &
                +vertices(1, j)*vertices(2, i) &
                -vertices(1, i)*vertices(2, j)
            j = i
        end do
        polygon_signed_area = 0.5_dp*polygon_signed_area
    end function polygon_signed_area

    logical function has_self_intersection(vertices, count, tolerance)
        real(dp), intent(in) :: vertices(:, :), tolerance
        integer, intent(in) :: count

        integer :: i, j, i_next, j_next

        has_self_intersection = .false.
        do i = 1, count
            i_next = merge(i + 1, 1, i < count)
            do j = i + 1, count
                j_next = merge(j + 1, 1, j < count)
                if (.not. segments_intersect(vertices(:, i), vertices(:, i_next), &
                        vertices(:, j), vertices(:, j_next), tolerance)) cycle
                if (adjacent_edges(i, j, count) .and. &
                        shared_endpoint_only(vertices(:, i), vertices(:, i_next), &
                        vertices(:, j), vertices(:, j_next), tolerance)) cycle
                has_self_intersection = .true.
                return
            end do
        end do
    end function has_self_intersection

    pure logical function adjacent_edges(first, second, count)
        integer, intent(in) :: first, second, count

        adjacent_edges = second == first + 1
        if (first == 1 .and. second == count) adjacent_edges = .true.
    end function adjacent_edges

    pure logical function shared_endpoint_only(a, b, c, d, tolerance)
        real(dp), intent(in) :: a(2), b(2), c(2), d(2), tolerance
        logical :: shared

        shared = .false.
        if (distance_squared(a, c) <= tolerance**2) then
            shared = .not. point_on_segment(b, c, d, tolerance)
        else if (distance_squared(a, d) <= tolerance**2) then
            shared = .not. point_on_segment(b, c, d, tolerance)
        else if (distance_squared(b, c) <= tolerance**2) then
            shared = .not. point_on_segment(a, c, d, tolerance)
        else if (distance_squared(b, d) <= tolerance**2) then
            shared = .not. point_on_segment(a, c, d, tolerance)
        end if
        shared_endpoint_only = shared
    end function shared_endpoint_only

    pure logical function segments_intersect(a, b, c, d, tolerance)
        real(dp), intent(in) :: a(2), b(2), c(2), d(2), tolerance
        real(dp) :: ab_c, ab_d, cd_a, cd_b

        ab_c = cross_product(a, b, c)
        ab_d = cross_product(a, b, d)
        cd_a = cross_product(c, d, a)
        cd_b = cross_product(c, d, b)
        segments_intersect = .false.
        if (opposite_or_zero(ab_c, ab_d, tolerance) .and. &
                opposite_or_zero(cd_a, cd_b, tolerance)) then
            if (ab_c*ab_d < 0.0_dp .or. cd_a*cd_b < 0.0_dp) then
                segments_intersect = .true.
            else if (point_on_segment(a, c, d, tolerance) .or. &
                    point_on_segment(b, c, d, tolerance) .or. &
                    point_on_segment(c, a, b, tolerance) .or. &
                    point_on_segment(d, a, b, tolerance)) then
                segments_intersect = .true.
            end if
        end if
    end function segments_intersect

    pure logical function opposite_or_zero(left, right, tolerance)
        real(dp), intent(in) :: left, right, tolerance

        opposite_or_zero = (left > tolerance .and. right < -tolerance) .or. &
            (left < -tolerance .and. right > tolerance) .or. &
            abs(left) <= tolerance .or. abs(right) <= tolerance
    end function opposite_or_zero

    pure logical function point_on_segment(point, first, second, tolerance)
        real(dp), intent(in) :: point(2), first(2), second(2), tolerance
        real(dp) :: lower(2), upper(2)

        lower = min(first, second) - tolerance
        upper = max(first, second) + tolerance
        point_on_segment = abs(cross_product(first, second, point)) <= tolerance &
            .and. all(point >= lower) .and. all(point <= upper)
    end function point_on_segment

    pure real(dp) function cross_product(first, second, point)
        real(dp), intent(in) :: first(2), second(2), point(2)

        cross_product = (second(1) - first(1))*(point(2) - first(2)) &
            -(second(2) - first(2))*(point(1) - first(1))
    end function cross_product

    subroutine parse_wall_line(line, r, z, n_tokens, status, message)
        character(len=*), intent(in) :: line
        real(dp), intent(out) :: r, z
        integer, intent(out) :: n_tokens, status
        character(len=*), intent(out) :: message

        character(len=256) :: token(3)
        character(len=:), allocatable :: work
        integer :: comment_hash, comment_bang, comment_pos
        integer :: first, last, nchar, ios

        r = 0.0_dp
        z = 0.0_dp
        n_tokens = 0
        status = WALL_IO_OK
        message = ''
        work = trim(adjustl(line))

        comment_hash = index(work, '#')
        comment_bang = index(work, '!')
        comment_pos = first_positive(comment_hash, comment_bang)
        if (comment_pos > 0) then
            if (comment_pos == 1) then
                work = ''
            else
                work = work(:comment_pos - 1)
            end if
        end if

        work = trim(adjustl(work))
        nchar = len_trim(work)
        first = 1
        do while (first <= nchar)
            do while (first <= nchar)
                if (.not. is_separator(work(first:first))) exit
                first = first + 1
            end do
            if (first > nchar) exit

            last = first
            do while (last <= nchar)
                if (is_separator(work(last:last))) exit
                last = last + 1
            end do
            n_tokens = n_tokens + 1
            if (n_tokens <= size(token)) then
                if (last - first > len(token(n_tokens))) then
                    status = WALL_IO_MALFORMED_ROW
                    message = 'wall row contains an oversized token'
                    return
                end if
                token(n_tokens) = ''
                token(n_tokens)(1:last - first) = work(first:last - 1)
            else
                status = WALL_IO_MALFORMED_ROW
                message = 'wall row contains more than two columns'
                return
            end if
            first = last + 1
        end do

        if (n_tokens == 0) return
        if (n_tokens /= 2) then
            status = WALL_IO_MALFORMED_ROW
            message = 'wall row must contain exactly two columns'
            return
        end if

        read (token(1), *, iostat=ios) r
        if (ios /= 0) then
            status = WALL_IO_MALFORMED_ROW
            message = 'wall R coordinate is not a real number'
            return
        end if
        read (token(2), *, iostat=ios) z
        if (ios /= 0) then
            status = WALL_IO_MALFORMED_ROW
            message = 'wall Z coordinate is not a real number'
            return
        end if
        if (.not. ieee_is_finite(r)) then
            status = WALL_IO_NONFINITE
            message = 'wall R coordinate is not finite'
            return
        end if
        if (.not. ieee_is_finite(z)) then
            status = WALL_IO_NONFINITE
            message = 'wall Z coordinate is not finite'
        end if
    end subroutine parse_wall_line

    integer function first_positive(first, second)
        integer, intent(in) :: first, second

        first_positive = 0
        if (first > 0) first_positive = first
        if (second > 0) then
            if (first_positive == 0) then
                first_positive = second
            else if (second < first_positive) then
                first_positive = second
            end if
        end if
    end function first_positive

    logical function is_separator(character_value)
        character(len=1), intent(in) :: character_value

        is_separator = character_value == ' ' .or. &
            character_value == achar(9) .or. character_value == achar(13)
    end function is_separator

    function wall_io_status_message(status) result(message)
        integer, intent(in) :: status
        character(len=128) :: message

        select case (status)
        case (WALL_IO_OK)
            message = 'ok'
        case (WALL_IO_MISSING_PATH)
            message = 'wall path is empty'
        case (WALL_IO_MISSING_FILE)
            message = 'wall file does not exist'
        case (WALL_IO_OPEN_ERROR)
            message = 'wall file could not be opened'
        case (WALL_IO_READ_ERROR)
            message = 'wall file could not be read'
        case (WALL_IO_MALFORMED_ROW)
            message = 'wall row is not exactly two real values'
        case (WALL_IO_NONFINITE)
            message = 'wall coordinate is not finite'
        case (WALL_IO_TOO_FEW_VERTICES)
            message = 'wall polygon has fewer than three vertices'
        case (WALL_IO_INVALID_UNITS)
            message = 'wall input units must be m or cm'
        case (WALL_IO_NONPOSITIVE_RADIUS)
            message = 'wall polygon contains a nonpositive R coordinate'
        case (WALL_IO_ZERO_AREA)
            message = 'wall polygon has zero signed area'
        case (WALL_IO_SELF_INTERSECTION)
            message = 'wall polygon self-intersects or has overlapping edges'
        case (WALL_IO_TOO_LONG_RECORD)
            message = 'wall record exceeds the supported maximum length'
        case (WALL_IO_CLOSE_ERROR)
            message = 'wall file could not be closed'
        case default
            message = 'unknown wall I/O error'
        end select
    end function wall_io_status_message

    subroutine set_message(message, value)
        character(len=*), intent(out), optional :: message
        character(len=*), intent(in) :: value

        if (present(message)) message = value
    end subroutine set_message

end module neort_wall_io
