module neort_gc_outward_interval
    !! Minimal outward-rounded interval arithmetic for generated certificates.
    !!
    !! This module contains no equilibrium or guiding-centre algebra.  Fortsym
    !! emits the physical expressions against this scalar type, while these
    !! operators provide an inclusion extension of ordinary real arithmetic.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_next_after
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    type, public :: gc_outward_interval_t
        real(dp) :: lo = 0.0_dp
        real(dp) :: hi = 0.0_dp
    end type gc_outward_interval_t

    interface operator(+)
        module procedure add_interval_interval
        module procedure add_interval_real
        module procedure add_real_interval
        module procedure add_interval_integer
        module procedure add_integer_interval
    end interface
    interface operator(-)
        module procedure subtract_interval_interval
        module procedure subtract_interval_real
        module procedure subtract_real_interval
        module procedure subtract_interval_integer
        module procedure subtract_integer_interval
        module procedure negate_interval
    end interface
    interface operator(*)
        module procedure multiply_interval_interval
        module procedure multiply_interval_real
        module procedure multiply_real_interval
        module procedure multiply_interval_integer
        module procedure multiply_integer_interval
    end interface
    interface operator(/)
        module procedure divide_interval_interval
        module procedure divide_interval_real
        module procedure divide_real_interval
        module procedure divide_interval_integer
        module procedure divide_integer_interval
    end interface
    interface operator(**)
        module procedure power_interval_integer
    end interface
    interface abs
        module procedure absolute_interval
    end interface
    interface sqrt
        module procedure square_root_interval
    end interface

    public :: gc_outward_interval
    public :: gc_outward_interval_is_valid
    public :: gc_outward_interval_contains_zero
    public :: gc_outward_interval_intersects
    public :: operator(+), operator(-), operator(*), operator(/), operator(**)
    public :: abs, sqrt

contains

    pure elemental function gc_outward_interval(lo, hi) result(value)
        real(dp), intent(in) :: lo, hi
        type(gc_outward_interval_t) :: value

        value%lo = lo
        value%hi = hi
    end function gc_outward_interval

    pure elemental logical function gc_outward_interval_is_valid(value)
        type(gc_outward_interval_t), intent(in) :: value

        gc_outward_interval_is_valid = ieee_is_finite(value%lo) .and. &
            ieee_is_finite(value%hi) .and. value%lo <= value%hi
    end function gc_outward_interval_is_valid

    pure elemental logical function gc_outward_interval_contains_zero(value)
        type(gc_outward_interval_t), intent(in) :: value

        gc_outward_interval_contains_zero = &
            gc_outward_interval_is_valid(value) .and. &
            value%lo <= 0.0_dp .and. value%hi >= 0.0_dp
    end function gc_outward_interval_contains_zero

    pure elemental logical function gc_outward_interval_intersects(value, lo, hi)
        type(gc_outward_interval_t), intent(in) :: value
        real(dp), intent(in) :: lo, hi

        gc_outward_interval_intersects = &
            gc_outward_interval_is_valid(value) .and. &
            ieee_is_finite(lo) .and. ieee_is_finite(hi) .and. lo <= hi .and. &
            value%hi >= lo .and. value%lo <= hi
    end function gc_outward_interval_intersects

    pure elemental function add_interval_interval(left, right) result(value)
        type(gc_outward_interval_t), intent(in) :: left, right
        type(gc_outward_interval_t) :: value

        value%lo = round_down(left%lo + right%lo)
        value%hi = round_up(left%hi + right%hi)
    end function add_interval_interval

    pure elemental function add_interval_real(left, right) result(value)
        type(gc_outward_interval_t), intent(in) :: left
        real(dp), intent(in) :: right
        type(gc_outward_interval_t) :: value

        value = left + point_interval(right)
    end function add_interval_real

    pure elemental function add_real_interval(left, right) result(value)
        real(dp), intent(in) :: left
        type(gc_outward_interval_t), intent(in) :: right
        type(gc_outward_interval_t) :: value

        value = point_interval(left) + right
    end function add_real_interval

    pure elemental function add_interval_integer(left, right) result(value)
        type(gc_outward_interval_t), intent(in) :: left
        integer, intent(in) :: right
        type(gc_outward_interval_t) :: value

        value = left + real(right, dp)
    end function add_interval_integer

    pure elemental function add_integer_interval(left, right) result(value)
        integer, intent(in) :: left
        type(gc_outward_interval_t), intent(in) :: right
        type(gc_outward_interval_t) :: value

        value = real(left, dp) + right
    end function add_integer_interval

    pure elemental function subtract_interval_interval(left, right) result(value)
        type(gc_outward_interval_t), intent(in) :: left, right
        type(gc_outward_interval_t) :: value

        value%lo = round_down(left%lo - right%hi)
        value%hi = round_up(left%hi - right%lo)
    end function subtract_interval_interval

    pure elemental function subtract_interval_real(left, right) result(value)
        type(gc_outward_interval_t), intent(in) :: left
        real(dp), intent(in) :: right
        type(gc_outward_interval_t) :: value

        value = left - point_interval(right)
    end function subtract_interval_real

    pure elemental function subtract_real_interval(left, right) result(value)
        real(dp), intent(in) :: left
        type(gc_outward_interval_t), intent(in) :: right
        type(gc_outward_interval_t) :: value

        value = point_interval(left) - right
    end function subtract_real_interval

    pure elemental function subtract_interval_integer(left, right) result(value)
        type(gc_outward_interval_t), intent(in) :: left
        integer, intent(in) :: right
        type(gc_outward_interval_t) :: value

        value = left - real(right, dp)
    end function subtract_interval_integer

    pure elemental function subtract_integer_interval(left, right) result(value)
        integer, intent(in) :: left
        type(gc_outward_interval_t), intent(in) :: right
        type(gc_outward_interval_t) :: value

        value = real(left, dp) - right
    end function subtract_integer_interval

    pure elemental function negate_interval(argument) result(value)
        type(gc_outward_interval_t), intent(in) :: argument
        type(gc_outward_interval_t) :: value

        value%lo = -argument%hi
        value%hi = -argument%lo
    end function negate_interval

    pure elemental function multiply_interval_interval(left, right) result(value)
        type(gc_outward_interval_t), intent(in) :: left, right
        type(gc_outward_interval_t) :: value
        real(dp) :: products_lo(4), products_hi(4)

        products_lo = [round_down(left%lo*right%lo), &
            round_down(left%lo*right%hi), &
            round_down(left%hi*right%lo), &
            round_down(left%hi*right%hi)]
        products_hi = [round_up(left%lo*right%lo), &
            round_up(left%lo*right%hi), &
            round_up(left%hi*right%lo), &
            round_up(left%hi*right%hi)]
        value%lo = minval(products_lo)
        value%hi = maxval(products_hi)
    end function multiply_interval_interval

    pure elemental function multiply_interval_real(left, right) result(value)
        type(gc_outward_interval_t), intent(in) :: left
        real(dp), intent(in) :: right
        type(gc_outward_interval_t) :: value

        value = left*point_interval(right)
    end function multiply_interval_real

    pure elemental function multiply_real_interval(left, right) result(value)
        real(dp), intent(in) :: left
        type(gc_outward_interval_t), intent(in) :: right
        type(gc_outward_interval_t) :: value

        value = point_interval(left)*right
    end function multiply_real_interval

    pure elemental function multiply_interval_integer(left, right) result(value)
        type(gc_outward_interval_t), intent(in) :: left
        integer, intent(in) :: right
        type(gc_outward_interval_t) :: value

        value = left*real(right, dp)
    end function multiply_interval_integer

    pure elemental function multiply_integer_interval(left, right) result(value)
        integer, intent(in) :: left
        type(gc_outward_interval_t), intent(in) :: right
        type(gc_outward_interval_t) :: value

        value = real(left, dp)*right
    end function multiply_integer_interval

    pure elemental function divide_interval_interval(left, right) result(value)
        type(gc_outward_interval_t), intent(in) :: left, right
        type(gc_outward_interval_t) :: value
        type(gc_outward_interval_t) :: reciprocal

        if (right%lo <= 0.0_dp .and. right%hi >= 0.0_dp) then
            value = invalid_interval()
            return
        end if
        reciprocal%lo = round_down(1.0_dp/right%hi)
        reciprocal%hi = round_up(1.0_dp/right%lo)
        if (reciprocal%lo > reciprocal%hi) then
            reciprocal = gc_outward_interval(reciprocal%hi, reciprocal%lo)
        end if
        value = left*reciprocal
    end function divide_interval_interval

    pure elemental function divide_interval_real(left, right) result(value)
        type(gc_outward_interval_t), intent(in) :: left
        real(dp), intent(in) :: right
        type(gc_outward_interval_t) :: value

        value = left/point_interval(right)
    end function divide_interval_real

    pure elemental function divide_real_interval(left, right) result(value)
        real(dp), intent(in) :: left
        type(gc_outward_interval_t), intent(in) :: right
        type(gc_outward_interval_t) :: value

        value = point_interval(left)/right
    end function divide_real_interval

    pure elemental function divide_interval_integer(left, right) result(value)
        type(gc_outward_interval_t), intent(in) :: left
        integer, intent(in) :: right
        type(gc_outward_interval_t) :: value

        value = left/real(right, dp)
    end function divide_interval_integer

    pure elemental function divide_integer_interval(left, right) result(value)
        integer, intent(in) :: left
        type(gc_outward_interval_t), intent(in) :: right
        type(gc_outward_interval_t) :: value

        value = real(left, dp)/right
    end function divide_integer_interval

    pure elemental function power_interval_integer(argument, exponent) result(value)
        type(gc_outward_interval_t), intent(in) :: argument
        integer, intent(in) :: exponent
        type(gc_outward_interval_t) :: value
        integer :: factor

        if (exponent < 0) then
            value = invalid_interval()
        else
            value = point_interval(1.0_dp)
            do factor = 1, exponent
                value = value*argument
            end do
        end if
    end function power_interval_integer

    pure elemental function absolute_interval(argument) result(value)
        type(gc_outward_interval_t), intent(in) :: argument
        type(gc_outward_interval_t) :: value

        if (.not. gc_outward_interval_is_valid(argument)) then
            value = invalid_interval()
        else if (argument%lo >= 0.0_dp) then
            value = argument
        else if (argument%hi <= 0.0_dp) then
            value%lo = round_down(-argument%hi)
            value%hi = round_up(-argument%lo)
        else
            value%lo = 0.0_dp
            value%hi = round_up(max(-argument%lo, argument%hi))
        end if
    end function absolute_interval

    pure elemental function square_root_interval(argument) result(value)
        type(gc_outward_interval_t), intent(in) :: argument
        type(gc_outward_interval_t) :: value

        if (.not. gc_outward_interval_is_valid(argument) .or. &
                argument%hi < 0.0_dp) then
            value = invalid_interval()
            return
        end if
        value%lo = 0.0_dp
        if (argument%lo > 0.0_dp) then
            value%lo = round_down(sqrt(argument%lo))
        end if
        value%hi = round_up(sqrt(max(0.0_dp, argument%hi)))
    end function square_root_interval

    pure elemental function point_interval(value_in) result(value)
        real(dp), intent(in) :: value_in
        type(gc_outward_interval_t) :: value

        value%lo = value_in
        value%hi = value_in
    end function point_interval

    pure elemental function invalid_interval() result(value)
        type(gc_outward_interval_t) :: value

        value%lo = 1.0_dp
        value%hi = -1.0_dp
    end function invalid_interval

    pure elemental real(dp) function round_down(value)
        real(dp), intent(in) :: value

        if (ieee_is_finite(value)) then
            round_down = ieee_next_after(value, -huge(value))
        else
            round_down = value
        end if
    end function round_down

    pure elemental real(dp) function round_up(value)
        real(dp), intent(in) :: value

        if (ieee_is_finite(value)) then
            round_up = ieee_next_after(value, huge(value))
        else
            round_up = value
        end if
    end function round_up

end module neort_gc_outward_interval
