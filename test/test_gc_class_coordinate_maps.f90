program test_gc_class_coordinate_maps
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none

    real(dp), parameter :: pi_value = 3.1415926535897932384626433832795_dp

    call test_numeric_inverses()
    call test_endpoint_truncation()
    call test_section_orientation()
    write (*, '(a)') 'test_gc_class_coordinate_maps OK'

contains

    subroutine test_numeric_inverses()
        integer :: left_type, right_type, sample
        real(dp) :: relmargin, widthclass, xbeg, xend, x, target, inverse

        relmargin = 1.0e-6_dp
        widthclass = 0.2_dp
        do left_type = 1, 4
            do right_type = 1, 4
                call independent_bounds(left_type, right_type, relmargin, &
                    widthclass, xbeg, xend)
                do sample = 1, 7
                    x = xbeg + (xend-xbeg)*real(sample, dp)/8.0_dp
                    target = map_value(left_type, right_type, x)
                    call require(map_derivative(left_type, right_type, x) > 0.0_dp, &
                        'class-coordinate monotonicity')
                    inverse = bisect_inverse(left_type, right_type, target, &
                        xbeg, xend)
                    call require(inverse >= xbeg .and. inverse <= xend, &
                        'numeric inverse bracket')
                    call require_relative(map_value(left_type, right_type, &
                        inverse), target, 5.0e-11_dp, 'numeric inverse residual')
                end do
            end do
        end do
    end subroutine test_numeric_inverses

    subroutine test_endpoint_truncation()
        integer :: left_type, right_type
        real(dp) :: relmargin, widthclass, ratio, xbeg, xend, coefficient
        real(dp) :: tail, expected_tail

        relmargin = 1.0e-7_dp
        widthclass = 0.3_dp
        ratio = relmargin/widthclass
        do left_type = 1, 4
            do right_type = 1, 4
                call independent_bounds(left_type, right_type, relmargin, &
                    widthclass, xbeg, xend)
                if (left_type >= 3) then
                    coefficient = logarithmic_coefficient(left_type)
                    tail = map_value(left_type, right_type, xbeg)
                    expected_tail = left_tail_prefactor(right_type) &
                        *ratio**(2.0_dp*coefficient)
                    call require_relative(tail, expected_tail, 3.0e-3_dp, &
                        'left endpoint asymptotic')
                end if
                if (right_type >= 3) then
                    coefficient = logarithmic_coefficient(right_type)
                    tail = 1.0_dp-map_value(left_type, right_type, xend)
                    expected_tail = right_tail_prefactor(left_type) &
                        *ratio**(2.0_dp*coefficient)
                    call require_relative(tail, expected_tail, 3.0e-3_dp, &
                        'right endpoint asymptotic')
                end if
            end do
        end do
    end subroutine test_endpoint_truncation

    subroutine test_section_orientation()
        integer :: left_type, right_type, sample
        real(dp) :: x, reverse_x, coordinate, reversed_coordinate
        real(dp) :: derivative, reversed_derivative

        do left_type = 1, 4
            do right_type = 1, 4
                do sample = 1, 4
                    x = orientation_sample(left_type, right_type, sample)
                    reverse_x = reverse_coordinate(left_type, right_type, x)
                    coordinate = section_coordinate(left_type, right_type, &
                        x, 1.2_dp, 1.8_dp)
                    reversed_coordinate = section_coordinate(right_type, &
                        left_type, reverse_x, 1.8_dp, 1.2_dp)
                    derivative = section_derivative(left_type, right_type, &
                        x, 1.2_dp, 1.8_dp)
                    reversed_derivative = section_derivative(right_type, &
                        left_type, reverse_x, 1.8_dp, 1.2_dp)
                    call require_close(coordinate, reversed_coordinate, &
                        2.0e-14_dp, 'section orientation')
                    call require_close(derivative, -reversed_derivative, &
                        2.0e-14_dp, 'section derivative orientation')
                end do
            end do
        end do
    end subroutine test_section_orientation

    pure real(dp) function map_value(left_type, right_type, x)
        integer, intent(in) :: left_type, right_type
        real(dp), intent(in) :: x

        select case (left_type)
        case (1)
            select case (right_type)
            case (1); map_value = x
            case (2); map_value = 1.0_dp-(1.0_dp-x)**2
            case default; map_value = tanh(x)
            end select
        case (2)
            select case (right_type)
            case (1); map_value = x**2
            case (2); map_value = 0.5_dp*(1.0_dp-cos(pi_value*x))
            case default; map_value = tanh(x)**2
            end select
        case (3, 4)
            select case (right_type)
            case (1); map_value = 1.0_dp+tanh(x)
            case (2); map_value = 1.0_dp-tanh(x)**2
            case default; map_value = 0.5_dp*(1.0_dp+tanh(x))
            end select
        end select
    end function map_value

    pure real(dp) function map_derivative(left_type, right_type, x)
        integer, intent(in) :: left_type, right_type
        real(dp), intent(in) :: x
        real(dp) :: hyperbolic

        hyperbolic = tanh(x)
        select case (left_type)
        case (1)
            select case (right_type)
            case (1); map_derivative = 1.0_dp
            case (2); map_derivative = 2.0_dp*(1.0_dp-x)
            case default; map_derivative = 1.0_dp/cosh(x)**2
            end select
        case (2)
            select case (right_type)
            case (1); map_derivative = 2.0_dp*x
            case (2); map_derivative = 0.5_dp*pi_value*sin(pi_value*x)
            case default; map_derivative = 2.0_dp*hyperbolic/cosh(x)**2
            end select
        case (3, 4)
            select case (right_type)
            case (1); map_derivative = 1.0_dp/cosh(x)**2
            case (2); map_derivative = -2.0_dp*hyperbolic/cosh(x)**2
            case default; map_derivative = 0.5_dp/cosh(x)**2
            end select
        end select
    end function map_derivative

    subroutine independent_bounds(left_type, right_type, relmargin, &
            widthclass, xbeg, xend)
        integer, intent(in) :: left_type, right_type
        real(dp), intent(in) :: relmargin, widthclass
        real(dp), intent(out) :: xbeg, xend
        real(dp) :: log_ratio

        log_ratio = log(relmargin/widthclass)
        xbeg = finite_left_bound(left_type, relmargin, log_ratio)
        xend = finite_right_bound(left_type, right_type, relmargin, log_ratio)
    end subroutine independent_bounds

    pure real(dp) function finite_left_bound(boundary_type, relmargin, log_ratio)
        integer, intent(in) :: boundary_type
        real(dp), intent(in) :: relmargin, log_ratio

        select case (boundary_type)
        case (1); finite_left_bound = 0.0_dp
        case (2); finite_left_bound = relmargin
        case (3); finite_left_bound = 0.5_dp*log_ratio
        case (4); finite_left_bound = 0.25_dp*log_ratio
        end select
    end function finite_left_bound

    pure real(dp) function finite_right_bound(left_type, boundary_type, &
            relmargin, log_ratio)
        integer, intent(in) :: left_type, boundary_type
        real(dp), intent(in) :: relmargin, log_ratio

        if (boundary_type >= 3) then
            select case (boundary_type)
            case (3); finite_right_bound = -0.5_dp*log_ratio
            case (4); finite_right_bound = -0.25_dp*log_ratio
            end select
        else if (left_type >= 3) then
            select case (boundary_type)
            case (1); finite_right_bound = 0.0_dp
            case (2); finite_right_bound = -relmargin
            end select
        else
            select case (boundary_type)
            case (1); finite_right_bound = 1.0_dp
            case (2); finite_right_bound = 1.0_dp-relmargin
            end select
        end if
    end function finite_right_bound

    pure real(dp) function logarithmic_coefficient(boundary_type)
        integer, intent(in) :: boundary_type

        if (boundary_type == 3) then
            logarithmic_coefficient = 0.5_dp
        else
            logarithmic_coefficient = 0.25_dp
        end if
    end function logarithmic_coefficient

    pure real(dp) function left_tail_prefactor(right_type)
        integer, intent(in) :: right_type

        if (right_type == 1) then
            left_tail_prefactor = 2.0_dp
        else if (right_type == 2) then
            left_tail_prefactor = 4.0_dp
        else
            left_tail_prefactor = 1.0_dp
        end if
    end function left_tail_prefactor

    pure real(dp) function right_tail_prefactor(left_type)
        integer, intent(in) :: left_type

        if (left_type == 1) then
            right_tail_prefactor = 2.0_dp
        else if (left_type == 2) then
            right_tail_prefactor = 4.0_dp
        else
            right_tail_prefactor = 1.0_dp
        end if
    end function right_tail_prefactor

    pure real(dp) function bisect_inverse(left_type, right_type, target, lo, hi)
        integer, intent(in) :: left_type, right_type
        real(dp), intent(in) :: target, lo, hi
        integer :: iteration
        real(dp) :: lower, upper, middle

        lower = lo
        upper = hi
        do iteration = 1, 100
            middle = 0.5_dp*(lower+upper)
            if (map_value(left_type, right_type, middle) < target) then
                lower = middle
            else
                upper = middle
            end if
        end do
        bisect_inverse = 0.5_dp*(lower+upper)
    end function bisect_inverse

    pure real(dp) function reverse_coordinate(left_type, right_type, x)
        integer, intent(in) :: left_type, right_type
        real(dp), intent(in) :: x

        if (left_type <= 2 .and. right_type <= 2) then
            reverse_coordinate = 1.0_dp-x
        else
            reverse_coordinate = -x
        end if
    end function reverse_coordinate

    pure real(dp) function section_coordinate(left_type, right_type, x, rlo, rhi)
        integer, intent(in) :: left_type, right_type
        real(dp), intent(in) :: x, rlo, rhi

        section_coordinate = rlo+(rhi-rlo)*map_value(left_type, right_type, x)
    end function section_coordinate

    pure real(dp) function section_derivative(left_type, right_type, x, rlo, rhi)
        integer, intent(in) :: left_type, right_type
        real(dp), intent(in) :: x, rlo, rhi

        section_derivative = (rhi-rlo)*map_derivative(left_type, right_type, x)
    end function section_derivative

    pure real(dp) function orientation_sample(left_type, right_type, sample)
        integer, intent(in) :: left_type, right_type, sample
        real(dp), parameter :: finite_samples(4) = [0.13_dp, 0.37_dp, &
            0.63_dp, 0.87_dp]
        real(dp), parameter :: singular_samples(4) = [-1.2_dp, -0.3_dp, &
            0.3_dp, 1.2_dp]

        if (left_type <= 2 .and. right_type <= 2) then
            orientation_sample = finite_samples(sample)
        else if (left_type <= 2) then
            orientation_sample = abs(singular_samples(sample))
        else if (right_type <= 2) then
            orientation_sample = -abs(singular_samples(sample))
        else
            orientation_sample = singular_samples(sample)
        end if
    end function orientation_sample

    subroutine require_close(actual, expected, tolerance, label)
        real(dp), intent(in) :: actual, expected, tolerance
        character(*), intent(in) :: label

        call require(abs(actual-expected) <= tolerance, label)
    end subroutine require_close

    subroutine require_relative(actual, expected, tolerance, label)
        real(dp), intent(in) :: actual, expected, tolerance
        character(*), intent(in) :: label
        real(dp) :: scale

        scale = max(abs(expected), tiny(1.0_dp))
        call require(abs(actual-expected)/scale <= tolerance, label)
    end subroutine require_relative

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message

        if (.not. condition) error stop message
    end subroutine require

end program test_gc_class_coordinate_maps
