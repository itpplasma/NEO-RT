program test_gc_potato_class_coordinate_map
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_positive_inf, &
        ieee_quiet_nan
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use gc_potato_class_coordinate_map, only: &
        evaluate_gc_potato_class_coordinate_map, &
        GC_POTATO_CLASS_MAP_INVALID_BOUNDS, &
        GC_POTATO_CLASS_MAP_INVALID_IFUNTYPE, &
        GC_POTATO_CLASS_MAP_INVALID_INPUT, GC_POTATO_CLASS_MAP_OK
    use util_for_test, only: pass_test
    implicit none

    real(dp), parameter :: pi_value = 3.1415926535897932384626433832795_dp

    call test_representative_oracles()
    call test_all_boundary_pairs()
    call test_endpoint_margin()
    call test_logarithmic_tendencies()
    call test_orientation_reversal()
    call test_invalid_inputs()
    write (*, '(a)') 'test_gc_potato_class_coordinate_map OK'
    call pass_test

contains

    subroutine test_representative_oracles()
        integer, parameter :: pair_count = 11
        integer, parameter :: pairs(pair_count) = &
            [11, 12, 21, 22, 13, 31, 14, 41, 34, 43, 44]
        real(dp), parameter :: samples(pair_count) = &
            [0.31_dp, 0.41_dp, 0.58_dp, 0.43_dp, 0.80_dp, -0.30_dp, &
            0.90_dp, -0.20_dp, 0.20_dp, -0.20_dp, 0.10_dp]
        integer :: i, status
        real(dp) :: xi, dxi, xbeg, xend, expected_xi, expected_dxi
        real(dp) :: expected_xbeg, expected_xend

        do i = 1, pair_count
            call evaluate_gc_potato_class_coordinate_map(pairs(i), samples(i), &
                0.03_dp, 0.20_dp, xi, dxi, xbeg, xend, status)
            call require(status == GC_POTATO_CLASS_MAP_OK, &
                'representative pair status')
            expected_xi = oracle_map(pairs(i), samples(i))
            expected_dxi = oracle_derivative(pairs(i), samples(i))
            call oracle_bounds(pairs(i), 0.03_dp, 0.20_dp, expected_xbeg, &
                expected_xend)
            call require_close(xi, expected_xi, 2.0e-13_dp, &
                'representative xi')
            call require_close(dxi, expected_dxi, 2.0e-13_dp, &
                'representative derivative')
            call require_close(xbeg, expected_xbeg, 2.0e-13_dp, &
                'representative xbeg')
            call require_close(xend, expected_xend, 2.0e-13_dp, &
                'representative xend')
        end do
    end subroutine test_representative_oracles

    subroutine test_all_boundary_pairs()
        integer :: left_kind, right_kind, ifuntype, status
        real(dp) :: xi, dxi, xbeg, xend, expected_xi, expected_dxi
        real(dp) :: expected_xbeg, expected_xend, x

        do left_kind = 1, 4
            do right_kind = 1, 4
                ifuntype = 10*left_kind+right_kind
                call oracle_bounds(ifuntype, 0.02_dp, 0.30_dp, &
                    expected_xbeg, expected_xend)
                x = expected_xbeg+0.37_dp*(expected_xend-expected_xbeg)
                call evaluate_gc_potato_class_coordinate_map(ifuntype, x, &
                    0.02_dp, 0.30_dp, xi, dxi, xbeg, xend, status)
                call require(status == GC_POTATO_CLASS_MAP_OK, &
                    'all pair status')
                expected_xi = oracle_map(ifuntype, x)
                expected_dxi = oracle_derivative(ifuntype, x)
                call require_close(xi, expected_xi, 3.0e-13_dp, &
                    'all pair xi')
                call require_close(dxi, expected_dxi, 3.0e-13_dp, &
                    'all pair derivative')
                call require_close(xbeg, expected_xbeg, 3.0e-13_dp, &
                    'all pair xbeg')
                call require_close(xend, expected_xend, 3.0e-13_dp, &
                    'all pair xend')
            end do
        end do
    end subroutine test_all_boundary_pairs

    subroutine test_endpoint_margin()
        integer :: status
        real(dp) :: xi, dxi, xbeg, xend

        call evaluate_gc_potato_class_coordinate_map(12, 0.41_dp, 0.07_dp, &
            0.30_dp, xi, dxi, xbeg, xend, status)
        call require(status == GC_POTATO_CLASS_MAP_OK, '12 margin status')
        call require_close(xend, 0.93_dp, 2.0e-14_dp, '12 right margin')

        call evaluate_gc_potato_class_coordinate_map(21, 0.41_dp, 0.07_dp, &
            0.30_dp, xi, dxi, xbeg, xend, status)
        call require(status == GC_POTATO_CLASS_MAP_OK, '21 margin status')
        call require_close(xbeg, 0.07_dp, 2.0e-14_dp, '21 left margin')

        call evaluate_gc_potato_class_coordinate_map(22, 0.41_dp, 0.07_dp, &
            0.30_dp, xi, dxi, xbeg, xend, status)
        call require(status == GC_POTATO_CLASS_MAP_OK, '22 margin status')
        call require_close(xbeg, 0.07_dp, 2.0e-14_dp, '22 left margin')
        call require_close(xend, 0.93_dp, 2.0e-14_dp, '22 right margin')
    end subroutine test_endpoint_margin

    subroutine test_logarithmic_tendencies()
        integer, parameter :: singular_pairs(6) = [31, 32, 34, 41, 43, 44]
        integer :: i, status
        real(dp) :: xi_small, dxi, xbeg_small, xend_small
        real(dp) :: xi_large, xbeg_large, xend_large, x_small, x_large
        real(dp) :: xi_right_small, xi_right_large

        do i = 1, size(singular_pairs)
            call evaluate_gc_potato_class_coordinate_map(singular_pairs(i), &
                0.0_dp, 1.0e-2_dp, 1.0_dp, xi_large, dxi, xbeg_large, &
                xend_large, status)
            call require(status == GC_POTATO_CLASS_MAP_OK, &
                'large-ratio logarithmic status')
            call evaluate_gc_potato_class_coordinate_map(singular_pairs(i), &
                0.0_dp, 1.0e-4_dp, 1.0_dp, xi_small, dxi, xbeg_small, &
                xend_small, status)
            call require(status == GC_POTATO_CLASS_MAP_OK, &
                'small-ratio logarithmic status')
            if (singular_pairs(i)/10 >= 3) then
                x_large = xbeg_large+1.0e-6_dp
                x_small = xbeg_small+1.0e-6_dp
                call evaluate_gc_potato_class_coordinate_map(singular_pairs(i), &
                    x_large, 1.0e-2_dp, 1.0_dp, xi_large, dxi, xbeg_large, &
                    xend_large, status)
                call evaluate_gc_potato_class_coordinate_map(singular_pairs(i), &
                    x_small, 1.0e-4_dp, 1.0_dp, xi_small, dxi, xbeg_small, &
                    xend_small, status)
                call require(xi_small < xi_large, &
                    'left logarithmic limit tends to left endpoint')
            end if
            if (mod(singular_pairs(i), 10) >= 3) then
                x_large = xend_large-1.0e-6_dp
                x_small = xend_small-1.0e-6_dp
                call evaluate_gc_potato_class_coordinate_map(singular_pairs(i), &
                    x_large, 1.0e-2_dp, 1.0_dp, xi_right_large, dxi, &
                    xbeg_large, xend_large, status)
                call evaluate_gc_potato_class_coordinate_map(singular_pairs(i), &
                    x_small, 1.0e-4_dp, 1.0_dp, xi_right_small, dxi, &
                    xbeg_small, xend_small, status)
                call require(1.0_dp-xi_right_small < 1.0_dp-xi_right_large, &
                    'right logarithmic limit tends to right endpoint')
            end if
        end do
    end subroutine test_logarithmic_tendencies

    subroutine test_orientation_reversal()
        integer :: left_kind, right_kind, ifuntype, reverse_type, status
        real(dp) :: x, reverse_x, xi, reverse_xi, dxi, reverse_dxi
        real(dp) :: xbeg, xend, reverse_xbeg, reverse_xend

        do left_kind = 1, 4
            do right_kind = 1, 4
                ifuntype = 10*left_kind+right_kind
                reverse_type = 10*right_kind+left_kind
                x = orientation_sample(left_kind, right_kind)
                if (left_kind <= 2 .and. right_kind <= 2) then
                    reverse_x = 1.0_dp-x
                else
                    reverse_x = -x
                end if
                call evaluate_gc_potato_class_coordinate_map(ifuntype, x, &
                    0.02_dp, 0.30_dp, xi, dxi, xbeg, xend, status)
                call require(status == GC_POTATO_CLASS_MAP_OK, &
                    'orientation status')
                call evaluate_gc_potato_class_coordinate_map(reverse_type, &
                    reverse_x, 0.02_dp, 0.30_dp, reverse_xi, reverse_dxi, &
                    reverse_xbeg, reverse_xend, status)
                call require(status == GC_POTATO_CLASS_MAP_OK, &
                    'reverse orientation status')
                call require_close(xi+reverse_xi, 1.0_dp, 4.0e-13_dp, &
                    'map orientation reversal')
                call require_close(dxi, reverse_dxi, 4.0e-13_dp, &
                    'derivative orientation reversal')
            end do
        end do
    end subroutine test_orientation_reversal

    subroutine test_invalid_inputs()
        integer :: status
        real(dp) :: xi, dxi, xbeg, xend, nan_value, inf_value

        nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
        inf_value = ieee_value(1.0_dp, ieee_positive_inf)
        call evaluate_gc_potato_class_coordinate_map(10, 0.2_dp, 0.1_dp, &
            0.3_dp, xi, dxi, xbeg, xend, status)
        call require(status == GC_POTATO_CLASS_MAP_INVALID_IFUNTYPE, &
            'unknown ifuntype rejected')
        call evaluate_gc_potato_class_coordinate_map(45, 0.2_dp, 0.1_dp, &
            0.3_dp, xi, dxi, xbeg, xend, status)
        call require(status == GC_POTATO_CLASS_MAP_INVALID_IFUNTYPE, &
            'unknown right boundary rejected')
        call evaluate_gc_potato_class_coordinate_map(11, nan_value, 0.1_dp, &
            0.3_dp, xi, dxi, xbeg, xend, status)
        call require(status == GC_POTATO_CLASS_MAP_INVALID_INPUT, &
            'nonfinite x rejected')
        call evaluate_gc_potato_class_coordinate_map(11, 0.2_dp, inf_value, &
            0.3_dp, xi, dxi, xbeg, xend, status)
        call require(status == GC_POTATO_CLASS_MAP_INVALID_INPUT, &
            'nonfinite margin rejected')
        call evaluate_gc_potato_class_coordinate_map(11, 0.2_dp, 0.1_dp, &
            0.0_dp, xi, dxi, xbeg, xend, status)
        call require(status == GC_POTATO_CLASS_MAP_INVALID_INPUT, &
            'nonpositive width rejected')
        call evaluate_gc_potato_class_coordinate_map(22, 0.2_dp, 0.6_dp, &
            1.0_dp, xi, dxi, xbeg, xend, status)
        call require(status == GC_POTATO_CLASS_MAP_INVALID_BOUNDS, &
            'unordered bounded interval rejected')
        call evaluate_gc_potato_class_coordinate_map(34, 0.2_dp, 1.0_dp, &
            1.0_dp, xi, dxi, xbeg, xend, status)
        call require(status == GC_POTATO_CLASS_MAP_INVALID_BOUNDS, &
            'zero logarithmic interval rejected')
    end subroutine test_invalid_inputs

    pure real(dp) function oracle_map(ifuntype, x)
        integer, intent(in) :: ifuntype
        real(dp), intent(in) :: x
        integer :: left_kind, right_kind

        left_kind = ifuntype/10
        right_kind = mod(ifuntype, 10)
        select case (left_kind)
        case (1)
            select case (right_kind)
            case (1); oracle_map = x
            case (2); oracle_map = 1.0_dp-(1.0_dp-x)**2
            case (3, 4); oracle_map = tanh(x)
            end select
        case (2)
            select case (right_kind)
            case (1); oracle_map = x**2
            case (2); oracle_map = 0.5_dp*(1.0_dp-cos(pi_value*x))
            case (3, 4); oracle_map = tanh(x)**2
            end select
        case (3, 4)
            select case (right_kind)
            case (1); oracle_map = 1.0_dp+tanh(x)
            case (2); oracle_map = 1.0_dp-tanh(x)**2
            case (3, 4); oracle_map = 0.5_dp*(1.0_dp+tanh(x))
            end select
        end select
    end function oracle_map

    pure real(dp) function oracle_derivative(ifuntype, x)
        integer, intent(in) :: ifuntype
        real(dp), intent(in) :: x
        integer :: left_kind, right_kind
        real(dp) :: hyperbolic

        left_kind = ifuntype/10
        right_kind = mod(ifuntype, 10)
        hyperbolic = tanh(x)
        select case (left_kind)
        case (1)
            select case (right_kind)
            case (1); oracle_derivative = 1.0_dp
            case (2); oracle_derivative = 2.0_dp*(1.0_dp-x)
            case (3, 4); oracle_derivative = 1.0_dp/cosh(x)**2
            end select
        case (2)
            select case (right_kind)
            case (1); oracle_derivative = 2.0_dp*x
            case (2); oracle_derivative = 0.5_dp*pi_value*sin(pi_value*x)
            case (3, 4); oracle_derivative = 2.0_dp*hyperbolic/cosh(x)**2
            end select
        case (3, 4)
            select case (right_kind)
            case (1); oracle_derivative = 1.0_dp/cosh(x)**2
            case (2); oracle_derivative = -2.0_dp*hyperbolic/cosh(x)**2
            case (3, 4); oracle_derivative = 0.5_dp/cosh(x)**2
            end select
        end select
    end function oracle_derivative

    subroutine oracle_bounds(ifuntype, relative_margin, relative_width, &
            xbeg, xend)
        integer, intent(in) :: ifuntype
        real(dp), intent(in) :: relative_margin, relative_width
        real(dp), intent(out) :: xbeg, xend
        integer :: left_kind, right_kind
        real(dp) :: log_ratio

        left_kind = ifuntype/10
        right_kind = mod(ifuntype, 10)
        log_ratio = log(relative_margin/relative_width)
        select case (left_kind)
        case (1); xbeg = 0.0_dp
        case (2); xbeg = relative_margin
        case (3); xbeg = 0.5_dp*log_ratio
        case (4); xbeg = 0.25_dp*log_ratio
        end select
        if (right_kind >= 3) then
            if (right_kind == 3) xend = -0.5_dp*log_ratio
            if (right_kind == 4) xend = -0.25_dp*log_ratio
        else if (left_kind >= 3) then
            if (right_kind == 1) xend = 0.0_dp
            if (right_kind == 2) xend = -relative_margin
        else
            if (right_kind == 1) xend = 1.0_dp
            if (right_kind == 2) xend = 1.0_dp-relative_margin
        end if
    end subroutine oracle_bounds

    pure real(dp) function orientation_sample(left_kind, right_kind)
        integer, intent(in) :: left_kind, right_kind

        if (left_kind <= 2 .and. right_kind <= 2) then
            orientation_sample = 0.37_dp
        else if (left_kind <= 2) then
            orientation_sample = 0.41_dp
        else if (right_kind <= 2) then
            orientation_sample = -0.41_dp
        else
            orientation_sample = 0.37_dp
        end if
    end function orientation_sample

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message

        if (.not. condition) error stop message
    end subroutine require

    subroutine require_close(actual, expected, tolerance, message)
        real(dp), intent(in) :: actual, expected, tolerance
        character(*), intent(in) :: message
        real(dp) :: scale

        scale = max(1.0_dp, max(abs(actual), abs(expected)))
        call require(abs(actual-expected) <= tolerance*scale, message)
    end subroutine require_close

end program test_gc_potato_class_coordinate_map
