module gc_potato_class_coordinate_map
    !! Runtime boundary/status seam for the generated POTATO class map.
    !!
    !! Boundary kind 1 is a regular rho boundary, 2 an inner/turning
    !! boundary, 3 a regular separatrix, and 4 an X-point.  ifuntype is
    !! 10*left_kind+right_kind.  The generated map is increasing in x; 34
    !! and 43 retain POTATO's same increasing logistic map and differ only in
    !! which physical endpoint is named left/right.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, &
        ieee_quiet_nan
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_gc_class_coordinate_map_generated, only: &
        evaluate_neort_gc_class_coordinate_map
    implicit none
    private

    integer, parameter, public :: GC_POTATO_CLASS_MAP_OK = 0
    integer, parameter, public :: GC_POTATO_CLASS_MAP_INVALID_IFUNTYPE = 1
    integer, parameter, public :: GC_POTATO_CLASS_MAP_INVALID_INPUT = 2
    integer, parameter, public :: GC_POTATO_CLASS_MAP_INVALID_BOUNDS = 3
    integer, parameter, public :: GC_POTATO_CLASS_MAP_NONFINITE_OUTPUT = 4

    public :: evaluate_gc_potato_class_coordinate_map
    public :: evaluate_gc_potato_class_bounds
    public :: evaluate_gc_potato_class_xpoint_limits

contains

    subroutine evaluate_gc_potato_class_coordinate_map(ifuntype, x, &
            relative_margin, relative_class_width, xi, dxi_dx, xbeg, xend, &
            status, xi_left_xpoint_limit, xi_right_xpoint_limit)
        integer, intent(in) :: ifuntype
        real(dp), intent(in) :: x, relative_margin, relative_class_width
        real(dp), intent(out) :: xi, dxi_dx, xbeg, xend
        real(dp), intent(out), optional :: xi_left_xpoint_limit, &
            xi_right_xpoint_limit
        integer, intent(out) :: status
        integer :: left_kind, right_kind
        real(dp) :: all_xi(4, 4), all_dxi(4, 4), all_xbeg(4, 4), all_xend(4, 4)
        real(dp) :: all_xpoint_left(4, 4), all_xpoint_right(4, 4)
        real(dp) :: nan_value

        nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
        xi = nan_value
        dxi_dx = nan_value
        xbeg = nan_value
        xend = nan_value
        status = GC_POTATO_CLASS_MAP_OK

        left_kind = ifuntype/10
        right_kind = mod(ifuntype, 10)
        if (left_kind < 1 .or. left_kind > 4) then
            status = GC_POTATO_CLASS_MAP_INVALID_IFUNTYPE
            return
        end if
        if (right_kind < 1 .or. right_kind > 4) then
            status = GC_POTATO_CLASS_MAP_INVALID_IFUNTYPE
            return
        end if
        if (.not. ieee_is_finite(x)) then
            status = GC_POTATO_CLASS_MAP_INVALID_INPUT
            return
        end if
        if (.not. ieee_is_finite(relative_margin)) then
            status = GC_POTATO_CLASS_MAP_INVALID_INPUT
            return
        end if
        if (.not. ieee_is_finite(relative_class_width)) then
            status = GC_POTATO_CLASS_MAP_INVALID_INPUT
            return
        end if
        if (relative_margin <= 0.0_dp) then
            status = GC_POTATO_CLASS_MAP_INVALID_INPUT
            return
        end if
        if (relative_class_width <= 0.0_dp) then
            status = GC_POTATO_CLASS_MAP_INVALID_INPUT
            return
        end if

        call evaluate_neort_gc_class_coordinate_map(x, relative_margin, &
            relative_class_width, &
            all_xi(1, 1), all_xi(1, 2), all_xi(1, 3), all_xi(1, 4), &
            all_xi(2, 1), all_xi(2, 2), all_xi(2, 3), all_xi(2, 4), &
            all_xi(3, 1), all_xi(3, 2), all_xi(3, 3), all_xi(3, 4), &
            all_xi(4, 1), all_xi(4, 2), all_xi(4, 3), all_xi(4, 4), &
            all_dxi(1, 1), all_dxi(1, 2), all_dxi(1, 3), all_dxi(1, 4), &
            all_dxi(2, 1), all_dxi(2, 2), all_dxi(2, 3), all_dxi(2, 4), &
            all_dxi(3, 1), all_dxi(3, 2), all_dxi(3, 3), all_dxi(3, 4), &
            all_dxi(4, 1), all_dxi(4, 2), all_dxi(4, 3), all_dxi(4, 4), &
            all_xbeg(1, 1), all_xbeg(1, 2), all_xbeg(1, 3), all_xbeg(1, 4), &
            all_xbeg(2, 1), all_xbeg(2, 2), all_xbeg(2, 3), all_xbeg(2, 4), &
            all_xbeg(3, 1), all_xbeg(3, 2), all_xbeg(3, 3), all_xbeg(3, 4), &
            all_xbeg(4, 1), all_xbeg(4, 2), all_xbeg(4, 3), all_xbeg(4, 4), &
            all_xend(1, 1), all_xend(1, 2), all_xend(1, 3), all_xend(1, 4), &
            all_xend(2, 1), all_xend(2, 2), all_xend(2, 3), all_xend(2, 4), &
            all_xend(3, 1), all_xend(3, 2), all_xend(3, 3), all_xend(3, 4), &
            all_xend(4, 1), all_xend(4, 2), all_xend(4, 3), all_xend(4, 4), &
            all_xpoint_left(1, 1), all_xpoint_left(1, 2), &
            all_xpoint_left(1, 3), all_xpoint_left(1, 4), &
            all_xpoint_left(2, 1), all_xpoint_left(2, 2), &
            all_xpoint_left(2, 3), all_xpoint_left(2, 4), &
            all_xpoint_left(3, 1), all_xpoint_left(3, 2), &
            all_xpoint_left(3, 3), all_xpoint_left(3, 4), &
            all_xpoint_left(4, 1), all_xpoint_left(4, 2), &
            all_xpoint_left(4, 3), all_xpoint_left(4, 4), &
            all_xpoint_right(1, 1), all_xpoint_right(1, 2), &
            all_xpoint_right(1, 3), all_xpoint_right(1, 4), &
            all_xpoint_right(2, 1), all_xpoint_right(2, 2), &
            all_xpoint_right(2, 3), all_xpoint_right(2, 4), &
            all_xpoint_right(3, 1), all_xpoint_right(3, 2), &
            all_xpoint_right(3, 3), all_xpoint_right(3, 4), &
            all_xpoint_right(4, 1), all_xpoint_right(4, 2), &
            all_xpoint_right(4, 3), all_xpoint_right(4, 4))

        xi = all_xi(left_kind, right_kind)
        dxi_dx = all_dxi(left_kind, right_kind)
        xbeg = all_xbeg(left_kind, right_kind)
        xend = all_xend(left_kind, right_kind)
        if (present(xi_left_xpoint_limit)) then
            xi_left_xpoint_limit = all_xpoint_left(left_kind, right_kind)
        end if
        if (present(xi_right_xpoint_limit)) then
            xi_right_xpoint_limit = all_xpoint_right(left_kind, right_kind)
        end if
        if (.not. ieee_is_finite(xi)) then
            status = GC_POTATO_CLASS_MAP_NONFINITE_OUTPUT
            return
        end if
        if (.not. ieee_is_finite(dxi_dx)) then
            status = GC_POTATO_CLASS_MAP_NONFINITE_OUTPUT
            return
        end if
        if (.not. ieee_is_finite(xbeg)) then
            status = GC_POTATO_CLASS_MAP_NONFINITE_OUTPUT
            return
        end if
        if (.not. ieee_is_finite(xend)) then
            status = GC_POTATO_CLASS_MAP_NONFINITE_OUTPUT
            return
        end if
        if (xbeg >= xend) then
            status = GC_POTATO_CLASS_MAP_INVALID_BOUNDS
            return
        end if
    end subroutine evaluate_gc_potato_class_coordinate_map

    subroutine evaluate_gc_potato_class_xpoint_limits(ifuntype, relative_margin, &
            relative_class_width, xi_left_xpoint_limit, xi_right_xpoint_limit, &
            status)
        integer, intent(in) :: ifuntype
        real(dp), intent(in) :: relative_margin, relative_class_width
        real(dp), intent(out) :: xi_left_xpoint_limit, xi_right_xpoint_limit
        integer, intent(out) :: status
        real(dp) :: xi_unused, dxi_dx_unused, xbeg_unused, xend_unused

        call evaluate_gc_potato_class_coordinate_map(ifuntype, 0.0_dp, &
            relative_margin, relative_class_width, xi_unused, dxi_dx_unused, &
            xbeg_unused, xend_unused, status, xi_left_xpoint_limit, &
            xi_right_xpoint_limit)
    end subroutine evaluate_gc_potato_class_xpoint_limits

    subroutine evaluate_gc_potato_class_bounds(ifuntype, relative_margin, &
            relative_class_width, xbeg, xend, status)
        integer, intent(in) :: ifuntype
        real(dp), intent(in) :: relative_margin, relative_class_width
        real(dp), intent(out) :: xbeg, xend
        integer, intent(out) :: status
        real(dp) :: xi_unused, dxi_dx_unused

        ! The generated table is shared with the coordinate evaluator.  Its
        ! bounds are independent of class_coordinate; this finite probe only
        ! avoids exposing that implementation detail to callers.
        call evaluate_gc_potato_class_coordinate_map(ifuntype, 0.0_dp, &
            relative_margin, relative_class_width, xi_unused, &
            dxi_dx_unused, xbeg, xend, status)
    end subroutine evaluate_gc_potato_class_bounds

end module gc_potato_class_coordinate_map
