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

contains

    subroutine evaluate_gc_potato_class_coordinate_map(ifuntype, x, &
            relative_margin, relative_class_width, xi, dxi_dx, xbeg, xend, &
            status)
        integer, intent(in) :: ifuntype
        real(dp), intent(in) :: x, relative_margin, relative_class_width
        real(dp), intent(out) :: xi, dxi_dx, xbeg, xend
        integer, intent(out) :: status
        integer :: left_kind, right_kind
        real(dp) :: all_xi(4, 4), all_dxi(4, 4), all_xbeg(4, 4), all_xend(4, 4)
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
            all_xend(4, 1), all_xend(4, 2), all_xend(4, 3), all_xend(4, 4))

        xi = all_xi(left_kind, right_kind)
        dxi_dx = all_dxi(left_kind, right_kind)
        xbeg = all_xbeg(left_kind, right_kind)
        xend = all_xend(left_kind, right_kind)
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

end module gc_potato_class_coordinate_map
