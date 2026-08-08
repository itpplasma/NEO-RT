module neort_gc_unit_quadrature
    !! Fixed Gauss-Legendre rules on [0,1].
    !!
    !! Fortnum owns the canonical [-1,1] rule.  Fortsym owns the affine
    !! coordinate and measure push-forward.  This module only composes those
    !! two providers and validates their finite, positive contract.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_quadrature, only: gauss_legendre
    use neort_gauss_legendre_interval_map_symbolic, only: &
        evaluate_neort_gauss_legendre_interval_map
    implicit none
    private

    integer, parameter, public :: GC_UNIT_QUADRATURE_SUCCESS = 0
    integer, parameter, public :: GC_UNIT_QUADRATURE_INVALID_ORDER = 1
    integer, parameter, public :: GC_UNIT_QUADRATURE_INVALID_RULE = 2

    public :: build_gc_unit_gauss_legendre

contains

    subroutine build_gc_unit_gauss_legendre(order, nodes, weights, status)
        integer, intent(in) :: order
        real(dp), allocatable, intent(out) :: nodes(:), weights(:)
        integer, intent(out) :: status

        real(dp) :: standard_node, standard_weight
        integer :: i

        status = GC_UNIT_QUADRATURE_INVALID_ORDER
        if (order < 1) return
        allocate(nodes(order), weights(order))
        call gauss_legendre(order, nodes, weights)
        do i = 1, order
            standard_node = nodes(i)
            standard_weight = weights(i)
            call evaluate_neort_gauss_legendre_interval_map(standard_node, &
                standard_weight, 0.0_dp, 1.0_dp, nodes(i), weights(i))
        end do
        status = GC_UNIT_QUADRATURE_INVALID_RULE
        if (.not. all(ieee_is_finite(nodes))) return
        if (.not. all(ieee_is_finite(weights))) return
        if (any(nodes <= 0.0_dp) .or. any(nodes >= 1.0_dp)) return
        if (any(weights <= 0.0_dp)) return
        if (order > 1) then
            if (any(nodes(2:) <= nodes(:order-1))) return
        end if
        status = GC_UNIT_QUADRATURE_SUCCESS
    end subroutine build_gc_unit_gauss_legendre

end module neort_gc_unit_quadrature
