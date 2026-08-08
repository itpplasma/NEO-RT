program test_gc_unit_quadrature
    !! Independent behavioral oracle: an n-point Gauss rule must integrate
    !! every monomial through degree 2n-1 on [0,1] exactly to 1/(p+1).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_gc_unit_quadrature, only: &
        GC_UNIT_QUADRATURE_INVALID_ORDER, GC_UNIT_QUADRATURE_SUCCESS, &
        build_gc_unit_gauss_legendre
    implicit none

    integer, parameter :: orders(4) = [2, 3, 5, 8]
    real(dp), allocatable :: nodes(:), weights(:)
    real(dp) :: actual, expected, tolerance
    integer :: case_index, order, power, status

    tolerance = 2.0e-13_dp
    do case_index = 1, size(orders)
        order = orders(case_index)
        call build_gc_unit_gauss_legendre(order, nodes, weights, status)
        call require(status == GC_UNIT_QUADRATURE_SUCCESS, &
            'unit Gauss rule construction failed')
        do power = 0, 2*order-1
            actual = sum(weights*nodes**power)
            expected = 1.0_dp/real(power+1, dp)
            call require(abs(actual-expected) <= tolerance, &
                'unit Gauss rule failed a polynomial moment')
        end do
    end do

    call build_gc_unit_gauss_legendre(0, nodes, weights, status)
    call require(status == GC_UNIT_QUADRATURE_INVALID_ORDER, &
        'invalid Gauss order was accepted')
    call require(.not. allocated(nodes) .and. .not. allocated(weights), &
        'invalid Gauss order allocated a rule')

    write (*, '(a)') 'test_gc_unit_quadrature: PASS'

contains

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message

        if (.not. condition) error stop message
    end subroutine require

end program test_gc_unit_quadrature
