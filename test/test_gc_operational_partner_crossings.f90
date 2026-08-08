program test_gc_operational_partner_crossings
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_gc_operational_fixed_points, only: &
        GC_FIXED_POINT_SUCCESS, find_gc_operational_fixed_points, &
        gc_operational_fixed_point_result_t
    use neort_gc_operational_partner_crossings, only: &
        GC_PARTNER_BOUNDARY_USUAL, GC_PARTNER_BOUNDARY_X, &
        GC_PARTNER_MISSING, GC_PARTNER_SUCCESS, &
        find_gc_operational_partner_crossings, &
        gc_operational_partner_options_t, gc_operational_partner_result_t
    use util_for_test, only: pass_test
    implicit none

    type(gc_operational_fixed_point_result_t) :: fixed_points
    type(gc_operational_partner_options_t) :: options
    type(gc_operational_partner_result_t) :: partners

    options%root%initial_intervals = 11
    options%x_exclusion_fraction = 1.0e-7_dp
    call find_gc_operational_fixed_points(partner_polynomial, &
        partner_classifier, 0.0_dp, 4.0_dp, options%root, fixed_points)
    call require(fixed_points%status == GC_FIXED_POINT_SUCCESS, &
        'manufactured fixed-point search failed')
    call require(fixed_points%n_x_points == 1, &
        'manufactured X point was not isolated')

    call find_gc_operational_partner_crossings(partner_polynomial, &
        fixed_points, 0.0_dp, 4.0_dp, options, partners)
    call require(partners%status == GC_PARTNER_SUCCESS, &
        'same-p_phi partner search failed')
    call require(partners%complete, 'partner set is incomplete')
    call require(partners%npairs == 1 .and. partners%nboundaries == 2, &
        'partner pair cardinality is wrong')
    call require(partners%boundaries(1)%kind == GC_PARTNER_BOUNDARY_X .and. &
        partners%boundaries(2)%kind == GC_PARTNER_BOUNDARY_USUAL, &
        'X and usual partner roles are wrong')
    call require(abs(partners%boundaries(1)%x-1.0_dp) < 1.0e-9_dp .and. &
        abs(partners%boundaries(2)%x-3.0_dp) < 1.0e-9_dp, &
        'same-p_phi partner locations are wrong')
    call require(abs(partners%boundaries(1)%canonical_momentum- &
        partners%boundaries(2)%canonical_momentum) < 1.0e-10_dp, &
        'paired boundaries do not share canonical momentum')
    call require(abs(partners%boundaries(2)%canonical_residual) < 1.0e-9_dp, &
        'regular partner residual did not converge')

    call find_gc_operational_fixed_points(no_partner_polynomial, &
        no_partner_classifier, 0.0_dp, 2.0_dp, options%root, fixed_points)
    call require(fixed_points%status == GC_FIXED_POINT_SUCCESS .and. &
        fixed_points%n_x_points == 1, 'no-partner X fixture is invalid')
    call find_gc_operational_partner_crossings(no_partner_polynomial, &
        fixed_points, 0.0_dp, 2.0_dp, options, partners)
    call require(partners%status == GC_PARTNER_MISSING, &
        'unpaired X point was accepted')
    call require(.not. partners%complete, &
        'unpaired X-point topology was marked complete')

    call pass_test

contains

    subroutine partner_polynomial(x, p_star, first, second, status)
        real(dp), intent(in) :: x
        real(dp), intent(out) :: p_star, first, second
        integer, intent(out) :: status

        p_star = (x-1.0_dp)**2*(x-3.0_dp)
        first = (x-1.0_dp)*(3.0_dp*x-7.0_dp)
        second = 6.0_dp*x-10.0_dp
        status = 0
    end subroutine partner_polynomial

    subroutine partner_classifier(x, discriminant, status)
        real(dp), intent(in) :: x
        real(dp), intent(out) :: discriminant
        integer, intent(out) :: status

        discriminant = 4.0_dp*(2.0_dp-x)
        status = 0
    end subroutine partner_classifier

    subroutine no_partner_polynomial(x, p_star, first, second, status)
        real(dp), intent(in) :: x
        real(dp), intent(out) :: p_star, first, second
        integer, intent(out) :: status

        p_star = (x-1.0_dp)**2
        first = 2.0_dp*(x-1.0_dp)
        second = 2.0_dp
        status = 0
    end subroutine no_partner_polynomial

    subroutine no_partner_classifier(x, discriminant, status)
        real(dp), intent(in) :: x
        real(dp), intent(out) :: discriminant
        integer, intent(out) :: status

        discriminant = 4.0_dp
        status = 0
        associate (unused_x => x)
        end associate
    end subroutine no_partner_classifier

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) error stop message
    end subroutine require

end program test_gc_operational_partner_crossings
