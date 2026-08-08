program test_gc_eqdsk_cut_endpoint_krawczyk
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_eqdsk_cut_endpoint_krawczyk_interval_symbolic, only: &
        evaluate_neort_eqdsk_cut_endpoint_krawczyk_interval
    use neort_eqdsk_cut_endpoint_newton_symbolic, only: &
        evaluate_neort_eqdsk_cut_endpoint_newton
    use neort_eqdsk_cut_endpoint_system_symbolic, only: &
        evaluate_neort_eqdsk_cut_endpoint_system
    use neort_gc_outward_interval, only: gc_outward_interval, &
        gc_outward_interval_t
    implicit none

    real(dp), parameter :: tolerance = 2.0e-14_dp
    real(dp) :: cut_residual, flux_residual
    real(dp) :: cut_jacobian_r, cut_jacobian_z
    real(dp) :: flux_jacobian_r, flux_jacobian_z
    real(dp) :: determinant, newton_r, newton_z
    real(dp) :: inverse_11, inverse_12, inverse_21, inverse_22
    type(gc_outward_interval_t) :: krawczyk_r, krawczyk_z
    type(gc_outward_interval_t) :: nonlinear_r, nonlinear_z
    type(gc_outward_interval_t) :: reversed_r, reversed_z

    ! Independent system oracle.  With psi=6, psi_sep=3, and target=1/2,
    ! the normalized-flux residual is 3/2 and its gradient is (4/3,-2/3).
    call evaluate_neort_eqdsk_cut_endpoint_system(0.1_dp, 2.0_dp, 1.0_dp, &
        6.0_dp, 4.0_dp, -2.0_dp, 3.0_dp, 0.5_dp, cut_residual, &
        flux_residual, cut_jacobian_r, cut_jacobian_z, flux_jacobian_r, &
        flux_jacobian_z)
    call require_close(cut_residual, 0.1_dp, 'cut residual')
    call require_close(flux_residual, 1.5_dp, 'flux residual')
    call require_close(cut_jacobian_r, 2.0_dp, 'cut R Jacobian')
    call require_close(cut_jacobian_z, 1.0_dp, 'cut Z Jacobian')
    call require_close(flux_jacobian_r, 4.0_dp/3.0_dp, &
        'flux R Jacobian')
    call require_close(flux_jacobian_z, -2.0_dp/3.0_dp, &
        'flux Z Jacobian')

    ! Independent affine-root oracle:
    !   N=2(R-1)+(Z-2), g=-(R-1)+3(Z-2).
    ! At (1.2,1.7), f=(0.1,-1.1), J=[[2,1],[-1,3]], det(J)=7,
    ! so one Newton step and the exact affine Krawczyk map both give (1,2).
    call evaluate_neort_eqdsk_cut_endpoint_newton(1.2_dp, 1.7_dp, 0.1_dp, &
        -1.1_dp, 2.0_dp, 1.0_dp, -1.0_dp, 3.0_dp, determinant, newton_r, &
        newton_z, inverse_11, inverse_12, inverse_21, inverse_22)
    call require_close(determinant, 7.0_dp, 'Jacobian determinant')
    call require_close(newton_r, 1.0_dp, 'affine Newton R')
    call require_close(newton_z, 2.0_dp, 'affine Newton Z')
    call require_close(inverse_11, 3.0_dp/7.0_dp, 'inverse 11')
    call require_close(inverse_12, -1.0_dp/7.0_dp, 'inverse 12')
    call require_close(inverse_21, 1.0_dp/7.0_dp, 'inverse 21')
    call require_close(inverse_22, 2.0_dp/7.0_dp, 'inverse 22')

    call evaluate_neort_eqdsk_cut_endpoint_krawczyk_interval( &
        point(1.2_dp), point(1.7_dp), point(0.1_dp), point(-1.1_dp), &
        point(2.0_dp), point(1.0_dp), point(-1.0_dp), point(3.0_dp), &
        point(2.0_dp), point(1.0_dp), point(-1.0_dp), point(3.0_dp), &
        gc_outward_interval(-0.2_dp, 0.2_dp), &
        gc_outward_interval(-0.3_dp, 0.3_dp), krawczyk_r, krawczyk_z)
    call require(krawczyk_r%lo <= 1.0_dp .and. krawczyk_r%hi >= 1.0_dp, &
        'affine Krawczyk R enclosure excluded the exact root')
    call require(krawczyk_z%lo <= 2.0_dp .and. krawczyk_z%hi >= 2.0_dp, &
        'affine Krawczyk Z enclosure excluded the exact root')
    call require(krawczyk_r%hi-krawczyk_r%lo < 1.0e-12_dp, &
        'affine Krawczyk R map did not contract to roundoff')
    call require(krawczyk_z%hi-krawczyk_z%lo < 1.0e-12_dp, &
        'affine Krawczyk Z map did not contract to roundoff')

    ! Nonlinear independent oracle with the same exact root (R,Z)=(1,2):
    !   N=2u+v+u^2, g=-u+3v+v^2, u=R-1, v=Z-2.
    ! At x0=(1.05,1.95), f=(0.0525,-0.1975) and
    ! J0=[[2.1,1],[-1,2.9]].  Over X=x0+[-0.1,0.1]^2, only the two
    ! diagonal derivatives vary, giving N_R=[1.9,2.3], g_Z=[2.7,3.1].
    call evaluate_neort_eqdsk_cut_endpoint_krawczyk_interval( &
        point(1.05_dp), point(1.95_dp), point(0.0525_dp), &
        point(-0.1975_dp), point(2.1_dp), point(1.0_dp), &
        point(-1.0_dp), point(2.9_dp), &
        gc_outward_interval(1.9_dp, 2.3_dp), point(1.0_dp), &
        point(-1.0_dp), gc_outward_interval(2.7_dp, 3.1_dp), &
        gc_outward_interval(-0.1_dp, 0.1_dp), &
        gc_outward_interval(-0.1_dp, 0.1_dp), nonlinear_r, nonlinear_z)
    call require(nonlinear_r%lo <= 1.0_dp .and. nonlinear_r%hi >= 1.0_dp, &
        'nonlinear Krawczyk R enclosure excluded the exact root')
    call require(nonlinear_z%lo <= 2.0_dp .and. nonlinear_z%hi >= 2.0_dp, &
        'nonlinear Krawczyk Z enclosure excluded the exact root')
    call require(nonlinear_r%lo > 0.95_dp .and. &
        nonlinear_r%hi < 1.15_dp, &
        'nonlinear Krawczyk R map is not strictly inside X')
    call require(nonlinear_z%lo > 1.85_dp .and. &
        nonlinear_z%hi < 2.05_dp, &
        'nonlinear Krawczyk Z map is not strictly inside X')

    ! Reversing N and its entire derivative row is a cut-orientation
    ! relabelling.  The endpoint and its Krawczyk enclosure must not move.
    call evaluate_neort_eqdsk_cut_endpoint_krawczyk_interval( &
        point(1.05_dp), point(1.95_dp), point(-0.0525_dp), &
        point(-0.1975_dp), point(-2.1_dp), point(-1.0_dp), &
        point(-1.0_dp), point(2.9_dp), &
        gc_outward_interval(-2.3_dp, -1.9_dp), point(-1.0_dp), &
        point(-1.0_dp), gc_outward_interval(2.7_dp, 3.1_dp), &
        gc_outward_interval(-0.1_dp, 0.1_dp), &
        gc_outward_interval(-0.1_dp, 0.1_dp), reversed_r, reversed_z)
    call require_close(reversed_r%lo, nonlinear_r%lo, &
        'cut reversal Krawczyk R lower bound')
    call require_close(reversed_r%hi, nonlinear_r%hi, &
        'cut reversal Krawczyk R upper bound')
    call require_close(reversed_z%lo, nonlinear_z%lo, &
        'cut reversal Krawczyk Z lower bound')
    call require_close(reversed_z%hi, nonlinear_z%hi, &
        'cut reversal Krawczyk Z upper bound')

    write (*, '(a)') 'test_gc_eqdsk_cut_endpoint_krawczyk OK'

contains

    pure function point(value) result(interval)
        real(dp), intent(in) :: value
        type(gc_outward_interval_t) :: interval

        interval = gc_outward_interval(value, value)
    end function point

    subroutine require_close(actual, expected, label)
        real(dp), intent(in) :: actual, expected
        character(*), intent(in) :: label

        call require(abs(actual-expected) <= tolerance*max(1.0_dp, &
            abs(expected)), trim(label)//' mismatch')
    end subroutine require_close

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message

        if (.not. condition) error stop message
    end subroutine require

end program test_gc_eqdsk_cut_endpoint_krawczyk
