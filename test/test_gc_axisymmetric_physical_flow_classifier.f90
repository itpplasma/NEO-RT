program test_gc_axisymmetric_physical_flow_classifier
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_gc_axisymmetric_physical_flow_classifier, only: &
        GC_PHYSICAL_FLOW_MISSING_STATIONARY_CERTIFICATE, &
        GC_PHYSICAL_FLOW_O, GC_PHYSICAL_FLOW_OUT_OF_DOMAIN, &
        GC_PHYSICAL_FLOW_SUCCESS, GC_PHYSICAL_FLOW_TURNING, &
        GC_PHYSICAL_FLOW_X, gc_axisymmetric_physical_flow_input_t, &
        gc_axisymmetric_physical_flow_result_t, &
        gc_axisymmetric_stationary_point_certificate_t, &
        classify_gc_axisymmetric_physical_fixed_point
    use util_for_test, only: pass_test
    implicit none

    type(gc_axisymmetric_physical_flow_input_t) :: input
    type(gc_axisymmetric_stationary_point_certificate_t) :: certificate
    type(gc_axisymmetric_physical_flow_result_t) :: result
    integer :: status
    real(dp), parameter :: sqrt_three = 1.7320508075688772935_dp
    real(dp), parameter :: tolerance = 2.0e-12_dp
    real(dp) :: hessian_rr, hessian_zz, littlejohn_factor
    real(dp) :: expected_discriminant

    certificate%certified = .true.
    certificate%certificate_id = 9101
    certificate%residual_r = 0.0_dp
    certificate%residual_z = 0.0_dp
    certificate%tolerance = 1.0e-12_dp

    ! Independent O-point oracle.  At this manufactured stationary point
    ! the physical psi-star Hessian is diagonal with the values below;
    ! Delta=-4*L**2*det(Hess psi-star) is therefore strictly negative.
    input = manufactured_input(0.0_dp)
    call classify_gc_axisymmetric_physical_fixed_point(input, certificate, &
        result, status)
    call require(status == GC_PHYSICAL_FLOW_SUCCESS, &
        'manufactured O point was rejected')
    call require(result%classification == GC_PHYSICAL_FLOW_O, &
        'negative physical-flow discriminant was not classified O')
    hessian_rr = -119.0_dp/12.0_dp
    hessian_zz = -35.0_dp/16.0_dp
    littlejohn_factor = sqrt_three/2.0_dp
    expected_discriminant = -4.0_dp*littlejohn_factor**2* &
        hessian_rr*hessian_zz
    call require_close(result%discriminant, expected_discriminant, &
        'independent O-point discriminant')
    call require(result%bparallel_star > 0.0_dp .and. &
        result%p_parallel > 0.0_dp, 'O-point physical scales are invalid')

    ! Independent X-point oracle: psi_Z=psi_RZ=0 makes the stationary
    ! residual independent of psi_ZZ at this point.  Setting psi_ZZ=-3
    ! reverses the physical psi-star Hessian determinant while preserving the
    ! same first-order stationary certificate.  Direct rational reduction of
    ! -4*L**2*det(Hess psi-star) gives 6.2475.
    input = manufactured_input(0.0_dp)
    input%psi_zz = -3.0_dp
    call classify_gc_axisymmetric_physical_fixed_point(input, certificate, &
        result, status)
    call require(status == GC_PHYSICAL_FLOW_SUCCESS, &
        'manufactured X point was rejected')
    call require(result%classification == GC_PHYSICAL_FLOW_X, &
        'positive physical-flow discriminant was not classified X')
    call require(result%discriminant > 0.0_dp, &
        'X-point discriminant sign was not preserved')
    call require_close(result%discriminant, 6.2475_dp, &
        'independent X-point discriminant')

    ! Sigma labels the signed parallel branch.  Only the branch sign and
    ! accepted-domain facts are asserted here; no sigma-reversal symmetry is
    ! assumed for the physical flow classification.
    input%sigma = -1.0_dp
    call classify_gc_axisymmetric_physical_fixed_point(input, certificate, &
        result, status)
    call require(status == GC_PHYSICAL_FLOW_SUCCESS .and. &
        result%p_parallel < 0.0_dp, 'negative sigma branch was not retained')

    ! A fixed point cannot be classified without the caller's certificate.
    certificate%certified = .false.
    input = manufactured_input(0.0_dp)
    call classify_gc_axisymmetric_physical_fixed_point(input, certificate, &
        result, status)
    call require(status == GC_PHYSICAL_FLOW_MISSING_STATIONARY_CERTIFICATE .and. &
        result%classification == 0, 'uncertified point was not unresolved')

    certificate%certified = .true.
    input = manufactured_input(0.0_dp)
    input%h0 = 0.5_dp
    call classify_gc_axisymmetric_physical_fixed_point(input, certificate, &
        result, status)
    call require(status == GC_PHYSICAL_FLOW_TURNING .and. &
        result%classification == 0, 'turning point was not rejected')

    input = manufactured_input(0.0_dp)
    input%radius = ieee_value(0.0_dp, ieee_quiet_nan)
    call classify_gc_axisymmetric_physical_fixed_point(input, certificate, &
        result, status)
    call require(status == 1 .and. result%classification == 0, &
        'nonfinite input was not rejected')

    input = manufactured_input(0.0_dp)
    input%sigma = 0.0_dp
    call classify_gc_axisymmetric_physical_fixed_point(input, certificate, &
        result, status)
    call require(status == GC_PHYSICAL_FLOW_OUT_OF_DOMAIN, &
        'invalid sigma was accepted')

    call pass_test

contains

    function manufactured_input(psi_rr) result(value)
        real(dp), intent(in) :: psi_rr
        type(gc_axisymmetric_physical_flow_input_t) :: value

        value%mass = 1.0_dp
        value%charge = 1.0_dp
        value%c_light = 1.0_dp
        value%radius = 1.0_dp
        value%h0 = 2.0_dp
        value%j_k = 0.25_dp
        value%sigma = 1.0_dp
        value%psi = 0.0_dp
        value%psi_r = 1.0_dp
        value%psi_z = 0.0_dp
        value%psi_rr = psi_rr
        value%psi_rz = 0.0_dp
        value%psi_zz = 1.0_dp
        value%psi_rrr = 0.0_dp
        value%psi_rrz = 0.0_dp
        value%psi_rzz = 0.0_dp
        value%psi_zzz = 0.0_dp
        value%psi_sep = 1.0_dp
        value%f = sqrt_three
        value%df_dpsihat = 0.0_dp
        value%d2f_dpsihat2 = 0.0_dp
        value%phi = 0.0_dp
        value%dphi_dpsi = 5.5_dp
        value%d2phi_dpsi2 = 0.5_dp
    end function manufactured_input

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) error stop message
    end subroutine require

    subroutine require_close(actual, expected, message)
        real(dp), intent(in) :: actual, expected
        character(len=*), intent(in) :: message
        real(dp) :: scale

        scale = max(1.0_dp, max(abs(actual), abs(expected)))
        call require(abs(actual-expected) <= tolerance*scale, message)
    end subroutine require_close

end program test_gc_axisymmetric_physical_flow_classifier
