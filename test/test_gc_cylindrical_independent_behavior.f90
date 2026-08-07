program test_gc_cylindrical_independent_behavior
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_cylindrical_littlejohn_symbolic, only: &
        evaluate_neort_cylindrical_littlejohn
    use neort_axisymmetric_pphi_symbolic, only: &
        evaluate_neort_axisymmetric_p_phi
    use util_for_test, only: pass_test
    implicit none

    real(dp), parameter :: mass = 1.0_dp, charge = 2.0_dp
    real(dp), parameter :: c_light = 3.0_dp, b0 = 5.0_dp
    real(dp), parameter :: radius = 2.0_dp, p_parallel = 0.7_dp
    real(dp), parameter :: mu = 0.2_dp, electric_field = 0.8_dp
    real(dp), parameter :: psi = 0.5_dp*b0*radius**2
    real(dp) :: output(22), expected_phi_dot, expected_pphi
    real(dp) :: pphi_output(3)

    ! Uniform straight B=B0 e_Z and Phi=E_R R.  The independent analytic
    ! E-cross-B result is v_phi=c E_R/B0 and precession dot(phi)=v_phi/R.
    call evaluate_neort_cylindrical_littlejohn(mass, charge, c_light, mu, b0, &
        electric_field*radius, radius, p_parallel, psi, 0.0_dp, 0.0_dp, b0, &
        0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
        0.0_dp, 0.0_dp, 0.0_dp, electric_field, 0.0_dp, 0.0_dp, &
        output(1), output(2), output(3), output(4), output(5), output(6), &
        output(7), output(8), output(9), output(10), output(11), output(12), &
        output(13), output(14), output(15), output(16), output(17), output(18), &
        output(19), output(20), output(21), output(22))

    expected_phi_dot = c_light*electric_field/(b0*radius)
    expected_pphi = charge/c_light*psi
    call require_close('analytic E cross B phi dot', output(16), expected_phi_dot)
    call require_close('analytic radial drift', output(14), 0.0_dp)
    call require_close('analytic parallel drift', output(15), 0.0_dp)
    call require_close('analytic parallel force', output(17), 0.0_dp)
    call require_close('analytic H conservation', output(20), 0.0_dp)
    call require_close('axisymmetric canonical p_phi', output(22), expected_pphi)
    call evaluate_neort_axisymmetric_p_phi(mass, charge, c_light, mu, radius, &
        p_parallel, psi, b0, 0.0_dp, 0.0_dp, 1.0_dp, &
        0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, electric_field/(b0*radius), &
        pphi_output(1), pphi_output(2), pphi_output(3))
    call require_close('generated axisymmetric P_phi', pphi_output(1), expected_pphi)
    call require_close('generated P_phi time derivative', pphi_output(2), 0.0_dp)
    call require_close('generated P_phi residual', pphi_output(3), 0.0_dp)

    write (*, '(A)') 'test_gc_cylindrical_independent_behavior OK'
    call pass_test

contains

    subroutine require_close(label, actual, expected)
        character(*), intent(in) :: label
        real(dp), intent(in) :: actual, expected
        real(dp) :: scale

        scale = max(1.0_dp, max(abs(actual), abs(expected)))
        if (abs(actual - expected) > 2.0e-12_dp*scale) then
            write (*, '(A,2(1X,ES24.16))') trim(label), actual, expected
            error stop 'independent cylindrical behavior oracle failed'
        end if
    end subroutine require_close

end program test_gc_cylindrical_independent_behavior
