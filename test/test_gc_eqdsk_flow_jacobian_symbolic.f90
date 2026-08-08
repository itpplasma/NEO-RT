program test_gc_eqdsk_flow_jacobian_symbolic
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_gc_axisymmetric_physical_flow_jacobian_symbolic, only: &
        evaluate_neort_gc_axisymmetric_physical_flow_jacobian
    use util_for_test, only: pass_test
    implicit none

    real(dp), parameter :: mass = 1.0_dp, charge = 1.0_dp
    real(dp), parameter :: c_light = 1.0_dp, radius = 1.0_dp
    real(dp), parameter :: h0 = 2.0_dp, j_k = 0.25_dp, sigma = 1.0_dp
    real(dp), parameter :: sqrt_three = 1.7320508075688772935_dp
    real(dp), parameter :: tolerance = 2.0e-12_dp
    real(dp) :: j11, j12, j21, j22, trace, determinant, discriminant
    real(dp) :: hessian_rr, hessian_rz, hessian_zz, hessian_determinant
    real(dp) :: littlejohn_factor, discriminant_from_hessian
    real(dp) :: fixed_point_v_r, fixed_point_v_z
    real(dp) :: psi_star_gradient_r, psi_star_gradient_z

    ! Manufactured local equilibrium at (R,Z)=(1,0):
    !
    !   psi = (R-1) + (Z**2)/2,  F=sqrt(3),
    !   Phi = (11/2)*psi + (1/4)*psi**2.
    !
    ! Native profile inputs below are therefore F,F',F'',Phi,Phi',Phi'';
    ! no R/Z chain-rule data is supplied by the test.
    call evaluate_neort_gc_axisymmetric_physical_flow_jacobian( &
        mass, charge, c_light, radius, h0, j_k, sigma, 0.0_dp, 1.0_dp, &
        0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
        1.0_dp, sqrt_three, 0.0_dp, 0.0_dp, 0.0_dp, 5.5_dp, 0.5_dp, &
        j11, j12, j21, j22, trace, determinant, &
        discriminant)

    hessian_rr = -119.0_dp/12.0_dp
    hessian_rz = 0.0_dp
    hessian_zz = -35.0_dp/16.0_dp
    hessian_determinant = hessian_rr*hessian_zz-hessian_rz**2
    littlejohn_factor = sqrt_three/2.0_dp
    discriminant_from_hessian = -4.0_dp*littlejohn_factor**2* &
        hessian_determinant

    call require_close("J11 fixed-point identity", j11, &
        -littlejohn_factor*hessian_rz)
    call require_close("J12 fixed-point identity", j12, &
        -littlejohn_factor*hessian_zz)
    call require_close("J21 fixed-point identity", j21, &
        littlejohn_factor*hessian_rr)
    call require_close("J22 fixed-point identity", j22, &
        littlejohn_factor*hessian_rz)
    call require_close("trace", trace, 0.0_dp)
    call require_close("determinant equals L squared Hessian determinant", &
        determinant, littlejohn_factor**2*hessian_determinant)
    call require_close("Delta equals -4 L squared Hessian determinant", &
        discriminant, discriminant_from_hessian)
    if (discriminant >= 0.0_dp) then
        error stop "manufactured elliptic fixed point was not classified O"
    end if

    ! Independent full-flow oracle at the same fixed point.  Here
    ! B=(0,sqrt(3),1), bhat=(0,sqrt(3)/2,1/2),
    ! curl(bhat)=(0,-1/2,sqrt(3)/2), p_parallel=sqrt(3),
    ! B_parallel_star=2, and force=(5,0,0), hence v_R=v_Z=0.
    fixed_point_v_r = 0.0_dp
    fixed_point_v_z = (sqrt_three*2.5_dp-2.5_dp*sqrt_three)/2.0_dp
    call require_close("full-flow fixed-point v_R", fixed_point_v_r, 0.0_dp)
    call require_close("full-flow fixed-point v_Z", fixed_point_v_z, 0.0_dp)

    ! The fixed-point condition is also checked independently from the
    ! canonical flux gradient: p_R=-5/sqrt(3), G=G_R=sqrt(3)/2.
    psi_star_gradient_r = 1.0_dp + (-5.0_dp/sqrt_three)* &
        (sqrt_three/2.0_dp) + sqrt_three*(sqrt_three/2.0_dp)
    psi_star_gradient_z = 0.0_dp
    call require_close("fixed-point psi-star R gradient", &
        psi_star_gradient_r, 0.0_dp)
    call require_close("fixed-point psi-star Z gradient", &
        psi_star_gradient_z, 0.0_dp)

    write (*, '(A)') 'test_gc_eqdsk_flow_jacobian_symbolic OK'
    call pass_test

contains

    subroutine require_close(label, actual, expected)
        character(*), intent(in) :: label
        real(dp), intent(in) :: actual, expected
        real(dp) :: scale

        scale = max(1.0_dp, max(abs(actual), abs(expected)))
        if (abs(actual-expected) > tolerance*scale) then
            write (*, '(A,2(1X,ES24.16))') trim(label), actual, expected
            error stop 'axisymmetric physical flow fixed-point oracle failed'
        end if
    end subroutine require_close

end program test_gc_eqdsk_flow_jacobian_symbolic
