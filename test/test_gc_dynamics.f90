program test_gc_dynamics
    !! Behavioral oracles for the coordinate-general guiding-center kernel.
    !! The expected uniform-field drifts are closed-form expressions, while
    !! conservation of H and J_perp follows from independent invariant
    !! definitions rather than from repository output.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_gc_dynamics, only: GC_SUCCESS, gc_field_sample_t, gc_rhs
    use util_for_test, only: pass_test

    implicit none

    type(gc_field_sample_t) :: field
    real(dp) :: grad_phi(3), xdot(3), pdot, xidot
    real(dp) :: p, xi, rho0, lambda, expected, hdot, jdot
    integer :: status

    call check_uniform_e_cross_b()
    call check_uniform_grad_b()
    call check_curvature_b_star()
    call check_invariants()
    call pass_test

contains

    subroutine check_uniform_e_cross_b()
        field = gc_field_sample_t()
        field%bmod = 2.0_dp
        field%sqrtg = 1.0_dp
        field%hcov = [0.0_dp, 0.0_dp, 1.0_dp]
        field%hcon = field%hcov
        grad_phi = [0.6_dp, 0.0_dp, 0.0_dp]
        p = 1.2_dp
        xi = 0.4_dp
        rho0 = 0.2_dp
        lambda = 0.5_dp

        call gc_rhs(field, grad_phi, rho0, lambda, p, xi, &
            xdot, pdot, xidot, status)
        if (status /= GC_SUCCESS) error stop 'uniform E x B RHS failed'

        expected = lambda*rho0*grad_phi(1)/(2.0_dp*field%bmod)
        call require_close('uniform E x B xdot(1)', xdot(1), 0.0_dp, 1.0e-14_dp)
        call require_close('uniform E x B xdot(2)', xdot(2), expected, 1.0e-14_dp)
        call require_close('uniform E x B parallel motion', xdot(3), p*xi, 1.0e-14_dp)
        call require_close('uniform E x B pdot', pdot, 0.0_dp, 1.0e-14_dp)
        call require_close('uniform E x B xidot', xidot, 0.0_dp, 1.0e-14_dp)
    end subroutine check_uniform_e_cross_b

    subroutine check_uniform_grad_b()
        field = gc_field_sample_t()
        field%bmod = 2.0_dp
        field%sqrtg = 1.0_dp
        field%grad_log_b = [0.3_dp, 0.0_dp, 0.0_dp]
        field%hcov = [0.0_dp, 0.0_dp, 1.0_dp]
        field%hcon = field%hcov
        grad_phi = 0.0_dp
        p = 1.2_dp
        xi = 0.4_dp
        rho0 = 0.2_dp
        lambda = 0.5_dp

        call gc_rhs(field, grad_phi, rho0, lambda, p, xi, &
            xdot, pdot, xidot, status)
        if (status /= GC_SUCCESS) error stop 'uniform grad-B RHS failed'

        expected = 0.5_dp*p**2*(1.0_dp - xi**2)/field%bmod &
            *lambda*rho0*field%grad_log_b(1)
        call require_close('uniform grad-B xdot(1)', xdot(1), 0.0_dp, 1.0e-14_dp)
        call require_close('uniform grad-B xdot(2)', xdot(2), expected, 1.0e-14_dp)
        call require_close('uniform grad-B parallel motion', xdot(3), p*xi, 1.0e-14_dp)
        call require_close('uniform grad-B pdot', pdot, 0.0_dp, 1.0e-14_dp)
    end subroutine check_uniform_grad_b

    subroutine check_curvature_b_star()
        !! A nonzero curl(h) must enter both h* and B_parallel*.
        field = gc_field_sample_t()
        field%bmod = 2.0_dp
        field%sqrtg = 1.0_dp
        field%hcov = [0.0_dp, 0.0_dp, 1.0_dp]
        field%hcon = field%hcov
        field%curl_h = [0.2_dp, -0.1_dp, 0.3_dp]
        grad_phi = 0.0_dp
        p = 1.2_dp
        xi = -0.4_dp
        rho0 = 0.2_dp
        lambda = 0.5_dp

        call gc_rhs(field, grad_phi, rho0, lambda, p, xi, &
            xdot, pdot, xidot, status)
        if (status /= GC_SUCCESS) error stop 'curvature B-star RHS failed'

        block
            real(dp) :: a_c(3), hstar(3), h_parallel_star
            a_c = field%curl_h*lambda*rho0/field%bmod
            hstar = field%hcon + p*xi*a_c
            h_parallel_star = 1.0_dp + p*xi*dot_product(a_c, field%hcov)
            call require_close('curvature B-star xdot(1)', xdot(1), &
                p*xi*hstar(1)/h_parallel_star, 1.0e-14_dp)
            call require_close('curvature B-star xdot(2)', xdot(2), &
                p*xi*hstar(2)/h_parallel_star, 1.0e-14_dp)
            call require_close('curvature B-star xdot(3)', xdot(3), &
                p*xi*hstar(3)/h_parallel_star, 1.0e-14_dp)
        end block
        call require_close('curvature B-star pdot', pdot, 0.0_dp, 1.0e-14_dp)
        call require_close('curvature B-star xidot', xidot, 0.0_dp, 1.0e-14_dp)
    end subroutine check_curvature_b_star

    subroutine check_invariants()
        !! Static-field guiding-center dynamics must conserve
        !! H=p**2+phihat and J_perp=p**2(1-xi**2)/B pointwise.
        field = gc_field_sample_t()
        field%bmod = 3.2_dp
        field%sqrtg = 1.7_dp
        field%grad_log_b = [0.13_dp, -0.07_dp, 0.19_dp]
        field%hcov = [0.2_dp, -0.3_dp, 0.7_dp]
        field%hcon = field%hcov/sum(field%hcov**2)
        field%curl_h = [-0.11_dp, 0.17_dp, 0.05_dp]
        grad_phi = [-0.23_dp, 0.31_dp, 0.09_dp]
        p = 1.37_dp
        xi = -0.42_dp
        rho0 = 0.08_dp
        lambda = -0.6_dp

        call gc_rhs(field, grad_phi, rho0, lambda, p, xi, &
            xdot, pdot, xidot, status)
        if (status /= GC_SUCCESS) error stop 'generic invariant RHS failed'

        hdot = 2.0_dp*p*pdot + dot_product(xdot, grad_phi)
        jdot = (2.0_dp*p*pdot*(1.0_dp - xi**2) &
            - 2.0_dp*p**2*xi*xidot &
            - p**2*(1.0_dp - xi**2)*dot_product(xdot, field%grad_log_b)) &
            /field%bmod
        call require_close('guiding-center energy conservation', hdot, 0.0_dp, 2.0e-13_dp)
        call require_close('guiding-center magnetic moment conservation', jdot, 0.0_dp, 2.0e-13_dp)
    end subroutine check_invariants

    subroutine require_close(label, actual, reference, tolerance)
        character(*), intent(in) :: label
        real(dp), intent(in) :: actual, reference, tolerance
        if (abs(actual - reference) > tolerance) then
            write(*, '(a,2(1x,es24.16),1x,a,1x,es12.4)') &
                trim(label), actual, reference, 'tolerance', tolerance
            error stop 'guiding-center behavioral oracle failed'
        end if
    end subroutine require_close

end program test_gc_dynamics
