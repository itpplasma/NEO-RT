module neort_gc_dynamics
    !! Coordinate-general, nonrelativistic guiding-center equations.
    !!
    !! Spatial vectors follow the repository's curvilinear convention:
    !! gradients and h are covariant, while velocities and curl(h) are
    !! contravariant.  The orbit-width scale multiplies every first-order rho
    !! occurrence in B*, B_parallel*, grad-B, and E x B terms.
    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none
    private

    integer, parameter, public :: GC_SUCCESS = 0
    integer, parameter, public :: GC_INVALID_INPUT = 1
    integer, parameter, public :: GC_SINGULAR_BSTAR = 2

    type, public :: gc_field_sample_t
        real(dp) :: bmod = 0.0_dp
        real(dp) :: sqrtg = 0.0_dp
        real(dp) :: grad_log_b(3) = 0.0_dp
        real(dp) :: hcov(3) = 0.0_dp
        real(dp) :: hcon(3) = 0.0_dp
        real(dp) :: curl_h(3) = 0.0_dp
        real(dp) :: psi = 0.0_dp
        real(dp) :: grad_psi(3) = 0.0_dp
    end type gc_field_sample_t

    public :: gc_rhs

contains

    pure subroutine gc_rhs(field, grad_phi, rho0, orbit_width_scale, p, xi, &
            xdot, pdot, xidot, status)
        type(gc_field_sample_t), intent(in) :: field
        real(dp), intent(in) :: grad_phi(3)
        real(dp), intent(in) :: rho0, orbit_width_scale, p, xi
        real(dp), intent(out) :: xdot(3), pdot, xidot
        integer, intent(out) :: status

        real(dp) :: rho, p_parallel, one_minus_xi2, mu_hat
        real(dp) :: a_phi(3), a_b(3), a_c(3), hstar(3), h_parallel_star

        xdot = 0.0_dp
        pdot = 0.0_dp
        xidot = 0.0_dp
        status = GC_INVALID_INPUT
        if (field%bmod <= 0.0_dp) return
        if (abs(field%sqrtg) <= tiny(field%sqrtg)) return
        if (p <= 0.0_dp) return
        if (rho0 == 0.0_dp) return
        if (abs(xi) > 1.0_dp + 100.0_dp*epsilon(xi)) return

        rho = orbit_width_scale*rho0
        p_parallel = p*xi
        one_minus_xi2 = max(0.0_dp, 1.0_dp - xi**2)
        mu_hat = 0.5_dp*p**2*one_minus_xi2/field%bmod

        a_phi = covariant_cross(field%hcov, grad_phi) &
            *rho/(2.0_dp*field%sqrtg*field%bmod)
        a_b = covariant_cross(field%hcov, field%grad_log_b) &
            *rho/field%sqrtg
        a_c = field%curl_h*rho/field%bmod
        hstar = field%hcon + p_parallel*a_c
        h_parallel_star = 1.0_dp + p_parallel*dot_product(a_c, field%hcov)
        if (abs(h_parallel_star) <= 100.0_dp*epsilon(h_parallel_star)) then
            status = GC_SINGULAR_BSTAR
            return
        end if

        xdot = (p_parallel*hstar + a_phi + mu_hat*a_b)/h_parallel_star
        pdot = -0.5_dp*dot_product(xdot, grad_phi)/p
        xidot = -0.5_dp*one_minus_xi2/h_parallel_star*( &
            dot_product(hstar, grad_phi)/p &
            + p*dot_product(hstar, field%grad_log_b) &
            + xi*dot_product(a_phi, field%grad_log_b))
        status = GC_SUCCESS
    end subroutine gc_rhs

    pure function covariant_cross(left, right) result(value)
        real(dp), intent(in) :: left(3), right(3)
        real(dp) :: value(3)

        value(1) = left(2)*right(3) - left(3)*right(2)
        value(2) = left(3)*right(1) - left(1)*right(3)
        value(3) = left(1)*right(2) - left(2)*right(1)
    end function covariant_cross

end module neort_gc_dynamics
