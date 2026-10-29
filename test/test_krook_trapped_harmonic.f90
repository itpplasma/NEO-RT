program test_krook_trapped_harmonic
    ! Circular large-aspect-ratio field, D at s=0.09, constant Krook rate,
    ! trapped mth=-40, third thermal speed node. Beyond eta/etatp-1 ~ 1e-3
    ! the bounce harmonic |H|/S falls to the orbit-integration floor ~3e-11,
    ! where the squared weight is noise. Without an orbit-error budget the
    ! pitch quadrature exhausted its pieces and NEO-RT stopped.
    ! Reference: fixed Gauss-Legendre panels in x = log(eta - etatp), which
    ! resolve the logarithmic layer at the trapped-passing boundary without
    ! the adaptive residual map, with an explicit Lorentzian.
    use iso_fortran_env, only: dp => real64
    use ieee_arithmetic, only: ieee_is_finite
    use fortnum_quadrature, only: gauss_legendre
    use neort_lib, only: neort_init, neort_prepare_splines, neort_setup_at_s
    use neort, only: prepare_transport, select_component, component_trapped
    use neort_transport, only: compute_krook_velocity, compute_transport_weight
    use neort_resonance, only: resonance_value
    use neort_mixture_collisions, only: effective_krook_rate
    use driftorbit, only: etamin, etamax, etatp, vth, mth
    use do_magfie_mod, only: s
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: ux = 0.0117197470703125_dp
    ! The orbit-error budget of this node is 1.8e-3 of its value.
    real(dp), parameter :: tolerance = 2.0e-3_dp
    real(dp) :: vmin, vmax, nu, values(3), coarse(3), fine(3), x(16), w(16)
    integer :: lower, upper
    call neort_init('driftorbit.in', 'in_file', 'in_file_pert')
    call neort_prepare_splines('plasma.in', 'profile.in')
    call neort_setup_at_s(s)
    call prepare_transport(lower, upper)
    call select_component(-40, component_trapped, vmin, vmax)
    call compute_krook_velocity(ux, values, -1.0_dp)
    if (.not. all(ieee_is_finite(values))) error stop 'nonfinite trapped sample'
    nu = effective_krook_rate(ux, real(mth, dp), 1.0_dp, ux*vth*0.01_dp)
    call gauss_legendre(size(x), x, w)
    call reference(64, coarse)
    call reference(128, fine)
    print '(a,3es24.15)', 'Krook quadrature:', values
    print '(a,3es24.15)', 'log-pitch panels:', fine
    print '(a,3es12.4)', 'relative deviation:', values/fine - 1.0_dp
    if (any(abs(coarse/fine - 1.0_dp) > 1.0e-4_dp)) &
        error stop 'log-pitch reference has not converged'
    if (any(abs(values/fine - 1.0_dp) > tolerance)) &
        error stop 'trapped high harmonic outside its orbit-error budget'
contains
    subroutine reference(panels, total)
        integer, intent(in) :: panels
        real(dp), intent(out) :: total(3)
        real(dp) :: xa, h, xm, eta, weight(3), res, dres
        integer :: k, j
        xa = log(etamin - etatp)
        h = (log(etamax - etatp) - xa)/real(panels, dp)
        total = 0.0_dp
        do k = 1, panels
            do j = 1, size(x)
                xm = xa + h*(real(k, dp) - 0.5_dp + 0.5_dp*x(j))
                eta = etatp + exp(xm)
                call compute_transport_weight(ux, eta, weight, -1.0_dp)
                call resonance_value(ux*vth, eta, res, dres)
                total = total + weight*nu/(pi*(res**2 + nu**2))*exp(xm)*0.5_dp*h*w(j)
            end do
        end do
    end subroutine reference
end program test_krook_trapped_harmonic
