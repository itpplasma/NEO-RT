module neort_thin_orbit_limit
    !! Numerical definition of the strict first-order toroidal precession.
    !!
    !! For return maps at fixed (H, mu, P_phi),
    !!
    !!   F(lambda) = (Delta phi(lambda) - 2*pi*W*q(s_*))/tau(lambda)
    !!   Omega_t^(1) = dF/dlambda at lambda=0.
    !!
    !! Three centered rungs provide an observed order and two Richardson
    !! estimates.  The module owns no field or integrator state.
    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none
    private

    integer, parameter, public :: THIN_LIMIT_SUCCESS = 0
    integer, parameter, public :: THIN_LIMIT_INVALID_INPUT = 1
    integer, parameter, public :: THIN_LIMIT_RETURN_ERROR = 2
    integer, parameter, public :: THIN_LIMIT_BASELINE_ERROR = 3
    integer, parameter, public :: THIN_LIMIT_TOPOLOGY_CHANGE = 4
    integer, parameter, public :: THIN_LIMIT_CONVERGENCE_ERROR = 5

    type, public :: orbit_return_t
        real(dp) :: period = 0.0_dp
        real(dp) :: delta_phi = 0.0_dp
        integer :: orbit_class = 0
        integer :: winding = 0
        integer :: status = 0
    end type orbit_return_t

    type, public :: thin_limit_result_t
        real(dp) :: omega = 0.0_dp
        real(dp) :: error_estimate = huge(1.0_dp)
        real(dp) :: observed_order = 0.0_dp
        real(dp) :: baseline_residual = huge(1.0_dp)
        real(dp) :: centered(3) = 0.0_dp
        real(dp) :: richardson_coarse = 0.0_dp
        real(dp) :: richardson_fine = 0.0_dp
        real(dp) :: lambda_used(3) = 0.0_dp
        integer :: refinement_count = 0
        integer :: status = THIN_LIMIT_INVALID_INPUT
    end type thin_limit_result_t

    public :: estimate_thin_limit

contains

    pure subroutine estimate_thin_limit(topological_delta_phi, base, lambda, &
            plus, minus, result, baseline_tolerance)
        real(dp), intent(in) :: topological_delta_phi
        type(orbit_return_t), intent(in) :: base
        real(dp), intent(in) :: lambda(3)
        type(orbit_return_t), intent(in) :: plus(3), minus(3)
        type(thin_limit_result_t), intent(out) :: result
        real(dp), intent(in), optional :: baseline_tolerance

        real(dp) :: tolerance, fplus, fminus, difference_coarse, difference_fine
        real(dp) :: ratio
        integer :: k

        result = thin_limit_result_t()
        result%lambda_used = lambda
        if (any(lambda <= 0.0_dp)) return
        if (abs(lambda(1)/lambda(2) - 2.0_dp) > 1.0e-12_dp) return
        if (abs(lambda(2)/lambda(3) - 2.0_dp) > 1.0e-12_dp) return
        if (base%status /= 0 .or. base%period <= 0.0_dp) then
            result%status = THIN_LIMIT_RETURN_ERROR
            return
        end if

        tolerance = 1.0e-10_dp*max(1.0_dp, abs(topological_delta_phi))
        if (present(baseline_tolerance)) tolerance = baseline_tolerance
        result%baseline_residual = base%delta_phi - topological_delta_phi
        if (abs(result%baseline_residual) > tolerance) then
            result%status = THIN_LIMIT_BASELINE_ERROR
            return
        end if

        do k = 1, 3
            if (plus(k)%status /= 0 .or. minus(k)%status /= 0 &
                .or. plus(k)%period <= 0.0_dp &
                .or. minus(k)%period <= 0.0_dp) then
                result%status = THIN_LIMIT_RETURN_ERROR
                return
            end if
            if (plus(k)%orbit_class /= base%orbit_class &
                .or. minus(k)%orbit_class /= base%orbit_class &
                .or. plus(k)%winding /= base%winding &
                .or. minus(k)%winding /= base%winding) then
                result%status = THIN_LIMIT_TOPOLOGY_CHANGE
                return
            end if
            fplus = (plus(k)%delta_phi - topological_delta_phi)/plus(k)%period
            fminus = (minus(k)%delta_phi - topological_delta_phi)/minus(k)%period
            result%centered(k) = (fplus - fminus)/(2.0_dp*lambda(k))
        end do

        result%richardson_coarse = &
            (4.0_dp*result%centered(2) - result%centered(1))/3.0_dp
        result%richardson_fine = &
            (4.0_dp*result%centered(3) - result%centered(2))/3.0_dp
        result%omega = result%richardson_fine
        result%error_estimate = abs(result%richardson_fine &
            - result%richardson_coarse)/15.0_dp

        difference_coarse = abs(result%centered(1) - result%centered(2))
        difference_fine = abs(result%centered(2) - result%centered(3))
        if (difference_fine > 100.0_dp*epsilon(1.0_dp) &
            *max(1.0_dp, maxval(abs(result%centered)))) then
            ratio = difference_coarse/difference_fine
            if (ratio > 0.0_dp) result%observed_order = log(ratio)/log(2.0_dp)
        else
            result%observed_order = huge(1.0_dp)
        end if
        result%status = THIN_LIMIT_SUCCESS
    end subroutine estimate_thin_limit

end module neort_thin_orbit_limit
