program test_collisional_layer
    ! Unit tests for the collisional boundary-layer factor:
    !   (1) complete elliptic integrals against known values,
    !   (2) R-table interpolation: grid nodes, clamping, far-limit 1,
    !   (3) physics: an interior resonance root approaches the collisionless
    !       (delta) result as collisionality -> 0 and is suppressed as
    !       collisionality -> infinity.
    use iso_fortran_env, only: dp => real64
    use neort_collisional_layer, only: layer_ratio, layer_delta_eta, &
        comp_ellint_k, comp_ellint_e

    implicit none

    real(dp), parameter :: PI_HALF = 1.5707963267948966_dp

    call test_elliptic_integrals
    call test_table_nodes
    call test_table_clamping
    call test_table_far_limit
    call test_interior_root_collisionless_limit
    call test_interior_root_suppression

    print *, "test_collisional_layer ... OK"

contains

    subroutine check(name, got, expected, tol)
        character(*), intent(in) :: name
        real(dp), intent(in) :: got, expected, tol

        if (abs(got - expected) > tol) then
            print *, "FAIL ", name, " got=", got, " expected=", expected
            error stop
        end if
    end subroutine check

    subroutine test_elliptic_integrals
        ! K(0) = E(0) = pi/2; reference values at m = 0.5.
        call check("K(0)", comp_ellint_k(0.0_dp), PI_HALF, 1.0e-7_dp)
        call check("E(0)", comp_ellint_e(0.0_dp), PI_HALF, 1.0e-7_dp)
        call check("K(0.5)", comp_ellint_k(0.5_dp), 1.8540746773_dp, 1.0e-7_dp)
        call check("E(0.5)", comp_ellint_e(0.5_dp), 1.3506438810_dp, 1.0e-7_dp)
    end subroutine test_elliptic_integrals

    subroutine test_table_nodes
        ! Exact grid nodes must return the tabulated entries.
        call check("R(0.5,0.5)", layer_ratio(0.5_dp, 0.5_dp), &
                   0.10539480507943665_dp, 1.0e-9_dp)
        call check("R(1,2)", layer_ratio(1.0_dp, 2.0_dp), &
                   1.7412945147047214_dp, 1.0e-9_dp)
        call check("R(64,0.5)", layer_ratio(64.0_dp, 0.5_dp), &
                   0.2628278437523768_dp, 1.0e-9_dp)
    end subroutine test_table_nodes

    subroutine test_table_clamping
        real(dp) :: r_lo, r_hi, r_mid

        ! Outside the grid clamps to the nearest node.
        call check("clamp-high == node(64,64)", &
                   layer_ratio(1.0e3_dp, 1.0e3_dp), &
                   layer_ratio(64.0_dp, 64.0_dp), 1.0e-12_dp)
        call check("clamp-low == node(0.5,0.5)", &
                   layer_ratio(1.0e-3_dp, 1.0e-3_dp), &
                   layer_ratio(0.5_dp, 0.5_dp), 1.0e-12_dp)

        ! Interpolation between two nodes stays between them (monotone edge).
        r_lo = layer_ratio(0.5_dp, 0.5_dp)
        r_hi = layer_ratio(1.0_dp, 0.5_dp)
        r_mid = layer_ratio(sqrt(2.0_dp)/2.0_dp, 0.5_dp)  ! log2 = -0.5
        if (r_mid < min(r_lo, r_hi) .or. r_mid > max(r_lo, r_hi)) then
            print *, "FAIL interpolation not bracketed", r_lo, r_mid, r_hi
            error stop
        end if
    end subroutine test_table_clamping

    subroutine test_table_far_limit
        ! An interior root many layer widths from both boundaries recovers the
        ! collisionless plateau R -> 1.
        call check("far-far R -> 1", layer_ratio(64.0_dp, 64.0_dp), &
                   1.0_dp, 0.03_dp)
    end subroutine test_table_far_limit

    subroutine test_interior_root_collisionless_limit
        real(dp) :: delta_eta, ratio
        real(dp), parameter :: ETATP = 0.9_dp, ETADT = 1.1_dp, ETA = 1.0_dp
        real(dp), parameter :: EPS = 0.1_dp, ABS_DOM = 1.0e5_dp
        real(dp), parameter :: NU_TINY = 1.0e-20_dp

        ! Vanishing collisionality: layer width -> 0, root sits many widths from
        ! both boundaries, so R must approach the delta result (1) within 2%.
        delta_eta = layer_delta_eta(NU_TINY, ETA, EPS, ETATP, ETADT, ABS_DOM)
        ratio = layer_ratio((ETADT - ETA)/delta_eta, (ETA - ETATP)/delta_eta)
        call check("interior collisionless limit R -> 1", ratio, 1.0_dp, 0.02_dp)
    end subroutine test_interior_root_collisionless_limit

    subroutine test_interior_root_suppression
        real(dp) :: delta_eta, ratio
        real(dp), parameter :: ETATP = 0.9_dp, ETADT = 1.1_dp, ETA = 1.0_dp
        real(dp), parameter :: EPS = 0.1_dp, ABS_DOM = 1.0e5_dp
        real(dp), parameter :: NU_HUGE = 1.0e30_dp

        ! Huge collisionality: layer width dwarfs the trapped domain, both scaled
        ! distances clamp to the near-boundary corner, so R is strongly
        ! suppressed below the delta result.
        delta_eta = layer_delta_eta(NU_HUGE, ETA, EPS, ETATP, ETADT, ABS_DOM)
        ratio = layer_ratio((ETADT - ETA)/delta_eta, (ETA - ETATP)/delta_eta)
        if (ratio >= 0.6_dp) then
            print *, "FAIL suppression: R =", ratio, " not < 0.6"
            error stop
        end if
    end subroutine test_interior_root_suppression

end program test_collisional_layer
