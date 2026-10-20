program test_krook_reference
    ! Retain the sensitive T co-passing m=4, velocity node 3 of the MARS L5/8
    ! reproducer. The fixture keeps 64 original equilibrium surfaces around
    ! s=.25 and a three-surface linear snapshot of all 641 perturbation modes.
    ! Reference: invariant-preserving passing samples, native Debug
    ! (GNU 16.2.1), rtol=1e-8.
    ! Qualified by independent angular period/drift integrals:
    ! sinh-graded trapezoidal quadrature, 16384/32768 nodes with Richardson
    ! extrapolation, fitted on x=1e-6*exp(0:.1:2). Independent 65536-node
    ! evaluations at x=1e-8..1e-6 agree within 1.6e-8. Transport at rtol=1e-8
    ! agrees with rtol=1e-7 within 2e-8. Main's old two-point extrapolation
    ! overestimates the near-separatrix torque by 5.2e-6. The independent
    ! field-line calculation, rather than main, qualifies this reference.
    use iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use neort_lib, only: neort_init, neort_setup_at_s
    use neort_species, only: selected_species
    use neort, only: set_to_passing_region
    use driftorbit, only: mth, sign_vpar, etamin, etamax
    use neort_transport, only: compute_krook_velocity
    use neort_krook, only: krook_rtol
    use neort_freq, only: Om_th, Om_tB
    use neort_orbit, only: Om_tB_local_v2
    use neort_profiles, only: vth
    use do_magfie_mod, only: do_magfie, s
    use driftorbit, only: etatp
    use fortnum_quadrature, only: gauss_legendre
    use util, only: pi
    implicit none
    real(dp), parameter :: expected(3) = [ &
        -2.5567993930219084e-9_dp, -8.9887544578422849e-9_dp, &
        3.5936183085097517e4_dp]
    real(dp) :: values(3), ux

    call neort_init('r0p50000.in', 'in_file', 'in_file_pert')
    selected_species = 2
    sign_vpar = 1
    call neort_setup_at_s(0.25_dp)
    mth = 4
    call set_to_passing_region(etamin, etamax)
    call check_angular_frequencies()
    ux = 1.0e-6_dp + 2.5_dp*(6.0_dp - 1.0e-6_dp)/8
    ! A 1e-5 integration request cannot qualify a 1e-6 reference comparison.
    krook_rtol = 1.0e-7_dp
    call compute_krook_velocity(ux, values, -1.0_dp)
    if (.not. all(ieee_is_finite(values))) error stop 'nonfinite Krook reference'
    print '(a,3es26.17e3)', 'Krook reference values:', values
    print *, 'Krook reference relative deviation:', abs(values/expected - 1.0_dp)
    if (any(abs(values - expected) > 1.0e-6_dp*abs(expected))) &
        error stop 'Krook transport exceeds the qualified reference bound'
contains
    subroutine check_angular_frequencies()
        real(dp) :: eta, coarse(2), fine(2), omega, drift, dv, de, error(2)
        integer :: k

        do k = 6, 8
            eta = etatp*(1.0_dp - 10.0_dp**(-k))
            call angular_integrals(eta, 2048, coarse)
            call angular_integrals(eta, 4096, fine)
            if (any(abs(coarse/fine - 1.0_dp) > 2.0e-8_dp)) &
                error stop 'MARS angular frequency oracle has not converged'
            call Om_th(vth, eta, omega, dv, de)
            call Om_tB(vth, eta, drift, dv, de)
            error = abs([omega*fine(1)/(2*pi*vth), &
                drift/(vth**2*fine(2))] - 1.0_dp)
            print *, 'MARS pitch gap and angular frequency errors:', k, error
            if (any(error > 1.0e-6_dp)) &
                error stop 'MARS separatrix frequencies exceed the angular oracle bound'
        end do
    end subroutine check_angular_frequencies

    subroutine angular_integrals(eta, panels, result_)
        ! Independent field-line quadrature; no time integrator or fitted model.
        real(dp), intent(in) :: eta
        integer, intent(in) :: panels
        real(dp), intent(out) :: result_(2)
        real(dp) :: nodes(32), weights(32), x(3), bmod, sqrtg
        real(dp) :: hder(3), hcovar(3), hctrvr(3), hcurl(3), dt, width
        integer :: panel, node

        call gauss_legendre(32, nodes, weights)
        result_ = 0.0_dp
        width = 2*pi/real(panels, dp)
        x(1:2) = [s, 0.0_dp]
        do panel = 1, panels
            do node = 1, size(nodes)
                x(3) = width*(real(panel, dp) - 0.5_dp + nodes(node)/2)
                call do_magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
                dt = width/2*weights(node)/ &
                    (sqrt(1.0_dp - eta*bmod)*abs(hctrvr(3)))
                result_ = result_ + dt*[1.0_dp, &
                    Om_tB_local_v2(eta, bmod, hder(1), hctrvr(3))]
            end do
        end do
        result_(2) = result_(2)/result_(1)
    end subroutine angular_integrals
end program test_krook_reference
