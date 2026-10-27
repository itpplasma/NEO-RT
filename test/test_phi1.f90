program test_phi1
    use iso_fortran_env, only: dp => real64
    use neort_lib, only: neort_init, neort_prepare_splines, neort_setup_at_s
    use neort, only: set_to_passing_region, set_to_trapped_region
    use neort_phi1, only: phi1, phi_mn, m0_phi
    use neort_transport, only: compute_transport_weight, transport_rates, Omth, &
                               perturbation_amplitude_bound
    use driftorbit, only: mth, etamin, etamax, sign_vpar, epsmn, m0, &
                   nopassing, comptorque, etatp, Bmin, drift_phase, fow_radial, pertfile
    use util, only: qi, mi
    implicit none
    real(dp) :: wb(3), we(3), twice(3), plus(3), minus(3), eta, u, cross(3)
    real(dp) :: y(7), rates(7), magnetic(7), b, bound
    complex(dp) :: expected
    integer :: region, charge, geometry

    call neort_init('driftorbit.in', 'in_file')
    nopassing = .false.
    comptorque = .true.
    call neort_prepare_splines('plasma.in', 'profile.in')
    call neort_setup_at_s(0.5_dp)
    phi1 = .true.
    m0_phi = m0
    mth = 0
    sign_vpar = 1.0_dp
    u = 0.8_dp
    do region = 1, 2
        if (region == 1) then
            call set_to_passing_region(etamin, etamax)
        else
            call set_to_trapped_region(etamin, etamax)
        end if
        eta = (etamin + etamax)/2
        epsmn = 1.0e-3_dp
        phi_mn = (0.0_dp, 0.0_dp)
        call compute_transport_weight(u, eta, wb)
        phi1 = .false.
        call compute_transport_weight(u, eta, twice)
        if (any(wb /= twice)) error stop 'zero Phi differs from switch off'
        phi1 = .true.
        epsmn = 0.0_dp
        phi_mn = cmplx(1.0e-3_dp*mi*(u*4.0e7_dp)**2/(2*qi), 0.0_dp, dp)
        call compute_transport_weight(u, eta, we)
        phi_mn = 2*phi_mn
        call compute_transport_weight(u, eta, twice)
        call check_weights(twice, 4*we, 'Phi-only quadratic scaling')
        phi_mn = phi_mn/2
        epsmn = 1.0e-3_dp
        call compute_transport_weight(u, eta, plus)
        phi_mn = -phi_mn
        call compute_transport_weight(u, eta, minus)
        call check_weights(plus + minus, 2*(wb + we), 'pi phase cross-term flip')
        cross = (plus - minus)/2
        if (abs(cross(1)) < 1.0e-3_dp*sqrt(wb(1)*we(1))) &
            error stop 'cross-term fixture is degenerate'
        print *, 'region, magnetic, electric, cross D11:', &
            region, wb(1), we(1), cross(1)
    end do

    ! Independent local oracle: no velocity/pitch weight on q*Phi, both signs,
    ! and exactly the magnetic orbit phase at t=0, trapped and passing t>0.
    drift_phase = .false.
    fow_radial = .false.
    pertfile = .false.
    b = Bmin
    y = 0.0_dp
    y(1) = 0.37_dp
    Omth = 0.8_dp
    do charge = -1, 1, 2
        qi = charge*abs(qi)
        do geometry = 0, 2
            eta = 0.5_dp*etatp
            if (geometry == 1) eta = 1.1_dp*etatp
            phi1 = .false.
            call transport_rates(2.0_dp, eta, 7, real(geometry, dp), y, b, &
                                 [0.0_dp, 0.0_dp, 0.0_dp], [0.0_dp, 1.0_dp, 0.0_dp], &
                                 [0.0_dp, 0.0_dp, 1.0_dp], magnetic)
            phi1 = .true.
            phi_mn = (0.03_dp, -0.02_dp)
            m0_phi = m0
            call transport_rates(2.0_dp, eta, 7, real(geometry, dp), y, b, &
                                 [0.0_dp, 0.0_dp, 0.0_dp], [0.0_dp, 1.0_dp, 0.0_dp], &
                                 [0.0_dp, 0.0_dp, 1.0_dp], rates)
            expected = qi*phi_mn/(mi*2.0_dp)* &
                       cmplx(magnetic(3), magnetic(4), dp)/((2 - eta*b)*epsmn)
            if (abs(cmplx(rates(3) - magnetic(3), rates(4) - magnetic(4), dp) &
             - expected) > 2.0e-13_dp*abs(expected)) error stop 'q Phi/E phase or units'
            call perturbation_amplitude_bound([0.5_dp, 0.0_dp, y(1)], &
                                              2.0_dp, eta, bound)
            if (bound < abs(cmplx(rates(3), rates(4), dp))) then
                print *, 'charge, geometry, bound, harmonic:', charge, geometry, &
                    bound, abs(cmplx(rates(3), rates(4), dp))
                error stop 'Hamiltonian amplitude bound omits Phi'
            end if
        end do
    end do
    print *, 'test_phi1 PASSED'
contains
    subroutine check_weights(actual, expected_weights, label)
        real(dp), intent(in) :: actual(3), expected_weights(3)
        character(*), intent(in) :: label
        if (any(abs(actual - expected_weights) > &
                2.0e-7_dp*max(abs(expected_weights), tiny(1.0_dp)))) then
            print *, label, actual, expected_weights
            error stop 'Phi transport oracle failed'
        end if
    end subroutine check_weights
end program test_phi1
