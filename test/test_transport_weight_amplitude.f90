program test_transport_weight_amplitude
    ! Linear response: the transport weight is proportional to |H_mn|**2 and so
    ! to epsmn**2 exactly. The bounce integral of H_n is strongly cancelling
    ! for slow passing particles; its accuracy must therefore not depend on
    ! the perturbation amplitude through a fixed absolute ODE tolerance.
    use iso_fortran_env, only: dp => real64
    use neort_lib, only: neort_init, neort_prepare_splines, neort_setup_at_s
    use neort, only: set_to_passing_region, set_to_trapped_region
    use neort_transport, only: compute_transport_weight
    use driftorbit, only: mth, etamin, etamax, sign_vpar, epsmn, nopassing, comptorque
    implicit none
    real(dp), parameter :: big = 1.0e-3_dp, small = 1.0e-10_dp
    real(dp) :: speeds(2) = [0.0625_dp, 1.0_dp], fractions(3) = [0.1_dp, 0.5_dp, 0.9_dp]
    real(dp) :: eta, reference(3), scaled(3), worst
    integer :: region, ks, kf, harmonic

    call neort_init('driftorbit.in', 'in_file')
    nopassing = .false.
    comptorque = .true.
    call neort_prepare_splines('plasma.in', 'profile.in')
    call neort_setup_at_s(0.5_dp)
    worst = 0.0_dp
    do region = 1, 2
        do harmonic = -2, 0
            mth = harmonic
            sign_vpar = 1.0_dp
            if (region == 1) then
                call set_to_passing_region(etamin, etamax)
            else
                call set_to_trapped_region(etamin, etamax)
            end if
            do ks = 1, size(speeds)
                do kf = 1, size(fractions)
                    eta = etamin + fractions(kf)*(etamax - etamin)
                    epsmn = big
                    call compute_transport_weight(speeds(ks), eta, reference)
                    epsmn = small
                    call compute_transport_weight(speeds(ks), eta, scaled)
                    scaled = scaled*(big/small)**2
                    if (any(reference <= 0.0_dp .and. [.true., .true., .false.])) &
                        error stop 'transport weight unexpectedly vanished'
                    worst = max(worst, maxval(abs(scaled - reference)/abs(reference), &
                                              mask=reference /= 0.0_dp))
                    if (any(abs(scaled - reference) > 1.0e-8_dp*abs(reference))) then
                        print *, 'region, mth, u, eta:', region, mth, speeds(ks), eta
                        print *, 'weight at epsmn=big:', reference
                        print *, 'weight at epsmn=small, scaled:', scaled
                        error stop 'transport weight is not quadratic in the amplitude'
                    end if
                end do
            end do
        end do
    end do
    print *, 'maximum relative deviation from epsmn**2 scaling =', worst
    print *, 'test_transport_weight_amplitude PASSED'
end program test_transport_weight_amplitude
