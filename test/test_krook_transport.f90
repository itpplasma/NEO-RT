program test_krook_transport
    use iso_fortran_env, only: dp => real64
    use fortnum_quadrature, only: gauss_legendre
    use neort_lib, only: neort_init, neort_prepare_splines, neort_setup_at_s
    use neort_transport, only: compute_transport_integral, compute_transport_weight
    use neort_freq, only: Om_ph, Om_th
    use neort_krook, only: krook_collisions, krook_nu, krook_model, krook_rtol
    use neort_profiles, only: vth, Om_tE
    use do_magfie_mod, only: q, iota, R0, a, s
    use driftorbit, only: mth, mph, etamin, etamax, B0, mi, &
                          comptorque, nopassing
    use util, only: pi, c, qi
    implicit none
    real(dp) :: baseline(3), disabled(3), zero_width(3), broadened(3), expected(3)
    real(dp) :: x(96), w(96), eta, omph, omth, dvdummy, detadummy
    real(dp) :: weight(3), residual, scale, D_plateau, dsdreff
    integer :: k

    call neort_init('driftorbit.in', 'in_file')
    nopassing = .true.
    call neort_prepare_splines('plasma.in', 'profile.in')
    call neort_setup_at_s(0.5_dp)
    mth = 0
    comptorque = .true.
    Om_tE = 1.0e6_dp
    krook_collisions = .false.
    krook_nu = 0.0_dp
    call compute_transport_integral(0.9_dp*vth, 1.1_dp*vth, 1, baseline(1:2), baseline(3))
    krook_nu = 1.0e7_dp
    call compute_transport_integral(0.9_dp*vth, 1.1_dp*vth, 1, disabled(1:2), disabled(3))
    krook_collisions = .true.
    krook_nu = 0.0_dp
    call compute_transport_integral(0.9_dp*vth, 1.1_dp*vth, 1, &
        zero_width(1:2), zero_width(3))
    if (any(disabled /= baseline)) &
        error stop 'disabled Krook path changed collisionless transport'
    if (any(zero_width /= baseline)) error stop 'zero width changed collisionless transport'
    krook_nu = 1.0e7_dp
    krook_model = 0
    krook_rtol = 1.0e-5_dp
    call compute_transport_integral(0.9_dp*vth, 1.1_dp*vth, 1, broadened(1:2), broadened(3))
    call gauss_legendre(96, x, w)
    expected = 0.0_dp
    do k = 1, size(x)
        eta = (etamin + etamax)/2.0_dp + (etamax - etamin)*x(k)/2.0_dp
        call Om_th(vth, eta, omth, dvdummy, detadummy)
        call Om_ph(vth, eta, omph, dvdummy, detadummy)
        residual = mth*omth + mph*omph
        call compute_transport_weight(1.0_dp, eta, weight)
        expected = expected + w(k)*weight*krook_nu/(pi*(residual**2 + krook_nu**2))
    end do
    expected = 0.2_dp*(etamax - etamin)*expected/2.0_dp
    D_plateau = pi*vth**3/(16.0_dp*R0*iota*(qi*B0/(mi*c))**2)
    dsdreff = 2.0_dp/a*sqrt(s)
    expected(1:2) = dsdreff**(-2)*expected(1:2)/D_plateau
    scale = maxval(abs(expected))
    if (scale <= 0.0_dp) error stop 'Krook transport oracle unexpectedly vanished'
    if (any(abs(broadened - expected) > 5.0e-4_dp*max(abs(expected), tiny(1.0_dp)))) then
        print *, 'broad-width transport:', broadened
        print *, 'independent eta quadrature:', expected
        error stop 'Krook transport normalization or dispatch mismatch'
    end if
    if (broadened(1) <= 0.0_dp) error stop 'finite width lost nonresonant transport'
    print *, 'Krook transport maximum relative error =', &
        maxval(abs(broadened-expected)/max(abs(expected), tiny(1.0_dp)))
    print *, 'test_krook_transport PASSED'
end program test_krook_transport
