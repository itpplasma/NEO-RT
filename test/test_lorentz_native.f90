program test_lorentz_native
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_quadrature, only: gauss_legendre
    use neort_config, only: config_t
    use neort_lib, only: neort_init, neort_compute_no_splines
    use neort_datatypes, only: transport_data_t
    use neort_lorentz_native, only: validate_lorentz_geometry, &
        sample_lorentz_orbit, compute_lorentz_velocity, compute_lorentz_transport
    use do_magfie_mod, only: do_magfie, s, q, dqds, psi_pr, sign_theta
    use neort_profiles, only: vth, Om_tE, ns, A1, A2
    use driftorbit, only: etatp, etadt, mph, nonlin, drift_phase, comptorque, epsmn
    use util, only: pi, mi, qi, c
    implicit none
    real(dp), parameter :: field_t = 2.0_dp, major_m = 5.0_dp, minor_m = 0.5_dp
    real(dp), parameter :: rate = 2.0e5_dp
    type(config_t) :: config
    type(transport_data_t) :: unused
    real(dp) :: mass, diffusion, omega, eta, expected(4), relative
    real(dp) :: x(40), w(40), rhs(3), torque(3), dcoef(2, 3), t_prefactor
    real(dp) :: old_q, old_dqds, saved_force(2), coarse_error, fine_error
    real(dp) :: old_rotation
    real(dp) :: junction_diffusion(3), max_orbit_error, boundary_distance, tp_error
    complex(dp) :: drive, work(3), reference(3)
    integer :: branch, k, info, refinement

    call write_chart('lorentz_circle.bc')
    config%s = 0.5_dp
    config%M_t = 0.01_dp
    config%qs = 1.0_dp
    config%ms = 2.014_dp
    config%vth = 4.0e7_dp
    config%epsmn = 1.0e-3_dp
    config%m0 = -6
    config%mph = 3
    config%comptorque = .false.
    config%magdrift = .false.
    config%noshear = .true.
    config%inp_swi = 8
    config%vsteps = 1
    config%mth_max_abs = 0
    call neort_init(config, 'lorentz_circle.bc')
    ! This initializes native geometry without evaluating any resonant roots:
    ! no magnetic drift and nonzero rotation put this harmonic off resonance.
    call neort_compute_no_splines(unused)
    call validate_lorentz_geometry(info)
    if (info /= 0) then
        print *, 'circular geometry status:', info
        error stop 'manufactured circle rejected'
    end if
    call gauss_legendre(40, x, w)
    max_orbit_error = 0.0_dp
    do branch = 1, 3
        if (branch == 1) then
            eta = (etatp + etadt)/2.0_dp
        else
            eta = 0.4_dp*etatp
        end if
        call independent_orbit(eta, expected)
        call sample_lorentz_orbit(1.0_dp, branch, eta, rate, 0.0_dp, mass, &
            diffusion, omega, drive)
        relative = max(abs(mass/expected(1) - 1.0_dp), &
            abs(real(drive, dp)/expected(2) - 1.0_dp), &
            abs(diffusion/expected(3) - 1.0_dp))
        max_orbit_error = max(max_orbit_error, relative)
        if (relative > 2.0e-7_dp) then
            print *, 'native orbit relative error:', branch, relative
            print *, mass, drive, diffusion, expected
            error stop 'native Lorentz orbit differs from independent quadrature'
        end if
        if (abs(aimag(drive)) > 1.0e-10_dp*abs(drive)) &
            error stop 'common physical phase failed to cancel'
        if (abs(omega/(mph*Om_tE) - 1.0_dp) > 1.0e-12_dp) &
            error stop 'passing streaming term did not cancel in common phase'
    end do
    do refinement = 1, 2
        boundary_distance = 10.0_dp**(-3 - refinement)
        do branch = 1, 3
            eta = etatp*(1.0_dp + &
                merge(boundary_distance, -boundary_distance, branch == 1))
            call sample_lorentz_orbit(1.0_dp, branch, eta, rate, 0.0_dp, mass, &
                junction_diffusion(branch), omega, drive)
        end do
        relative = abs(junction_diffusion(1)/sum(junction_diffusion(2:3)) - 1.0_dp)
        if (refinement == 1) tp_error = relative
    end do
    ! The approach contains delta*log(delta); a decade refinement reduces it
    ! by about 0.12 here, rather than requiring equality at finite separation.
    if (relative > 8.0e-4_dp .or. relative > 0.2_dp*tp_error) &
        error stop 'native TP limiting flux weights disagree'
    print *, 'native orbit max relative error and TP flux limit:', &
        max_orbit_error, relative
    ! Actual-field zero-diffusion limit: compare complex reaction response with
    ! a separate nested Gauss quadrature, not the graph's assembly or callback.
    reference = (0.0_dp, 0.0_dp)
    do branch = 1, 3
        do k = 1, size(x)
            if (branch == 1) then
                eta = etatp + (etadt - etatp)*(1.0_dp + x(k))/2.0_dp
                rhs(1) = (etadt - etatp)/2.0_dp
            else
                eta = etatp*(1.0_dp + x(k))/2.0_dp
                rhs(1) = etatp/2.0_dp
            end if
            call independent_orbit(eta, expected)
            reference(branch) = reference(branch) + rhs(1)*w(k)*expected(1)* &
                expected(2)**2/cmplx(rate, mph*Om_tE, dp)
        end do
    end do
    do refinement = 1, 2
        call compute_lorentz_velocity(1.0_dp, 4*refinement, 0.0_dp, rate, work, info)
        if (info /= 0) error stop 'native complex reaction limit failed'
        relative = maxval(abs(work - reference)/abs(reference))
        if (refinement == 1) coarse_error = relative
        if (refinement == 2) fine_error = relative
    end do
    print *, 'native reaction-limit relative errors:', coarse_error, fine_error
    if (fine_error > 0.015_dp .or. fine_error >= coarse_error) &
        error stop 'native complex reaction limit did not refine'

    ! A no-torque call must not consult forces initialized by an earlier run.
    old_rotation = Om_tE
    Om_tE = 0.0_dp
    A1 = 0.0_dp
    A2 = 0.0_dp
    call compute_lorentz_transport(0.9_dp*vth, 1.1_dp*vth, 1, 3, &
        deflection_rate, dcoef, torque, info)
    if (info /= 2) then
        print *, 'no-torque first singular status:', info, comptorque, Om_tE
        error stop 'no-torque singular response read unused forces'
    end if
    Om_tE = old_rotation

    ! Torque is direct dissipative susceptibility work with the verified legacy
    ! causal normalization. Zero physical thermodynamic forcing gives zero.
    ns = 2.0e13_dp
    comptorque = .true.
    A1 = -1.4_dp
    A2 = 0.3_dp
    call compute_lorentz_transport(0.9_dp*vth, 1.1_dp*vth, 1, 3, &
        deflection_rate, dcoef, torque, info)
    if (info /= 0) error stop 'native Lorentz transport failed'
    call compute_lorentz_velocity(1.0_dp, 3, rate, 0.0_dp, work, info)
    t_prefactor = sign(1.0_dp, psi_pr*q*sign_theta)*pi**1.5_dp*mph**2* &
        ns*c*vth/qi*exp(-1.0_dp)*(A1 + A2)
    rhs = 0.2_dp*t_prefactor*real(work, dp)/pi
    if (maxval(abs(torque - rhs)) > 1.0e-10_dp*maxval(abs(rhs))) &
        error stop 'native Lorentz torque normalization differs from direct work'
    if (maxval(abs(torque)) <= 0.0_dp) error stop 'actual Lorentz torque vanished'
    saved_force = [A1, A2]
    comptorque = .false.
    Om_tE = 0.0_dp
    call compute_lorentz_transport(0.9_dp*vth, 1.1_dp*vth, 1, 3, &
        deflection_rate, dcoef, torque, info)
    if (info /= 2) error stop 'previous forcing changed no-torque singular status'
    comptorque = .true.
    Om_tE = old_rotation
    A1 = 0.0_dp
    A2 = 0.0_dp
    call compute_lorentz_transport(0.9_dp*vth, 1.1_dp*vth, 1, 3, &
        deflection_rate, dcoef, torque, info)
    if (any(torque /= 0.0_dp)) error stop 'zero drive produced Lorentz torque'
    A1 = saved_force(1)
    A2 = saved_force(2)
    qi = -qi
    mi = mi/3670.0_dp
    call compute_lorentz_transport(0.9_dp*vth, 1.1_dp*vth, 1, 3, &
        deflection_rate, dcoef, torque, info)
    if (info /= 0) error stop 'negative charge rejected by real Lorentz path'
    epsmn = 0.0_dp
    call compute_lorentz_transport(0.9_dp*vth, 1.1_dp*vth, 1, 3, &
        deflection_rate, dcoef, torque, info)
    if (info /= 0 .or. any(torque /= 0.0_dp)) &
        error stop 'zero Hamiltonian drive produced Lorentz torque'
    epsmn = config%epsmn

    old_q = q
    old_dqds = dqds
    q = 2.1_dp
    call validate_lorentz_geometry(info)
    if (info /= -13) error stop 'noninteger nq accepted'
    q = old_q
    dqds = 0.01_dp
    call validate_lorentz_geometry(info)
    if (info /= -12) error stop 'physical q shear accepted'
    dqds = old_dqds
    drift_phase = .true.
    call validate_lorentz_geometry(info)
    if (info /= -10) error stop 'drift phase accepted'
    drift_phase = .false.
    nonlin = .true.
    call validate_lorentz_geometry(info)
    if (info /= -10) error stop 'nonlinear response accepted'
    print *, 'test_lorentz_native PASSED'

contains

    function deflection_rate(speed) result(frequency)
        real(dp), intent(in) :: speed
        real(dp) :: frequency
        frequency = rate
        if (speed <= 0.0_dp) error stop 'nonpositive rate sample speed'
    end function deflection_rate

    subroutine independent_orbit(eta, result_)
        real(dp), intent(in) :: eta
        real(dp), intent(out) :: result_(4)
        real(dp) :: theta, phi, kappa, factor, bmod, sqrtg, hder(3), hcovar(3)
        real(dp) :: hctrvr(3), hcurl(3), point(3), bmin, bmax, dl_dt, jacobian
        real(dp) :: hamiltonian, inverse_b, orbit_weight
        integer :: point_index
        bmin = 1.0_dp/etadt
        bmax = 1.0_dp/etatp
        result_ = 0.0_dp
        do point_index = 1, size(x)
            if (eta > etatp) then
                kappa = sqrt((1.0_dp/eta - bmin)/(bmax - bmin))
                phi = pi*x(point_index)/2.0_dp
                theta = 2.0_dp*asin(kappa*sin(phi))
                jacobian = 2.0_dp*kappa/ &
                    sqrt((1.0_dp - eta*bmin)*(1.0_dp - kappa**2*sin(phi)**2))
                factor = pi
            else
                theta = pi*(1.0_dp + x(point_index))
                point = [s, 0.0_dp, theta]
                call do_magfie(point, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
                jacobian = 1.0_dp/sqrt(1.0_dp - eta*bmod)
                factor = pi
            end if
            point = [s, 0.0_dp, theta]
            call do_magfie(point, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
            dl_dt = vth*abs(hctrvr(3))
            orbit_weight = factor*w(point_index)*jacobian/dl_dt
            hamiltonian = mi*vth**2/2.0_dp*config%epsmn*(2.0_dp - eta*bmod)
            inverse_b = 1.0_dp/bmod - eta
            result_(1) = result_(1) + orbit_weight
            result_(2) = result_(2) + orbit_weight*hamiltonian
            result_(3) = result_(3) + orbit_weight*inverse_b
        end do
        result_(2) = result_(2)/result_(1)
        result_(3) = 2.0_dp*rate*eta*result_(3)
    end subroutine independent_orbit

    subroutine write_chart(path)
        character(*), intent(in) :: path
        real(dp), parameter :: surfaces(3) = [0.25_dp, 0.5_dp, 0.75_dp]
        real(dp) :: radius, bph, bth
        integer :: unit, surface
        open(newunit=unit, file=path, status='replace')
        write(unit, '(a)') 'CC manufactured exact cosine circular field'
        write(unit, '(a)') 'CC constant q=2, one physical magnetic well'
        write(unit, '(a)') 'CC m0b n0b nsurf nper flux a R'
        write(unit, '(a)') 'CC'
        write(unit, '(a)') 'CC'
        write(unit, *) 1, 0, 3, 1, -pi*minor_m**2*field_t, minor_m, major_m
        do surface = 1, size(surfaces)
            radius = minor_m*sqrt(surfaces(surface))
            bph = -field_t*major_m*1.0e6_dp/0.2_dp
            bth = -field_t*radius**2*0.5_dp/major_m*1.0e6_dp/0.2_dp
            write(unit, '(a)') 'CC s iota Jpol Itor pprime sqrtg00'
            write(unit, '(a)') 'CC units A A Pa m3'
            write(unit, *) surfaces(surface), 0.5_dp, bph, bth, 0.0_dp, 0.0_dp
            write(unit, '(a)') 'CC m n rmnc zmnc vmnc bmnc'
            write(unit, *) 0, 0, major_m, 0.0_dp, 0.0_dp, field_t
            write(unit, *) 1, 0, radius, radius, 0.0_dp, -field_t*radius/major_m
        end do
        close(unit)
    end subroutine write_chart

end program test_lorentz_native
