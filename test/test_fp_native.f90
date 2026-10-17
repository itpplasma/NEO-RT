program test_fp_native
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_quadrature, only: gauss_legendre
    use neort_lib, only: config_t, mixture_t, neort_init, &
        neort_prepare_mixture, neort_setup_at_s
    use neort_species, only: selected_species
    use neort_fp_native, only: compute_fp_transport, fp_physical_rates
    use neort_fp_mode, only: tp_fp_collisions
    use neort_lorentz_native, only: sample_lorentz_orbit, &
        compute_lorentz_transport
    use neort_profiles, only: vth, ns, A1, A2, Om_tE
    use driftorbit, only: etatp, etadt, mph
    use do_magfie_mod, only: psi_pr, q, sign_theta, do_magfie, s, &
        dBthcovds, dBphcovds
    use util, only: pi, mi, qi, c
    implicit none
    real(dp), parameter :: elementary_charge = 1.602176634e-19_dp
    real(dp), parameter :: epsilon0 = 8.8541878188e-12_dp
    real(dp), parameter :: mass_kg(3) = &
        [3.345228e-27_dp, 5.017842e-27_dp, 9.1091e-31_dp]
    real(dp), parameter :: charge(3) = [1.0_dp, 1.0_dp, -1.0_dp]
    real(dp), parameter :: density(3) = [5.0e18_dp, 5.0e18_dp, 1.0e19_dp]
    real(dp), parameter :: temperature = 2000.0_dp, rate_scale = 0.7_dp
    real(dp), parameter :: field_t = 2.0_dp, major_m = 5.0_dp, minor_m = 0.5_dp
    character(len=1), parameter :: name(3) = ['D', 'T', 'e']
    type(config_t) :: config
    type(mixture_t) :: mixture
    real(dp) :: torque(3), lorentz_torque(3), diffusion(2, 3), expected_torque(3)
    real(dp) :: expected_speed, expected_force(2), nud, duu, expected_rate(2)
    real(dp) :: eta, speed, tau(2), pitch, omega, drift(2), inverse_b(2)
    real(dp) :: drift_reference, drift_coarse, drift_relative
    real(dp) :: coarse_error, fine_error, scale, snapshot(5), direct_prefactor
    real(dp) :: tail_torque(3, 4), tail_coarse, tail_fine, tail_final, umax
    complex(dp) :: work(2, 3), hamiltonian(2)
    integer :: species, branch, point, refinement, info
    logical :: physical_response

    call write_chart('fp_circle.bc')
    mixture%reference_s = 0.5_dp
    mixture%omega_e = -700.0_dp
    mixture%coulomb_log = 20.0_dp
    allocate(mixture%species(3))
    do species = 1, 3
        mixture%species(species)%name = name(species)
        mixture%species(species)%mass_kg = mass_kg(species)
        mixture%species(species)%charge_e = charge(species)
        mixture%species(species)%density_m3 = density(species)
        mixture%species(species)%temperature_ev = temperature
        mixture%species(species)%density_gradient_m3 = -0.4_dp*density(species)
        mixture%species(species)%temperature_gradient_ev = 0.2_dp*temperature
    end do
    config%s = mixture%reference_s
    config%epsmn = 1.0e-3_dp
    config%m0 = -6
    config%mph = 3
    config%comptorque = .true.
    config%magdrift = .true.
    config%magdrift_passing = 1
    config%noshear = .true.
    config%inp_swi = 8
    config%vsteps = 1
    config%mth_max_abs = 0
    config%tp_fp_collisions = .true.
    config%fp_rate_scale = rate_scale
    call neort_init(config, 'fp_circle.bc')
    call neort_prepare_mixture(mixture)
    selected_species = 1
    call neort_setup_at_s(config%s)

    ! Cache scaling must hold for actual orbits, including a direct magnetic
    ! drift observation. Subtracting two large total frequencies is avoided.
    scale = 0.0_dp
    do branch = 1, 3
        eta = 0.4_dp*etatp
        if (branch == 1) eta = (etatp + etadt)/2.0_dp
        do point = 1, 2
            speed = 0.65_dp*point
            call sample_lorentz_orbit(speed, branch, eta, 1.0_dp, 0.0_dp, &
                tau(point), pitch, omega, hamiltonian(point), &
                drift_frequency=drift(point), inverse_b_average=inverse_b(point))
            call independent_drift(speed, eta, 48, drift_coarse)
            call independent_drift(speed, eta, 96, drift_reference)
            call close(drift_coarse, drift_reference, 1.0e-9_dp, &
                'independent theta quadrature convergence')
            drift_relative = abs(drift(point) - drift_reference)/abs(drift_reference)
            print *, 'native drift/reference error:', branch, speed, drift_relative
            ! Fixed native absolute ODE tolerances act on the small tau*Omega
            ! accumulator. Bound each sample against angle quadrature, rather
            ! than treating their difference as exact u-scaling information.
            call close(drift(point), drift_reference, 1.0e-6_dp, &
                'native magnetic drift versus theta quadrature')
            tau(point) = speed*tau(point)
            hamiltonian(point) = hamiltonian(point)/speed**2
            drift(point) = drift(point)/speed**2
        end do
        call close(tau(1), tau(2), 3.0e-7_dp, 'u times native bounce time')
        if (abs(hamiltonian(2)) == 0.0_dp) error stop 'native H vanished'
        if (abs(hamiltonian(1) - hamiltonian(2)) > &
                3.0e-7_dp*abs(hamiltonian(2))) &
            error stop 'native H/u squared is not speed independent'
        call close(inverse_b(1), inverse_b(2), 3.0e-7_dp, 'native inverse B')
        scale = max(scale, abs(drift(1)))
    end do
    if (scale == 0.0_dp) error stop 'drift scaling tested a disabled drift'
    call compute_fp_transport(4, 3, 2.0_dp, fp_physical_rates, &
        torque, work, physical_response, info)
    if (info /= 0 .or. physical_response) &
        error stop 'magnetic-drift native full-FP solve failed'
    if (.not. all(ieee_is_finite(torque)) .or. &
            .not. all(ieee_is_finite(real(work, dp))) .or. &
            .not. all(ieee_is_finite(aimag(work)))) &
        error stop 'nonfinite magnetic-drift native FP response'

    ! Keep electric detuning nonzero for ordinary two-column responses. The
    ! physical mixture initializes exact SI quantities for D, T and electrons.
    config%magdrift = .false.
    call neort_init(config, 'fp_circle.bc')
    call neort_prepare_mixture(mixture)
    do species = 1, 3
        selected_species = species
        call neort_setup_at_s(config%s)
        expected_speed = 100.0_dp*sqrt(2*elementary_charge*temperature/ &
            mass_kg(species))
        call close(vth, expected_speed, 3.0e-15_dp, 'exact SI thermal speed')
        call close(mi, mass_kg(species)*1.0e3_dp, 3.0e-15_dp, 'SI mass')
        call close(qi/c, charge(species)*elementary_charge/10.0_dp, &
            3.0e-15_dp, 'exact SI charge over c')
        call close(ns, density(species)*1.0e-6_dp, 3.0e-15_dp, 'SI density')
        call close(Om_tE, mixture%omega_e, 3.0e-15_dp, 'common rotation')
        expected_force(2) = 0.2_dp
        expected_force(1) = -0.7_dp - charge(species)*elementary_charge/10.0_dp/ &
            (temperature*elementary_charge*1.0e7_dp)* &
            sign_theta*psi_pr/q*mixture%omega_e
        call close(A1, expected_force(1), 4.0e-15_dp, 'exact SI A1')
        call close(A2, expected_force(2), 3.0e-15_dp, 'temperature force')
        call independent_rates(0.8_dp*expected_speed/100.0_dp, expected_rate)
        call fp_physical_rates(0.8_dp*vth, nud, duu)
        call close(nud, expected_rate(1), 2.0e-9_dp, 'all-background nuD')
        call close(duu, expected_rate(2), 2.0e-9_dp, 'all-background Dvv/vth2')
        call compute_fp_transport(4, 3, 2.0_dp, fp_physical_rates, &
            torque, work, physical_response, info)
        if (info /= 0 .or. physical_response) &
            error stop 'native full-FP two-column solve failed'
        if (.not. all(ieee_is_finite(torque))) &
            error stop 'nonfinite native full-FP torque'
        if (.not. all(ieee_is_finite(real(work, dp))) .or. &
                .not. all(ieee_is_finite(aimag(work)))) &
            error stop 'nonfinite native full-FP work'
        if (maxval(abs(torque)) == 0.0_dp) error stop 'full-FP torque vanished'
        ! Direct electromagnetic torque uses q_C/10, not legacy Gaussian qe.
        ! Eliminating beta and the imaginary source fixes sign and pi factors.
        direct_prefactor = sign(1.0_dp, psi_pr*q*sign_theta)*sqrt(pi)* &
            real(mph, dp)**2*density(species)*1.0e-6_dp*expected_speed/ &
            (charge(species)*elementary_charge/10.0_dp)
        expected_torque = direct_prefactor*real( &
            expected_force(1)*work(1, :) + expected_force(2)*work(2, :), dp)
        scale = maxval(abs(expected_torque))
        if (maxval(abs(torque - expected_torque)) > 2.0e-13_dp*scale) &
            error stop 'native FP direct torque normalization'
    end do

    ! Disable only energy diffusion, without changing the active FP-normalized
    ! physical species state. Compare against an independent speed quadrature
    ! of the native Lorentz solver on exactly the same finite interval.
    selected_species = 1
    call neort_setup_at_s(config%s)
    snapshot = [vth, qi/c, mi, A1, A2]
    ! Finite speed boundaries reflect. Hold h=0.25 while moving that wall, so
    ! this tests truncation independently of speed-grid refinement.
    do point = 1, 4
        umax = real(point + 1, dp)
        call compute_fp_transport(4*(point + 1), 3, umax, fp_physical_rates, &
            tail_torque(:, point), work, physical_response, info)
        if (info /= 0 .or. physical_response) &
            error stop 'finite-reflecting full-FP tail solve failed'
    end do
    if (.not. all(ieee_is_finite(tail_torque))) &
        error stop 'nonfinite full-FP tail torque'
    scale = maxval(abs(tail_torque(:, 4)))
    if (scale == 0.0_dp) error stop 'full-FP tail torque vanished'
    tail_coarse = maxval(abs(tail_torque(:, 2) - tail_torque(:, 1)))/scale
    tail_fine = maxval(abs(tail_torque(:, 3) - tail_torque(:, 2)))/scale
    tail_final = maxval(abs(tail_torque(:, 4) - tail_torque(:, 3)))/scale
    print *, 'full-FP reflecting-wall relative changes:', &
        tail_coarse, tail_fine, tail_final
    if (tail_fine >= tail_coarse .or. tail_final >= tail_fine .or. &
            tail_final > 0.01_dp) error stop 'FP speed truncation did not converge'
    call compute_lorentz_transport(0.0_dp, 2.0_dp*vth, 32, 3, &
        deflection_rate, diffusion, lorentz_torque, info)
    if (info /= 0) error stop 'same-state Lorentz reference failed'
    scale = maxval(abs(lorentz_torque))
    if (scale == 0.0_dp) error stop 'same-state Lorentz torque vanished'
    do refinement = 1, 2
        call compute_fp_transport(8*refinement, 3, 2.0_dp, fp_physical_rates, &
            torque, work, physical_response, info, energy_multiplier=0.0_dp)
        if (info /= 0 .or. physical_response) error stop 'energy-off FP failed'
        fine_error = maxval(abs(torque - lorentz_torque))/scale
        if (refinement == 1) coarse_error = fine_error
    end do
    print *, 'energy-off FP/Lorentz torque errors:', coarse_error, fine_error
    if (fine_error > 0.01_dp .or. fine_error >= 0.9_dp*coarse_error) &
        error stop 'energy-off FP torque did not converge to Lorentz'
    if (any(snapshot /= [vth, qi/c, mi, A1, A2])) &
        error stop 'cross-mode comparison changed physical species state'

    ! At zero electric and magnetic detuning the separate density sources need
    ! not be compatible. Zero PHYSICAL forcing must still give the zero state.
    mixture%omega_e = 0.0_dp
    mixture%species%density_gradient_m3 = 0.0_dp
    mixture%species%temperature_gradient_ev = 0.0_dp
    call neort_prepare_mixture(mixture)
    selected_species = 1
    call neort_setup_at_s(config%s)
    if (A1 /= 0.0_dp .or. A2 /= 0.0_dp .or. Om_tE /= 0.0_dp) &
        error stop 'zero physical drive fixture is not zero'
    call compute_fp_transport(3, 3, 2.0_dp, fp_physical_rates, &
        torque, work, physical_response, info)
    if (info /= 0 .or. .not. physical_response) &
        error stop 'zero physical drive at zero detuning rejected'
    if (any(torque /= 0.0_dp) .or. any(abs(work) /= 0.0_dp)) &
        error stop 'zero physical drive produced FP response'
    tp_fp_collisions = .false.
    call compute_fp_transport(3, 3, 2.0_dp, fp_physical_rates, &
        torque, work, physical_response, info)
    if (info /= -20) error stop 'inactive native FP state was accepted'
    print *, 'test_fp_native PASSED'

contains

    subroutine independent_drift(ux, eta_, order, frequency)
        real(dp), intent(in) :: ux, eta_
        integer, intent(in) :: order
        real(dp), intent(out) :: frequency
        real(dp) :: x(order), weights(order), theta, phi, kappa, jacobian
        real(dp) :: bmin, bmax, bmod, sqrtg, hder(3), hcovar(3), hctrvr(3)
        real(dp) :: hcurl(3), position(3), orbit_weight, total_weight, local_drift
        integer :: node
        call gauss_legendre(order, x, weights)
        bmin = 1.0_dp/etadt
        bmax = 1.0_dp/etatp
        frequency = 0.0_dp
        total_weight = 0.0_dp
        do node = 1, order
            if (eta_ > etatp) then
                kappa = sqrt((1.0_dp/eta_ - bmin)/(bmax - bmin))
                phi = pi*x(node)/2.0_dp
                theta = 2.0_dp*asin(kappa*sin(phi))
                jacobian = 2.0_dp*kappa/sqrt((1.0_dp - eta_*bmin)* &
                    (1.0_dp - kappa**2*sin(phi)**2))
            else
                theta = pi*(1.0_dp + x(node))
            end if
            position = [s, 0.0_dp, theta]
            call do_magfie(position, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
            if (eta_ < etatp) jacobian = 1.0_dp/sqrt(1.0_dp - eta_*bmod)
            orbit_weight = pi*weights(node)*jacobian/(ux*vth*abs(hctrvr(3)))
            ! Explicit local guiding-centre drift, with zero physical shear.
            ! Read the native field, but do not call its derived frequency.
            local_drift = mi*q*(ux*vth)**2/ &
                (2.0_dp*(qi/c)*sign_theta*psi_pr*bmod)* &
                (-(2.0_dp - eta_*bmod)*bmod*hder(1) + &
                2.0_dp*(1.0_dp - eta_*bmod)*hctrvr(3)* &
                (dBthcovds + q*dBphcovds))
            frequency = frequency + orbit_weight*real(mph, dp)*local_drift
            total_weight = total_weight + orbit_weight
        end do
        frequency = frequency/total_weight
    end subroutine independent_drift

    subroutine close(actual, expected, tolerance, label)
        real(dp), intent(in) :: actual, expected, tolerance
        character(*), intent(in) :: label
        if (.not. ieee_is_finite(actual)) error stop 'nonfinite native observation'
        if (abs(actual - expected) > tolerance*max(abs(actual), abs(expected))) then
            print *, label, actual, expected
            error stop 'native FP oracle mismatch'
        end if
    end subroutine close

    function deflection_rate(speed_cm_s) result(rate)
        real(dp), intent(in) :: speed_cm_s
        real(dp) :: rate, ignored
        call fp_physical_rates(speed_cm_s, rate, ignored)
    end function deflection_rate

    subroutine independent_rates(speed_m_s, result_)
        real(dp), intent(in) :: speed_m_s
        real(dp), intent(out) :: result_(2)
        real(dp) :: xb, chandrasekhar, gamma, thermal_speed
        integer :: background
        result_ = 0.0_dp
        do background = 1, 3
            thermal_speed = sqrt(2*elementary_charge*temperature/mass_kg(background))
            xb = speed_m_s/thermal_speed
            ! For x >= 8 the omitted relative correction to G is below
            ! 1.5e-27, far below binary64 roundoff; avoid a trapped exp tail.
            if (xb >= 8.0_dp) then
                chandrasekhar = 1.0_dp/(2*xb**2)
            else
                chandrasekhar = (erf(xb) - &
                    2*xb*exp(-xb**2)/sqrt(pi))/(2*xb**2)
            end if
            gamma = density(background)*mixture%coulomb_log* &
                (charge(selected_species)*charge(background)*elementary_charge**2)**2/ &
                (4*pi*epsilon0**2*mass_kg(selected_species)**2)
            result_(1) = result_(1) + gamma* &
                (erf(xb) - chandrasekhar)/speed_m_s**3
            result_(2) = result_(2) + gamma*chandrasekhar/speed_m_s/ &
                (2*elementary_charge*temperature/mass_kg(selected_species))
        end do
        result_ = rate_scale*result_
    end subroutine independent_rates

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
end program test_fp_native
