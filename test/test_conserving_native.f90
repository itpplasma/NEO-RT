program test_conserving_native
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, &
                                             ieee_is_finite
    use neort_lib, only: config_t, mixture_t, neort_init, &
                         neort_prepare_mixture, neort_setup_at_s
    use neort_species, only: selected_species
    use neort_conserving_native, only: conserving_species_t, &
                           prepare_conserving_species, assemble_conserving_collisions, &
                                       conserving_native_point
    use neort_landau_pair_quadrature, only: landau_quadrature_t, &
                                            initialize_landau_quadrature
    use neort_lorentz_native, only: sample_lorentz_orbit
    use neort_transport, only: transport_omth => Omth, timestep_transport
    use neort_orbit, only: nvar
    use driftorbit, only: native_mth => mth, native_sign => sign_vpar, etatp
    use fortnum_quadrature, only: gauss_legendre
    use do_magfie_mod, only: psi_pr, q
    use neort_profiles, only: vth, Om_tE
    use util, only: pi
    implicit none
    character(len=1), parameter :: names(3) = ['D', 'T', 'e']
    real(dp), parameter :: elementary_charge = 1.602176634e-19_dp
    real(dp), parameter :: mass_kg(3) = &
                           [3.345228e-27_dp, 5.017842e-27_dp, 9.1091e-31_dp]
    real(dp), parameter :: charge(3) = [1.0_dp, 1.0_dp, -1.0_dp]
    real(dp), parameter :: density(3) = [5.0e18_dp, 5.0e18_dp, 1.0e19_dp]
    real(dp), parameter :: temperature = 2000.0_dp
    real(dp), parameter :: field_t = 2.0_dp, major_m = 5.0_dp, minor_m = 0.5_dp
    type(config_t) :: config
    type(mixture_t) :: mixture
    type(conserving_species_t) :: cache(3)
    type(landau_quadrature_t) :: quadrature
    complex(dp) :: collision(51, 51), vector(51), work
    real(dp) :: omega
    real(dp) :: scale, residual, density_mass, expected_mass, theta, xi, u
    real(dp) :: x(128), weights(128), orb_mass, orb_precession, local_mass
    real(dp) :: mass,pitch,frequency,eta,lo,hi,gamma,precession,error(3),coarse_error(3)
    real(dp) :: kpair, energy_j, vrel, friction, heat, volume, flow_coeff(3), observed
    complex(dp) :: drive, orb_h, local_h
    integer :: species, info, n, i, offset, point, arm, order, refinement
    real(dp), parameter :: sample_pitch(3) = [0.05_dp, 0.9_dp, -0.9_dp]

    call write_chart('conserving_circle.bc')
    mixture%reference_s = 0.5_dp
    mixture%omega_e = -700.0_dp
    mixture%coulomb_log = 20.0_dp
    allocate (mixture%species(3))
    do species = 1, 3
        mixture%species(species)%name = names(species)
        mixture%species(species)%mass_kg = mass_kg(species)
        mixture%species(species)%charge_e = charge(species)
        mixture%species(species)%density_m3 = density(species)
        mixture%species(species)%temperature_ev = temperature
        mixture%species(species)%density_gradient_m3 = -0.4_dp*density(species)
        mixture%species(species)%temperature_gradient_ev = 0.2_dp*temperature
    end do
    config%s = mixture%reference_s
    config%epsmn = 1.0e-3_dp
    config%m0 = -5
    config%mph = 3
    config%comptorque = .true.
    config%magdrift = .true.
    config%magdrift_passing = 1
    config%noshear = .true.
    config%inp_swi = 8
    config%vsteps = 1
    config%mth_max_abs = 0
    config%tp_fp_collisions = .true.
    call neort_init(config, 'conserving_circle.bc')
    call neort_prepare_mixture(mixture)
    ! A local Hamiltonian at t=0 has no bounce-time phase. It must not depend
    ! on a frequency installed by an earlier orbit integration.
    transport_omth = ieee_value(0.0_dp, ieee_quiet_nan)
    do species = 1, 3
        selected_species = species
        call neort_setup_at_s(config%s)
        call prepare_conserving_species(2, 2, 3.0_dp, 48, 8, cache(species), info)
        if (info /= 0) then
            print *, 'native preparation status:', info
            if (info == -32) print *, 'finite M/S/P/ReH/ImH/ReH2/ImH2:', &
                all(ieee_is_finite(cache(species)%mass)), &
                all(ieee_is_finite(cache(species)%streaming)), &
                all(ieee_is_finite(cache(species)%precession)), &
                all(ieee_is_finite(real(cache(species)%hamiltonian, dp))), &
                all(ieee_is_finite(aimag(cache(species)%hamiltonian))), &
                all(ieee_is_finite(real(cache(species)%hamiltonian2, dp))), &
                all(ieee_is_finite(aimag(cache(species)%hamiltonian2)))
            error stop 'native preparation rejected'
        end if
        n = cache(species)%basis%ngraph
        density_mass = sum(cache(species)%mass(:n, :n))
        expected_mass = density(species)*sum(cache(species)%theta_weight* &
                                             cache(species)%jacobian_si)
       call close (density_mass, expected_mass, 2.0e-9_dp, 'full physical density mass')
        scale = maxval(abs(cache(species)%mass))*vth/major_m/100.0_dp
        residual = maxval(abs(cache(species)%streaming + &
                              transpose(cache(species)%streaming)))/scale
        if (residual > 2.0e-12_dp) error stop 'native streaming is not skew'
        if (maxval(abs(cache(species)%precession - &
                       transpose(cache(species)%precession))) > &
            2.0e-15_dp*maxval(abs(cache(species)%precession))) &
            error stop 'native precession symmetry'
        if (maxval(abs(cache(species)%precession(:n, 16))) == 0) &
            error stop 'passing/flow precession cross blocks omitted'
        do point = 1, 8
            theta = 0.13_dp + 0.7_dp*point
            u = 0.4_dp + 0.13_dp*point
            do arm = 1, 3
                xi = sample_pitch(arm)
                call check_local_phase(theta, u, xi)
            end do
        end do
       print *, 'physical mass and streaming residual:', species, density_mass, residual
    end do
    ! Independent orbit-time quadrature: constant graph restriction of M,P,H.
    ! The theta/xi assembly above has neither bounce time nor orbit-averaged H.
    selected_species = 1
    call neort_setup_at_s(config%s)
    do refinement = 1, 2
        order = 64*refinement
        call gauss_legendre(order, x(:order), weights(:order))
        orb_mass = 0
        orb_precession = 0
        orb_h = (0.0_dp, 0.0_dp)
        gamma = density(1)*1.0e-6_dp*vth*abs(psi_pr/q)/sqrt(pi)
        do arm = 1, 3
            if (arm == 1) then
                lo = cache(1)%basis%etatp
                hi = cache(1)%basis%etadt
            else
                lo = 0
                hi = cache(1)%basis%etatp
            end if
            do point = 1, order
                eta = ((1 - x(point))*lo + (1 + x(point))*hi)/2
                call sample_lorentz_orbit(1.0_dp, arm, eta, 1.0_dp, 0.0_dp, &
                              mass, pitch, frequency, drive, drift_frequency=precession)
                local_mass = gamma*(hi - lo)/2*weights(point)*mass
                orb_mass = orb_mass + local_mass*sqrt(pi)/4
                orb_precession = orb_precession + local_mass*sqrt(pi)* &
                                 (real(config%mph, dp)*Om_tE/4 + 3*precession/8)
                orb_h = orb_h + local_mass*3*sqrt(pi)/8*drive
            end do
        end do
        n = cache(1)%basis%ngraph
        local_mass = sum(cache(1)%mass(:n, :n))
        local_h = sum(cache(1)%hamiltonian(:n, :))
        error(1) = abs(local_mass - orb_mass)/abs(orb_mass)
   error(2) = abs(sum(cache(1)%precession(:n, :n)) - orb_precession)/abs(orb_precession)
        error(3) = abs(local_h - orb_h)/abs(orb_h)
        print *, 'native graph M/P/H versus actual bounce:', order, error
        if (refinement == 1) coarse_error = error
    end do
    if (any(error > 0.4_dp*coarse_error)) &
        error stop 'independent orbit quadrature did not converge'
    if (maxval(error) > 3.0e-4_dp) error stop 'native graph restriction mismatch'
    ! Full physical-theta coupled Landau blocks retain all admitted invariants.
    call initialize_landau_quadrature(2, 2, 3, 6, quadrature)
    call assemble_conserving_collisions(cache, 48, quadrature, collision, info)
    if (info /= 0) error stop 'native collision assembly rejected'
    scale = maxval(abs(collision))
    if (scale <= 0) error stop 'empty native collision matrix'
    if (maxval(abs(collision - transpose(conjg(collision)))) > 3.0e-15_dp*scale) &
        error stop 'native Landau Hermiticity'
    residual = 0
    do species = 1, 3
        vector = (0.0_dp, 0.0_dp)
        offset = (species - 1)*17
        vector(offset + 1:offset + 15) = 1
        residual = max(residual, maxval(abs(matmul(collision, vector)))/scale)
    end do
    vector = (0.0_dp, 0.0_dp)
    do species = 1, 3
        vector((species - 1)*17 + 16) = mass_kg(species)*cache(species)%vth_m_s/ &
                                        (mass_kg(1)*cache(1)%vth_m_s)
    end do
    residual = max(residual, maxval(abs(matmul(collision, vector)))/scale)
    vector = (0.0_dp, 0.0_dp)
    vector([17, 34, 51]) = 1
    residual = max(residual, maxval(abs(matmul(collision, vector)))/scale)
    if (residual > 1.0e-13_dp) error stop 'projected conservation residual'
    do i = 1, 51
        vector(i) = cmplx(sin(real(i, dp)), cos(0.7_dp*i), dp)
    end do
    work = dot_product(vector, matmul(collision, vector))
    if (real(work, dp) <= 0) error stop 'native entropy positivity'
    if (maxval(abs(collision(:17, 18:34))) == 0) &
        error stop 'cross species response blocks omitted'
    ! Independent common-T Maxwellian relative-flow and heat exchange. Their
    ! analytic Gaussian relative-velocity integrals contain no graph basis.
    energy_j = elementary_charge*temperature
    volume = sum(cache(1)%theta_weight*cache(1)%jacobian_si)
    flow_coeff = mass_kg*[cache(1)%vth_m_s, cache(2)%vth_m_s, &
                          cache(3)%vth_m_s]/energy_j
    do species = 2, 3
        kpair = 2*pi*(charge(1)*charge(species)*elementary_charge**2/ &
                      (4*pi*8.8541878188e-12_dp))**2*mixture%coulomb_log
        vrel = sqrt(cache(1)%vth_m_s**2 + cache(species)%vth_m_s**2)
        friction = 4*kpair*density(1)*density(species)/(3*sqrt(pi)*energy_j*vrel)
        heat = 3*friction/(mass_kg(1) + mass_kg(species))
        observed = -real(collision(16, (species - 1)*17 + 16), dp)* &
                   flow_coeff(1)*flow_coeff(species)
        call close (observed, volume*friction/energy_j, 3.0e-13_dp, &
                    'absolute native pair friction')
        observed = -real(collision(17, species*17), dp)
        call close (observed, volume*heat, 3.0e-13_dp, 'absolute native pair heat')
    end do
    print *, 'native collision null residual and entropy:', residual, real(work, dp)
    print *, 'test_conserving_native PASSED'
contains
    subroutine check_local_phase(theta_, u_, xi_)
        real(dp), intent(in) :: theta_, u_, xi_
        real(dp) :: state(nvar), derivative(nvar), field, eta_, time, harmonic
        real(dp) :: saved_sign, saved_frequency
        complex(dp) :: local_h, expected, finite_time_h
        integer :: saved_mth, direction
        integer(int64) :: poisoned_frequency

        saved_mth = native_mth
        saved_sign = native_sign
        saved_frequency = transport_omth
        transport_omth = ieee_value(0.0_dp, ieee_quiet_nan)
        poisoned_frequency = transfer(transport_omth, poisoned_frequency)
        call conserving_native_point(theta_, u_, xi_, omega, local_h)
        expected = temperature*elementary_charge*1.0e7_dp*u_**2* &
                   (1 + xi_**2)*config%epsmn*exp(cmplx(0.0_dp, theta_, dp))
        if (.not. all(ieee_is_finite([real(local_h, dp), aimag(local_h)]))) &
            error stop 'local Hamiltonian depends on unset bounce frequency'
        if (abs(local_h - expected) > 5.0e-14_dp*abs(expected)) &
            error stop 'native local Hamiltonian phase/normalization'
        if (native_mth /= saved_mth .or. &
            transfer(native_sign, 0_int64) /= transfer(saved_sign, 0_int64)) &
            error stop 'local Hamiltonian changed orbit context'
        if (transfer(transport_omth, poisoned_frequency) /= poisoned_frequency) &
            error stop 'local Hamiltonian changed bounce frequency'

        ! Independent finite-time partner: preserve the actual orbit phase for
        ! nonzero harmonics on trapped and both passing branches.
        field = field_t*1.0e4_dp*(1.0_dp - &
                    minor_m*sqrt(config%s)/major_m*cos(theta_))
        eta_ = (1.0_dp - xi_**2)/field
        native_mth = 2
        if (eta_ <= etatp) native_mth = -4
        harmonic = real(native_mth, dp)
        if (eta_ <= etatp) harmonic = harmonic + q*config%mph
        state = 0.0_dp
        state(1) = theta_
        state(2) = u_*vth*xi_
        call timestep_transport(u_*vth, eta_, nvar, 0.0_dp, state, derivative)
        finite_time_h = temperature*elementary_charge*1.0e7_dp*u_**2* &
                        cmplx(derivative(3), derivative(4), dp)
        if (.not. all(ieee_is_finite([real(finite_time_h, dp), &
                                      aimag(finite_time_h)]))) &
            error stop 'zero-time nonzero harmonic read bounce frequency'
        if (abs(finite_time_h - expected) > 5.0e-14_dp*abs(expected)) &
            error stop 'native zero-time Hamiltonian phase changed'
        transport_omth = 1700.0_dp
        do direction = -1, 1, 2
            time = direction*0.00073_dp
            call timestep_transport(u_*vth, eta_, nvar, time, state, derivative)
            finite_time_h = temperature*elementary_charge*1.0e7_dp*u_**2* &
                            cmplx(derivative(3), derivative(4), dp)
            finite_time_h = finite_time_h/ &
                exp(cmplx(0.0_dp, -harmonic*time*transport_omth, dp))
            if (.not. all(ieee_is_finite([real(finite_time_h, dp), &
                                          aimag(finite_time_h)]))) &
                error stop 'native finite-time Hamiltonian is not finite'
            if (abs(finite_time_h - expected) > 5.0e-14_dp*abs(expected)) &
                error stop 'native finite-time Hamiltonian phase changed'
        end do
        native_mth = saved_mth
        native_sign = saved_sign
        transport_omth = saved_frequency
    end subroutine check_local_phase

    subroutine close (actual, expected, tolerance, label)
        real(dp), intent(in) :: actual, expected, tolerance
        character(*), intent(in) :: label
        if (abs(actual - expected) > tolerance*abs(expected)) then
            print *, label, actual, expected
            error stop 'independent native normalization'
        end if
    end subroutine close
    subroutine write_chart(path)
        character(*), intent(in) :: path
        real(dp), parameter :: surfaces(3) = [0.25_dp, 0.5_dp, 0.75_dp]
        real(dp) :: radius, bph, bth
        integer :: unit, surface
        open (newunit=unit, file=path, status='replace')
        write (unit, '(a)') 'CC manufactured exact cosine circular field'
        write (unit, '(a)') 'CC constant q=2, one physical magnetic well'
        write (unit, '(a)') 'CC m0b n0b nsurf nper flux a R'
        write (unit, '(a)') 'CC'
        write (unit, '(a)') 'CC'
        write (unit, *) 1, 0, 3, 1, -pi*minor_m**2*field_t, minor_m, major_m
        do surface = 1, size(surfaces)
            radius = minor_m*sqrt(surfaces(surface))
            bph = -field_t*major_m*1.0e6_dp/0.2_dp
            bth = -field_t*radius**2*0.5_dp/major_m*1.0e6_dp/0.2_dp
            write (unit, '(a)') 'CC s iota Jpol Itor pprime sqrtg00'
            write (unit, '(a)') 'CC units A A Pa m3'
            write (unit, *) surfaces(surface), 0.5_dp, bph, bth, 0.0_dp, 0.0_dp
            write (unit, '(a)') 'CC m n rmnc zmnc vmnc bmnc'
            write (unit, *) 0, 0, major_m, 0.0_dp, 0.0_dp, field_t
            write (unit, *) 1, 0, radius, radius, 0.0_dp, -field_t*radius/major_m
        end do
        close (unit)
    end subroutine write_chart
end program test_conserving_native
