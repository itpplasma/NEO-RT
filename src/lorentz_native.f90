module neort_lorentz_native
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use neort_lorentz_graph, only: lorentz_branch_t, solve_lorentz_graph, &
        neort_lorentz_diffusion
    implicit none
    private

    abstract interface
        function lorentz_rate_i(speed_cm_s) result(rate)
            import :: dp
            real(dp), intent(in) :: speed_cm_s
            real(dp) :: rate
        end function lorentz_rate_i
    end interface

    public :: lorentz_rate_i, validate_lorentz_geometry, sample_lorentz_orbit
    public :: compute_lorentz_velocity, compute_lorentz_transport

contains

    subroutine validate_lorentz_geometry(info)
        use do_magfie_mod, only: do_magfie, inp_swi, s, q, dqds
#ifdef NEORT_STANDALONE
        use do_magfie_mod, only: m0b, n0b, nmode, modes0
#endif
        use driftorbit, only: mph, drift_phase, fow_radial, nonlin, supban, &
            nopassing, Bmin, Bmax, magdrift, magdrift_passing
        use neort_orbit, only: th0
        use util, only: pi
        integer, intent(out) :: info
        real(dp) :: field(32), bmod, sqrtg, hder(3), hcovar(3), hctrvr(3)
        real(dp) :: hcurl(3), x(3), mean_b, harmonic_b, theta
        integer :: k

        info = 0
        if (drift_phase .or. fow_radial .or. nonlin .or. supban .or. nopassing) then
            info = -10
            return
        end if
        if (magdrift .and. magdrift_passing /= 1) then
            info = -10
            return
        end if
        if (inp_swi /= 8 .or. abs(th0) > 1.0e-13_dp) then
            info = -11
            return
        end if
#ifndef NEORT_STANDALONE
        info = -11
        return
#else
        ! Sampling alone cannot certify a spectrum: m=33 aliases m=1 on a
        ! 32-point grid. Restrict this starter to the two actual input modes.
        if (m0b /= 1 .or. n0b /= 0 .or. nmode /= 2) then
            info = -14
            return
        end if
        if (.not. allocated(modes0)) then
            info = -14
            return
        end if
        if (any(modes0(1, :, 1) /= [0.0_dp, 1.0_dp]) .or. &
                any(modes0(1, :, 2) /= 0.0_dp)) then
            info = -14
            return
        end if
#endif
        if (.not. all(ieee_is_finite([q, dqds, Bmin, Bmax]))) then
            info = -12
            return
        end if
        if (abs(dqds) > 1.0e-11_dp*max(1.0_dp, abs(q))) then
            info = -12
            return
        end if
        if (abs(mph*q - real(nint(mph*q), dp)) > 1.0e-11_dp) then
            info = -13
            return
        end if
        if (Bmin <= 0.0_dp .or. Bmax <= Bmin) then
            info = -14
            return
        end if
        ! First supported geometry is a circular cosine well, with its minimum
        ! at the native orbit launch point theta=0. Reject additional wells and
        ! shaping rather than assume their bounce phases glue identically.
        do k = 1, size(field)
            theta = 2.0_dp*pi*real(k - 1, dp)/size(field)
            x = [s, 0.0_dp, theta]
            call do_magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
            field(k) = bmod
        end do
        mean_b = sum(field)/size(field)
        harmonic_b = 0.0_dp
        do k = 1, size(field)
            theta = 2.0_dp*pi*real(k - 1, dp)/size(field)
            harmonic_b = harmonic_b + 2.0_dp*field(k)*cos(theta)/size(field)
        end do
        if (harmonic_b >= 0.0_dp) then
            info = -14
            return
        end if
        do k = 1, size(field)
            theta = 2.0_dp*pi*real(k - 1, dp)/size(field)
            if (abs(field(k) - mean_b - harmonic_b*cos(theta)) > &
                1.0e-10_dp*mean_b) then
                info = -14
                return
            end if
        end do
    end subroutine validate_lorentz_geometry

    subroutine sample_lorentz_orbit(ux, branch, eta, nud, nu, mass, diffusion, &
            omega, drive, drift_frequency, inverse_b_average)
        use util, only: pi, mi
        use do_magfie_mod, only: q
        use driftorbit, only: mth, mph, sign_vpar, etatp, etadt
        use neort_profiles, only: vth
        use neort_orbit, only: bounce_fast, bounce_time, nvar
        use driftorbit, only: sign_vpar_htheta
        use neort_profiles, only: Om_tE
        use neort_transport, only: transport_omth => Omth
        real(dp), intent(in) :: ux, eta, nud, nu
        integer, intent(in) :: branch
        real(dp), intent(out) :: mass, diffusion, omega
        complex(dp), intent(out) :: drive
        real(dp), intent(out), optional :: drift_frequency
        real(dp), intent(out), optional :: inverse_b_average
        real(dp) :: omth, average(nvar)
        real(dp) :: old_sign, inverse_b
        integer :: old_mth, status

        if (ux <= 0.0_dp) error stop 'Lorentz orbit needs positive speed'
        if (eta <= 0.0_dp .or. eta >= etadt .or. eta == etatp) &
            error stop 'Lorentz coefficients must be sampled inside physical edges'
        if (nud < 0.0_dp .or. nu < 0.0_dp) &
            error stop 'Lorentz and reaction rates must be nonnegative'
        old_sign = sign_vpar
        old_mth = mth
        select case (branch)
        case (1)
            if (eta <= etatp) error stop 'trapped Lorentz sample outside trapped edge'
            sign_vpar = 1.0_dp
            mth = 0
        case (2, 3)
            if (eta >= etatp) error stop 'passing Lorentz sample outside passing edge'
            sign_vpar = merge(1.0_dp, -1.0_dp, branch == 2)
            mth = -nint(mph*q)
        case default
            error stop 'Lorentz graph requires trapped and both passing branches'
        end select
        ! Integrate the actual orbit at every quadrature point. The graph
        ! reaches the physical endpoints beyond the legacy spline knots.
        mass = bounce_time(ux*vth, eta)
        omth = sign_vpar_htheta*2.0_dp*pi/mass
        transport_omth = omth
        call bounce_fast(ux*vth, eta, mass, average, timestep_lorentz, status, &
            zero_integrals=.true.)
        if (status /= 2) error stop 'Lorentz native bounce integral failed'
        omega = real(mph, dp)*(Om_tE + average(7))
        if (present(drift_frequency)) drift_frequency = real(mph, dp)*average(7)
        drive = (mi*(ux*vth)**2/2.0_dp)*cmplx(average(3), average(4), dp)
        inverse_b = average(5) + eta
        if (present(inverse_b_average)) inverse_b_average = inverse_b
        call neort_lorentz_diffusion(nud, eta, mass, inverse_b, diffusion)
        sign_vpar = old_sign
        mth = old_mth
    end subroutine sample_lorentz_orbit

    subroutine timestep_lorentz(v, eta, neq, time, state, derivative)
        use neort_transport, only: timestep_transport
        use do_magfie_mod, only: do_magfie, s
        use neort_orbit, only: Om_tB_local_v2
        use driftorbit, only: magdrift, magdrift_passing, etatp
        real(dp), intent(in) :: v, eta, time
        integer, intent(in) :: neq
        real(dp), intent(in) :: state(neq)
        real(dp), intent(out) :: derivative(neq)
        real(dp) :: x(3), bmod, sqrtg, hder(3), hcovar(3), hctrvr(3), hcurl(3)

        call timestep_transport(v, eta, neq, time, state, derivative)
        x = [s, 0.0_dp, state(1)]
        call do_magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
        derivative(5) = 1.0_dp/bmod - eta
        derivative(7) = 0.0_dp
        if (magdrift) then
            if (eta > etatp .or. magdrift_passing > 0) then
                derivative(7) = Om_tB_local_v2(eta, bmod, hder(1), &
                    hctrvr(3))*v**2
            end if
        end if
    end subroutine timestep_lorentz

    subroutine compute_lorentz_velocity(ux, pitch_steps, nud, nu, work, info)
        use driftorbit, only: etatp, etadt
        real(dp), intent(in) :: ux, nud, nu
        integer, intent(in) :: pitch_steps
        complex(dp), intent(out) :: work(3)
        integer, intent(out) :: info
        type(lorentz_branch_t) :: branches(3)
        real(dp) :: coordinate
        integer :: branch, k

        work = (0.0_dp, 0.0_dp)
        if (pitch_steps < 2) then
            info = -15
            return
        end if
        do branch = 1, 3
            allocate(branches(branch)%eta(pitch_steps + 1))
            do k = 1, pitch_steps + 1
                coordinate = real(k - 1, dp)/pitch_steps
                if (branch == 1) then
                    branches(branch)%eta(k) = etatp + (etadt - etatp)*coordinate
                else
                    branches(branch)%eta(k) = etatp*(1.0_dp - coordinate)
                end if
            end do
        end do
        ! The native state belongs to the caller's surface/species thread.
        ! Avoid creating new uninitialized threadprivate native orbit contexts.
        call solve_lorentz_graph(branches, sample, info, parallel_samples=.false.)
        if (info /= 0) return
        work = branches%work
    contains
        subroutine sample(branch, eta, mass, diffusion, omega, reaction, drive)
            integer, intent(in) :: branch
            real(dp), intent(in) :: eta
            real(dp), intent(out) :: mass, diffusion, omega, reaction
            complex(dp), intent(out) :: drive
            reaction = nu
            call sample_lorentz_orbit(ux, branch, eta, nud, nu, mass, diffusion, &
                omega, drive)
        end subroutine sample
    end subroutine compute_lorentz_velocity

    subroutine compute_lorentz_transport(vmin, vmax, vsteps, pitch_steps, &
            deflection_rate, diffusion, torque, info)
        use util, only: pi, mi, qi, c
        use do_magfie_mod, only: R0, iota, a, s
        use neort_profiles, only: vth, A1, A2
        use driftorbit, only: B0, comptorque
        use neort_transport, only: D11int, D12int, Tphi_int
        use neort_lorentz_controls, only: lorentz_susceptibility_valid
        real(dp), intent(in) :: vmin, vmax
        integer, intent(in) :: vsteps, pitch_steps
        procedure(lorentz_rate_i) :: deflection_rate
        real(dp), intent(out) :: diffusion(2, 3), torque(3)
        integer, intent(out) :: info
        complex(dp) :: work(3)
        real(dp) :: ux, du, nodal_rate, integral(3)
        real(dp) :: plateau, coordinate_factor
        integer :: velocity
        logical :: zero_physical_force

        diffusion = 0.0_dp
        torque = 0.0_dp
        lorentz_susceptibility_valid = .true.
        zero_physical_force = .false.
        if (comptorque) zero_physical_force = A1 == 0.0_dp .and. A2 == 0.0_dp
        call validate_lorentz_geometry(info)
        if (info /= 0) return
        if (pitch_steps < 2) then
            info = -15
            return
        end if
        if (vsteps < 1 .or. vmin < 0.0_dp .or. vmax <= vmin) then
            info = -16
            return
        end if
        du = (vmax - vmin)/(vsteps*vth)
        do velocity = 1, vsteps
            ux = vmin/vth + (velocity - 0.5_dp)*du
            nodal_rate = deflection_rate(ux*vth)
            if (nodal_rate == 0.0_dp) then
                ! A finite pitch mesh has no certified retarded delta limit.
                ! Only the exactly unforced physical response is defined here.
                info = -17
                if (zero_physical_force) then
                    diffusion = 0.0_dp
                    torque = 0.0_dp
                    lorentz_susceptibility_valid = .false.
                    info = 0
                end if
                return
            end if
            call compute_lorentz_velocity(ux, pitch_steps, nodal_rate, 0.0_dp, &
                work, info)
            if ((info == 2 .or. info == -5) .and. zero_physical_force) then
                ! Physical g has exactly zero source at every speed. The
                ! auxiliary susceptibility z has no stationary solution here;
                ! its zero diagnostic placeholders must be marked undefined.
                diffusion = 0.0_dp
                torque = 0.0_dp
                lorentz_susceptibility_valid = .false.
                info = 0
                return
            end if
            if (info /= 0) return
            ! Real susceptibility work replaces tau*|H|^2*delta(R).
            ! The 1/pi is required by the causal resolvent's delta limit.
            integral = real(work, dp)/pi
            diffusion(1, :) = diffusion(1, :) + &
                du*D11int(ux, 1.0_dp, 1.0_dp)*integral
            diffusion(2, :) = diffusion(2, :) + &
                du*D12int(ux, 1.0_dp, 1.0_dp)*integral
            if (comptorque) then
                torque = torque + du*Tphi_int(ux, 1.0_dp, 1.0_dp)*integral
            end if
        end do
        plateau = pi*vth**3/(16.0_dp*R0*iota*(qi*B0/(mi*c))**2)
        coordinate_factor = 2.0_dp/a*sqrt(s)
        diffusion = coordinate_factor**(-2)*diffusion/plateau
    end subroutine compute_lorentz_transport

end module neort_lorentz_native
