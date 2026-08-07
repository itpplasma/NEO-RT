module neort_transport
    use iso_fortran_env, only: dp => real64
    use util, only: imun, pi, c, qi
    use logger, only: trace, debug, warning, error
    use do_magfie_mod, only: do_magfie, s, a, R0, iota, q, psi_pr, eps, inp_swi, &
        bphcov, dbthcovds, dbphcovds, q, dqds, sign_theta, Bthcov
    use do_magfie_pert_mod, only: do_magfie_pert_amp
    use neort_magfie, only: dVds, B0
    use neort_profiles, only: ni1, Om_tE
    use neort_nonlin, only: nonlinear_attenuation
    use neort_freq, only: Om_th, Om_ph
    use neort_gc_frequency_splines, only: &
        evaluate_gc_phase_average_surface, gc_spline_q, &
        evaluate_gc_full_orbit_frequency_surface, &
        evaluate_gc_full_orbit_phase_average_surface
    use neort_gc_frequency_provider, only: gc_full_orbit_frequency_result_t, &
        GC_FREQUENCY_SUCCESS
    use neort_gc_full_resonance, only: GC_RESONANCE_SUCCESS, &
        find_gc_resonances
    use neort_gc_orbit_integrator, only: GC_ORBIT_TRAPPED, GC_ORBIT_PASSING, &
        gc_orbit_average_t
    use neort_orbit, only: bounce_fast, nvar, noshear, poloidal_velocity
    use neort_resonance, only: driftorbit_coarse, driftorbit_root
    use driftorbit, only: vth, mth, mph, mi, B0, Bmin, Bmax, comptorque, epsmn, &
        etamin, etamax, A1, A2, nlev, pertfile, nonlin, m0, etatp, etadt, &
        sign_vpar_htheta, sign_vpar, frequency_model, &
        FREQUENCY_MODEL_GC_THIN, FREQUENCY_MODEL_GC_FULL

    implicit none

    real(dp) :: Omth, dOmthdv, dOmthdeta

contains

    pure function fmt_dbg(msg1, v1, msg2, v2, msg3, v3, msg4, v4) result(s)
        ! Helper to compose a short debug line
        character(*), intent(in) :: msg1, msg2
        character(*), intent(in), optional :: msg3, msg4
        real(dp), intent(in) :: v1, v2
        real(dp), intent(in), optional :: v3, v4
        character(len=256) :: s
        character(len=64) :: a1, a2, a3, a4
        a3 = ''; a4 = ''
        write(a1,'(ES12.5)') v1
        write(a2,'(ES12.5)') v2
        if (present(v3)) write(a3,'(ES12.5)') v3
        if (present(v4)) write(a4,'(ES12.5)') v4
        if (present(msg4)) then
            s = trim(msg1)//trim(a1)//' '//trim(msg2)//trim(a2)//' '//trim(msg3)//trim(a3)//' '//trim(msg4)//trim(a4)
        else if (present(msg3)) then
            s = trim(msg1)//trim(a1)//' '//trim(msg2)//trim(a2)//' '//trim(msg3)//trim(a3)
        else
            s = trim(msg1)//trim(a1)//' '//trim(msg2)//trim(a2)
        end if
    end function fmt_dbg

    ! original contains follows

    pure function D11int(ux, taub, Hmn2)
        real(dp) :: D11int
        real(dp), intent(in) :: ux, taub, Hmn2

        D11int = pi**(3.0_dp / 2.0_dp) * mph**2 * c**2 * q * vth &
            / (qi**2 * dVds * abs(psi_pr)) * ux**3 * exp(-ux**2) * taub * Hmn2
    end function D11int

    pure function D12int(ux, taub, Hmn2)
        real(dp) :: D12int
        real(dp), intent(in) :: ux, taub, Hmn2

        D12int = D11int(ux, taub, Hmn2) * ux**2
    end function D12int

    pure function Tphi_int(ux, taub, Hmn2)
        real(dp) :: Tphi_int
        real(dp), intent(in) :: ux, taub, Hmn2

        Tphi_int = sign(1.0_dp, psi_pr * q * sign_theta) * pi**(3.0_dp / 2.0_dp) * mph**2 * ni1 * &
            c * vth / qi &
            * ux**3 * exp(-ux**2) * taub * Hmn2 * (A1 + A2 * ux**2)
    end function Tphi_int

    subroutine compute_transport_integral(vmin, vmax, vsteps, D, T)
        ! compute transport integral via midpoint rule
        real(dp), intent(in) :: vmin, vmax
        integer, intent(in) :: vsteps
        real(dp), intent(out) :: D(2), T ! Transport coefficients D and torque density T
        real(dp) :: D_plateau, dsdreff ! Plateau diffusion coefficient and ds/dreff=<|grad s|>
        real(dp) :: ux, du, dD11, dD12, dT, v, eta
        real(dp) :: eta_res(2)
        real(dp) :: taub, bounceavg(nvar), Omph, dOmphdv, dOmphdeta
        real(dp) :: q_fieldline
        integer :: istate_dv
        integer :: direct_status, orbit_class
        type(gc_orbit_average_t) :: direct_average
        type(gc_full_orbit_frequency_result_t) :: full_frequency
        real(dp) :: Hmn2, attenuation_factor
        real(dp) :: roots(nlev, 3)
        integer :: nroots, kr, ku
        real(dp) :: full_root_values(nlev), full_root_derivatives(nlev)

        call debug(fmt_dbg('compute_transport_integral: vmin=', vmin, ' vmax=', vmax, ' vsteps=', dble(vsteps)))

        D = 0.0_dp
        T = 0.0_dp
        du = (vmax - vmin) / (vsteps * vth)
        ux = vmin / vth + du / 2.0_dp

        do ku = 1, vsteps
            v = ux * vth
            if (frequency_model == FREQUENCY_MODEL_GC_FULL) then
                call collect_full_orbit_roots(v, full_root_values, &
                    full_root_derivatives, nroots, direct_status)
                if (direct_status /= GC_RESONANCE_SUCCESS) then
                    call warning(fmt_dbg('full GC resonance search failed: v=', &
                        v, ' status=', dble(direct_status)))
                    nroots = 0
                end if
            else
                call driftorbit_coarse(v, etamin, etamax, roots, nroots)
            end if
            ! No explicit nroots==0 guard: the do-loop below is empty when
            ! nroots==0, and a bare `cycle` here would skip the trailing
            ! `ux = ux + du` velocity-grid increment and stall the sweep.
            do kr = 1, nroots
                if (frequency_model == FREQUENCY_MODEL_GC_FULL) then
                    eta_res = [full_root_values(kr), full_root_derivatives(kr)]
                else
                    eta_res = driftorbit_root(v, 1.0e-8_dp * abs(Om_tE), &
                        roots(kr, 1), roots(kr, 2))
                end if
                if (eta_res(1) < 0.0_dp) cycle
                if (abs(eta_res(2)) <= sqrt(epsilon(eta_res(2))) &
                        *max(1.0_dp, abs(Om_tE)) &
                        /max(etamax - etamin, tiny(eta_res(2)))) then
                    call warning(fmt_dbg('ill-conditioned resonance skipped: eta=', &
                        eta_res(1), ' derivative=', eta_res(2)))
                    cycle
                end if
                eta = eta_res(1)

                call Om_th(v, eta, Omth, dOmthdv, dOmthdeta)

                taub = 2.0_dp * pi / abs(Omth)
                if (frequency_model == FREQUENCY_MODEL_GC_FULL) then
                    orbit_class = merge(GC_ORBIT_TRAPPED, GC_ORBIT_PASSING, &
                        eta > etatp)
                    call evaluate_gc_full_orbit_frequency_surface(v, eta, &
                        int(sign_vpar), orbit_class, taub, full_frequency, &
                        direct_status)
                    if (direct_status /= GC_FREQUENCY_SUCCESS) cycle
                    Omth = full_frequency%omega_b
                    Omph = full_frequency%omega_phi
                    taub = full_frequency%period
                    call evaluate_gc_full_orbit_phase_average_surface(v, eta, &
                        int(sign_vpar), orbit_class, taub, mth, mph, &
                        evaluate_direct_perturbation, direct_average, direct_status)
                    if (direct_status /= 0) cycle
                    bounceavg = 0.0_dp
                    bounceavg(3) = real(direct_average%perturbation_average)
                    bounceavg(4) = aimag(direct_average%perturbation_average)
                    bounceavg(5) = direct_average%inverse_b_average
                    bounceavg(6) = direct_average%b_average
                    istate_dv = 2
                else if (frequency_model == FREQUENCY_MODEL_GC_THIN) then
                    orbit_class = merge(GC_ORBIT_TRAPPED, GC_ORBIT_PASSING, &
                        eta > etatp)
                    call Om_ph(v, eta, Omph, dOmphdv, dOmphdeta)
                    q_fieldline = gc_spline_q()
                    call evaluate_gc_phase_average_surface(v, eta, &
                        int(sign_vpar), orbit_class, taub, Omth, Omph, &
                        q_fieldline, mth, mph, evaluate_direct_perturbation, &
                        direct_average, direct_status)
                    if (direct_status /= 0) then
                        call error(fmt_dbg('direct GC orbit average failed: mth=', &
                            dble(mth), ' eta=', eta))
                    end if
                    bounceavg = 0.0_dp
                    bounceavg(3) = real(direct_average%perturbation_average)
                    bounceavg(4) = aimag(direct_average%perturbation_average)
                    bounceavg(5) = direct_average%inverse_b_average
                    bounceavg(6) = direct_average%b_average
                    istate_dv = 2
                else
                    call bounce_fast(v, eta, taub, bounceavg, &
                        timestep_transport, istate_dv)
                end if
                if (istate_dv == -1) then
                    call error(fmt_dbg('VODE MXSTEP: mth=', dble(mth), ' ux=', ux, ' eta=', eta, ' taub=', taub))
                else if (istate_dv /= 2) then
                    call warning(fmt_dbg('dvode istate=', dble(istate_dv), ' at mth=', dble(mth), ' ux=', ux, ' eta=', eta))
                else
                    if (abs(eta - etatp) < 1.0e-8_dp*etatp) then
                        call trace(fmt_dbg('near etatp: mth=', dble(mth), ' ux=', ux, ' eta=', eta, ' taub=', taub))
                    end if
                end if
                Hmn2 = (bounceavg(3)**2 + bounceavg(4)**2) * (mi * (ux * vth)**2 / 2.0_dp)**2
                attenuation_factor = nonlinear_attenuation(ux, eta, bounceavg, Omth, &
                    dOmthdv, dOmthdeta, Hmn2)

                dD11 = du * D11int(ux, taub, Hmn2) / abs(eta_res(2))
                dD12 = du * D12int(ux, taub, Hmn2) / abs(eta_res(2))
                D(1) = D(1) + dD11 * attenuation_factor
                D(2) = D(2) + dD12 * attenuation_factor

                if (comptorque) then
                    dT = du * Tphi_int(ux, taub, Hmn2) / abs(eta_res(2))
                    T = T + dT * attenuation_factor
                end if
            end do
            ux = ux + du
        end do

        D_plateau = pi * vth**3 / (16.0_dp * R0 * iota * (qi * B0 / (mi * c))**2)
        dsdreff = 2.0_dp / a * sqrt(s) ! TODO: Use exact value instead of this approximation
        D = dsdreff**(-2) * D / D_plateau

        call debug("compute_transport_integral complete")

    contains

        subroutine collect_full_orbit_roots(velocity, root_values, &
                root_derivatives, root_count, root_status)
            real(dp), intent(in) :: velocity
            real(dp), intent(out) :: root_values(:), root_derivatives(:)
            integer, intent(out) :: root_count, root_status

            real(dp) :: region_roots(size(root_values))
            real(dp) :: region_derivatives(size(root_values))
            real(dp) :: lower, upper
            integer :: region_count, region_status

            root_values = 0.0_dp
            root_derivatives = 0.0_dp
            root_count = 0
            root_status = GC_RESONANCE_SUCCESS

            lower = etamin
            upper = min(etamax, etatp*(1.0_dp - 1.0e-8_dp))
            if (upper > lower) then
                call find_gc_resonances(full_residual, lower, upper, nlev, &
                    max(1.0e-8_dp*abs(Om_tE), 1.0e-6_dp), 1.0e-10_dp, &
                    region_roots, region_derivatives, region_count, region_status)
                if (region_status /= GC_RESONANCE_SUCCESS) then
                    root_status = region_status
                    return
                end if
                root_values(1:region_count) = region_roots(1:region_count)
                root_derivatives(1:region_count) = &
                    region_derivatives(1:region_count)
                root_count = region_count
            end if

            lower = max(etamin, etatp*(1.0_dp + 1.0e-8_dp))
            upper = etamax
            if (upper > lower .and. root_count < size(root_values)) then
                call find_gc_resonances(full_residual, lower, upper, nlev, &
                    max(1.0e-8_dp*abs(Om_tE), 1.0e-6_dp), 1.0e-10_dp, &
                    region_roots, region_derivatives, region_count, region_status)
                if (region_status /= GC_RESONANCE_SUCCESS) then
                    root_status = region_status
                    return
                end if
                region_count = min(region_count, size(root_values) - root_count)
                root_values(root_count + 1:root_count + region_count) = &
                    region_roots(1:region_count)
                root_derivatives(root_count + 1:root_count + region_count) = &
                    region_derivatives(1:region_count)
                root_count = root_count + region_count
            end if
        end subroutine collect_full_orbit_roots

        subroutine full_residual(pitch, residual, residual_status)
            real(dp), intent(in) :: pitch
            real(dp), intent(out) :: residual
            integer, intent(out) :: residual_status

            real(dp) :: thin_omega_b, thin_dv, thin_deta, period_estimate
            integer :: local_class
            type(gc_full_orbit_frequency_result_t) :: local_frequency

            residual = 0.0_dp
            call Om_th(v, pitch, thin_omega_b, thin_dv, thin_deta)
            if (thin_omega_b == 0.0_dp) then
                residual_status = 1
                return
            end if
            period_estimate = 2.0_dp*pi/abs(thin_omega_b)
            local_class = merge(GC_ORBIT_TRAPPED, GC_ORBIT_PASSING, &
                pitch > etatp)
            call evaluate_gc_full_orbit_frequency_surface(v, pitch, &
                int(sign_vpar), local_class, period_estimate, local_frequency, &
                residual_status)
            if (residual_status /= GC_FREQUENCY_SUCCESS) return
            residual = real(mth, dp)*local_frequency%omega_b &
                +real(mph, dp)*local_frequency%omega_phi
        end subroutine full_residual

    end subroutine compute_transport_integral

    subroutine evaluate_direct_perturbation(position, bmod, amplitude, status)
        real(dp), intent(in) :: position(3), bmod
        complex(dp), intent(out) :: amplitude
        integer, intent(out) :: status
        complex(dp) :: field_amplitude

        amplitude = (0.0_dp, 0.0_dp)
        status = 1
        if (bmod <= 0.0_dp) return
        if (pertfile) then
            call do_magfie_pert_amp(position, field_amplitude)
            amplitude = epsmn*field_amplitude/bmod
        else
            amplitude = epsmn*exp(imun*real(m0, dp)*position(3))
        end if
        status = 0
    end subroutine evaluate_direct_perturbation

    subroutine timestep_transport(v, eta, neq, t, y, ydot)
        !
        !  Timestep function for orbit integration.
        !  Includes poloidal angle theta and parallel velocity.
        !  More integrands may be added starting from y(3)
        !

        real(dp), intent(in) :: v, eta
        integer, intent(in) :: neq
        real(dp), intent(in) :: t
        real(dp), intent(in) :: y(neq)
        real(dp), intent(out) :: ydot(neq)

        ! BEGIN TODO: remove all of this after refactoring and re-use routine in orbit
        ! for y(1:3)
        real(dp) :: bmod, sqrtg, x(3), hder(3), hcovar(3), hctrvr(3), hcurl(3), Om_tB_v
        real(dp) :: t0, parallel_phase
        complex(dp) :: epsn, Hn ! relative amplitude of perturbation field epsn=Bn/B0
        ! and Hamiltonian Hn = (H - H0)_n

        x(1) = s
        x(2) = 0.0_dp
        x(3) = y(1)
        call do_magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
        call poloidal_velocity(v, eta, bmod, hctrvr(3), hder(3), y(2), ydot)
        if (inp_swi == 11) then
            ydot(7) = y(2)*hctrvr(2)
            parallel_phase = mph*y(7)
        else
            ydot(7) = 0.0_dp
            parallel_phase = q*mph*y(1)
        end if

        ! evaluate orbit averages of Hamiltonian perturbation
        if (pertfile) then
            call do_magfie_pert_amp(x, epsn)
            epsn = epsmn * epsn / bmod
        else
            epsn = epsmn * exp(imun * m0 * y(1))
        end if

        if (eta > etatp) then
            !t0 = 0.25*2*pi/Omth ! Different starting position in orbit
            t0 = 0.0_dp
            Hn = (2.0_dp - eta * bmod) * epsn * exp(imun * (parallel_phase - mth * (t - &
                t0) * Omth))
        else
            Hn = (2.0_dp - eta * bmod) * epsn * exp(imun * (parallel_phase - (mth + q * mph) &
                * t * Omth))
        end if
        ydot(3) = real(Hn)
        ydot(4) = aimag(Hn)

        ! evaluate orbit averages for nonlinear attenuation
        if (nonlin) then
            ydot(5) = 1.0_dp / bmod
            ydot(6) = bmod
        else
            ydot(5:6) = 0.0_dp
        end if

        ! Zero any trailing integrands this routine does not compute (the unused
        ! abs(B) slot, y(7)). The DOP853 solver carries every component through
        ! its stage combinations, so an uninitialised derivative leaves a
        ! denormal in the state that trips the underflow FPE trap.
        if (inp_swi /= 11) ydot(7:) = 0.0_dp
    end subroutine timestep_transport

end module neort_transport
