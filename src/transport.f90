module neort_transport
    use neort_mixture_collisions, only: finite_krook_active, effective_krook_rate
    use iso_fortran_env, only: dp => real64, error_unit
    use ieee_arithmetic, only: ieee_is_finite
    use util, only: imun, pi, c, qi
    use logger, only: trace, debug, warning, error
    use do_magfie_mod, only: do_magfie, s, a, R0, iota, q, psi_pr, eps, &
        bphcov, dbthcovds, dbphcovds, q, dqds, sign_theta, Bthcov
    use do_magfie_pert_mod, only: do_magfie_pert_amp, &
        magnetic_amplitude_bound => perturbation_amplitude_bound, &
        do_magfie_pert_kappa_amp, curvature_drive
    use neort_phi1, only: phi1, do_elec_pert_amp, phi1_amplitude_bound
    use neort_magfie, only: dVds, B0
    use neort_profiles, only: ns, Om_tE
    use neort_species, only: selected_species
    use neort_collisional_layer, only: collisional_layer, collisional_layer_factor, &
        deflection_frequency
    use neort_nonlin, only: nonlinear_attenuation
    use neort_freq, only: Om_th, Om_ph, Om_tB
    use neort_orbit, only: bounce_fast, nvar, poloidal_velocity, Om_tB_local_v2, &
        bounce_rtol, bounce_atol, th0, vpar
    use neort_resonance, only: driftorbit_coarse, driftorbit_root, resonance_value
    use neort_resonance_cutout, only: resonance_cutout, resonance_cutout_width, &
        CUTOUT_OFF, bounce_pitch_diffusion, cutout_weight, cutout_curvature
    use neort_krook, only: krook_collisions, krook_nu, krook_model, krook_harmonic, &
        krook_geometry_trapped, krook_geometry_passing, krook_eta_steps, &
        krook_max_depth, krook_rtol, krook_atol, krook_rate, integrate_krook
    use kinetic_energy_selection, only: maxwellian_speed_density, &
        energy_sample_count, energy_sample_speed, energy_sample_weight
    use driftorbit, only: vth, mth, mph, mi, B0, Bmin, Bmax, comptorque, epsmn, &
        pertfile_scale, &
        etamin, etamax, A1, A2, nlev, pertfile, nonlin, m0, etatp, etadt, &
        sign_vpar_htheta, sign_vpar, drift_phase, magdrift, magdrift_passing, &
        fow_radial

  implicit none

  ! Relative tolerance of the adaptive velocity integral of the collisionless
  ! resonance components; zero selects the deprecated vsteps-node midpoint
  ! rule. Folds of a resonance line make the integrand singular like
  ! 1/sqrt(u - u_fold), where the midpoint rule converges like sqrt(du) and
  ! erratically.
  real(dp) :: vtol = 1.0e-3_dp
  integer, parameter :: vtol_initial_panels = 16, vtol_max_panels = 4000
  ! Error budget of the serial driver (compute_transport), which runs the
  ! components one by one: coarse sizing pass and per-component share.
  logical :: serial_coarse = .false.
  real(dp) :: serial_atol(3) = 0.0_dp, serial_sum(3) = 0.0_dp

  ! These values remain live while bounce_fast invokes timestep_transport.
  ! Each transport worker needs an independent callback state.
  real(dp) :: Omth, dOmthdv, dOmthdeta
  !$omp threadprivate (Omth, dOmthdv, dOmthdeta)

  abstract interface
    function drift_osc_i(t) result(Om_osc)
      import :: dp
      real(dp), intent(in) :: t
      real(dp) :: Om_osc
    end function drift_osc_i
  end interface

  ! Test seam: when associated (and drift_phase is set), replaces the local
  ! oscillating toroidal drift Om_tB(theta) - <Om_tB> [rad/s] by a prescribed
  ! function of orbit time, so tests can use analytic Bessel sidebands.
  procedure(drift_osc_i), pointer :: drift_osc_override => null()

  ! Contributions of one velocity node to the D11, D12 and torque integrals:
  ! one term per resonant root (collisionless) or the finite-width pitch
  ! integral (Krook). assemble_velocity_samples sums them in node and root
  ! order, independently of where and in which order they were computed.
  type :: velocity_sample_t
      integer :: nterm = 0
      real(dp), allocatable :: term(:, :)  ! (3, capacity); only 1:nterm is used
      ! Adaptive integrals: coarse = initial panels only; atol is the share of
      ! the error budget of the summed harmonics granted to this component.
      logical :: coarse = .false.
      real(dp) :: atol(3) = 0.0_dp
  end type velocity_sample_t

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

        D11int = maxwellian_speed_density(ux)*D11mono(ux, taub, Hmn2)
    end function D11int

    pure function D11mono(ux, taub, Hmn2)
        real(dp) :: D11mono
        real(dp), intent(in) :: ux, taub, Hmn2

        D11mono = pi**2/4.0_dp*mph**2*c**2*q*vth &
                  /(qi**2*dVds*abs(psi_pr))*ux*taub*Hmn2
    end function D11mono

    pure function D12int(ux, taub, Hmn2)
        real(dp) :: D12int
        real(dp), intent(in) :: ux, taub, Hmn2

        D12int = D11int(ux, taub, Hmn2) * ux**2
    end function D12int

    pure function D12mono(ux, taub, Hmn2)
        real(dp) :: D12mono
        real(dp), intent(in) :: ux, taub, Hmn2

        D12mono = D11mono(ux, taub, Hmn2)*ux**2
    end function D12mono

    pure function Tphi_int(ux, taub, Hmn2)
        real(dp) :: Tphi_int
        real(dp), intent(in) :: ux, taub, Hmn2

        Tphi_int = maxwellian_speed_density(ux)*Tphi_mono(ux, taub, Hmn2)
    end function Tphi_int

    pure function Tphi_mono(ux, taub, Hmn2)
        real(dp) :: Tphi_mono
        real(dp), intent(in) :: ux, taub, Hmn2

        Tphi_mono = sign(1.0_dp, psi_pr*q*sign_theta)*pi**2/4.0_dp*mph**2 &
                    *ns*c*vth/qi*ux*taub*Hmn2*(A1 + A2*ux**2)
    end function Tphi_mono

    subroutine compute_transport_integral(vmin, vmax, vsteps, energy_ratio, D, T)
        ! Velocity integral of one resonance component, node by node in order.
        real(dp), intent(in) :: vmin, vmax
        integer, intent(in) :: vsteps
        real(dp), intent(in) :: energy_ratio
        real(dp), intent(out) :: D(2), T  ! Transport coefficients D and torque density T
        type(velocity_sample_t), allocatable :: samples(:)
        integer :: ku

        call debug(fmt_dbg('compute_transport_integral: vmin=', vmin, ' vmax=', vmax, ' vsteps=', dble(vsteps)))

        allocate(samples(velocity_sample_count(vsteps, energy_ratio)))
        call prepare_velocity_samples(samples)
        samples%coarse = serial_coarse
        do ku = 1, size(samples)
            samples(ku)%atol = serial_atol
        end do
        do ku = 1, size(samples)
            call compute_velocity_sample(ku, vmin, vmax, vsteps, energy_ratio, samples(ku))
        end do
        if (serial_coarse) then
            do ku = 1, size(samples)
                if (samples(ku)%nterm > 0) serial_sum = serial_sum &
                    + abs(sum(samples(ku)%term(:, 1:samples(ku)%nterm), dim=2))
            end do
        end if
        call assemble_velocity_samples(vmin, vmax, vsteps, energy_ratio, samples, D, T)

        call debug("compute_transport_integral complete")
    end subroutine compute_transport_integral

    integer function velocity_sample_count(vsteps, energy_ratio)
        integer, intent(in) :: vsteps
        real(dp), intent(in) :: energy_ratio

        velocity_sample_count = energy_sample_count(energy_ratio, vsteps)
        ! The adaptive integral of a component is one sample.
        if (adaptive_velocity(energy_ratio)) velocity_sample_count = 1
    end function velocity_sample_count

    logical function adaptive_velocity(energy_ratio)
        ! Adaptive integration applies to the collisionless thermal integral.
        real(dp), intent(in) :: energy_ratio

        adaptive_velocity = vtol > 0.0_dp .and. energy_ratio <= 0.0_dp .and. &
            .not. finite_krook_active()
    end function adaptive_velocity

    subroutine compute_velocity_sample(ku, vmin, vmax, vsteps, energy_ratio, sample)
        ! Contributions of velocity node ku. The nodes are independent: each
        ! depends only on the per-species state and the resonance component
        ! (mth, sign_vpar, etamin, etamax) held by the calling thread.
        integer, intent(in) :: ku
        real(dp), intent(in) :: vmin, vmax
        integer, intent(in) :: vsteps
        real(dp), intent(in) :: energy_ratio
        type(velocity_sample_t), intent(inout) :: sample
        logical :: ok

        if (collisional_layer .and. finite_krook_active()) &
            error stop 'collisional_layer and finite Krook broadening cannot be combined'
        if (finite_krook_active()) then
            block
                character(len=4096) :: ledger_path
                call resonance_ledger_path(ledger_path)
            end block
            call compute_krook_sample(ku, vmin, vmax, vsteps, energy_ratio, sample)
        else if (adaptive_velocity(energy_ratio)) then
            ok = .false.
            if (resonance_cutout == CUTOUT_OFF .and. .not. ledger_active()) &
                call compute_pitch_sample(vmin, vmax, sample, ok)
            if (.not. ok) call compute_adaptive_sample(vmin, vmax, sample)
        else
            call compute_collisionless_sample(ku, vmin, vmax, vsteps, energy_ratio, sample)
        end if
    end subroutine compute_velocity_sample

    subroutine assemble_velocity_samples(vmin, vmax, vsteps, energy_ratio, samples, D, T)
        ! Sums the node contributions in node and root order and normalizes.
        real(dp), intent(in) :: vmin, vmax
        integer, intent(in) :: vsteps
        real(dp), intent(in) :: energy_ratio
        type(velocity_sample_t), intent(in) :: samples(:)
        real(dp), intent(out) :: D(2), T
        real(dp) :: du, integral(3), D_plateau, dsdreff
        integer :: ku, kr

        if (energy_ratio > 0.0_dp) then
            du = 0.0_dp
        else
            du = (vmax - vmin)/(vsteps*vth)
        end if
        D_plateau = pi * vth**3 / (16.0_dp * R0 * iota * (qi * B0 / (mi * c))**2)
        dsdreff = 2.0_dp / a * sqrt(s)  ! TODO: Use exact value instead of this approximation

        if (finite_krook_active()) then
            integral = 0.0_dp
            do ku = 1, size(samples)
                ! thermal weights from compute_transport_weight already carry the
                ! Maxwellian speed density; monoenergetic ones are weight 1.
                if (energy_ratio > 0.0_dp) then
                    integral = integral + samples(ku)%term(:, 1)
                else
                    integral = integral + du*samples(ku)%term(:, 1)
                end if
            end do
            D = dsdreff**(-2)*integral(1:2)/D_plateau
            T = integral(3)
        else
            D = 0.0_dp
            T = 0.0_dp
            do ku = 1, size(samples)
                do kr = 1, samples(ku)%nterm
                    D(1) = D(1) + samples(ku)%term(1, kr)
                    D(2) = D(2) + samples(ku)%term(2, kr)
                    if (comptorque) T = T + samples(ku)%term(3, kr)
                end do
            end do
            D = dsdreff**(-2) * D / D_plateau
        end if
    end subroutine assemble_velocity_samples

    impure elemental subroutine prepare_velocity_samples(sample)
        ! Reserve result storage before dispatch, never in a velocity worker.
        ! nlev bounds the root search; Krook needs only its pitch integral.
        type(velocity_sample_t), intent(inout) :: sample
        integer :: capacity

        capacity = nlev
        if (finite_krook_active()) capacity = 1
        if (allocated(sample%term)) deallocate(sample%term)
        allocate(sample%term(3, capacity))
        sample%nterm = 0
    end subroutine prepare_velocity_samples

    subroutine store_sample_terms(terms, sample)
        real(dp), intent(in) :: terms(:, :)
        type(velocity_sample_t), intent(inout) :: sample

        sample%nterm = size(terms, 2)
        sample%term(:, 1:sample%nterm) = terms
    end subroutine store_sample_terms

    subroutine compute_collisionless_sample(ku, vmin, vmax, vsteps, energy_ratio, sample)
        ! Midpoint-rule node ku: one term per resonant root.
        integer, intent(in) :: ku
        real(dp), intent(in) :: vmin, vmax
        integer, intent(in) :: vsteps
        real(dp), intent(in) :: energy_ratio
        type(velocity_sample_t), intent(inout) :: sample
        real(dp) :: ux, du, terms(3, nlev)
        integer :: nterm

        if (energy_ratio > 0.0_dp) then
            du = 0.0_dp
        else
            du = (vmax - vmin)/(vsteps*vth)
        end if
        ux = energy_sample_speed(energy_ratio, ku, vmin/vth, vmax/vth, vsteps)
        call collisionless_point(ux, energy_sample_weight(energy_ratio, ux, du), terms, nterm)
        call store_sample_terms(terms(:, 1:nterm), sample)
    end subroutine compute_collisionless_sample

    subroutine compute_adaptive_sample(vmin, vmax, sample)
        ! The velocity integral of one collisionless component as a single
        ! term: adaptive 7-15 point Gauss-Kronrod in u = v/vth on [vmin, vmax],
        ! bisecting the panel with the largest relative error estimate until
        ! the summed estimate of each of D11, D12 and the torque is below vtol
        ! of its integral. Panels start uniform and are refined in a fixed
        ! order, so the result does not depend on threads.
        real(dp), intent(in) :: vmin, vmax
        type(velocity_sample_t), intent(inout) :: sample
        real(dp), allocatable :: lo(:), hi(:), val(:, :), err(:, :)
        real(dp) :: total(3), scale(3), weight(3), a, b
        integer :: np, k, worst

        allocate (lo(vtol_max_panels), hi(vtol_max_panels), &
            val(3, vtol_max_panels), err(3, vtol_max_panels))
        np = vtol_initial_panels
        do k = 1, np
            lo(k) = vmin/vth + (vmax - vmin)/vth*real(k - 1, dp)/np
            hi(k) = vmin/vth + (vmax - vmin)/vth*real(k, dp)/np
            call kronrod_panel(lo(k), hi(k), val(:, k), err(:, k))
        end do
        do
            total = sum(val(:, 1:np), dim=2)
            scale = max(vtol*abs(total), sample%atol)
            if (sample%coarse .or. all(sum(err(:, 1:np), dim=2) <= scale)) exit
            ! A vanishing integral, such as the torque when it is not
            ! computed, has no tolerance and does not steer refinement.
            weight = 0.0_dp
            do k = 1, 3
                if (scale(k) > 0.0_dp) weight(k) = vtol/scale(k)
            end do
            if (np + 1 > vtol_max_panels) then
                call warning(fmt_dbg('adaptive velocity integral: panel limit, mth=', &
                    dble(mth), ' relative error=', &
                    maxval(sum(err(:, 1:np), dim=2)*weight)))
                exit
            end if
            worst = maxloc(err(1, 1:np)*weight(1) + err(2, 1:np)*weight(2) &
                + err(3, 1:np)*weight(3), dim=1)
            a = lo(worst)
            b = hi(worst)
            hi(worst) = 0.5_dp*(a + b)
            np = np + 1
            lo(np) = 0.5_dp*(a + b)
            hi(np) = b
            call kronrod_panel(lo(worst), hi(worst), val(:, worst), err(:, worst))
            call kronrod_panel(lo(np), hi(np), val(:, np), err(:, np))
        end do
        sample%nterm = 1
        sample%term(:, 1) = total
        ! The ledger gets the roots of the final panels with their quadrature
        ! weights, so that its rows sum to the native integral.
        block
            character(len=4096) :: ledger_path
            call resonance_ledger_path(ledger_path)
            if (len_trim(ledger_path) > 0 .and. .not. sample%coarse) then
                do k = 1, np
                    call record_panel(lo(k), hi(k))
                end do
            end if
        end block
    end subroutine compute_adaptive_sample

    subroutine record_panel(a, b)
        ! Ledger rows of the 15 Kronrod nodes on [a, b].
        real(dp), intent(in) :: a, b
        real(dp) :: c, h, x(15), w(15), wg(15), terms(3, nlev)
        integer :: k, nterm

        call kronrod_rule(x, w, wg)
        c = 0.5_dp*(a + b)
        h = 0.5_dp*(b - a)
        do k = 1, 15
            if (c + h*x(k) <= 0.0_dp) cycle
            call collisionless_point(c + h*x(k), h*w(k)*maxwellian_speed_density(c + h*x(k)), &
                terms, nterm)
        end do
    end subroutine record_panel

    subroutine compute_pitch_sample(vmin, vmax, sample, ok)
        ! The velocity integral of one collisionless component, taken along
        ! the resonance line in pitch. At fixed eta the resonance
        ! R = mph Om_ph + mth Om_th is a quadratic in u = v/vth, so
        ! int du sum_eta-roots F/|R_eta| = int deta sum_u-roots F/|R_u| with
        ! closed-form roots. On segments between the double roots and the
        ! points where a root crosses vmin or vmax the integrand is smooth up
        ! to 1/sqrt at the double roots, which eta = m - h cos(theta) removes,
        ! so adaptive Gauss-Kronrod in theta converges fast. ok = .false. when
        ! R is not quadratic in u (e.g. supban); the caller then integrates in u.
        real(dp), intent(in) :: vmin, vmax
        type(velocity_sample_t), intent(inout) :: sample
        logical, intent(out) :: ok
        integer, parameter :: nscan = 1024, nend = 24, max_panels = 2000, &
            initial_panels = 4
        real(dp), allocatable :: grid(:), fn(:, :), brk(:), mid(:), half(:), &
            lo(:), hi(:), val(:, :), err(:, :)
        integer, allocatable :: seg(:)
        real(dp) :: ulo, uhi, coef(3), total(3), scale(3), weight(3), a, b, t, f1(3)
        integer :: ng, nb, ns, np, k, j, worst
        logical :: fit_ok

        ok = .false.
        ulo = vmin/vth
        uhi = vmax/vth
        ! Scan in a parameter clustered at both ends, where the separatrix
        ! and the pitch limits concentrate the structure.
        ng = nscan - 1 + 2*nend
        allocate (grid(ng), fn(3, ng))
        do k = 1, nscan - 1
            grid(k) = 0.5_dp*(1.0_dp - cos(pi*real(k, dp)/nscan))
        end do
        do k = 1, nend
            t = 0.5_dp*(1.0_dp - cos(pi/nscan))*10.0_dp**(-10.0_dp*real(k, dp)/nend)
            grid(nscan - 1 + k) = t
            grid(nscan - 1 + nend + k) = 1.0_dp - t
        end do
        call sort_values(grid)
        grid = etamin + (etamax - etamin)*grid
        do k = 1, ng
            call resonance_quadratic(grid(k), coef, fit_ok)
            if (.not. fit_ok) return
            fn(:, k) = event_functions(coef, ulo, uhi)
        end do
        ! Segment ends: the pitch limits and the sign changes of the double-root
        ! discriminant and of R at ulo and uhi, refined by bisection.
        allocate (brk(2 + 3*ng))
        nb = 1
        brk(1) = etamin
        do k = 1, ng - 1
            do j = 1, 3
                if (fn(j, k)*fn(j, k + 1) < 0.0_dp) then
                    nb = nb + 1
                    brk(nb) = bisect_event(j, grid(k), grid(k + 1), fn(j, k), ulo, uhi)
                end if
            end do
        end do
        nb = nb + 1
        brk(nb) = etamax
        call sort_values(brk(1:nb))
        ns = nb - 1
        allocate (mid(ns), half(ns), lo(max_panels), hi(max_panels), &
            val(3, max_panels), err(3, max_panels), seg(max_panels))
        np = 0
        do k = 1, ns
            mid(k) = 0.5_dp*(brk(k) + brk(k + 1))
            half(k) = 0.5_dp*(brk(k + 1) - brk(k))
            if (half(k) <= 0.0_dp) cycle
            ! Several initial panels, so that one Gauss-Kronrod pair cannot
            ! miss structure that both rules sample alike.
            do j = 1, initial_panels
                if (np >= max_panels) exit
                np = np + 1
                seg(np) = k
                lo(np) = pi*real(j - 1, dp)/initial_panels
                hi(np) = pi*real(j, dp)/initial_panels
                call pitch_panel(mid(k), half(k), lo(np), hi(np), ulo, uhi, &
                    val(:, np), err(:, np))
            end do
        end do
        do
            total = sum(val(:, 1:np), dim=2)
            scale = max(vtol*abs(total), sample%atol)
            if (sample%coarse .or. all(sum(err(:, 1:np), dim=2) <= scale)) exit
            weight = 0.0_dp
            do k = 1, 3
                if (scale(k) > 0.0_dp) weight(k) = vtol/scale(k)
            end do
            if (np + 1 > max_panels) then
                call warning(fmt_dbg('pitch resonance integral: panel limit, mth=', &
                    dble(mth), ' relative error=', &
                    maxval(sum(err(:, 1:np), dim=2)*weight)))
                exit
            end if
            worst = maxloc(err(1, 1:np)*weight(1) + err(2, 1:np)*weight(2) &
                + err(3, 1:np)*weight(3), dim=1)
            a = lo(worst)
            b = hi(worst)
            hi(worst) = 0.5_dp*(a + b)
            np = np + 1
            seg(np) = seg(worst)
            lo(np) = 0.5_dp*(a + b)
            hi(np) = b
            k = seg(worst)
            call pitch_panel(mid(k), half(k), lo(worst), hi(worst), ulo, uhi, &
                val(:, worst), err(:, worst))
            call pitch_panel(mid(k), half(k), lo(np), hi(np), ulo, uhi, val(:, np), &
                err(:, np))
        end do
        f1 = 0.0_dp
        if (np > 0) f1 = sum(val(:, 1:np), dim=2)
        sample%nterm = 1
        sample%term(:, 1) = f1
        ok = .true.
    end subroutine compute_pitch_sample

    logical function ledger_active()
        ! The pitch integral keeps no per-root ledger rows.
        character(len=4096) :: ledger_path

        call resonance_ledger_path(ledger_path)
        ledger_active = len_trim(ledger_path) > 0
    end function ledger_active

    subroutine resonance_quadratic(eta, coef, ok)
        ! R(u vth, eta) = coef(1) + coef(2) u + coef(3) u^2, fitted at u = 1, 2, 3
        ! and checked at u = 4 and 1/2.
        real(dp), intent(in) :: eta
        real(dp), intent(out) :: coef(3)
        logical, intent(out) :: ok
        real(dp) :: r(5), dr, size_scale
        real(dp), parameter :: u(5) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 0.5_dp]
        integer :: k

        do k = 1, 5
            call resonance_value(u(k)*vth, eta, r(k), dr)
        end do
        coef(3) = 0.5_dp*(r(3) - 2.0_dp*r(2) + r(1))
        coef(2) = r(2) - r(1) - 3.0_dp*coef(3)
        coef(1) = r(1) - coef(2) - coef(3)
        size_scale = abs(coef(1)) + 4.0_dp*abs(coef(2)) + 16.0_dp*abs(coef(3))
        ok = all(ieee_is_finite(r)) .and. size_scale > 0.0_dp
        if (.not. ok) return
        ok = abs(coef(1) + 4.0_dp*coef(2) + 16.0_dp*coef(3) - r(4)) <= 1.0e-9_dp*size_scale &
            .and. abs(coef(1) + 0.5_dp*coef(2) + 0.25_dp*coef(3) - r(5)) <= 1.0e-9_dp*size_scale
    end subroutine resonance_quadratic

    pure function event_functions(coef, ulo, uhi) result(f)
        ! Discriminant of R(u) = 0 and R at the speed limits: the number of
        ! roots in (ulo, uhi) changes only where one of them changes sign.
        real(dp), intent(in) :: coef(3), ulo, uhi
        real(dp) :: f(3)

        f(1) = coef(2)**2 - 4.0_dp*coef(1)*coef(3)
        f(2) = coef(1) + coef(2)*ulo + coef(3)*ulo**2
        f(3) = coef(1) + coef(2)*uhi + coef(3)*uhi**2
    end function event_functions

    function bisect_event(j, a0, b0, fa0, ulo, uhi) result(x)
        integer, intent(in) :: j
        real(dp), intent(in) :: a0, b0, fa0, ulo, uhi
        real(dp) :: x, a, b, fa, fx(3), coef(3)
        integer :: it
        logical :: fit_ok

        a = a0
        b = b0
        fa = fa0
        do it = 1, 100
            x = 0.5_dp*(a + b)
            if (x <= a .or. x >= b) exit
            call resonance_quadratic(x, coef, fit_ok)
            fx = event_functions(coef, ulo, uhi)
            if (fx(j) == 0.0_dp) exit
            if (fa*fx(j) < 0.0_dp) then
                b = x
            else
                a = x
                fa = fx(j)
            end if
        end do
    end function bisect_event

    subroutine pitch_panel(m, h, a, b, ulo, uhi, val, err)
        ! 15-point Kronrod integral in theta on [a, b] of the pitch integrand on
        ! the segment eta = m - h cos(theta), and its difference to the embedded
        ! 7-point Gauss rule.
        real(dp), intent(in) :: m, h, a, b, ulo, uhi
        real(dp), intent(out) :: val(3), err(3)
        real(dp) :: x(15), w(15), wg(15), th, f(3)
        integer :: k

        call kronrod_rule(x, w, wg)
        val = 0.0_dp
        err = 0.0_dp
        do k = 1, 15
            th = 0.5_dp*(a + b) + 0.5_dp*(b - a)*x(k)
            f = pitch_sum(m - h*cos(th), ulo, uhi)*h*sin(th)
            val = val + w(k)*f
            err = err + wg(k)*f
        end do
        val = 0.5_dp*(b - a)*val
        err = abs(val - 0.5_dp*(b - a)*err)
    end subroutine pitch_panel

    function pitch_sum(eta, ulo, uhi) result(f)
        ! Integrand in eta: the resonant speeds at this pitch, each weighted
        ! by the Maxwellian speed density over |dR/du|.
        real(dp), intent(in) :: eta, ulo, uhi
        real(dp) :: f(3), coef(3), disc, qr, u(2), jac
        integer :: nu, k
        logical :: fit_ok

        f = 0.0_dp
        call resonance_quadratic(eta, coef, fit_ok)
        disc = coef(2)**2 - 4.0_dp*coef(1)*coef(3)
        if (disc < 0.0_dp) return
        nu = 0
        if (coef(3) == 0.0_dp) then
            if (coef(2) == 0.0_dp) return
            nu = 1
            u(1) = -coef(1)/coef(2)
        else
            qr = -0.5_dp*(coef(2) + sign(sqrt(disc), coef(2)))
            nu = 1
            u(1) = qr/coef(3)
            if (qr /= 0.0_dp) then
                nu = 2
                u(2) = coef(1)/qr
            end if
        end if
        do k = 1, nu
            if (u(k) <= ulo .or. u(k) >= uhi) cycle
            jac = abs(coef(2) + 2.0_dp*coef(3)*u(k))
            if (jac <= 0.0_dp) cycle
            f = f + pitch_root_terms(u(k), eta, jac)
        end do
    end function pitch_sum

    function pitch_root_terms(ux, eta, jac) result(term)
        ! D11, D12 and torque integrands of one resonant root over |dR/du|,
        ! with the same orbit physics as collisionless_point.
        real(dp), intent(in) :: ux, eta, jac
        real(dp) :: term(3), v, taub, bounceavg(nvar), Hmn2, attenuation_factor, &
            wt, layer_factor, res, dres_deta
        integer :: istate_dv

        v = ux*vth
        call Om_th(v, eta, Omth, dOmthdv, dOmthdeta)
        taub = 2.0_dp*pi/abs(Omth)
        call transport_bounce(v, eta, taub, bounceavg, istate_dv)
        if (istate_dv == -1) then
            call error(fmt_dbg('VODE MXSTEP: mth=', dble(mth), ' ux=', ux, ' eta=', eta, &
                ' taub=', taub))
        else if (istate_dv /= 2) then
            call warning(fmt_dbg('dvode istate=', dble(istate_dv), ' at mth=', dble(mth), &
                ' ux=', ux, ' eta=', eta))
        end if
        Hmn2 = (bounceavg(3)**2 + bounceavg(4)**2)*(mi*(ux*vth)**2/2.0_dp)**2
        attenuation_factor = nonlinear_attenuation(ux, eta, bounceavg, Omth, &
            dOmthdv, dOmthdeta, Hmn2)
        wt = maxwellian_speed_density(ux)/jac
        term(1) = wt*D11mono(ux, taub, Hmn2)
        term(2) = wt*D12mono(ux, taub, Hmn2)
        term(3) = 0.0_dp
        if (comptorque) term(3) = wt*Tphi_mono(ux, taub, Hmn2)
        if (collisional_layer .and. eta > etatp) then
            call resonance_value(v, eta, res, dres_deta)
            layer_factor = collisional_layer_factor(deflection_frequency(v), eta, eps, &
                etatp, etadt, abs(dres_deta))
            term = term*layer_factor
        end if
        term = term*attenuation_factor
    end function pitch_root_terms

    pure subroutine sort_values(x)
        real(dp), intent(inout) :: x(:)
        real(dp) :: t
        integer :: i, j

        do i = 2, size(x)
            t = x(i)
            j = i - 1
            do while (j >= 1)
                if (x(j) <= t) exit
                x(j + 1) = x(j)
                j = j - 1
            end do
            x(j + 1) = t
        end do
    end subroutine sort_values

    subroutine kronrod_panel(a, b, val, err)
        ! 15-point Kronrod integral of the root sums on [a, b] and its
        ! difference to the embedded 7-point Gauss rule, per quantity.
        real(dp), intent(in) :: a, b
        real(dp), intent(out) :: val(3), err(3)
        real(dp) :: x(15), w(15), wg(15), f(3), gauss(3)
        integer :: k

        call kronrod_rule(x, w, wg)
        val = 0.0_dp
        gauss = 0.0_dp
        do k = 1, 15
            f = root_sum(0.5_dp*(a + b) + 0.5_dp*(b - a)*x(k))
            val = val + w(k)*f
            gauss = gauss + wg(k)*f
        end do
        val = 0.5_dp*(b - a)*val
        err = abs(val - 0.5_dp*(b - a)*gauss)
    end subroutine kronrod_panel

    pure subroutine kronrod_rule(x, w, wg)
        ! Nodes of the 15-point Kronrod rule on [-1, 1], its weights, and the
        ! weights of the embedded 7-point Gauss rule (zero at Kronrod-only nodes).
        real(dp), intent(out) :: x(15), w(15), wg(15)
        real(dp), parameter :: xk(8) = [0.991455371120812639206854697526329_dp, &
            0.949107912342758524526189684047851_dp, 0.864864423359769072789712788640926_dp, &
            0.741531185599394439863864773280788_dp, 0.586087235467691130294144845693013_dp, &
            0.405845151377397166906606412076961_dp, 0.207784955007898467600689403773245_dp, &
            0.0_dp]
        real(dp), parameter :: wk(8) = [0.022935322010529224963732008058970_dp, &
            0.063092092629978553290700663189204_dp, 0.104790010322250183839876322541518_dp, &
            0.140653259715525918745189590510238_dp, 0.169004726639267902826583426598550_dp, &
            0.190350578064785409913256402421014_dp, 0.204432940075298892414161999234649_dp, &
            0.209482141084727828012999174891714_dp]
        real(dp), parameter :: wgp(4) = [0.129484966168869693270611432679082_dp, &
            0.279705391489276667901467771423780_dp, 0.381830050505118944950369775488975_dp, &
            0.417959183673469387755102040816327_dp]
        integer :: k

        wg = 0.0_dp
        do k = 1, 7
            x(k) = -xk(k)
            x(16 - k) = xk(k)
            w(k) = wk(k)
            w(16 - k) = wk(k)
            if (mod(k, 2) == 0) then
                wg(k) = wgp(k/2)
                wg(16 - k) = wgp(k/2)
            end if
        end do
        x(8) = 0.0_dp
        w(8) = wk(8)
        wg(8) = wgp(4)
    end subroutine kronrod_rule

    function root_sum(ux) result(f)
        ! D11, D12 and torque integrands at speed ux, summed over the roots,
        ! with the Maxwellian speed density.
        real(dp), intent(in) :: ux
        real(dp) :: f(3), terms(3, nlev)
        integer :: nterm

        f = 0.0_dp
        if (ux <= 0.0_dp) return
        call collisionless_point(ux, maxwellian_speed_density(ux), terms, nterm, &
            record=.false.)
        f = sum(terms(:, 1:nterm), dim=2)
    end function root_sum

    subroutine collisionless_point(ux, sample_weight, terms, nterm, record)
        ! Contributions of the resonant roots at normalised speed ux, each
        ! weighted by sample_weight: one term per root. record = .false.
        ! keeps the evaluation out of the resonance ledger.
        real(dp), intent(in) :: ux, sample_weight
        real(dp), intent(out) :: terms(3, nlev)
        integer, intent(out) :: nterm
        logical, intent(in), optional :: record
        real(dp) :: dD11, dD12, dT, v, eta
        real(dp) :: eta_res(2)
        real(dp) :: taub, bounceavg(nvar)
        integer :: istate_dv
        real(dp) :: Hmn2, attenuation_factor, layer_factor, layer_nu
        real(dp) :: cutout_nu, cutout_resonance_weight
        logical :: cutout_inside_layer
        real(dp) :: roots(nlev, 3)
        integer :: nroots, kr
        character(len=4096) :: ledger_path
        character(len=4) :: orbit_branch
        character(len=1024), allocatable :: ledger_rows(:)
        integer :: ledger_count
        logical :: ledger_enabled
        real(dp) :: Omph, dOmphdv, dOmphdeta, resonance_residual

        call resonance_ledger_path(ledger_path)
        ledger_enabled = len_trim(ledger_path) > 0
        if (present(record)) ledger_enabled = ledger_enabled .and. record
        ledger_count = 0
        if (ledger_enabled) then
            allocate(ledger_rows(nlev))
            if (etamax < etatp) then
                if (sign_vpar > 0.0_dp) then
                    orbit_branch = 'cop'
                else
                    orbit_branch = 'ctr'
                end if
            else
                orbit_branch = 'trap'
            end if
        end if

        nterm = 0
        v = ux * vth
        if (collisional_layer .and. etamax > etatp) &
            layer_nu = deflection_frequency(v)
        if (resonance_cutout /= CUTOUT_OFF) cutout_nu = deflection_frequency(v)
        call driftorbit_coarse(v, etamin, etamax, roots, nroots)
        do kr = 1, nroots
            eta_res = driftorbit_root(v, 1.0e-8_dp * abs(Om_tE), roots(kr, 1), roots(kr, 2))
            if (eta_res(1) < 0.0_dp) cycle  ! bracket-failure sentinel
            eta = eta_res(1)

            call Om_th(v, eta, Omth, dOmthdv, dOmthdeta)

            taub = 2.0_dp * pi / abs(Omth)
            call transport_bounce(v, eta, taub, bounceavg, istate_dv)
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

            cutout_inside_layer = .false.
            if (resonance_cutout /= CUTOUT_OFF) &
                call resonance_cutout_weight(v, eta, abs(eta_res(2)), cutout_nu, &
                                            cutout_resonance_weight, cutout_inside_layer)
            if (cutout_inside_layer) then
                dD11 = sample_weight*D11mono(ux, taub, Hmn2)*cutout_resonance_weight
                dD12 = sample_weight*D12mono(ux, taub, Hmn2)*cutout_resonance_weight
            else
                dD11 = sample_weight*D11mono(ux, taub, Hmn2)/abs(eta_res(2))
                dD12 = sample_weight*D12mono(ux, taub, Hmn2)/abs(eta_res(2))
            end if
            if (collisional_layer) then
                layer_factor = 1.0_dp
                if (eta > etatp) layer_factor = collisional_layer_factor( &
                    layer_nu, eta, eps, etatp, etadt, abs(eta_res(2)))
                dD11 = dD11 * layer_factor
                dD12 = dD12 * layer_factor
            end if
            nterm = nterm + 1
            terms(1, nterm) = dD11 * attenuation_factor
            terms(2, nterm) = dD12 * attenuation_factor
            terms(3, nterm) = 0.0_dp
            if (comptorque) then
                if (cutout_inside_layer) then
                    dT = sample_weight*Tphi_mono(ux, taub, Hmn2)*cutout_resonance_weight
                else
                    dT = sample_weight*Tphi_mono(ux, taub, Hmn2)/abs(eta_res(2))
                end if
                if (collisional_layer) dT = dT * layer_factor
                terms(3, nterm) = dT * attenuation_factor
            end if
            if (ledger_enabled) then
                call Om_ph(v, eta, Omph, dOmphdv, dOmphdeta)
                resonance_residual = real(mth, dp) * Omth + real(mph, dp) * Omph
                ledger_count = ledger_count + 1
                write(ledger_rows(ledger_count), *) &
                    s, selected_species, mth, mph, trim(orbit_branch), ux, eta, &
                    Omth, Omph, resonance_residual, eta_res(2), taub, Hmn2, &
                    attenuation_factor, sample_weight, dT*attenuation_factor, istate_dv
            end if
        end do
        if (ledger_enabled) &
            call append_resonance_ledger(ledger_path, ledger_rows(:ledger_count))
    end subroutine collisionless_point

    subroutine resonance_cutout_weight(v, eta, abs_R_eta, nu_d, weight, inside_layer)
        ! Finite replacement of 1/|R_eta| (doc/resonance_cutout.md).
        real(dp), intent(in) :: v, eta, abs_R_eta, nu_d
        real(dp), intent(out) :: weight
        logical, intent(out) :: inside_layer

        call cutout_weight(resonance_cutout, resonance_cutout_width, abs_R_eta, &
            resonance_curvature(v, eta), bounce_pitch_diffusion(nu_d, eta, etatp, etadt), &
            weight, inside_layer)
    end subroutine resonance_cutout_weight

    function resonance_curvature(v, eta) result(abs_R_eta_eta)
        ! |d2R/deta2| by central differences of the analytic dR/deta. The
        ! step stays inside the orbit-class interval and below a quarter of
        ! the distance to its ends, so it resolves the logarithmic curvature
        ! next to eta_tp and never crosses it. A root on an interval end gets
        ! zero (no layer); fold_layer_width maps a nonfinite curvature to no
        ! layer as well.
        real(dp), parameter :: RELATIVE_STEP = 1.0e-4_dp, END_FRACTION = 0.25_dp
        real(dp), intent(in) :: v, eta
        real(dp) :: abs_R_eta_eta
        real(dp) :: h, res, dlo, dhi

        abs_R_eta_eta = 0.0_dp
        h = min(RELATIVE_STEP*(etamax - etamin), &
                END_FRACTION*min(eta - etamin, etamax - eta))
        if (.not. h > 0.0_dp) return
        call resonance_value(v, eta - h, res, dlo)
        call resonance_value(v, eta + h, res, dhi)
        call cutout_curvature(dlo, dhi, h, abs_R_eta_eta)
        abs_R_eta_eta = abs(abs_R_eta_eta)
    end function resonance_curvature

    subroutine resonance_ledger_path(path)
        character(len=*), intent(out) :: path
        integer :: env_length, env_status
        path = ''
        call get_environment_variable('NEORT_RESONANCE_LEDGER', path, &
            length=env_length, status=env_status)
        if (env_status == -1) error stop 'NEORT_RESONANCE_LEDGER path is too long'
        if (env_status /= 0 .or. env_length == 0) then
            path = ''
            return
        end if
        if (.not. comptorque) error stop 'NEORT_RESONANCE_LEDGER requires comptorque=.true.'
        if (finite_krook_active()) &
            error stop 'NEORT_RESONANCE_LEDGER requires collisionless transport'
    end subroutine resonance_ledger_path

    subroutine append_resonance_ledger(path, rows)
        character(len=*), intent(in) :: path, rows(:)
        integer :: unit, k, file_size
        logical :: exists
        ! Each worker buffers one speed sample. Separate append units must never
        ! overlap: locking just WRITE would leave independent file positions stale.
        !$omp critical(neort_resonance_ledger_io)
        inquire(file=trim(path), exist=exists, size=file_size)
        open(newunit=unit, file=trim(path), status='unknown', &
            position='append', action='write')
        if (.not. exists .or. file_size == 0) then
            write(unit, '(A)') '# s species mth mph branch ux eta Omth Omph residual ' &
                //'dres_deta taub Hmn2 attenuation sample_weight dTds istate'
        end if
        do k = 1, size(rows)
            write(unit, '(A)') trim(rows(k))
        end do
        close(unit)
        !$omp end critical(neort_resonance_ledger_io)
    end subroutine append_resonance_ledger

    subroutine compute_krook_sample(ku, vmin, vmax, vsteps, energy_ratio, sample)
        ! Node ku: its finite-width pitch integral as the single term.
        integer, intent(in) :: ku
        real(dp), intent(in) :: vmin, vmax
        integer, intent(in) :: vsteps
        real(dp), intent(in) :: energy_ratio
        type(velocity_sample_t), intent(inout) :: sample
        real(dp) :: ux, values(3, 1)

        ux = energy_sample_speed(energy_ratio, ku, vmin/vth, vmax/vth, vsteps)
        call compute_krook_velocity(ux, values(:, 1), energy_ratio)
        call store_sample_terms(values, sample)
    end subroutine compute_krook_sample

    subroutine compute_krook_velocity(ux, values, energy_ratio)
        real(dp), intent(in) :: ux
        real(dp), intent(in), optional :: energy_ratio
        real(dp), intent(out) :: values(3)
        real(dp) :: effective_harmonic, geometry, nu
        real(dp) :: natural(3), cross(3), absolute_tolerance(3), amplitude_scale, amplitude
        real(dp) :: amplitude_error
        integer :: status
        real(dp), parameter :: cross_rtol = 1.0e-2_dp
        effective_harmonic = real(mth, dp)
        geometry = krook_geometry_trapped
        if (etamax < etatp) then
            geometry = krook_geometry_passing
            if (krook_harmonic == 1) effective_harmonic = effective_harmonic + q*mph
        end if
        nu = effective_krook_rate(ux, effective_harmonic, geometry, ux*vth*0.01_dp)
        if (.not. fow_radial) then
            ! A nearly vanishing Fourier harmonic cannot be resolved relative
            ! to itself by an ODE controlled relative to the unaveraged drive.
            ! Bound that drive by the Fourier triangle inequality, and carry
            ! its orbit error e*S into the squared weight, which errs by at
            ! most 2e*S*|H| + e^2*S^2. The cross term is integrated directly:
            ! Young's bound e^2*S^2/rtol can exceed the value itself where the
            ! harmonic is resolved, and would hide quadrature error there.
            ! The phase sensitivity is 2*pi*|mth+nq| times the period error
            ! on passing orbits. Trapped orbits carry the secular phase
            ! mth*Omth*t and the bounded nq*theta, hence 2*pi*(|mth|+|nq|);
            ! their high harmonics fall to this noise floor in the deep well.
            if (pertfile) then
                call magnetic_amplitude_bound([s, 0.0_dp, th0], amplitude)
                amplitude_scale = abs(pertfile_scale)*amplitude/Bmin
            else
                amplitude_scale = abs(epsmn)
            end if
            if (etamax < etatp) then
                amplitude_error = (1.0_dp + &
                    2.0_dp*pi*abs(real(mth, dp) + q*mph))*bounce_rtol
            else
                amplitude_error = (1.0_dp + &
                    2.0_dp*pi*(abs(real(mth, dp)) + abs(q*mph)))*bounce_rtol
            end if
            call integrate_krook(residual, envelope, etamin, etamax, nu, &
                krook_eta_steps, krook_rtol, 0.0_dp, natural, status, krook_max_depth)
            if (status /= 0) call report_krook_failure(status, ux, nu, &
                effective_harmonic, geometry)
            ! An error estimate needs only two digits; its own floor is e*S.
            call integrate_krook(residual, cross_weight, etamin, etamax, nu, &
                krook_eta_steps, cross_rtol, 0.0_dp, cross, status, krook_max_depth, &
                amplitude_error*abs(natural))
            if (status /= 0) call report_krook_failure(status, ux, nu, &
                effective_harmonic, geometry)
            absolute_tolerance = 2.0_dp*amplitude_error*((1.0_dp + cross_rtol)* &
                abs(cross) + amplitude_error*abs(natural)) + &
                amplitude_error**2*abs(natural)
            call integrate_krook(residual, weight, etamin, etamax, nu, krook_eta_steps, &
                krook_rtol, krook_atol, values, status, krook_max_depth, &
                absolute_tolerance)
        else
            call integrate_krook(residual, weight, etamin, etamax, nu, krook_eta_steps, &
                krook_rtol, krook_atol, values, status, krook_max_depth)
        end if
        if (status /= 0) call report_krook_failure(status, ux, nu, &
            effective_harmonic, geometry)
    contains
        subroutine residual(eta, res, derivative)
            real(dp), intent(in) :: eta
            real(dp), intent(out) :: res, derivative
            call resonance_value(ux*vth, eta, res, derivative)
        end subroutine residual
        subroutine envelope(eta, result_)
            real(dp), intent(in) :: eta
            real(dp), intent(out) :: result_(:)
            real(dp) :: taub, drive_bound, energy_squared
            logical :: mono
            call Om_th(ux*vth, eta, Omth, dOmthdv, dOmthdeta)
            taub = 2.0_dp*pi/abs(Omth)
            drive_bound = (2.0_dp - eta*Bmin)*amplitude_scale
            if (phi1) call perturbation_amplitude_bound( &
                [s, 0.0_dp, th0], ux*vth, eta, drive_bound)
            energy_squared = (mi*(ux*vth)**2/2.0_dp*drive_bound)**2
            mono = .false.
            if (present(energy_ratio)) mono = energy_ratio > 0.0_dp
            result_(3) = 0.0_dp
            if (mono) then
                result_(1) = D11mono(ux, taub, energy_squared)
                result_(2) = D12mono(ux, taub, energy_squared)
                if (comptorque) result_(3) = Tphi_mono(ux, taub, energy_squared)
            else
                result_(1) = D11int(ux, taub, energy_squared)
                result_(2) = D12int(ux, taub, energy_squared)
                if (comptorque) result_(3) = Tphi_int(ux, taub, energy_squared)
            end if
        end subroutine envelope
        subroutine weight(eta, result_)
            real(dp), intent(in) :: eta
            real(dp), intent(out) :: result_(:)
            call compute_transport_weight(ux, eta, result_, energy_ratio)
        end subroutine weight
        subroutine cross_weight(eta, result_)
            ! Prefactor*S*|H|: both factors carry the same transport prefactor.
            real(dp), intent(in) :: eta
            real(dp), intent(out) :: result_(:)
            real(dp) :: squared(size(result_)), bound(size(result_))
            call compute_transport_weight(ux, eta, squared, energy_ratio)
            call envelope(eta, bound)
            result_ = sqrt(abs(squared*bound))
        end subroutine cross_weight
    end subroutine compute_krook_velocity

    subroutine perturbation_amplitude_bound(x, v, eta, amplitude)
        real(dp), intent(in) :: x(3), v, eta
        real(dp), intent(out) :: amplitude
        real(dp) :: magnetic, electric
        magnetic = abs(epsmn)
        if (pertfile) then
            call magnetic_amplitude_bound(x, magnetic)
            magnetic = abs(pertfile_scale)*magnetic/Bmin
        end if
        amplitude = (2.0_dp - eta*Bmin)*magnetic
        if (.not. phi1) return
        call phi1_amplitude_bound(x, electric)
        amplitude = amplitude + abs(qi)*electric/(mi*v**2/2.0_dp)
        ! Round outward so the bound also covers a rounded complex phase/norm.
        if (amplitude > 0.0_dp) amplitude = nearest(amplitude, 1.0_dp)
    end subroutine perturbation_amplitude_bound

    subroutine report_krook_failure(status, ux, nu, effective_harmonic, geometry)
        integer, intent(in) :: status
        real(dp), intent(in) :: ux, nu, effective_harmonic, geometry
        ! Preserve the returned status and the consuming thread's physical state.
        ! The eta values are integration bounds, not a guessed failing node.
        !$omp critical(krook_failure_context)
        write(error_unit, '(a,6(a,i0),*(a,es25.17e3))') &
            'KROOK_QUADRATURE_FAILURE', ' status=', status, &
            ' species_index=', selected_species, ' mth=', mth, ' mph=', mph, &
            ' eta_steps=', krook_eta_steps, ' max_depth=', krook_max_depth, &
            ' s=', s, ' mass_g=', mi, ' charge_statc=', qi, &
            ' density_cm-3=', ns, ' omega_E_s-1=', Om_tE, ' ux=', ux, &
            ' vth_cm_s=', vth, ' nu_s-1=', nu, ' eta_min=', etamin, &
            ' eta_max=', etamax, ' eta_tp=', etatp, &
            ' sign_vpar=', sign_vpar, &
            ' effective_harmonic=', effective_harmonic, ' geometry=', geometry, &
            ' rtol=', krook_rtol, ' atol=', krook_atol
        flush(error_unit)
        !$omp end critical(krook_failure_context)
        error stop 'finite-width Krook quadrature did not converge'
    end subroutine report_krook_failure

    subroutine compute_transport_weight(ux, eta, values, energy_ratio)
        ! With energy_ratio > 0 the weights omit the Maxwellian speed density
        ! (monoenergetic mode, sample weight 1); otherwise thermal weights.
        real(dp), intent(in) :: ux, eta
        real(dp), intent(in), optional :: energy_ratio
        logical :: mono
        real(dp), intent(out) :: values(:)
        real(dp) :: taub, bounceavg(nvar), Hmn2
        integer :: istate
        call Om_th(ux*vth, eta, Omth, dOmthdv, dOmthdeta)
        taub = 2.0_dp*pi/abs(Omth)
        call transport_bounce(ux*vth, eta, taub, bounceavg, istate)
        if (istate /= 2) error stop 'orbit integration failed in Krook quadrature'
        Hmn2 = (bounceavg(3)**2 + bounceavg(4)**2)*(mi*(ux*vth)**2/2.0_dp)**2
        mono = .false.
        if (present(energy_ratio)) mono = energy_ratio > 0.0_dp
        values(3) = 0.0_dp
        if (mono) then
            values(1) = D11mono(ux, taub, Hmn2)
            values(2) = D12mono(ux, taub, Hmn2)
            if (comptorque) values(3) = Tphi_mono(ux, taub, Hmn2)
        else
            values(1) = D11int(ux, taub, Hmn2)
            values(2) = D12int(ux, taub, Hmn2)
            if (comptorque) values(3) = Tphi_int(ux, taub, Hmn2)
        end if
    end subroutine compute_transport_weight

    subroutine transport_bounce(v, eta, taub, bounceavg, istate)
        ! Bounce averages of the transport weight, H_mn**2 = |<H_n>|**2.
        ! int H_n dt, with the dimensionless H_n = drive*exp(i*phase) (transport_rates),
        ! is an oscillatory, often strongly cancelling integral of size up to
        ! taub*max|H_n|. The fixed bounce_atol = 1e-10 can exceed that scale
        ! (small eps_n, short taub) and then leaves these components without
        ! error control: H_mn**2 jumped by O(1) between neighbouring eta as the
        ! DVODE step sequence changed. Their absolute tolerance is therefore
        ! bounce_rtol relative to taub*max|H_n|, which makes the weight (to rounding)
        ! quadratic in the perturbation amplitude.
        integer, parameter :: samples = 8
        real(dp), intent(in) :: v, eta, taub
        real(dp), intent(out) :: bounceavg(nvar)
        integer, intent(out) :: istate
        real(dp) :: atol(nvar), y(nvar), ydot(nvar), scale
        integer :: k
        atol = bounce_atol
        scale = 0.0_dp
        y = 0.0_dp
        do k = 0, samples - 1
            y(1) = th0 + 2.0_dp*pi*real(k, dp)/real(samples, dp)
            call timestep_transport(v, eta, nvar, 0.0_dp, y, ydot)
            scale = max(scale, abs(cmplx(ydot(3), ydot(4), dp)))
        end do
        scale = bounce_rtol*taub*scale
        if (scale > 0.0_dp) atol(3:4) = min(bounce_atol, scale)
        if (eta > etatp) then
            call trapped_bounce(v, eta, taub, bounceavg, istate, atol)
        else
            call passing_bounce(v, eta, taub, bounceavg, istate, atol)
        end if
    end subroutine transport_bounce

    subroutine passing_bounce(v, eta, taub, bounceavg, istate, atol_in)
        ! Bounce averages of a passing orbit over its closed transit.
        ! The supplied taub and the module Omth come from the interpolated
        ! frequency Om_th, which differs from the orbit period by the spline
        ! error. A time integration over [0, taub] with the phase
        ! -(mth + q*mph)*Omth*t then neither closes theta nor the phase, and
        ! the strongly cancelling harmonic picks up an error of that relative
        ! size times max|H_n|; for ITER C12 at m = -11 this changed its squared
        ! weight by a factor of about 200, and its noise in eta stopped
        ! the Krook pitch quadrature. Here the orbit is parametrised by the
        ! poloidal angle travelled, theta = th0 + direction*phi with
        ! dt/dphi = 1/|v_parallel*h^theta|, over exactly one transit. A first
        ! pass gives the period, whose frequency replaces Omth in the phase
        ! during the second pass, as in trapped_bounce. The two passes take
        ! different steps, and their periods differ by the ODE error (about
        ! 1e-11 for ITER C12 at s = 0.09). The phase then did not close on the
        ! second pass by (mth + q*mph)*2*pi times that, which for harmonics
        ! cancelling to |<H_n>| ~ 1e-7*max|H_n| made the weight noisy by 1e-4
        ! in eta and exhausted the Krook piece budget. The second pass
        ! therefore also integrates d/dOmth of the harmonic and moves it to
        ! the frequency of its own period, to first order in the difference.
        use fortnum_ode_vode, only: vode_state_t, vode_init, vode_integrate_to
        use fortnum_status, only: fortnum_status_t, FORTNUM_OK, &
            FORTNUM_CONVERGENCE_ERROR
        real(dp), intent(in) :: v, eta, taub, atol_in(nvar)
        real(dp), intent(out) :: bounceavg(nvar)
        integer, intent(out) :: istate
        integer, parameter :: neq = nvar + 2
        real(dp) :: direction, sign_vparallel, saved_omth, y0(neq), atol(neq)
        real(dp) :: position(3), bmod, sqrtg, hder(3), hcovar(3), hctrvr(3), hcurl(3)
        real(dp) :: harmonic, period_omth
        real(dp), allocatable :: yend(:)
        integer :: pass
        logical :: averages
        type(vode_state_t) :: vstate
        type(fortnum_status_t) :: status

        ! Same orientation as bounce_fast: v_parallel = sign(h^theta)*sign_vpar*|v_par|
        ! at th0, so theta moves in the direction of sign_vpar.
        position = [s, 0.0_dp, th0]
        call do_magfie(position, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
        sign_vparallel = sign(1.0_dp, hctrvr(3))*sign(1.0_dp, sign_vpar)
        direction = sign(1.0_dp, sign_vparallel*hctrvr(3))
        atol(:nvar) = atol_in
        atol(1) = bounce_rtol*taub
        atol(nvar + 1:) = atol(3:4)*taub
        harmonic = real(mth, dp) + q*mph
        bounceavg = 0.0_dp
        saved_omth = Omth
        do pass = 1, 2
            averages = pass == 2
            y0 = 0.0_dp
            call vode_init(vstate, neq, 0.0_dp, y0)
            call vode_integrate_to(phi_rhs, vstate, 2.0_dp*pi, bounce_rtol, atol, &
                                   yend, status)
            if (status%code == FORTNUM_OK) then
                istate = 2
            else if (status%code == FORTNUM_CONVERGENCE_ERROR) then
                istate = -1
            else
                istate = 0
            end if
            if (istate /= 2) exit
            period_omth = sign(2.0_dp*pi/yend(1), saved_omth)
            if (pass == 2) then
                yend(3) = yend(3) + (period_omth - Omth)*yend(nvar + 1)
                yend(4) = yend(4) + (period_omth - Omth)*yend(nvar + 2)
            end if
            Omth = period_omth
        end do
        Omth = saved_omth
        if (istate /= 2) return
        bounceavg(3:) = yend(3:nvar)/yend(1)

    contains

        subroutine phi_rhs(phi, y, dydphi, rhs_ctx)
            real(dp), intent(in) :: phi
            real(dp), intent(in) :: y(:)
            real(dp), intent(out) :: dydphi(:)
            class(*), intent(in), optional :: rhs_ctx
            real(dp) :: dtdphi, vparallel, orbit(nvar), rate(nvar), x(3)
            real(dp) :: bmod_, sqrtg_, hder_(3), hcovar_(3), hctrvr_(3), hcurl_(3)
            associate (dummy => rhs_ctx)
            end associate
            x = [s, 0.0_dp, th0 + direction*phi]
            call do_magfie(x, bmod_, sqrtg_, hder_, hcovar_, hctrvr_, hcurl_)
            vparallel = vpar(v, eta, bmod_)
            dtdphi = 1.0_dp/(vparallel*abs(hctrvr_(3)))
            dydphi = 0.0_dp
            dydphi(1) = dtdphi
            if (.not. averages) return
            orbit = 0.0_dp
            orbit(1) = x(3)
            orbit(2) = sign(vparallel, sign_vparallel)
            orbit(3:) = y(3:nvar)
            call transport_rates(v, eta, nvar, y(1), orbit, bmod_, hder_, hcovar_, &
                                 hctrvr_, rate)
            dydphi(3:nvar) = rate(3:)*dtdphi
            ! d/dOmth of the harmonic exp(-i*harmonic*Omth*t)*H_n.
            dydphi(nvar + 1) = harmonic*y(1)*rate(4)*dtdphi
            dydphi(nvar + 2) = -harmonic*y(1)*rate(3)*dtdphi
        end subroutine phi_rhs
    end subroutine passing_bounce

    subroutine trapped_bounce(v, eta, taub, bounceavg, istate, atol_in)
        ! Bounce averages of a trapped orbit over its closed orbit.
        ! Time integration of (theta, v_parallel) accumulates an energy error
        ! which near the trapped-passing boundary is amplified by
        ! 1/(eta*B_max - 1) through the time spent near the X-point; the weight
        ! then became noisy at the 1e-4 level and the Krook pitch quadrature
        ! did not converge. Here the orbit is parametrised by the angle chi,
        !   theta = theta_m - direction*Delta*cos(chi),
        ! between the turning points theta_m -+ Delta, and t(chi) is integrated
        ! together with the bounce integrands. v_parallel follows from the
        ! energy invariant, and dt/dchi is regular at the turning points.
        use fortnum_quadrature, only: gauss_legendre
        use fortnum_ode_vode, only: vode_state_t, vode_init, vode_integrate_to
        use fortnum_status, only: fortnum_status_t, FORTNUM_OK, &
            FORTNUM_CONVERGENCE_ERROR
        ! Below this 1 - eta*B is integrated from the turning point: the direct
        ! difference has relative rounding error eps/(1 - eta*B), which would
        ! break the ODE step control next to the turning points.
        real(dp), parameter :: gap_switch = 1.0e-6_dp
        real(dp), intent(in) :: v, eta, taub, atol_in(nvar)
        real(dp), intent(out) :: bounceavg(nvar)
        integer, intent(out) :: istate
        real(dp) :: lower, upper, middle, half, direction, sign_htheta, chi0, saved_omth
        integer :: pass
        logical :: averages
        real(dp) :: nodes(8), weights(8), y0(nvar), atol(nvar), position(3)
        real(dp) :: bmod, sqrtg, hder(3), hcovar(3), hctrvr(3), hcurl(3)
        real(dp), allocatable :: yend(:)
        type(vode_state_t) :: vstate
        type(fortnum_status_t) :: status

        call gauss_legendre(size(nodes), nodes, weights)
        direction = sign(1.0_dp, sign_vpar)
        position = [s, 0.0_dp, th0]
        call do_magfie(position, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
        sign_htheta = sign(1.0_dp, hctrvr(3))
        lower = turning_point(-1.0_dp)
        upper = turning_point(1.0_dp)
        middle = (lower + upper)/2.0_dp
        half = (upper - lower)/2.0_dp
        chi0 = acos(max(-1.0_dp, min(1.0_dp, direction*(middle - th0)/half)))
        atol = atol_in
        atol(1) = bounce_rtol*taub
        bounceavg = 0.0_dp
        ! The bounce harmonics exp(-i*mth*Omth*t) must close on the orbit:
        ! a first pass gives its period, which replaces the interpolated
        ! bounce frequency in the phase during the second pass.
        saved_omth = Omth
        do pass = 1, 2
            averages = pass == 2
            y0 = 0.0_dp
            call vode_init(vstate, nvar, chi0, y0)
            call vode_integrate_to(chi_rhs, vstate, chi0 + 2.0_dp*pi, bounce_rtol, &
                                   atol, yend, status)
            if (status%code == FORTNUM_OK) then
                istate = 2
            else if (status%code == FORTNUM_CONVERGENCE_ERROR) then
                istate = -1
            else
                istate = 0
            end if
            if (istate /= 2) exit
            Omth = sign(2.0_dp*pi/yend(1), saved_omth)
        end do
        Omth = saved_omth
        if (istate /= 2) return
        bounceavg(3:) = yend(3:)/yend(1)

    contains

        subroutine chi_rhs(chi, y, dydchi, rhs_ctx)
            real(dp), intent(in) :: chi
            real(dp), intent(in) :: y(:)
            real(dp), intent(out) :: dydchi(:)
            class(*), intent(in), optional :: rhs_ctx
            real(dp) :: theta, gap, dtdchi, distance, slope, orbit(nvar), rate(nvar)
            real(dp) :: bmod_, sqrtg_, hder_(3), hcovar_(3), hctrvr_(3), hcurl_(3), x(3)
            associate (dummy => rhs_ctx)
            end associate
            theta = middle - direction*half*cos(chi)
            x = [s, 0.0_dp, theta]
            call do_magfie(x, bmod_, sqrtg_, hder_, hcovar_, hctrvr_, hcurl_)
            gap = 1.0_dp - eta*bmod_
            if (gap > gap_switch) then
                dtdchi = half*abs(sin(chi))/(v*sqrt(gap)*abs(hctrvr_(3)))
            else if (cos(chi) > 0.0_dp) then
                ! Turning point at chi = 0, distance 2*Delta*sin(chi/2)**2.
                distance = 2.0_dp*half*sin(chi/2.0_dp)**2
                slope = mean_slope(theta, middle - direction*half)
                gap = eta*slope*distance
                dtdchi = sqrt(2.0_dp*half/(eta*slope))*abs(cos(chi/2.0_dp))/ &
                    (v*abs(hctrvr_(3)))
            else
                distance = 2.0_dp*half*cos(chi/2.0_dp)**2
                slope = mean_slope(theta, middle + direction*half)
                gap = eta*slope*distance
                dtdchi = sqrt(2.0_dp*half/(eta*slope))*abs(sin(chi/2.0_dp))/ &
                    (v*abs(hctrvr_(3)))
            end if
            dydchi = 0.0_dp
            dydchi(1) = dtdchi
            if (.not. averages) return
            orbit = 0.0_dp
            orbit(1) = theta
            orbit(2) = direction*sign(1.0_dp, sin(chi))*sign_htheta*v*sqrt(gap)
            orbit(3:) = y(3:)
            call transport_rates(v, eta, nvar, y(1), orbit, bmod_, hder_, hcovar_, &
                                 hctrvr_, rate)
            dydchi(3:) = rate(3:)*dtdchi
        end subroutine chi_rhs

        real(dp) function mean_slope(theta, turning)
            ! |dB/dtheta| averaged between theta and the turning point, so
            ! that 1 - eta*B(theta) = eta*mean_slope*|turning - theta|.
            real(dp), intent(in) :: theta, turning
            real(dp) :: x(3), bmod_, sqrtg_, hder_(3), hcovar_(3), hctrvr_(3), hcurl_(3)
            integer :: k
            mean_slope = 0.0_dp
            x(1:2) = [s, 0.0_dp]
            do k = 1, size(nodes)
                x(3) = theta + (turning - theta)*(nodes(k) + 1.0_dp)/2.0_dp
                call do_magfie(x, bmod_, sqrtg_, hder_, hcovar_, hctrvr_, hcurl_)
                mean_slope = mean_slope + weights(k)*bmod_*hder_(3)
            end do
            mean_slope = abs(mean_slope)/2.0_dp
        end function mean_slope

        real(dp) function turning_point(side)
            ! First root of 1 - eta*B from th0 towards side. A grid scan finds
            ! a sign change or a local maximum of B that reflects the orbit;
            ! bisection then resolves the root to rounding.
            integer, parameter :: steps = 256
            real(dp), intent(in) :: side
            real(dp) :: inside, outside, theta, gap, slope, previous_slope
            integer :: k
            inside = th0
            call gap_and_slope(inside, side, gap, previous_slope)
            do k = 1, steps
                theta = th0 + side*2.0_dp*pi*real(k, dp)/real(steps, dp)
                call gap_and_slope(theta, side, gap, slope)
                if (gap < 0.0_dp) then
                    outside = theta
                    exit
                end if
                if (previous_slope > 0.0_dp .and. slope <= 0.0_dp) then
                    outside = local_maximum(inside, theta, side)
                    call gap_and_slope(outside, side, gap, slope)
                    if (gap < 0.0_dp) exit
                end if
                inside = theta
                previous_slope = slope
            end do
            if (k > steps) error stop 'trapped orbit without turning point'
            do
                theta = (inside + outside)/2.0_dp
                if (theta == inside .or. theta == outside) exit
                call gap_and_slope(theta, side, gap, slope)
                if (gap < 0.0_dp) then
                    outside = theta
                else
                    inside = theta
                end if
            end do
            turning_point = inside
        end function turning_point

        real(dp) function local_maximum(a, b, side)
            real(dp), intent(in) :: a, b, side
            real(dp) :: left, right, theta, gap, slope
            left = a
            right = b
            do
                theta = (left + right)/2.0_dp
                if (theta == left .or. theta == right) exit
                call gap_and_slope(theta, side, gap, slope)
                if (slope > 0.0_dp) then
                    left = theta
                else
                    right = theta
                end if
            end do
            local_maximum = left
        end function local_maximum

        subroutine gap_and_slope(theta, side, gap, slope)
            ! 1 - eta*B and dB/dtheta along side at theta.
            real(dp), intent(in) :: theta, side
            real(dp), intent(out) :: gap, slope
            real(dp) :: x(3), bmod_, sqrtg_, hder_(3), hcovar_(3), hctrvr_(3), hcurl_(3)
            x = [s, 0.0_dp, theta]
            call do_magfie(x, bmod_, sqrtg_, hder_, hcovar_, hctrvr_, hcurl_)
            gap = 1.0_dp - eta*bmod_
            slope = side*bmod_*hder_(3)
        end subroutine gap_and_slope
    end subroutine trapped_bounce

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
        real(dp) :: bmod, sqrtg, x(3), hder(3), hcovar(3), hctrvr(3), hcurl(3)

        x(1) = s
        x(2) = 0.0_dp
        x(3) = y(1)
        call do_magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
        call transport_rates(v, eta, neq, t, y, bmod, hder, hcovar, hctrvr, ydot)
    end subroutine timestep_transport

    subroutine transport_rates(v, eta, neq, t, y, bmod, hder, hcovar, hctrvr, ydot)
        ! timestep_transport for the field at theta = y(1).
        real(dp), intent(in) :: v, eta
        integer, intent(in) :: neq
        real(dp), intent(in) :: t
        real(dp), intent(in) :: y(neq), bmod, hder(3), hcovar(3), hctrvr(3)
        real(dp), intent(out) :: ydot(neq)

        ! BEGIN TODO: remove all of this after refactoring and re-use routine in orbit
        ! for y(1:3)
        real(dp) :: x(3), Om_tB_v
        real(dp) :: t0, Om_osc, OmtB_avg, dOmtBdv, dOmtBdeta, Omph_avg
        real(dp) :: dOmphdv, dOmphdeta
        real(dp) :: phase_drift, bmod1, vparallel, x1(3), sqrtg1, hder1(3), hcovar1(3), &
                    hctrvr1(3), hcurl1(3)
        complex(dp) :: epsn, Hn  ! relative amplitude of perturbation field epsn=Bn/B0
        ! and Hamiltonian Hn = (H - H0)_n

        x(1) = s
        x(2) = 0.0_dp
        x(3) = y(1)
        vparallel = y(2)
        if (eta < etatp) then
            ! A passing orbit never reverses. Enforce its energy invariant:
            ! evolving v_parallel accumulates an ODE energy error which is
            ! amplified by 1/(1-eta*B) near the trapped-passing boundary.
            vparallel = sign(vpar(v, eta, bmod), y(2))
        end if
        call poloidal_velocity(v, eta, bmod, hctrvr(3), hder(3), vparallel, ydot)

        ! evaluate orbit averages of Hamiltonian perturbation, optionally at the
        ! first-order guiding-centre radius (finite radial orbit width). Only the
        ! perturbation amplitude and bmod are taken there; everything else below
        ! (Hn phase, drift phase) stays on the label surface s.
        x1 = x
        bmod1 = bmod
        if (fow_radial) then
            x1(1) = s + fow_displacement(vparallel, hcovar(2))
            call do_magfie_displaced(x1, bmod1, sqrtg1, hder1, hcovar1, hctrvr1, hcurl1)
        end if
        if (pertfile) then
            call do_magfie_pert_amp(x1, epsn)
            epsn = pertfile_scale * epsn / bmod1
        else
            epsn = epsmn * exp(imun * m0 * y(1))
        end if

        ! Oscillating part of the toroidal drift phase,
        ! int (Om_tB(theta) - <Om_tB>) dt, accumulated in y(7).  Its toroidal
        ! phase is mph times this integral; <Om_tB> is the value that enters
        ! the resonance condition (Om_ph - Om_tE for trapped orbits).
        Om_osc = 0.0_dp
        phase_drift = 0.0_dp
        if (drift_phase) then
            if (associated(drift_osc_override)) then
                Om_osc = drift_osc_override(t)
            else if (magdrift .and. (eta > etatp .or. magdrift_passing > 0)) then
                if (eta > etatp) then
                    call Om_ph(v, eta, Omph_avg, dOmphdv, dOmphdeta)
                    OmtB_avg = Omph_avg - Om_tE
                else
                    call Om_tB(v, eta, OmtB_avg, dOmtBdv, dOmtBdeta)
                end if
                Om_tB_v = Om_tB_local_v2(eta, bmod, hder(1), hctrvr(3))
                Om_osc = Om_tB_v * v**2 - OmtB_avg
            end if
            phase_drift = mph * y(7)
        end if

        if (abs(t) <= 0.0_dp) then
            ! Local Hamiltonian samples have no bounce-time phase and need no
            ! frequency installed by a previous orbit integration.
            Hn = (2.0_dp - eta * bmod1) * epsn * &
                 exp(imun * (q * mph * y(1) + phase_drift))
        else if (eta > etatp) then
            !t0 = 0.25*2*pi/Omth ! Different starting position in orbit
            t0 = 0.0_dp
            Hn = (2.0_dp - eta * bmod1) * epsn * exp(imun * (q * mph * (y(1)) - mth * (t - &
                                                                                       t0) * Omth + phase_drift))
        else
            Hn = (2.0_dp - eta * bmod1) * epsn * exp(imun * (q * mph * (y(1)) - (mth + q * mph) &
                                                             * t * Omth + phase_drift))
        end if
        ! Preserve the default arithmetic above for bitwise Fast-build compatibility.
        if (curvature_drive) &
            call general_drive(eta, t, y(1), phase_drift, x1, bmod1, epsn, Hn)
        if (phi1) call add_electric_drive(v, eta, t, y(1), phase_drift, x1, Hn)
        ydot(3) = real(Hn)
        ydot(4) = aimag(Hn)

        ! evaluate orbit averages for nonlinear attenuation
        if (nonlin) then
            ydot(5) = 1.0_dp / bmod
            ydot(6) = bmod
        else
            ydot(5:6) = 0.0_dp
        end if

        ! y(7) integrates the oscillating drift phase (zero unless drift_phase).
        ! Zero any trailing integrands this routine does not compute. The DOP853 solver carries every component through
        ! its stage combinations, so an uninitialised derivative leaves a
        ! denormal in the state that trips the underflow FPE trap.
        ydot(7) = Om_osc
        ydot(8:) = 0.0_dp
    end subroutine transport_rates

    subroutine general_drive(eta, t, theta, phase_drift, x, bmod, epsn, Hn)
        ! General ideal drive in the u=0 representation (doc/hamiltonian_drive.md).
        real(dp), intent(in) :: eta, t, theta, phase_drift, x(3), bmod
        complex(dp), intent(in) :: epsn
        complex(dp), intent(out) :: Hn
        complex(dp) :: kappan, phase_factor
        real(dp) :: w_lagrangian, w_curvature

        phase_factor = orbit_phase(eta, t, theta, phase_drift)
        call do_magfie_pert_kappa_amp(x, kappan)
        kappan = pertfile_scale*kappan/bmod
        call drive_weight_leaf(eta*bmod, w_lagrangian, w_curvature)
        Hn = (w_lagrangian*epsn + w_curvature*kappan)*phase_factor
    end subroutine general_drive

    subroutine add_electric_drive(v, eta, t, theta, phase_drift, x, Hn)
        real(dp), intent(in) :: v, eta, t, theta, phase_drift, x(3)
        complex(dp), intent(inout) :: Hn
        complex(dp) :: electric

        call do_elec_pert_amp(x, electric)
        Hn = Hn + qi*electric/(mi*v**2/2.0_dp)*orbit_phase(eta, t, theta, phase_drift)
    end subroutine add_electric_drive

    function orbit_phase(eta, t, theta, phase_drift) result(phase_factor)
        real(dp), intent(in) :: eta, t, theta, phase_drift
        complex(dp) :: phase_factor

        if (abs(t) <= 0.0_dp) then
            phase_factor = exp(imun*(q*mph*theta + phase_drift))
        else if (eta > etatp) then
            phase_factor = exp(imun*(q*mph*theta - mth*t*Omth + phase_drift))
        else
            phase_factor = exp(imun*(q*mph*theta - (mth + q*mph)*t*Omth + phase_drift))
        end if
    end function orbit_phase

    function fow_displacement(v_par, hphi) result(ds)
        ! First-order radial guiding-centre displacement s(t) - s from the thin
        ! orbit on s. Conservation of p_phi = m v_par h_phi + (e/c) A_phi(s) with
        ! A_phi' = -sign_theta*iota*psi_pr (x = (s, phi, theta) has Jacobian
        ! -sqrtg, test_fow_gc_orbit) gives
        !   Delta s = m c q v_par h_phi / (e sign_theta psi_pr).
        ! s is the p_phi label, (e/c) A_phi(s) = p_phi, for trapped and passing
        ! orbits alike. This is the label of NEO-RT's frequencies: the shear term
        ! of Om_tB is q' <Delta s thetadot> with this Delta s. The orbit mean is
        ! therefore not removed; for passing orbits it is <Delta s> ~ +-rho_pol.
        real(dp), intent(in) :: v_par, hphi
        real(dp) :: ds

        ds = mi * c * q * v_par * hphi / (qi * sign_theta * psi_pr)
    end function fow_displacement

    subroutine do_magfie_displaced(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
        ! do_magfie at a displaced radius without side effects: do_magfie
        ! overwrites the flux-surface module quantities of the label surface s,
        ! which the caller and drift_phase (Om_tB_local_v2) still need. (The
        ! standalone B0h is also overwritten; nothing reads it after init.)
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: bmod, sqrtg, hder(3), hcovar(3), hctrvr(3), hcurl(3)
        real(dp) :: saved(7)

        saved = [Bthcov, Bphcov, dBthcovds, dBphcovds, q, dqds, iota]
        call do_magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
        Bthcov = saved(1)
        Bphcov = saved(2)
        dBthcovds = saved(3)
        dBphcovds = saved(4)
        q = saved(5)
        dqds = saved(6)
        iota = saved(7)
    end subroutine do_magfie_displaced

    include 'generated/drive_leaf.inc'

end module neort_transport
