module neort_fp_native
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_quadrature, only: gauss_legendre
    use neort_fp_tensor, only: fp_branch_t, solve_fp_tensor, &
        neort_fp_tensor_coeff, neort_fp_phase_scale
    implicit none
    private
    public :: compute_fp_transport, fp_physical_rates, fp_rate_i

    abstract interface
        subroutine fp_rate_i(speed_cm_s, nud, duu)
            import :: dp
            real(dp), intent(in) :: speed_cm_s
            real(dp), intent(out) :: nud, duu
        end subroutine fp_rate_i
    end interface

contains

    subroutine fp_physical_rates(speed_cm_s, nud, duu)
        use neort_mixture_collisions, only: selected_fp_coefficients
        use neort_fp_controls, only: fp_rate_scale
        real(dp), intent(in) :: speed_cm_s
        real(dp), intent(out) :: nud, duu
        call selected_fp_coefficients(speed_cm_s/100.0_dp, nud, duu)
        nud = fp_rate_scale*nud
        duu = fp_rate_scale*duu
    end subroutine fp_physical_rates

    subroutine compute_fp_transport(speed_steps, pitch_steps, umax, rates, &
            torque, work, physical_response, info, energy_multiplier)
        use neort_lorentz_native, only: validate_lorentz_geometry, &
            sample_lorentz_orbit
        use neort_profiles, only: vth, ns, A1, A2, Om_tE
        use neort_fp_controls, only: validate_fp_controls, validate_fp_background
        use neort_fp_mode, only: tp_fp_collisions
        use driftorbit, only: etatp, etadt, mph, comptorque
        use do_magfie_mod, only: psi_pr, q, sign_theta
        use util, only: pi, qi, c
        integer, intent(in) :: speed_steps, pitch_steps
        real(dp), intent(in) :: umax
        procedure(fp_rate_i) :: rates
        real(dp), intent(out) :: torque(3)
        complex(dp), intent(out) :: work(2, 3)
        logical, intent(out) :: physical_response
        integer, intent(out) :: info
        real(dp), intent(in), optional :: energy_multiplier
        type(fp_branch_t) :: branches(3)
        real(dp), allocatable :: u(:), tau(:, :, :), inverse_b(:, :, :)
        real(dp), allocatable :: drift(:, :, :), nud(:, :), duu(:, :)
        complex(dp), allocatable :: hamiltonian(:, :, :)
        real(dp) :: nodes(3), weights(3), eta, unit_diffusion, omega
        real(dp) :: beta, prefactor, energy_scale, coordinate, speed, electric_frequency
        complex(dp) :: coefficients(2), physical_work
        integer :: branch, cell, point, node

        torque = 0.0_dp
        work = (0.0_dp, 0.0_dp)
        physical_response = .false.
        ! The native profile/source constants are SI-consistent only while
        ! this mode is active. A direct caller must honor that state contract.
        if (.not. tp_fp_collisions) then
            info = -20
            return
        end if
        call validate_fp_controls()
        call validate_fp_background()
        call validate_lorentz_geometry(info)
        if (info /= 0) return
        if (pitch_steps < 2 .or. speed_steps < 2 .or. &
                .not. ieee_is_finite(umax) .or. umax <= 0.0_dp) then
            info = -18
            return
        end if
        energy_scale = 1.0_dp
        if (present(energy_multiplier)) energy_scale = energy_multiplier
        if (.not. ieee_is_finite(energy_scale) .or. energy_scale < 0.0_dp) then
            info = -18
            return
        end if
        call gauss_legendre(3, nodes, weights)
        allocate(u(speed_steps + 1), nud(3, speed_steps), duu(3, speed_steps))
        do node = 1, size(u)
            u(node) = umax*real(node - 1, dp)/speed_steps
        end do
        do cell = 1, speed_steps
            do point = 1, 3
                speed = ((1.0_dp - nodes(point))*u(cell) + &
                    (1.0_dp + nodes(point))*u(cell + 1))/2.0_dp
                call rates(speed*vth, nud(point, cell), duu(point, cell))
            end do
        end do
        if (.not. all(ieee_is_finite(nud)) .or. &
                .not. all(ieee_is_finite(duu))) then
            info = -19
            return
        end if
        if (any(nud <= 0.0_dp) .or. any(duu < 0.0_dp)) then
            info = -19
            return
        end if
        allocate(tau(3, pitch_steps, 3), inverse_b(3, pitch_steps, 3))
        allocate(drift(3, pitch_steps, 3), hamiltonian(3, pitch_steps, 3))
        do branch = 1, 3
            allocate(branches(branch)%eta(pitch_steps + 1))
            do node = 1, pitch_steps + 1
                coordinate = real(node - 1, dp)/pitch_steps
                if (branch == 1) then
                    branches(branch)%eta(node) = etatp + &
                        (etadt - etatp)*coordinate
                else
                    branches(branch)%eta(node) = etatp*(1.0_dp - coordinate)
                end if
            end do
            ! Cache the actual native geometry on the caller's species thread.
            ! In this guarded thin-orbit phase, u*tau, H/u² and OmegaB/u²
            ! are independent of u. Later callback sampling is read-only.
            do cell = 1, pitch_steps
                do point = 1, 3
                    eta = ((1.0_dp - nodes(point))*branches(branch)%eta(cell) + &
                        (1.0_dp + nodes(point))*branches(branch)%eta(cell + 1))/2
                    call sample_lorentz_orbit(1.0_dp, branch, eta, 1.0_dp, &
                        0.0_dp, tau(point, cell, branch), unit_diffusion, omega, &
                        hamiltonian(point, cell, branch), &
                        drift_frequency=drift(point, cell, branch), &
                        inverse_b_average=inverse_b(point, cell, branch))
                end do
            end do
        end do
        beta = real(mph, dp)*q/((qi/c)*sign_theta*psi_pr)
        electric_frequency = real(mph, dp)*Om_tE
        coefficients = cmplx(0.0_dp, -beta, dp)*[A1, A2]
        call solve_fp_tensor(u, branches, sample, info, &
            parallel_samples=.true., physical_coefficients=coefficients, &
            physical_response=physical_response)
        if (info /= 0) return
        prefactor = -real(mph, dp)*sqrt(pi)*ns*vth*abs(psi_pr/q)
        do branch = 1, 3
            work(:, branch) = branches(branch)%work
            if (physical_response) then
                physical_work = work(1, branch)
            else
                physical_work = cmplx(0.0_dp, -beta, dp)* &
                    (A1*work(1, branch) + A2*work(2, branch))
            end if
            if (comptorque) torque(branch) = prefactor*aimag(physical_work)
        end do

    contains

        subroutine sample(branch_, ux, eta_, w, p, e, frequency, drive)
            integer, intent(in) :: branch_
            real(dp), intent(in) :: ux, eta_
            real(dp), intent(out) :: w, p, e, frequency
            complex(dp), intent(out) :: drive(2)
            real(dp) :: fraction, local, hreal, himag
            integer :: cu, cp, qu, qp
            cu = min(speed_steps, 1 + int(ux/umax*speed_steps))
            local = 2.0_dp*(ux - u(cu))/(u(cu + 1) - u(cu)) - 1.0_dp
            qu = minloc(abs(nodes - local), dim=1)
            fraction = (eta_ - branches(branch_)%eta(1))/ &
                (branches(branch_)%eta(pitch_steps + 1) - &
                    branches(branch_)%eta(1))*pitch_steps
            cp = min(pitch_steps, 1 + int(fraction))
            local = 2.0_dp*(fraction - real(cp - 1, dp)) - 1.0_dp
            qp = minloc(abs(nodes - local), dim=1)
            call neort_fp_tensor_coeff(ux, tau(qp, cp, branch_), nud(qu, cu), &
                energy_scale*duu(qu, cu), eta_, inverse_b(qp, cp, branch_), w, p, e)
            call neort_fp_phase_scale(ux, &
                real(hamiltonian(qp, cp, branch_), dp), &
                aimag(hamiltonian(qp, cp, branch_)), electric_frequency, &
                drift(qp, cp, branch_), hreal, himag, frequency)
            drive(1) = cmplx(hreal, himag, dp)
            drive(2) = ux**2*drive(1)
        end subroutine sample
    end subroutine compute_fp_transport
end module neort_fp_native
