program test_omage_prime_prog
    use iso_fortran_env, only: dp => real64
    use do_magfie_mod, only: do_magfie_init, params0, sign_theta, s
    use do_magfie_pert_mod, only: do_magfie_pert_init
    use neort, only: check_magfie, write_magfie_data_to_files, init_profiles, &
                     init, set_to_trapped_region
    use neort_config, only: read_and_set_config
    use neort_main, only: runname
    use neort_datatypes, only: magfie_data_t
    use neort_profiles, only: read_and_init_plasma_input, read_and_init_profile_input, M_t
    use neort_orbit, only: th0, nvar, bounce_time, vpar, vperp, bounce_integral, &
                            bounce_fast, poloidal_velocity, evaluate_bfield_local
    use neort_freq, only: Om_th, Om_ph, d_Om_ds
    use neort_resonance, only: driftorbit_coarse, driftorbit_root
    use neort_transport, only: timestep_transport
    use neort_nonlin, only: omega_prime
    use driftorbit
    implicit none

    real(dp) :: DELTA = 1.0e-8_dp

    type :: freq_data_t
        real(dp) :: s, v, eta
        real(dp) :: J(3), Jbar(3)
        real(dp) :: Om, Ompr_old, Ompr_new
        real(dp) :: dOmds, dOmdv, dOmdeta, dOmdpph
        real(dp) :: Omth, dOmthds, dOmthdv, dOmthdeta
        real(dp) :: Omph, dOmphds, dOmphdv, dOmphdeta
        real(dp) :: Om_tE, dOm_tEds
    end type

    real(dp), allocatable :: spl_psi_pol(:, :)

    call main
contains

    subroutine main
        integer, parameter :: NUM_SAMPLES = 32
        real(dp) :: s0, v0, eta0
        real(dp) :: eta_res(2)

        type(freq_data_t) :: freq_data(NUM_SAMPLES)

        call get_command_argument(1, runname)
        call read_and_set_config(runname)
        call setup
        s0 = s
        v0 = 0.9_dp * vth
        eta_res = first_resonance(v0)
        eta0 = eta_res(1)
        !eta0 = 0.5_dp*(etamin + etamax)

        call spline_psi_pol
        call sample_values(s0, v0, eta0, freq_data)
        call show_results(freq_data)
    end subroutine main

    subroutine sample_values(s0, v0, eta0, freq_data)
        real(dp), intent(in) :: s0, v0, eta0
        type(freq_data_t), intent(out) :: freq_data(:)

        integer :: i
        real(dp) :: xi

        freq_data(1)%s = s0
        s = freq_data(1)%s
        freq_data(1)%v = v0
        freq_data(1)%eta = eta0
        call test_omega_prime(freq_data(1))
        do i = 2, size(freq_data)
            call random_number(xi)
            freq_data(i)%s = s0 * (1.0_dp + DELTA * (xi - 0.5_dp))
            s = freq_data(i)%s
            call random_number(xi)
            freq_data(i)%v = v0 * (1.0_dp + DELTA * (xi - 0.5_dp))
            call random_number(xi)
            freq_data(i)%eta = eta0 + (etamax - etamin) * DELTA * (xi - 0.5_dp)
            call test_omega_prime(freq_data(i))
        end do
    end subroutine sample_values

    subroutine show_results(freq_data)
        type(freq_data_t), intent(in) :: freq_data(:)

        character(1024) :: temp_file
        integer :: funit, i

        temp_file = "/tmp/test_omega_prime.dat"

        print *, "Writing results to ", trim(adjustl(temp_file))

        open (newunit=funit, file=temp_file, recl=8192)
        write (funit, *) "s v eta J1 J2 J3 Jbar1 Jbar2 Jbar3 Om " &
            //"Ompr_old Ompr_new dOmds dOmdv dOmdeta dOmdpph Omth dOmthds " &
            //"dOmthdv dOmthdeta Omph dOmphds dOmphdv dOmphdeta Om_tE dOm_tEds"

        do i = 1, size(freq_data)
            write (funit, *) freq_data(i)
        end do

        close (funit)

    end subroutine show_results

    subroutine setup
        type(magfie_data_t) :: magfie_data
        call do_magfie_init("in_file")
        if (pertfile) call do_magfie_pert_init("in_file_pert")
        call init_profiles(R0)
        call read_and_init_plasma_input("plasma.in", s)
        call read_and_init_profile_input("profile.in", s, R0, efac, bfac)
        call init
        call check_magfie(magfie_data)
        call write_magfie_data_to_files(magfie_data, runname)

        mth = -3
        sign_vpar = 1
        call set_to_trapped_region(etamin, etamax)
    end subroutine setup

    function first_resonance(v) result(eta_res)
        real(dp), intent(in) :: v
        real(dp) :: eta_res(2)
        real(dp) :: roots(nlev, 3)
        integer :: nroots
        call driftorbit_coarse(v, etamin, etamax, roots, nroots)
        eta_res = driftorbit_root(v, 1.0e-8_dp * abs(Om_tE), roots(1, 1), roots(1, 2))
    end function first_resonance

    subroutine test_omega_prime(f)
        type(freq_data_t), intent(inout) :: f
        real(dp) :: bounceavg(nvar)

        call setup

        call compute_frequencies(f)

        print *, "dOmthds_test: ", dOmthds_test(f%v, f%eta)

        call bounce_fast(f%v, f%eta, 2.0_dp * pi / abs(f%Omth), bounceavg, timestep_transport)

        f%Ompr_old = omega_prime_old(f%v / vth, f%eta, f%Omth, f%dOmdv, f%dOmdeta, f%dOmdpph)
        f%Ompr_new = omega_prime(f%v / vth, f%eta, bounceavg, f%Omth, f%dOmdv, f%dOmdeta, f%dOmdpph)

        call compute_invariants(f%v, f%eta, f%J)

        ! Transformed actions for Galileo transformation near resonance
        f%Jbar(1) = f%J(1)
        f%Jbar(2) = f%J(2) - mth / mph * f%J(3)
        f%Jbar(3) = f%J(3) / mph
    end subroutine test_omega_prime

    subroutine compute_frequencies(f)
        type(freq_data_t), intent(inout) :: f

        real(dp) :: taub, bounceavg(nvar)

        call Om_th(f%v, f%eta, f%Omth, f%dOmthdv, f%dOmthdeta)
        taub = 2.0_dp * pi / abs(f%Omth)
        call bounce_fast(f%v, f%eta, taub, bounceavg, timestep_transport)
        call Om_ph(f%v, f%eta, f%Omph, f%dOmphdv, f%dOmphdeta)
        call d_Om_ds(f%v, f%eta, taub, f%dOmthds, f%dOmphds)
        f%Om = mth * f%Omth + mph * f%Omph
        f%dOmdv = mth * f%dOmthdv + mph * f%dOmphdv
        f%dOmdeta = mth * f%dOmthdeta + mph * f%dOmphdeta
        f%dOmds = mth * f%dOmthds + mph * f%dOmphds
        f%dOmdpph = -(qi / c * iota * sign_theta * psi_pr)**(-1) * f%dOmds
        f%Om_tE = Om_tE
        f%dOm_tEds = dOm_tEds
    end subroutine compute_frequencies

    subroutine compute_invariants(v, eta, J)
        real(dp), intent(in) :: v, eta
        real(dp), intent(out) :: J(3)

        integer, parameter :: neq = 3
        real(dp) :: taub, bounceavg(neq)
        real(dp) :: bounceint(neq + 1)
        real(dp) :: bmod, htheta
        real(dp) :: y0(neq)

        ! Initialize bounce-averated quantities y0. Their meaning
        ! is defined inside subroutine timestep (thin orbit integration)
        call evaluate_bfield_local(bmod, htheta)
        y0(1) = th0
        y0(2) = sign(1.0_dp, htheta) * sign_vpar * vpar(v, eta, bmod)
        y0(3) = 1.0e-15_dp

        taub = 2.0_dp * pi / (vperp(v, eta, bmod) * iota / R0 * sqrt(eps / 2.0_dp))
        bounceint = bounce_integral(v, eta, neq, y0, taub / 5.0_dp, timestep_invariants)
        taub = bounceint(1)
        bounceavg = bounceint(2:) / taub

        print *, "taub: ", taub
        print *, "bounceavg: ", bounceavg

        J(1) = Jperp(v, eta)
        J(2) = bounceint(1 + 3)
        J(3) = pphi()

        print *, "s: ", s
        print *, "psi_pol: ", psi_pol()

        print *, "J: ", J
    end subroutine compute_invariants

    subroutine timestep_invariants(v, eta, neq, t, y, ydot)
        real(dp), intent(in) :: v, eta, t
        integer, intent(in) :: neq
        real(dp), intent(in) :: y(neq)
        real(dp), intent(out) :: ydot(neq)

        real(dp) :: bmod, sqrtg, x(3), hder(3), hcovar(3), hctrvr(3), hcurl(3)

        x(1) = s
        x(2) = 0.0_dp
        x(3) = y(1)

        call do_magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
        call poloidal_velocity(v, eta, bmod, hctrvr(3), hder(3), y(2), ydot(:2))

        ydot(3) = mi * ydot(2)**2  ! Jpar
    end subroutine timestep_invariants

    pure function Jperp(v, eta)
        real(dp) :: Jperp
        real(dp), intent(in) :: v, eta
        Jperp = 0.5_dp * mi * v**2 * mi * c / qi * eta
    end function Jperp

    function pphi()
        real(dp) :: pphi
        pphi = -qi / c * psi_pol()
    end function pphi

    function psi_pol()
        real(dp) :: psi_pol, psi_pol_result(3)
        psi_pol_result = spline_val_0(spl_psi_pol, s)
        psi_pol = psi_pol_result(1)
    end function psi_pol

    subroutine psi_pol_grid(psi)
        real(dp), allocatable, intent(out) :: psi(:)

        integer :: i

        associate (sx => params0(:, 1), iota => params0(:, 2))
            allocate (psi(size(sx)))
            psi(1) = 0.0_dp
            do i = 2, size(sx)
                psi(i) = psi(i - 1) &
                         + psi_pr*trapz_step(sx(i - 1), sx(i), iota(i - 1), iota(i))
            end do
        end associate
    end subroutine psi_pol_grid

    subroutine spline_psi_pol
        real(dp), allocatable :: psi(:)

        associate (sx => params0(:, 1))
            call psi_pol_grid(psi)
            allocate (spl_psi_pol(size(sx) - 1, 5))
            spl_psi_pol = spline_coeff(sx, psi)
        end associate
    end subroutine spline_psi_pol

    function trapz(x, y)
        real(dp), intent(in) :: x(:), y(:)
        real(dp) :: trapz
        integer :: i

        trapz = 0.0_dp
        do i = 2, size(x)
            trapz = trapz + trapz_step(x(i - 1), x(i), y(i - 1), y(i))
        end do
    end function trapz

    function trapz_step(x1, x2, y1, y2)
        real(dp), intent(in) :: x1, x2, y1, y2
        real(dp) :: trapz_step

        trapz_step = 0.5_dp * (x2 - x1) * (y1 + y2)
    end function trapz_step

    function dOmthds_test(v, eta)
        real(dp) :: dOmthds_test
        real(dp), intent(in) :: v, eta

        real(dp) :: s0, ds
        real(dp) :: Omth1, Omth2

        ds = DELTA

        s0 = s
        s = s0 + ds
        Omth1 = 2.0_dp * pi / bounce_time(v, eta)

        s = s0 - ds
        Omth2 = 2.0_dp * pi / bounce_time(v, eta)

        dOmthds_test = (Omth1 - Omth2) / (2.0_dp * ds)

        s = s0
    end function dOmthds_test

    function omega_prime_old(ux, eta, Omth, dOmdv, dOmdeta, dOmdpph)
        real(dp), intent(in) :: ux, eta, Omth, dOmdv, dOmdeta, dOmdpph
        real(dp) :: omega_prime_old
        omega_prime_old = mth * (eta * dOmdeta - ux * vth / 2 * dOmdv) / (mi * (ux * vth)**2 / &
                                                                          (2.0_dp * Omth)) + dOmdpph
    end function omega_prime_old

end program test_omage_prime_prog
