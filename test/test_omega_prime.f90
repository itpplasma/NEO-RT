program test_omage_prime_prog
    use neort, only: read_control, check_magfie, init_profile, init_plasma, init_run, &
                     compute_transport_harmonic, runname, s, M_t
    use driftorbit
    implicit none

    real(8) :: DELTA = 1d-6

    type :: freq_data_t
        real(8) :: s, v, eta
        real(8) :: J(3), Jbar(3)
        real(8) :: Om, Ompr_old, Ompr_new
        real(8) :: dOmds, dOmdv, dOmdeta, dOmdpph
        real(8) :: Omth, dOmthds, dOmthdv, dOmthdeta
        real(8) :: Omph, dOmphds, dOmphdv, dOmphdeta
    end type

    real(8), allocatable :: spl_psi_pol(:, :)

    call main
contains

    subroutine main
        integer, parameter :: NUM_SAMPLES = 32
        character(1024), parameter :: TEST_RUN = "driftorbit64.new"
        real(8) :: s0, v0, eta0
        real(8) :: eta_res(2)

        type(freq_data_t) :: freq_data(NUM_SAMPLES)

        runname = TEST_RUN

        call read_control
        call setup
        s0 = s
        v0 = 0.9d0*vth
        eta_res = first_resonance(v0)
        !eta0 = eta_res(1)
        eta0 = 0.5d0*(etamin + etamax)

        call spline_psi_pol
        call sample_values(s0, v0, eta0, freq_data)
        call show_results(freq_data)
    end subroutine main

    subroutine sample_values(s0, v0, eta0, freq_data)
        real(8), intent(in) :: s0, v0, eta0
        type(freq_data_t), intent(out) :: freq_data(:)

        integer :: i
        real(8) :: xi

        freq_data(1)%s = s0
        freq_data(1)%v = v0
        freq_data(1)%eta = eta0
        call test_omega_prime(freq_data(1))
        do i = 2, size(freq_data)
            call random_number(xi)
            freq_data(i)%s = s0*(1d0 + DELTA*(xi-0.5d0))
            call random_number(xi)
            freq_data(i)%v = v0*(1d0 + DELTA*(xi-0.5d0))
            call random_number(xi)
            freq_data(i)%eta = eta0 + (etamax - etamin)*DELTA*(xi-0.5d0)
            call test_omega_prime(freq_data(i))
        end do
    end subroutine sample_values

    subroutine show_results(freq_data)
        type(freq_data_t), intent(in) :: freq_data(:)

        character(1024) :: temp_file
        integer :: funit, i

        temp_file = "/tmp/test_omega_prime.dat"

        print *, "Writing results to ", trim(adjustl(temp_file))

        open(newunit=funit, file=temp_file, recl=8192)
        write(funit, *) "s v eta J1 J2 J3 Jbar1 Jbar2 Jbar3 Om " &
            // "Ompr_old Ompr_new dOmds dOmdv dOmdeta dOmdpph Omth dOmthds " &
            // "dOmthdv dOmthdeta Omph dOmphds dOmphdv dOmphdeta"

        do i = 1, size(freq_data)
            write(funit, *) freq_data(i)
        end do

        close(funit)

    end subroutine show_results

    subroutine setup
        call do_magfie_init
        call init_profile
        call init_plasma
        call init_run(use_thermodynamic_profiles=.true.)
        call check_magfie

        mth = -3
        sigv = 1
        call get_trapped_region(etamin, etamax)
    end subroutine setup

    function first_resonance(v) result(eta_res)
        real(8), intent(in) :: v
        real(8) :: eta_res(2)
        real(8) :: roots(nlev, 3)
        integer :: nroots
        call driftorbit_coarse(v, etamin, etamax, roots, nroots)
        eta_res = driftorbit_root(v, 1d-8*abs(Om_tE), roots(1, 1), roots(1, 2))
    end function first_resonance

    subroutine test_omega_prime(f)
        type(freq_data_t), intent(inout) :: f
        real(8) :: bounceavg(nvar)

        call setup

        call compute_frequencies(f)

        print *, "dOmthds_test: ", dOmthds_test(f%v, f%eta)

        call bounce_fast(f%v, f%eta, 2d0*pi/abs(f%Omth), bounceavg)

        f%Ompr_old = omega_prime(f%v/vth, f%eta, f%Omth, f%dOmdv, f%dOmdeta, f%dOmdpph)
        f%Ompr_new = omega_prime_new(f%v/vth, f%eta, bounceavg, f%Omth, f%dOmdv, f%dOmdeta, f%dOmdpph)

        call compute_invariants(f%v, f%eta, f%J)

        ! Transformed actions for Galileo transformation near resonance
        f%Jbar(1) = f%J(1)
        f%Jbar(2) = f%J(2) - mth/mph*f%J(3)
        f%Jbar(3) = f%J(3)/mph
    end subroutine test_omega_prime

    subroutine compute_frequencies(f)
        type(freq_data_t), intent(inout) :: f

        real(8) :: taub, bounceavg(nvar)

        call Om_th(f%v, f%eta, f%Omth, f%dOmthdv, f%dOmthdeta)
        taub = 2d0*pi/abs(f%Omth)
        call bounce_fast(f%v, f%eta, taub, bounceavg)
        call Om_ph(f%v, f%eta, f%Omph, f%dOmphdv, f%dOmphdeta)
        call d_Om_ds(f%v, f%eta, f%dOmthds, f%dOmphds)
        f%Om = mth*f%Omth + mph*f%Omph
        f%dOmdv = mth*f%dOmthdv + mph*f%dOmphdv
        f%dOmdeta = mth*f%dOmthdeta + mph*f%dOmphdeta
        f%dOmds = mth*f%dOmthds + mph*f%dOmphds
        f%dOmdpph = -(qi/c*iota*psi_pr)**(-1)*f%dOmds
    end subroutine compute_frequencies

    subroutine compute_invariants(v, eta, J)
        real(8), intent(in) :: v, eta
        real(8), intent(out) :: J(3)

        integer, parameter :: neq = 3
        real(8) :: taub, bounceavg(neq)
        real(8) :: bounceint(neq + 1)
        real(8) :: bmod
        real(8) :: y0(neq)

        ! Initialize bounce-averated quantities y0. Their meaning
        ! is defined inside subroutine timestep (thin orbit integration)
        call evaluate_bfield_local(bmod)
        y0(1) = th0
        y0(2) = vpar(v, eta, bmod)
        y0(3) = 1d-15

        taub = 2.0*pi/(vperp(v, eta, bmod)*iota/R0*sqrt(eps/2d0))
        bounceint = bounce_integral(v, eta, neq, y0, taub/5d0, timestep_invariants)
        taub = bounceint(1)
        bounceavg = bounceint(2:)/taub

        print *, "taub: ", taub
        print *, "bounceavg: ", bounceavg

        J(1) = Jperp(v, eta)
        J(2) = bounceint(1+3)
        J(3) = pphi()

        print *, "psi_pol: ", psi_pol()

        print *, "J: ", J
    end subroutine compute_invariants

    subroutine timestep_invariants(v, eta, neq, t, y, ydot)
        real(8), intent(in) :: v, eta, t
        integer, intent(in) :: neq
        real(8), intent(in) :: y(neq)
        real(8), intent(out) :: ydot(neq)

        real(8) :: bmod, sqrtg, x(3), hder(3), hcovar(3), hctrvr(3), hcurl(3)

        x(1) = s
        x(2) = 0d0
        x(3) = y(1)

        call do_magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
        call timestep_poloidal_internal(v, eta, bmod, hctrvr(3), hder(3), 2, t, y(:2), ydot(:2))

        ydot(3) = mi*ydot(2)**2  ! Jpar
    end subroutine timestep_invariants

    function pphi()
        real(8) :: pphi
        pphi = -qi/c*psi_pol()
    end function pphi

    function psi_pol()
        real(8) :: psi_pol, psi_pol_result(3)
        psi_pol_result = spline_val_0(spl_psi_pol, s)
        psi_pol = psi_pol_result(1)
    end function psi_pol

    subroutine psi_pol_grid(psi_pol)
        real(8), allocatable, intent(out) :: psi_pol(:)

        integer :: i

        associate(sx => params0(:, 1), iota => params0(:, 2))
            allocate(psi_pol(size(sx)))
            psi_pol(1) = 0d0
            do i = 2, size(sx)
                psi_pol(i) = psi_pol(i-1) &
                    + psi_pr*trapz_step(sx(i-1), sx(i), iota(i-1), iota(i))
            end do
            psi_pol = psi_pol
        end associate
    end subroutine psi_pol_grid

    subroutine spline_psi_pol
        real(8), allocatable :: psi(:)

        associate(sx => params0(:, 1))
            call psi_pol_grid(psi)
            allocate(spl_psi_pol(size(sx) - 1, 5))
            spl_psi_pol = spline_coeff(sx, psi)
        end associate
    end subroutine spline_psi_pol

    function trapz(x, y)
        real(8), intent(in) :: x(:), y(:)
        real(8) :: trapz
        integer :: i

        trapz = 0d0
        do i = 2, size(x)
            trapz = trapz + trapz_step(x(i-1), x(i), y(i-1), y(i))
        end do
    end function trapz

    function trapz_step(x1, x2, y1, y2)
        real(8), intent(in) :: x1, x2, y1, y2
        real(8) :: trapz_step

        trapz_step = 0.5d0*(x2 - x1)*(y1 + y2)
    end function trapz_step

    function dOmthds_test(v, eta)
        real(8) :: dOmthds_test
        real(8), intent(in) :: v, eta

        real(8) :: s0, ds
        real(8) :: Omth1, Omth2

        ds = DELTA

        s0 = s
        s = s0 + ds
        Omth1 = 1d0 / bounce_time(v, eta)


        s = s0 - ds
        Omth2 = 1d0 / bounce_time(v, eta)

        dOmthds_test = (Omth1 - Omth2) / (2d0*ds)

        s = s0
    end function dOmthds_test

end program test_omage_prime_prog
