program test_omage_prime_prog
    use neort, only: read_control, check_magfie, init_profile, init_plasma, init_run, &
                     compute_transport_harmonic, runname, s, M_t
    use driftorbit
    implicit none

    type :: freq_data_t
        real(8) :: s, u, eta
        real(8) :: J(3), Jbar(3), Om, Ompr_old, Ompr_new
    end type

    call main
contains

    subroutine main
        integer, parameter :: NUM_SAMPLES = 16
        character(1024), parameter :: TEST_RUN = "driftorbit64.new"
        real(8) :: s0, ux0, eta0
        real(8) :: delta = 1d-10
        real(8) :: eta_res(2)

        type(freq_data_t) :: freq_data(NUM_SAMPLES)

        runname = TEST_RUN

        call read_control
        call setup
        s0 = s
        ux0 = v/vth
        eta_res = first_resonance()
        eta0 = eta_res(1)

        call sample_values(s0, ux0, eta0, delta, freq_data)
        call show_results(freq_data)
    end subroutine main

    subroutine sample_values(s0, ux0, eta0, delta, freq_data)
        real(8), intent(in) :: s0, ux0, eta0, delta
        type(freq_data_t), intent(out) :: freq_data(:)

        integer :: i
        real(8) :: xi

        freq_data(1)%s = s0
        freq_data(1)%u = ux0
        freq_data(1)%eta = eta0
        call test_omega_prime(freq_data(1))
        do i = 2, size(freq_data)
            call random_number(xi)
            freq_data(i)%s = s0*(1d0 + delta*(xi-0.5d0))
            call random_number(xi)
            freq_data(i)%u = ux0*(1d0 + delta*(xi-0.5d0))
            call random_number(xi)
            freq_data(i)%eta = eta0*(1d0 + delta*(xi-0.5d0))
            call test_omega_prime(freq_data(i))
        end do
    end subroutine sample_values

    subroutine show_results(freq_data)
        type(freq_data_t), intent(in) :: freq_data(:)

        character(1024) :: temp_file
        integer :: funit, i

        temp_file = "/tmp/test_omega_prime.dat"

        print *, "Writing results to ", trim(adjustl(temp_file))

        open(newunit=funit, file=temp_file, status='unknown')
        write(funit, *) "# s, u, eta, J(1), J(2), J(3), Jbar(1), Jbar(2), Jbar(3), Om, Ompr_old, Ompr_new"

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
        etamin = (1 + epst)*etatp
        etamax = (1 - epst)*etadt
    end subroutine setup

    function first_resonance() result(eta_res)
        real(8) :: eta_res(2)
        real(8) :: roots(nlev, 3)
        integer :: nroots
        call driftorbit_coarse(etamin, etamax, roots, nroots)
        eta_res = driftorbit_root(1d-8*abs(Om_tE), roots(1, 1), roots(1, 2))
    end function first_resonance

    subroutine test_omega_prime(freq_data)
        type(freq_data_t), intent(inout) :: freq_data
        real(8) :: dOmdv, dOmdeta, Ompr_old, Ompr_new, dOmdpph, Omth, dOmthdeta, vpar
        real(8) :: bounceavg(nvar)

        s = freq_data%s
        call setup

        freq_data%J(1) = Jperp()
        ! TODO: J(2) = Jpar, J(3) = pphi

        call compute_frequencies(freq_data%u, freq_data%eta, Omth, freq_data%Om, dOmdv, dOmdeta, dOmdpph)
        call bounce_fast(2d0*pi/abs(Omth), bounceavg)

        freq_data%Ompr_old = omega_prime(freq_data%u, Omth, dOmdv, dOmdeta, dOmdpph)
        freq_data%Ompr_new = omega_prime_new(freq_data%u, bounceavg,Omth, dOmdv, dOmdeta, dOmdpph)
    end subroutine test_omega_prime

    subroutine compute_frequencies(ux, etax, Omth, Om, dOmdv, dOmdeta, dOmdpph)
        real(8), intent(in) :: ux, etax
        real(8), intent(out) :: Omth, Om, dOmdv, dOmdeta, dOmdpph

        real(8) :: Omph, dOmphdv, dOmphdeta, dOmphds, dOmthds, dOmthdv, dOmthdeta

        real(8) :: taub, bounceavg(nvar)

        v = ux*vth
        eta = etax

        call Om_th(Omth, dOmthdv, dOmthdeta)
        taub = 2d0*pi/abs(Omth)
        call bounce_fast(taub, bounceavg)
        call Om_ph(Omph, dOmphdv, dOmphdeta)
        call d_Om_ds(taub, dOmthds, dOmphds)
        Om = mth*Omth + mph*Omph
        dOmdv = mth*dOmthdv + mph*dOmphdv
        dOmdeta = mth*dOmthdeta + mph*dOmphdeta
        dOmdpph = -(qi/c*iota*psi_pr)**(-1)*(mth*dOmthds + mph*dOmphds)
    end subroutine compute_frequencies

end program test_omage_prime_prog
