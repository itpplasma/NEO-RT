program test_parallel
    use iso_fortran_env, only: dp => real64
    use omp_lib
    use neort_datatypes, only: transport_data_t
    use do_magfie_mod, only: read_boozer_file, init_magfie_at_s, do_magfie, s, R0, bfac, &
                             psi_pr, Bthcov, Bphcov, q, iota, eps, B0h
    use do_magfie_pert_mod, only: init_mph_from_shared
    use neort_profiles, only: read_plasma_input, prepare_plasma_splines, init_plasma_at_s, &
                              prepare_profile_splines, init_profile_at_s, vth, dvthds, &
                              ni1, Ti1, M_t, Om_tE
    use neort, only: init, compute_transport
    use neort_config, only: read_and_set_config
    use driftorbit, only: efac, B0, Bmin, Bmax, dVds, etatp, etadt
    use neort_magfie, only: init_flux_surface_average
    use util, only: readdata
    use logger, only: set_log_level

    implicit none

    integer, parameter :: N_THREADS = 8
    real(dp), parameter :: TOL_INIT = 1.0d-12
    real(dp), parameter :: TOL_TRANSPORT = 1.0d-10

    type :: init_results_t
        real(dp) :: s, psi_pr, Bthcov, Bphcov, q, iota, eps, B0h
        real(dp) :: vth, dvthds, ni1, Ti1, M_t, Om_tE
        real(dp) :: B0, Bmin, Bmax, dVds, etatp, etadt
        real(dp) :: bmod, sqrtg
    end type

    call run_test()

contains

    subroutine run_test()
        type(init_results_t) :: init_results(N_THREADS)
        type(transport_data_t) :: transport_results(N_THREADS)
        integer :: nplasma
        real(dp) :: am1, am2, Z1, Z2
        real(dp), allocatable :: plasma(:, :), profile_data(:, :)
        real(dp) :: test_s
        integer :: i, thread_id
        logical :: all_match
        real(dp) :: x(3), bmod, sqrtg, bder(3), hcovar(3), hctrvr(3), hcurl(3)

        print *, "============================================"
        print *, "Parallel test with", N_THREADS, "threads"
        print *, "============================================"
        print *

        call omp_set_num_threads(N_THREADS)
        call set_log_level(-1)

        print *, "Reading input files..."
        call read_and_set_config("driftorbit")
        call read_boozer_file("in_file")
        call read_plasma_input("plasma.in", nplasma, am1, am2, Z1, Z2, plasma)
        call prepare_plasma_splines(nplasma, am1, am2, Z1, Z2, plasma)
        deallocate (plasma)
        call readdata("profile.in", 2, profile_data)
        call prepare_profile_splines(profile_data)
        deallocate (profile_data)
        print *, "Input files read."
        print *

        test_s = 0.5d0

        !$omp parallel private(thread_id, x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
        thread_id = omp_get_thread_num() + 1

        s = test_s
        call init_magfie_at_s()
        call init_plasma_at_s()
        call init_profile_at_s(R0, efac, bfac)
        call init_flux_surface_average(s)
        call init_mph_from_shared()

        x = [s, 0.3d0, 0.7d0]
        call do_magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)

        init_results(thread_id)%s = s
        init_results(thread_id)%psi_pr = psi_pr
        init_results(thread_id)%Bthcov = Bthcov
        init_results(thread_id)%Bphcov = Bphcov
        init_results(thread_id)%q = q
        init_results(thread_id)%iota = iota
        init_results(thread_id)%eps = eps
        init_results(thread_id)%B0h = B0h
        init_results(thread_id)%vth = vth
        init_results(thread_id)%dvthds = dvthds
        init_results(thread_id)%ni1 = ni1
        init_results(thread_id)%Ti1 = Ti1
        init_results(thread_id)%M_t = M_t
        init_results(thread_id)%Om_tE = Om_tE
        init_results(thread_id)%B0 = B0
        init_results(thread_id)%Bmin = Bmin
        init_results(thread_id)%Bmax = Bmax
        init_results(thread_id)%dVds = dVds
        init_results(thread_id)%etatp = etatp
        init_results(thread_id)%etadt = etadt
        init_results(thread_id)%bmod = bmod
        init_results(thread_id)%sqrtg = sqrtg

        call init()
        call compute_transport(transport_results(thread_id))

        !$omp end parallel

        print *, "All threads completed. Comparing results..."
        print *

        all_match = .true.

        print *, "=== Initialization results ==="
        do i = 2, N_THREADS
            call compare_init_results(init_results(1), init_results(i), i, all_match)
        end do

        print *
        print *, "=== Transport results ==="
        do i = 2, N_THREADS
            call compare_transport_results(transport_results(1), transport_results(i), i, all_match)
        end do

        print *
        if (all_match) then
            print *, "============================================"
            print *, "ALL TESTS PASSED"
            print *, "============================================"
        else
            print *, "============================================"
            print *, "TEST FAILED"
            print *, "============================================"
            error stop "Parallel test failed"
        end if

    end subroutine run_test

    subroutine compare_init_results(ref, other, thread_num, all_match)
        type(init_results_t), intent(in) :: ref, other
        integer, intent(in) :: thread_num
        logical, intent(inout) :: all_match

        call check_value("s", ref%s, other%s, thread_num, TOL_INIT, all_match)
        call check_value("psi_pr", ref%psi_pr, other%psi_pr, thread_num, TOL_INIT, all_match)
        call check_value("Bthcov", ref%Bthcov, other%Bthcov, thread_num, TOL_INIT, all_match)
        call check_value("Bphcov", ref%Bphcov, other%Bphcov, thread_num, TOL_INIT, all_match)
        call check_value("q", ref%q, other%q, thread_num, TOL_INIT, all_match)
        call check_value("iota", ref%iota, other%iota, thread_num, TOL_INIT, all_match)
        call check_value("eps", ref%eps, other%eps, thread_num, TOL_INIT, all_match)
        call check_value("B0h", ref%B0h, other%B0h, thread_num, TOL_INIT, all_match)
        call check_value("vth", ref%vth, other%vth, thread_num, TOL_INIT, all_match)
        call check_value("dvthds", ref%dvthds, other%dvthds, thread_num, TOL_INIT, all_match)
        call check_value("ni1", ref%ni1, other%ni1, thread_num, TOL_INIT, all_match)
        call check_value("Ti1", ref%Ti1, other%Ti1, thread_num, TOL_INIT, all_match)
        call check_value("M_t", ref%M_t, other%M_t, thread_num, TOL_INIT, all_match)
        call check_value("Om_tE", ref%Om_tE, other%Om_tE, thread_num, TOL_INIT, all_match)
        call check_value("B0", ref%B0, other%B0, thread_num, TOL_INIT, all_match)
        call check_value("Bmin", ref%Bmin, other%Bmin, thread_num, TOL_INIT, all_match)
        call check_value("Bmax", ref%Bmax, other%Bmax, thread_num, TOL_INIT, all_match)
        call check_value("dVds", ref%dVds, other%dVds, thread_num, TOL_INIT, all_match)
        call check_value("etatp", ref%etatp, other%etatp, thread_num, TOL_INIT, all_match)
        call check_value("etadt", ref%etadt, other%etadt, thread_num, TOL_INIT, all_match)
        call check_value("bmod", ref%bmod, other%bmod, thread_num, TOL_INIT, all_match)
        call check_value("sqrtg", ref%sqrtg, other%sqrtg, thread_num, TOL_INIT, all_match)
    end subroutine compare_init_results

    subroutine compare_transport_results(ref, other, thread_num, all_match)
        type(transport_data_t), intent(in) :: ref, other
        integer, intent(in) :: thread_num
        logical, intent(inout) :: all_match
        integer :: k

        call check_value("M_t", ref%summary%M_t, other%summary%M_t, thread_num, &
                         TOL_TRANSPORT, all_match)
        call check_value("Dco(1)", ref%summary%Dco(1), other%summary%Dco(1), thread_num, &
                         TOL_TRANSPORT, all_match)
        call check_value("Dco(2)", ref%summary%Dco(2), other%summary%Dco(2), thread_num, &
                         TOL_TRANSPORT, all_match)
        call check_value("Dctr(1)", ref%summary%Dctr(1), other%summary%Dctr(1), thread_num, &
                         TOL_TRANSPORT, all_match)
        call check_value("Dctr(2)", ref%summary%Dctr(2), other%summary%Dctr(2), thread_num, &
                         TOL_TRANSPORT, all_match)
        call check_value("Dt(1)", ref%summary%Dt(1), other%summary%Dt(1), thread_num, &
                         TOL_TRANSPORT, all_match)
        call check_value("Dt(2)", ref%summary%Dt(2), other%summary%Dt(2), thread_num, &
                         TOL_TRANSPORT, all_match)

        if (size(ref%harmonics) == size(other%harmonics)) then
            do k = 1, min(3, size(ref%harmonics))
                call check_value("Dresco(1)", ref%harmonics(k)%Dresco(1), &
                                 other%harmonics(k)%Dresco(1), thread_num, TOL_TRANSPORT, all_match)
                call check_value("Drest(1)", ref%harmonics(k)%Drest(1), &
                                 other%harmonics(k)%Drest(1), thread_num, TOL_TRANSPORT, all_match)
            end do
        else
            print *, "  ERROR: harmonics size mismatch:", size(ref%harmonics), "vs", &
                size(other%harmonics)
            all_match = .false.
        end if
    end subroutine compare_transport_results

    subroutine check_value(name, expected, actual, thread_num, tol, all_match)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: expected, actual, tol
        integer, intent(in) :: thread_num
        logical, intent(inout) :: all_match
        real(dp) :: diff, rel_diff

        diff = abs(expected - actual)
        if (abs(expected) > 1d-30) then
            rel_diff = diff / abs(expected)
        else
            rel_diff = diff
        end if

        if (diff > tol .and. rel_diff > tol) then
            print *, "  ERROR:", name, "mismatch in thread", thread_num
            print *, "    Expected:", expected
            print *, "    Actual:  ", actual
            print *, "    Diff:", diff, "Rel:", rel_diff
            all_match = .false.
        end if
    end subroutine check_value

end program test_parallel
