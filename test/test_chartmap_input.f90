program test_chartmap_input
    ! Behavioral test: inp_swi=10 chartmap path (read_boozer_chartmap_file +
    ! do_magfie_chartmap) cross-checked against the inp_swi=8 .bc path.
    !
    ! Phase 1 (always): verify that read_boozer_file dispatches to
    !   read_boozer_chartmap_file and that do_magfie returns finite, physically
    !   sane values (bmod>0, sqrtg>0, iota>0) at a reference point.
    !
    ! Phase 2 (when both CHARTMAP_FILE and BC_FILE env vars are set): compare
    !   bmod, iota, sqrtg, Bthcov, Bphcov between the chartmap and .bc paths at
    !   the same point.  Tolerances account for the interpolation difference
    !   between the Fourier-series .bc representation and the chartmap rho grid.
    !
    ! The CMake fixture gen_circ_chartmap generates circ_chartmap.nc from circ.bc
    ! and injects CHARTMAP_FILE/BC_FILE so both phases run unconditionally in CI.

    use iso_fortran_env, only: dp => real64
    use do_magfie_mod, only: inp_swi, bfac, read_boozer_file, set_s, &
        init_magfie_at_s, magfie_thread_init, do_magfie, &
        iota, Bthcov, Bphcov, psi_pr, a, dBthcovds, dBphcovds
    use util, only: pi
    use logger, only: set_log_level

    implicit none

    character(len=512) :: chartmap_path, bc_path_env
    integer :: stat
    real(dp) :: x(3), bmod_cm, sqrtg_cm
    real(dp) :: bder_cm(3), hcovar_cm(3), hctrvr_cm(3), hcurl_cm(3)
    real(dp) :: iota_cm, Bth_cm, Bph_cm, psip_cm, a_cm
    real(dp) :: bmod_bc, sqrtg_bc
    real(dp) :: bder_bc(3), hcovar_bc(3), hctrvr_bc(3), hcurl_bc(3)
    real(dp) :: iota_bc, Bth_bc, Bph_bc, psip_bc
    real(dp), parameter :: tol_rel = 1.0e-4_dp ! tightened: FP-accurate fixture from libneo #347
    logical :: run_crosscheck

    call set_log_level(-1)

    ! --- Phase 1: chartmap field evaluation ---

    call get_environment_variable("CHARTMAP_FILE", chartmap_path, status=stat)
    if (stat /= 0 .or. len_trim(chartmap_path) == 0) then
        print *, "FAIL: CHARTMAP_FILE not set."
        print *, "Set CHARTMAP_FILE=/path/to/circ_chartmap.nc and BC_FILE=/path/to/circ.bc"
        print *, "or run via CTest which injects them via the gen_circ_chartmap fixture."
        error stop "test_chartmap_input: CHARTMAP_FILE required"
    end if

    inp_swi = 10
    bfac = 1.0_dp

    call read_boozer_file(trim(chartmap_path))
    call set_s(0.5_dp)
    call init_magfie_at_s()

    x(1) = 0.5_dp
    x(2) = 0.0_dp
    x(3) = 0.0_dp

    call do_magfie(x, bmod_cm, sqrtg_cm, bder_cm, hcovar_cm, hctrvr_cm, hcurl_cm)

    iota_cm = iota
    Bth_cm  = Bthcov
    Bph_cm  = Bphcov
    psip_cm = psi_pr
    a_cm = a

    print *, "Phase 1: chartmap bmod   =", bmod_cm
    print *, "         chartmap sqrtg  =", sqrtg_cm
    print *, "         chartmap iota   =", iota_cm
    print *, "         chartmap Bthcov =", Bth_cm
    print *, "         chartmap Bphcov =", Bph_cm
    print *, "         chartmap psi_pr =", psip_cm
    print *, "         chartmap a      =", a_cm
    print *, "         chartmap hctrvr =", hctrvr_cm

    if (.not. (bmod_cm > 0.0_dp)) then
        print *, "ERROR: chartmap bmod <= 0, got", bmod_cm
        error stop "test_chartmap_input phase 1 failed: bmod"
    end if
    if (.not. (a_cm > 0.0_dp)) then
        print *, "ERROR: chartmap minor radius <= 0, got", a_cm
        error stop "test_chartmap_input phase 1 failed: minor radius"
    end if
    if (.not. (sqrtg_cm > 0.0_dp)) then
        print *, "ERROR: chartmap sqrtg <= 0, got", sqrtg_cm
        error stop "test_chartmap_input phase 1 failed: sqrtg"
    end if
    if (.not. (iota_cm > 0.0_dp)) then
        print *, "ERROR: chartmap iota <= 0, got", iota_cm
        error stop "test_chartmap_input phase 1 failed: iota"
    end if
    if (.not. (Bth_cm > 0.0_dp)) then
        print *, "ERROR: chartmap Bthcov <= 0 (sign bug), got", Bth_cm
        error stop "test_chartmap_input phase 1 failed: Bthcov sign"
    end if
    if (.not. (Bph_cm > 0.0_dp)) then
        print *, "ERROR: chartmap Bphcov <= 0 (sign bug), got", Bph_cm
        error stop "test_chartmap_input phase 1 failed: Bphcov sign"
    end if
    if (.not. (hctrvr_cm(3) > 0.0_dp)) then
        print *, "ERROR: chartmap hctrvr(3) <= 0 (sqrtg sign bug), got", hctrvr_cm(3)
        error stop "test_chartmap_input phase 1 failed: hctrvr sign"
    end if

    print *, "Phase 1 OK: chartmap bmod=", bmod_cm, " iota=", iota_cm

    ! --- Phase 2: cross-check against .bc if available ---

    call get_environment_variable("BC_FILE", bc_path_env, status=stat)
    run_crosscheck = (stat == 0 .and. len_trim(bc_path_env) > 0)

    if (.not. run_crosscheck) then
        print *, "Phase 2 SKIP: BC_FILE not set."
    else
        inp_swi = 8
        bfac = 1.0_dp

        ! Reset per-thread array allocation flags; inp_swi changed since last call.
        call magfie_thread_init()

        call read_boozer_file(trim(bc_path_env))
        call set_s(x(1))
        call init_magfie_at_s()
        call do_magfie(x, bmod_bc, sqrtg_bc, bder_bc, hcovar_bc, hctrvr_bc, hcurl_bc)

        iota_bc = iota
        Bth_bc  = Bthcov
        Bph_bc  = Bphcov
        psip_bc = psi_pr

        print *, "Phase 2: bc bmod   =", bmod_bc, " chartmap bmod  =", bmod_cm
        print *, "         bc iota   =", iota_bc, " chartmap iota  =", iota_cm
        print *, "         bc sqrtg  =", sqrtg_bc, " chartmap sqrtg =", sqrtg_cm
        print *, "         bc Bthcov =", Bth_bc,  " chartmap Bth   =", Bth_cm
        print *, "         bc Bphcov =", Bph_bc,  " chartmap Bph   =", Bph_cm
        print *, "         bc psi_pr =", psip_bc, " chartmap psi_pr=", psip_cm

        call assert_close("bmod",   bmod_cm,   bmod_bc,   tol_rel)
        call assert_close("iota",   iota_cm,   iota_bc,   tol_rel)
        call assert_close("sqrtg",  sqrtg_cm,  sqrtg_bc,  tol_rel)
        call assert_close("Bthcov", Bth_cm,    Bth_bc,    tol_rel)
        call assert_close("Bphcov", Bph_cm,    Bph_bc,    tol_rel)
        call assert_close("psi_pr", psip_cm,   psip_bc,   tol_rel)
        ! Legacy chartmaps have no aminor_m metadata; libneo then derives the
        ! equal-area edge radius, which differs slightly from the .bc header.
        call assert_close("minor radius", a_cm, a, 2.0e-3_dp)

        print *, "Phase 2 OK: all quantities agree to relative tol", tol_rel
    end if

    print *, "============================================"
    print *, "test_chartmap_input: PASSED"
    print *, "============================================"

contains

    subroutine assert_close(label, got, ref, tol)
        character(len=*), intent(in) :: label
        real(dp), intent(in) :: got, ref, tol

        real(dp) :: rel_err

        rel_err = abs(got - ref) / (abs(ref) + tiny(1.0_dp))
        if (rel_err > tol) then
            print *, "ERROR: ", trim(label), " mismatch: bc=", ref, " chartmap=", got, &
                " rel_err=", rel_err, " tol=", tol
            error stop "test_chartmap_input cross-check failed"
        end if
    end subroutine assert_close

end program test_chartmap_input
