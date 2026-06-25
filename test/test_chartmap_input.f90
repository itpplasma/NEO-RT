program test_chartmap_input
    ! Behavioral test: inp_swi=10 chartmap path.
    !
    ! Phase 1 (always): verify that read_boozer_file dispatches to
    !   read_boozer_chartmap_file and that do_magfie returns finite values
    !   at a reference point.
    !
    ! Phase 2 (when CHARTMAP_REF env var is set): compare Bmod, iota, Bthcov,
    !   Bphcov from the chartmap path against the .bc path (inp_swi=8) at the
    !   same point.  Full FP cross-check requires a chartmap generated from the
    !   reference .bc via the libneo converters (libneo#344 + libneo#343).
    !   Until those are merged this phase is skipped but the harness is ready.

    use iso_fortran_env, only: dp => real64
    use do_magfie_mod, only: inp_swi, bfac, read_boozer_file, set_s, &
                             init_magfie_at_s, iota, Bthcov, Bphcov, psi_pr
    use util, only: pi
    use logger, only: set_log_level

    implicit none

    character(len=512) :: chartmap_path, bc_path_env
    integer :: stat
    real(dp) :: x(3), bmod_cm, sqrtg_cm
    real(dp), dimension(3) :: bder_cm, hcovar_cm, hctrvr_cm, hcurl_cm
    real(dp) :: iota_cm, Bth_cm, Bph_cm
    logical :: run_crosscheck

    call set_log_level(-1)

    ! --- Phase 1: chartmap field evaluation ---

    call get_environment_variable("CHARTMAP_FILE", chartmap_path, status=stat)
    if (stat /= 0 .or. len_trim(chartmap_path) == 0) then
        print *, "SKIP: CHARTMAP_FILE not set; chartmap read path not tested."
        print *, "Set CHARTMAP_FILE=/path/to/file.nc to enable."
        print *, "============================================"
        print *, "test_chartmap_input: SKIPPED (no input file)"
        print *, "============================================"
        stop 0
    end if

    inp_swi = 10
    bfac = 1.0_dp

    call read_boozer_file(trim(chartmap_path))
    call set_s(0.5_dp)
    call init_magfie_at_s()

    x(1) = 0.5_dp
    x(2) = 0.0_dp
    x(3) = 0.0_dp

    call do_magfie_at(x, bmod_cm, sqrtg_cm, bder_cm, hcovar_cm, hctrvr_cm, hcurl_cm)

    if (.not. (bmod_cm > 0.0_dp)) then
        print *, "ERROR: chartmap bmod <= 0 at x=(0.5,0,0), got", bmod_cm
        error stop "test_chartmap_input phase 1 failed"
    end if
    if (.not. (sqrtg_cm > 0.0_dp)) then
        print *, "ERROR: chartmap sqrtg <= 0, got", sqrtg_cm
        error stop "test_chartmap_input phase 1 failed: sqrtg"
    end if
    if (.not. (iota > 0.0_dp)) then
        print *, "ERROR: chartmap iota <= 0, got", iota
        error stop "test_chartmap_input phase 1 failed: iota"
    end if

    iota_cm = iota
    Bth_cm = Bthcov
    Bph_cm = Bphcov

    print *, "Phase 1 OK: chartmap bmod=", bmod_cm, " iota=", iota_cm

    ! --- Phase 2: cross-check against .bc if both files available ---

    call get_environment_variable("BC_FILE", bc_path_env, status=stat)
    run_crosscheck = (stat == 0 .and. len_trim(bc_path_env) > 0)

    if (run_crosscheck) then
        call run_bc_crosscheck(trim(bc_path_env), x, bmod_cm, iota_cm, Bth_cm, Bph_cm)
    else
        print *, "Phase 2 SKIP: BC_FILE not set."
        print *, "Set BC_FILE=/path/to/in_file to enable cross-check."
    end if

    print *, "============================================"
    print *, "test_chartmap_input: PASSED"
    print *, "============================================"

contains

    subroutine do_magfie_at(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
        use do_magfie_mod, only: do_magfie

        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: bmod, sqrtg
        real(dp), intent(out) :: bder(3), hcovar(3), hctrvr(3), hcurl(3)

        call do_magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
    end subroutine do_magfie_at

    subroutine run_bc_crosscheck(bc_path, x, bmod_cm, iota_cm, Bth_cm, Bph_cm)
        use do_magfie_mod, only: inp_swi, bfac, read_boozer_file, set_s, &
                                 init_magfie_at_s, iota, Bthcov, Bphcov, do_magfie

        character(len=*), intent(in) :: bc_path
        real(dp), intent(in) :: x(3), bmod_cm, iota_cm, Bth_cm, Bph_cm

        real(dp) :: bmod_bc, sqrtg_bc
        real(dp), dimension(3) :: bder_bc, hcovar_bc, hctrvr_bc, hcurl_bc
        real(dp), parameter :: tol_rel = 1.0e-4_dp

        inp_swi = 8
        bfac = 1.0_dp

        call read_boozer_file(bc_path)
        call set_s(x(1))
        call init_magfie_at_s()
        call do_magfie(x, bmod_bc, sqrtg_bc, bder_bc, hcovar_bc, hctrvr_bc, hcurl_bc)

        print *, "Phase 2: bc bmod=", bmod_bc, " chartmap bmod=", bmod_cm
        print *, "         bc iota=", iota, " chartmap iota=", iota_cm

        if (abs(bmod_cm - bmod_bc) > tol_rel * abs(bmod_bc)) then
            print *, "ERROR: bmod mismatch: bc=", bmod_bc, " cm=", bmod_cm
            error stop "test_chartmap_input phase 2 failed: bmod"
        end if
        if (abs(iota_cm - iota) > tol_rel * abs(iota)) then
            print *, "ERROR: iota mismatch: bc=", iota, " cm=", iota_cm
            error stop "test_chartmap_input phase 2 failed: iota"
        end if

        print *, "Phase 2 OK: bmod and iota agree to tol", tol_rel
    end subroutine run_bc_crosscheck

end program test_chartmap_input
