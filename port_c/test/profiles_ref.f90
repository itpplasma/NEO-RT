! Reference dump of profile + collision state from the real Fortran modules,
! for cross-checking the C profiles/collis port. Order matches test_profiles.c.
program profiles_ref
    use do_magfie_mod, only: set_s, read_boozer_file, init_magfie_at_s, do_magfie, &
                             inp_swi, bfac, R0, psi_pr, q
    use neort_profiles, only: read_and_init_plasma_input, read_and_init_profile_input, &
                              init_thermodynamic_forces, vth, M_t, Om_tE, A1, A2, &
                              ni1, ni2, Ti1, Ti2, Te
    use util, only: qi, mi
    use collis_alp, only: efcolf, velrat, enrat
    implicit none
    integer, parameter :: dp = kind(1.0d0)
    real(dp) :: x(3), bmod, sqrtg, bder(3), hcovar(3), hctrvr(3), hcurl(3), s
    character(len=256) :: bpath, ppath, rpath

    call get_command_argument(1, bpath)  ! in_file
    call get_command_argument(2, ppath)  ! plasma.in
    call get_command_argument(3, rpath)  ! profile.in
    s = 0.5d0
    inp_swi = 9
    bfac = 1.0d0
    call read_boozer_file(trim(bpath))
    call set_s(s)
    call init_magfie_at_s()
    x = [s, 0.0d0, 0.0d0]
    call do_magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)

    call read_and_init_plasma_input(trim(ppath), s)
    call read_and_init_profile_input(trim(rpath), s, R0, 1.0d0, 1.0d0)
    call init_thermodynamic_forces(psi_pr, q)

    write (*, '(ES24.16E3)') vth
    write (*, '(ES24.16E3)') M_t
    write (*, '(ES24.16E3)') Om_tE
    write (*, '(ES24.16E3)') A1
    write (*, '(ES24.16E3)') A2
    write (*, '(ES24.16E3)') ni1
    write (*, '(ES24.16E3)') ni2
    write (*, '(ES24.16E3)') Ti1
    write (*, '(ES24.16E3)') Ti2
    write (*, '(ES24.16E3)') Te
    write (*, '(ES24.16E3)') qi
    write (*, '(ES24.16E3)') mi
    write (*, '(ES24.16E3)') efcolf(1)
    write (*, '(ES24.16E3)') efcolf(2)
    write (*, '(ES24.16E3)') efcolf(3)
    write (*, '(ES24.16E3)') velrat(1)
    write (*, '(ES24.16E3)') velrat(2)
    write (*, '(ES24.16E3)') velrat(3)
    write (*, '(ES24.16E3)') enrat(1)
    write (*, '(ES24.16E3)') enrat(2)
    write (*, '(ES24.16E3)') enrat(3)
end program profiles_ref
