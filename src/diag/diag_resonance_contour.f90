module diag_resonance_contour
    use iso_fortran_env, only: dp => real64
    use neort, only: init
    use neort_config, only: read_and_set_config
    use neort_main, only: runname
    use neort_profiles, only: read_and_init_profile_input, read_and_init_plasma_input, init_profiles, &
        vth, M_t, Om_tE
    use neort_freq, only: Om_th, Om_ph
    use driftorbit, only: mth, mph, etatp, etadt, epsst_spl, sign_vpar, magdrift, pertfile
    use do_magfie_mod, only: R0, s, inp_swi, do_magfie_init
    use do_magfie_pert_mod, only: do_magfie_pert_init
    implicit none

contains

    subroutine run_resonance_contour_diag(arg_runname)
        character(*), intent(in) :: arg_runname
        logical :: file_exists
        integer, parameter :: neta = 160
        integer :: u, i, k, nsurf
        real(dp), allocatable :: profile(:, :)
        real(dp) :: eta, eta0, eta1, ux, v, Omth, dOmthdv, dOmthdeta
        real(dp) :: Omph, dOmphdv, dOmphdeta, residual, rho

        runname = trim(adjustl(arg_runname))
        call read_and_set_config(trim(runname)//".in")
        inp_swi = 9
        call do_magfie_init("in_file")
        if (pertfile) call do_magfie_pert_init("in_file_pert")
        call init_profiles(R0)

        inquire(file="plasma.in", exist=file_exists)
        if (file_exists) call read_and_init_plasma_input("plasma.in", s)
        inquire(file="profile.in", exist=file_exists)
        if (file_exists) call read_and_init_profile_input("profile.in", s, R0, 1.0_dp, 1.0_dp)

        call read_profile_table("profile.in", profile)
        nsurf = size(profile, 1)

        ux = 1.5_dp
        mth = 0
        mph = 2
        magdrift = .true.
        sign_vpar = 1.0_dp

        open(newunit=u, file=trim(runname)//"_resonance_contour_neort.dat", status="replace", action="write")
        write(u, '(A)') "# rho_pol s eta ux mth mph residual"
        do k = 1, nsurf
            s = profile(k, 1)
            M_t = profile(k, 2)
            vth = profile(k, 3)
            Om_tE = 0.0_dp
            v = ux*vth
            call init
            eta0 = etatp*(1.0_dp + 1.0e-5_dp)
            eta1 = etatp + (etadt - etatp)*(1.0_dp - epsst_spl)
            rho = sqrt(max(0.0_dp, min(1.0_dp, s)))
            do i = 1, neta
                eta = eta0 + (eta1 - eta0)*real(i - 1, dp)/real(neta - 1, dp)
                call Om_th(v, eta, Omth, dOmthdv, dOmthdeta)
                call Om_ph(v, eta, Omph, dOmphdv, dOmphdeta)
                residual = real(mth, dp)*Omth + real(mph, dp)*Omph
                write(u, '(7ES18.9)') rho, s, eta, ux, real(mth, dp), real(mph, dp), residual
            end do
            write(u, *)
        end do
        close(u)
    end subroutine run_resonance_contour_diag

    subroutine read_profile_table(path, profile)
        character(*), intent(in) :: path
        real(dp), allocatable, intent(out) :: profile(:, :)
        integer :: u, n, ios
        real(dp) :: row(3)

        n = 0
        open(newunit=u, file=path, status="old", action="read")
        do
            read(u, *, iostat=ios) row
            if (ios /= 0) exit
            n = n + 1
        end do
        rewind(u)
        allocate(profile(n, 3))
        do n = 1, size(profile, 1)
            read(u, *) profile(n, :)
        end do
        close(u)
    end subroutine read_profile_table

end module diag_resonance_contour
