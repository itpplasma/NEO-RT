! Direct frequency/resonance scan using the same NEO-RT setup as production.
!
! The diagnostic deliberately calls neort_setup_at_s, Om_th, and Om_ph from
! the linked library.  It therefore traces the installed field/profile
! splines and the production orbit-frequency path instead of rebuilding an
! approximate frequency model in a post-processing script.
module diag_resonance_scan
    use iso_fortran_env, only: dp => real64
    use neort_lib, only: neort_init, neort_prepare_splines, neort_setup_at_s
    use neort, only: set_to_passing_region, set_to_trapped_region
    use neort_profiles, only: vth
    use neort_freq, only: Om_th, Om_ph
    use driftorbit, only: mth, mph, etatp, sign_vpar
    use do_magfie_mod, only: s, q
    implicit none

contains

    subroutine run_resonance_scan_diag(arg_runname, ux_arg)
        character(*), intent(in) :: arg_runname
        real(dp), intent(in), optional :: ux_arg
        character(len=16) :: branch_name
        integer, parameter :: nmth = 11, neta = 180
        integer :: u, ur, k, j, i, branch
        real(dp), allocatable :: profile(:, :)
        real(dp) :: ux, v, eta0, eta1, eta, eta_prev, res_prev
        real(dp) :: omth, domthdv, domthdeta
        real(dp) :: omph, domphdv, domphdeta
        real(dp) :: residual, dresdeta, rho, eta_root
        real(dp) :: mth_real, mph_real

        ux = 1.5_dp
        if (present(ux_arg)) ux = ux_arg
        if (ux <= 0.0_dp) error stop "UX must be positive"

        call neort_init(trim(arg_runname)//".in", "in_file", "in_file_pert")
        call neort_prepare_splines("plasma.in", "profile.in")
        call read_profile_table("profile.in", profile)

        open(newunit=u, file=trim(arg_runname)//"_resonance_scan.dat", status="replace", action="write")
        open(newunit=ur, file=trim(arg_runname)//"_resonance_roots.dat", status="replace", action="write")
        write(u, '(A)') "# s_tor rho_tor branch ux mth mph eta eta_over_etatp Omth_s-1 Omph_s-1 residual_s-1 dres_deta_s-1 q"
        write(ur, '(A)') "# s_tor rho_tor branch ux mth mph eta eta_over_etatp Omth_s-1 Omph_s-1 dres_deta_s-1 q"

        do k = 1, size(profile, 1)
            call neort_setup_at_s(profile(k, 1))
            rho = sqrt(max(0.0_dp, s))
            v = ux * vth

            do branch = 1, 3
                select case (branch)
                case (1)
                    branch_name = "passing_co"
                    sign_vpar = 1
                    call set_to_passing_region(eta0, eta1)
                case (2)
                    branch_name = "passing_ctr"
                    sign_vpar = -1
                    call set_to_passing_region(eta0, eta1)
                case (3)
                    branch_name = "trapped"
                    sign_vpar = 1
                    call set_to_trapped_region(eta0, eta1)
                end select

                do j = -nmth / 2, nmth / 2
                    mth = j
                    mth_real = real(mth, dp)
                    mph_real = real(mph, dp)
                    res_prev = 0.0_dp
                    eta_prev = eta0

                    do i = 1, neta
                        eta = eta0 + (eta1 - eta0) * real(i - 1, dp) / real(neta - 1, dp)
                        call Om_th(v, eta, omth, domthdv, domthdeta)
                        call Om_ph(v, eta, omph, domphdv, domphdeta)
                        residual = mth_real * omth + mph_real * omph
                        dresdeta = mth_real * domthdeta + mph_real * domphdeta
                        write(u, '(2ES18.9,1X,A,1X,ES18.9,1X,2I5,1X,7ES18.9)') &
                            s, rho, trim(branch_name), ux, mth, mph, eta, eta / etatp, &
                            omth, omph, residual, dresdeta, q

                        if (i > 1 .and. res_prev /= 0.0_dp) then
                            if ((residual > 0.0_dp .and. res_prev < 0.0_dp) .or. &
                                (residual < 0.0_dp .and. res_prev > 0.0_dp)) then
                                eta_root = 0.5_dp * (eta_prev + eta)
                                call Om_th(v, eta_root, omth, domthdv, domthdeta)
                                call Om_ph(v, eta_root, omph, domphdv, domphdeta)
                                dresdeta = mth_real * domthdeta + mph_real * domphdeta
                                write(ur, '(2ES18.9,1X,A,1X,ES18.9,1X,2I5,1X,6ES18.9)') &
                                    s, rho, trim(branch_name), ux, mth, mph, &
                                    eta_root, eta_root / etatp, omth, omph, dresdeta, q
                            end if
                        end if
                        res_prev = residual
                        eta_prev = eta
                    end do
                end do
            end do
        end do

        close(u)
        close(ur)
    end subroutine run_resonance_scan_diag


    subroutine read_profile_table(path, data)
        character(*), intent(in) :: path
        real(dp), allocatable, intent(out) :: data(:, :)
        integer :: n, ios, u
        real(dp) :: row(3)

        n = 0
        open(newunit=u, file=path, status="old", action="read")
        do
            read(u, *, iostat=ios) row
            if (ios /= 0) exit
            n = n + 1
        end do
        rewind(u)
        allocate(data(n, 3))
        do n = 1, size(data, 1)
            read(u, *) data(n, :)
        end do
        close(u)
    end subroutine read_profile_table

end module diag_resonance_scan
