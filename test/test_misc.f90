program test_misc
    use iso_fortran_env, only: dp => real64
    use util, only: readdata, disp, c, qi, mi, pi
    use driftorbit, only: etamin, etamax, etatp, etadt, epsst_spl, epst_spl, epssp_spl, &
        epsp_spl, sign_vpar, mth, nlev, nopassing, epsp, epst, dVds
    use do_magfie_mod, only: do_magfie, psi_pr, q, s, eps, R0
    use do_magfie_pert_mod, only: mph
    use neort_profiles, only: vth, Om_tE, read_and_init_plasma_input
    use neort_main, only: main, runname
    use neort_orbit, only: bounce, nvar
    use neort_resonance, only: driftorbit_root, driftorbit_coarse
    use neort_magfie, only: init_flux_surface_average
    use neort_freq, only: Om_th, Om_ph, Om_tB, d_Om_ds

    implicit none

    call main

    call test_profile
    call test_torfreq
    call test_resline

contains

    subroutine test_torfreq
        integer, parameter :: n = 100
        integer :: k
        real(dp) :: Omph, dOmphdv, dOmphdeta, dOmphds
        real(dp) :: Omth, dOmthdv, dOmthdeta, dOmthds
        real(dp) :: OmtB, dOmtBdv, dOmtBdeta
        real(dp) :: aa, b
        real(dp) :: taub, bounceavg(nvar)
        real(dp) :: v, eta

        v = vth

        etamin = etatp
        etamax = etatp + (etadt - etatp) * (1.0_dp - epsst_spl)

        b = log(epst_spl)
        aa = 1.0_dp / (n - 1.0_dp) * (log(etamax / etamin - 1.0_dp) - b)
        eta = etamax

        sign_vpar = 1

        call disp("test_torfreq: vth        = ", vth)
        call disp("test_torfreq: v/vth      = ", v/vth)
        call disp("test_torfreq: mph        = ", 1.0_dp*mph)
        call disp("test_torfreq: mth        = ", 1.0_dp*mth)
        call disp("test_torfreq: Om_tE      = ", Om_tE)
        call disp("test_torfreq: Om_tB_ref  = ", c*mi*vth**2/(2*qi*psi_pr))
        call bounce(v, eta, taub, bounceavg)
        call disp("test_torfreq: Om_tB_ba   = ", vth**2*bounceavg(3))
        call disp("test_torfreq: etamin = ", etamin)
        call disp("test_torfreq: etamax = ", etamax)
        call Om_th(v, eta, Omth, dOmthdv, dOmthdeta)

        call disp("test_torfreq: Om_th_approx    = ", v/(q*R0*sqrt(2.0_dp/eps)))
        call disp("test_torfreq: Om_th_deeptrap  = ", Omth)

        open (unit=9, file=trim(adjustl(runname))//"_torfreq.out", recl=1024)
        write (9, *) "!1:eta                   "// &
            " 2:etatp                  "// &
            " 3:etadt                  "// &
            " 4:Om_tE                  "// &
            " 5:OmtB                   "// &
            " 6:dOmtbdv                "// &
            " 7:dOmtbdeta              "// &
            " 8:Omth                   "// &
            " 9:dOmthdv                "// &
            "10:dOmthdeta              "// &
            "11:dOmthds                "// &
            "12:Omph                   "// &
            "13:dOmphdv                "// &
            "14:dOmphdeta              "// &
            "15:dOmphds                "
        do k = 0, n - 1
            eta = etamin*(1.0_dp + exp(aa*k + b))
            call Om_ph(v, eta, Omph, dOmphdv, dOmphdeta)
            call Om_th(v, eta, Omth, dOmthdv, dOmthdeta)
            call Om_tB(v, eta, OmtB, dOmtBdv, dOmtBdeta)
            call d_Om_ds(v, eta, 2.0_dp*pi/Omth, dOmthds, dOmphds)
            write (9, *) eta, etatp, etadt, &
                Om_tE, OmtB, dOmtbdv, dOmtbdeta, &
                Omth, dOmthdv, dOmthdeta, dOmthds, &
                Omph, dOmphdv, dOmphdeta, dOmphds
        end do
        close (unit=9)
    end subroutine test_torfreq

    subroutine test_resline
        integer, parameter :: n = 500
        integer :: k
        real(dp) :: vmin, vmax, v
        real(dp) :: etarest(2), etaresp(2)
        real(dp) :: roots(nlev, 3)
        integer :: nroots, kr

        print *, "test_resline"

        vmin = 1.0e-6_dp * vth
        vmax = 10.0_dp * vth

        etaresp = etatp
        etarest = etatp
        sign_vpar = 1

        open (unit=9, file=trim(adjustl(runname))//"_resline_pco.out", recl=1024)
        open (unit=10, file=trim(adjustl(runname))//"_resline_pct.out", recl=1024)
        open (unit=11, file=trim(adjustl(runname))//"_resline_t.out", recl=1024)
        do k = 0, n - 1
            v = vmin + k / (n - 1.0_dp) * (vmax - vmin)

            if (.not. nopassing) then
                ! resonance (passing)
                sign_vpar = 1
                call driftorbit_coarse(v, etatp*epsp, etatp*(1 - epsp), roots, nroots)
                do kr = 1, nroots
                    etaresp = driftorbit_root(v, 1.0e-8_dp*abs(Om_tE), roots(kr, 1), roots(kr, 2))
                    write (9, *) v/vth, kr, etaresp(1), 0.0_dp, etatp
                end do
                sign_vpar = -1
                call driftorbit_coarse(v, etatp*epsp, etatp*(1 - epsp), roots, nroots)
                do kr = 1, nroots
                    etaresp = driftorbit_root(v, 1.0e-8_dp*abs(Om_tE), roots(kr, 1), roots(kr, 2))
                    write (10, *) v/vth, kr, etaresp(1), 0.0_dp, etatp
                end do
            end if

            ! resonance (trapped)
            sign_vpar = 1
            call driftorbit_coarse(v, etatp*(1 + epst), etadt*(1 - epst), roots, nroots)
            do kr = 1, nroots
                etarest = driftorbit_root(v, 1.0e-8_dp*abs(Om_tE), roots(kr, 1), roots(kr, 2))
                write (11, *) v/vth, kr, etarest(1), etatp, etadt
            end do
        end do
        close (unit=11)
        close (unit=10)
        close (unit=9)
    end subroutine test_resline

    subroutine test_profile
        integer :: k
        real(dp), allocatable :: data(:, :)

        print *, "test_profile"

        call readdata("profile.in", 3, data)

        open (unit=9, file=trim(adjustl(runname))//"_profile.out", recl=1024)
        write (9, *) "#s deltas=c*mi*vth*hcovar(2)*q/(qi*psi_pr) dVds q psi_pr"
        do k = 1, size(data, 1)
            s = data(k, 1)
            call output_flux_surface_data(9)
        end do
        s = 0.96_dp
        call output_flux_surface_data(9)
        s = 0.97_dp
        call output_flux_surface_data(9)
        s = 0.98_dp
        call output_flux_surface_data(9)
        s = 0.99_dp
        call output_flux_surface_data(9)
        s = 1.0_dp
        call output_flux_surface_data(9)

        close (unit=9)

        deallocate (data)
    end subroutine test_profile

    subroutine output_flux_surface_data(unit)
        integer, intent(in) :: unit
        real(dp) :: bmod, sqrtg, x(3), hder(3), hcovar(3), hctrvr(3), hcurl(3)
        x(1) = s
        x(2) = 0.0_dp
        x(3) = 0.0_dp
        call read_and_init_plasma_input("plasma.in", s)
        call init_flux_surface_average(s)
        call do_magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
        write (unit, *) x(1), c*mi*vth*hcovar(2)*q/(qi*psi_pr), dVds, q, psi_pr, vth
    end subroutine output_flux_surface_data
end program test_misc
