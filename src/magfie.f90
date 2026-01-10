module neort_magfie
    use iso_fortran_env, only: dp => real64
    use util, only: disp, pi
    use do_magfie_mod, only: do_magfie, eps, iota
    use neort_orbit, only: th0
    use driftorbit, only: B0, Bmin, Bmax, etadt, etatp, dVds
    use logger, only: log_result

    implicit none

contains

    subroutine init_flux_surface_average(s)
        ! Calculate the flux surface areas for normalization
        ! Note: thread initialization is now handled by auto-init pattern in each module
        real(dp), intent(in) :: s
        integer, parameter :: nth = 1000
        integer :: k
        real(dp) :: thrange(nth), dth
        real(dp) :: bmod, sqrtg, x(3), hder(3), hcovar(3), hctrvr(3), hcurl(3)
        character(len=256) :: buffer

        write(buffer, "(A,ES12.5)") "        s: ", s
        call log_result(buffer)

        thrange = -pi + (/(k*2*pi/nth, k=1, nth)/)

        dth = thrange(2) - thrange(1)
        x(1) = s
        x(2) = 0.0_dp
        x(3) = 0.0_dp

        dVds = 0.0_dp
        B0 = 0.0_dp
        write(buffer, "(A,ES12.5)") " eps orig: ", eps
        call log_result(buffer)
        eps = 0.0_dp

        Bmin = -1.0_dp
        Bmax = 0.0_dp

        do k = 1, nth
            x(3) = thrange(k)
            call do_magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
            dVds = dVds + abs(sqrtg)*dth
            B0 = B0 + bmod*dth
            eps = eps - cos(x(3))*bmod*dth

            ! TODO: do fine search for minima and maxima
            if ((Bmin < 0) .or. (bmod < Bmin)) then
                Bmin = bmod
                th0 = x(3)
            end if
            if (bmod > Bmax) Bmax = bmod
        end do

        dVds = 2.0_dp * pi * dVds
        B0 = B0 / (2.0_dp * pi)
        eps = eps / (B0 * pi)

        etatp = 1.0_dp / Bmax
        etadt = 1.0_dp / Bmin

        write(buffer, "(A,ES12.5)") " eps calc: ", eps
        call log_result(buffer)
        write(buffer, "(A,ES12.5)") "      th0: ", th0
        call log_result(buffer)
        write(buffer, "(A,ES12.5)") "     dVds: ", dVds
        call log_result(buffer)
        write(buffer, "(A,ES12.5)") "       B0: ", B0
        call log_result(buffer)
        write(buffer, "(A,2ES12.5)") "Bmin,Bmax: ", Bmin, Bmax
        call log_result(buffer)
        write(buffer, "(A,3ES12.5)") "        x: ", x
        call log_result(buffer)
        write(buffer, "(A,ES12.5)") "     bmod: ", bmod
        call log_result(buffer)
        write(buffer, "(A,ES12.5)") "    sqrtg: ", sqrtg
        call log_result(buffer)
        write(buffer, "(A,3ES12.5)") "     hder: ", hder
        call log_result(buffer)
        write(buffer, "(A,3ES12.5)") "   hcovar: ", hcovar
        call log_result(buffer)
        write(buffer, "(A,3ES12.5)") "   hctrvr: ", hctrvr
        call log_result(buffer)
        write(buffer, "(A,3ES12.5)") "    hcurl: ", hcurl
        call log_result(buffer)

        write(buffer, "(A,ES12.5)") "init_flux_surface_average: iota =", iota
        call log_result(buffer)

    end subroutine init_flux_surface_average

end module neort_magfie
