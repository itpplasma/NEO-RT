module neort_magfie
    use util, only: disp, pi
    use do_magfie_mod, only: s, do_magfie, eps, iota
    use neort_orbit, only: th0
    use driftorbit, only: B0, Bmin, Bmax, etadt, etatp, dVds

    implicit none

contains

    subroutine init_fsa
        ! Calculate the flux surface areas for normalization
        integer, parameter :: nth = 1000
        integer :: k
        real(8) :: thrange(nth), dth
        real(8) :: bmod, sqrtg, x(3), hder(3), hcovar(3), hctrvr(3), hcurl(3)

        print *, "       s: ", s

        thrange = -pi + (/(k*2*pi/nth, k=1, nth)/)

        dth = thrange(2) - thrange(1)
        x(1) = s
        x(2) = 0d0
        x(3) = 0d0

        dVds = 0d0
        B0 = 0d0
        print *, " eps orig: ", eps
        eps = 0d0

        Bmin = -1d0
        Bmax = 0d0

        do k = 1, nth
            x(3) = thrange(k)
            call do_magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
            dVds = dVds + sqrtg*dth
            B0 = B0 + bmod*dth
            eps = eps - cos(x(3))*bmod*dth

            ! TODO: do fine search
            if ((Bmin < 0) .or. (bmod < Bmin)) then
                Bmin = bmod
                th0 = x(3)
            end if
            if (bmod > Bmax) Bmax = bmod
        end do

        dVds = 2d0*pi*dVds
        B0 = B0/(2d0*pi)
        eps = eps/(B0*pi)

        print *, " eps calc: ", eps
        print *, "      th0: ", th0
        print *, "     dVds: ", dVds
        print *, "       B0: ", B0
        print *, "Bmin,Bmax: ", Bmin, Bmax
        print *, "        x: ", x
        print *, "     bmod: ", bmod
        print *, "    sqrtg: ", sqrtg
        print *, "     hder: ", hder
        print *, "   hcovar: ", hcovar
        print *, "   hctrvr: ", hctrvr
        print *, "    hcurl: ", hcurl

        call disp('init_fsa: iota       = ', iota)

    end subroutine init_fsa

    subroutine init_misc
        ! TODO: fine search for minima and maxima
        etatp = 1d0/Bmax
        etadt = 1d0/Bmin
    end subroutine init_misc


end module neort_magfie
