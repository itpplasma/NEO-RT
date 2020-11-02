module orbit_thin
    implicit none

    contains

    subroutine timestep(t, y, ydot)
        !
        !  Timestep function for thin orbit integration.
        !  Includes poloidal angle theta and parallel velocity.
        !  More integrands may be added starting from y(3)
        !

        real(8), intent (in) :: t
        real(8), intent (in) :: y(*)
        real(8), intent (out) :: ydot(*)

        real(8) :: bmod, sqrtg, x(3), hder(3), hcovar(3), hctrvr(3), hcurl(3)
        real(8) :: eta

        ! TODO
        ! call magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)

        eta = bmod*(1.0d0 - (y(5)/y(4))**2)

        ydot(1) = 0.0d0                                             ! r
        ydot(2) = 0.0d0                                             ! varphi
        ydot(3) = y(5)*hctrvr(3)                                    ! theta
        ydot(4) = 0.0d0                                             ! p
        ydot(5) = -y(4)**2*eta/2d0*hctrvr(3)*hder(3)*bmod           ! v_par
    end subroutine timestep
end module orbit_thin
