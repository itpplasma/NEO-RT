module orbit
    use do_magfie_mod, only: do_magfie, q, iota, psi_pr, Bphcov, Bthcov
    use driftorbit, only: v, vth, eta, mi, c, qi, Om_tE

    implicit none
    save

contains

    function vecprod(u, v)
        real(8) :: vecprod(3)
        real(8), intent(in) :: u(3), v(3)

        ! this is for right-handed system x(1)-x(3)-x(2) and not 123 !!
        vecprod(1) = u(3)*v(2) - u(2)*v(3)
        vecprod(3) = u(2)*v(1) - u(1)*v(2)
        vecprod(2) = u(1)*v(3) - u(3)*v(1)
    end function vecprod

    subroutine step_zeroorder(neq, t, y, ydot)
        integer, intent(in) :: neq
        real(8), intent(in) :: t
        real(8), intent(in) :: y(neq)
        real(8), intent(out) :: ydot(neq)

        real(8) :: bmod, sqrtg, hder(3), hcovar(3), hctrvr(3), hcurl(3)

        call do_magfie(y(1:3), bmod, sqrtg, hder, hcovar, hctrvr, hcurl)

        ydot(1) = 0d0                                               ! r
        ydot(2) = y(4)*hctrvr(2)                                    ! phi
        ydot(3) = y(4)*hctrvr(3)                                    ! theta
        ydot(4) = -v**2*eta/2d0*hctrvr(3)*hder(3)*bmod              ! v_par

    end subroutine step_zeroorder

    subroutine step_full(neq, t, y, ydot)
        ! TODO: add gradPhi for fast rotation

        integer, intent(in) :: neq
        real(8), intent(in) :: t
        real(8), intent(in) :: y(neq)
        real(8), intent(out) :: ydot(neq)

        real(8) :: bmod, sqrtg, hder(3), hcovar(3), hctrvr(3), hcurl(3)

        real(8) :: hstar(3)       ! hstar = h + vpar/(omc)*curl(h)
        real(8) :: gradH(3)       ! gradH = (mu*gradB + e*gradPhi)
        real(8) :: hstarpar       ! hstar_par = h*hstar
        real(8) :: Er             ! Er = dphi/dr = psi_pr/(q*c)*Om_tE, no fast rotation yet
        real(8) :: sqgbmod2

        call do_magfie(y(1:3), bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
        sqgbmod2 = psi_pr*(Bphcov + iota*Bthcov)

        hstar = hctrvr + mi*c/(qi*bmod)*y(4)*hcurl
        hstarpar = 1d0 + mi*c/(qi*bmod)*y(4)*sum(hcovar*hcurl)

        !Er = psi_pr/(q*c)*Om_tE
        Er = 0.0
        gradH = mi*v**2/2d0*eta*bmod*hder + (/qi*Er, 0d0, 0d0/)

        ydot(1:3) = 1d0/hstarpar*(y(4)*hstar + c/(qi*sqgbmod2)*bmod*vecprod(hcovar, gradH))
        ydot(4) = -1d0/(mi*hstarpar)*sum(ydot(1:3)/y(4)*gradH)

        !print *, eta, (1d0-y(4)**2/v**2)/bmod

    end subroutine step_full

    subroutine step_rela(neq, t, y, ydot)
        integer, intent(in) :: neq
        real(8), intent(in) :: t
        real(8), intent(in) :: y(neq)
        real(8), intent(out) :: ydot(neq)

        real(8) :: tau, z(neq + 1), zdot(neq + 1)

        tau = vth*t
        z(1:3) = y(1:3)
        z(4) = v/vth
        z(5) = y(4)/v

        call velo(tau, z, zdot)

        ydot(1:3) = zdot(1:3)*vth
        ydot(4) = zdot(5)*v*vth

    end subroutine step_rela

    subroutine velo(tau, z, vz)
        !
        !
        !  Computes the components of the 5-D velocity          -  vz
        !  for given set of phase space coordinates             -  z.
        !
        !  Warning !!!  The dimensionless time is used (see below)
        !
        !  Order of the coordinates is the following:
        !  z(i) = x(i)  for i=1,2,3     - spatial coordinates with real
        !                                 dimension of the general covariant
        !                                 space coordinate system
        !  z(4) = p                     - momentum  module normalized to
        !                                 thermal momentum and sqrt(2);
        !  z(5) = alambd                - cosine of the pitch-angle
        !
        !  Input parameters:
        !            formal:  tau    -   dimensionless time: tau=sqrt(2*T/m)*t
        !                     z      -   see above
        !            common:  rmu    -   inverse relativistic temperature
        !                     ro0    -   Larmor radius for the reference
        !                                magnetic field and temperature:
        !                                ro0=sqrt(2*T/m)/(e*B_ref/m*c)
        !                     p0     -   dimensionless momentum module in the
        !                                initial point
        !                     alamb0 -   cos of pitch-angle in the initial point
        !  Output parameters:
        !            formal:  vz     -   see above
        !
        !  Called routines: magfie, elefie

        implicit none

        integer :: i

        double precision tau, z, vz
        double precision x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl
        double precision derphi
        double precision p, alambd, p2, ovmu, gamma2, gamma, ppar, vpa, coala
        double precision rmumag, rovsqg, rosqgb, rovbm
        double precision a_phi, a_b, a_c, hstar
        double precision s_hc, hpstar, phidot, blodot, bra
        double precision pardeb
        double precision ro0
        double precision rmu

        dimension z(5), vz(5)
        dimension x(3), bder(3), hcovar(3), hctrvr(3), hcurl(3)
        dimension derphi(3)
        dimension a_phi(3), a_b(3), a_c(3), hstar(3)

        ro0 = vth*mi*c/qi ! Larmor radius scaling
        rmu = 1d5         ! non-relativistic limit

        do i = 1, 3
            x(i) = z(i)
        end do
        !
        ! in magfie: x(i)   - set of 3 curvilinear space coordinates (input)
        !            bmod   - dimensionless magnetic field module: bmod=B/B_ref
        !            sqrtg  - Jacobian of space coordinates (square root of
        !                     metric tensor
        !            bder   - derivatives of logarithm of bmod over space coords
        !                     (covariant vector)
        !            hcovar - covariant components of the unit vector along
        !                     the magnetic field
        !            hctrvr - contravariant components of the unit vector along
        !                     the magnetic field
        !            hcurl  - contravariant components of the curl of this vector
        !
        call do_magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
        ! in elefie: x(i)   - space coords (input, see above)
        !            derphi - derivatives of the dimensionless electric potential
        !                     phihat=e*phi/T over space coords (covar. vector)
        !
        derphi = 0d0
        !      derphi(1) = psi_pr/(q*c)*Om_tE*qi/(mi*vth**2/2d0)

        p = z(4)
        alambd = z(5)

        p2 = p*p
        ovmu = 2.d0/rmu
        gamma2 = p2*ovmu + 1.d0
        gamma = dsqrt(gamma2)
        ppar = p*alambd
        ! vpa - dimensionless parallel velocity: vpa=v_parallel/sqrt(2*T/m)
        vpa = ppar/gamma
        coala = (1.d0 - alambd**2)
        ! rmumag - magnetic moment
        rmumag = .5d0*p2*coala/bmod

        rovsqg = ro0/sqrtg
        rosqgb = .5d0*rovsqg/bmod
        rovbm = ro0/bmod

        ! minus signs due to right-handed system r,theta,phi = x(1,3,2)
        a_phi(1) = -(hcovar(2)*derphi(3) - hcovar(3)*derphi(2))*rosqgb
        a_b(1) = -(hcovar(2)*bder(3) - hcovar(3)*bder(2))*rovsqg
        a_phi(2) = -(hcovar(3)*derphi(1) - hcovar(1)*derphi(3))*rosqgb
        a_b(2) = -(hcovar(3)*bder(1) - hcovar(1)*bder(3))*rovsqg
        a_phi(3) = -(hcovar(1)*derphi(2) - hcovar(2)*derphi(1))*rosqgb
        a_b(3) = -(hcovar(1)*bder(2) - hcovar(2)*bder(1))*rovsqg

        s_hc = 0.d0
        do i = 1, 3
            a_c(i) = hcurl(i)*rovbm
            s_hc = s_hc + a_c(i)*hcovar(i)
            hstar(i) = hctrvr(i) + ppar*a_c(i)
        end do
        hpstar = 1.d0 + ppar*s_hc

        ! velocities in the coordinate space
        !
        ! phidot - derivative of the dmls el. potential over dmls time
        ! blodot - derivative of the logarith of the mag. field module over dmls time
        phidot = 0.d0
        blodot = 0.d0
        do i = 1, 3
            bra = vpa*hstar(i) + a_phi(i) + a_b(i)*rmumag/gamma
            vz(i) = bra/hpstar
            phidot = phidot + vz(i)*derphi(i)
            blodot = blodot + vz(i)*bder(i)
        end do
        !
        ! velocities in the phase space
        !
        vz(4) = -0.5d0*gamma*phidot/p
        vz(5) = -(0.5d0*coala/hpstar)*(sum(hstar*derphi)/p &
                                       + p*sum(hstar*bder)/gamma + alambd*sum(a_phi*bder))
    end subroutine velo

end module orbit
