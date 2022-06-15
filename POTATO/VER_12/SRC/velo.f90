!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine elefie(x, phi_elec, derphi)
!
  use field_eq_mod,      only : psif,dpsidr,dpsidz
!
  implicit none
!
  double precision :: phi_elec,dPhi_dpsi
  double precision, dimension(3) :: x,derphi
!
  call phielec_of_psi(psif,phi_elec,dPhi_dpsi)
!
  derphi(1) = dPhi_dpsi*dpsidr
  derphi(2) = 0d0
  derphi(3) = dPhi_dpsi*dpsidz
!
  end subroutine elefie
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      subroutine velo(tau,z,vz)
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
!
!NEW in version 4 =>      use parmot_mod, only : rmu,ro0
      use parmot_mod, only : rmu,ro0,gradpsiast, & !<=NEW in version 4
                             dpsiast_dR,dpsiast_dZ !<=NEW in version 4
      use field_eq_mod, only : ierrfield
!
      implicit none
!
      integer :: i
!
      double precision, intent(in)  :: tau, z
      double precision, intent(out) :: vz
      double precision x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl
      double precision derphi,phi_elec
      double precision p,alambd,p2,ovmu,gamma2,gamma,ppar,vpa,coala
      double precision rmumag,rovsqg,rosqgb,rovbm
      double precision a_phi,a_b,a_c,hstar
      double precision s_hc,hpstar,phidot,blodot,bra
      double precision pardeb
!
      dimension z(5), vz(5)
      dimension x(3),bder(3),hcovar(3),hctrvr(3),hcurl(3)
      dimension derphi(3)
      dimension a_phi(3),a_b(3),a_c(3),hstar(3)
!
      do 1 i=1,3
        x(i)=z(i)
 1    continue
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
      call magfie(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
!
! TODO: error handling magfie
      if(ierrfield.ne.0) then
!        vz=0.d0
!        return
         print *,'velo: magfie error',ierrfield
         stop
      endif
! in elefie: x(i)   - space coords (input, see above)
!            derphi - derivatives of the dimensionless electric potential
!                     phihat=e*phi/T over space coords (covar. vector)
!
      call elefie(x,phi_elec,derphi)
!
      p=z(4)
      alambd=z(5)
!
      p2=p*p
      ovmu=2.d0/rmu
      gamma2=p2*ovmu+1.d0
      gamma=dsqrt(gamma2)
      ppar=p*alambd
! vpa - dimensionless parallel velocity: vpa=v_parallel/sqrt(2*T/m)
      vpa=ppar/gamma
      coala=(1.d0-alambd**2)
! rmumag - magnetic moment
      rmumag=.5d0*p2*coala/bmod
!
      rovsqg=ro0/sqrtg
      rosqgb=.5d0*rovsqg/bmod
      rovbm=ro0/bmod
!
      a_phi(1)=(hcovar(2)*derphi(3)-hcovar(3)*derphi(2))*rosqgb
      a_b(1)=(hcovar(2)*bder(3)-hcovar(3)*bder(2))*rovsqg
      a_phi(2)=(hcovar(3)*derphi(1)-hcovar(1)*derphi(3))*rosqgb
      a_b(2)=(hcovar(3)*bder(1)-hcovar(1)*bder(3))*rovsqg
      a_phi(3)=(hcovar(1)*derphi(2)-hcovar(2)*derphi(1))*rosqgb
      a_b(3)=(hcovar(1)*bder(2)-hcovar(2)*bder(1))*rovsqg
!
      s_hc=0.d0
      do i=1,3
        a_c(i)=hcurl(i)*rovbm
        s_hc=s_hc+a_c(i)*hcovar(i)
        hstar(i)=hctrvr(i)+ppar*a_c(i)
      enddo
      hpstar=1.d0+ppar*s_hc
!
! velocities in the coordinate space
!
! phidot - derivative of the dmls el. potential over dmls time
! blodot - derivative of the logarith of the mag. field module over dmls time
      phidot=0.d0
      blodot=0.d0
      do i=1,3
        bra=vpa*hstar(i)+a_phi(i)+a_b(i)*rmumag/gamma
        vz(i)=bra/hpstar
        phidot=phidot+vz(i)*derphi(i)
        blodot=blodot+vz(i)*bder(i)
      enddo
!
! velocities in the phase space
!
      vz(4)=-0.5d0*gamma*phidot/p
!      vz(5)=coala/(alambd+dsign(1.d-32,alambd))*(vz(4)/p-0.5d0*blodot)
      vz(5)=-(0.5d0*coala/hpstar)*(sum(hstar*derphi)/p                 &
            + p*sum(hstar*bder)/gamma+alambd*sum(a_phi*bder))
!
      if(gradpsiast) then
        dpsiast_dR=bmod*hpstar*x(1)/vpa
        dpsiast_dZ=dpsiast_dR
        dpsiast_dR=dpsiast_dR*vz(3)
        dpsiast_dZ=-dpsiast_dZ*vz(1)
      endif
!
      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
