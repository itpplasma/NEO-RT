module transport
  use common
  implicit none

  contains

  function Hm(alpha, mth, mph)
  ! Computes Hamiltonian perturbation harmonic Hmn
    complex(8) :: Hm                 ! Complex harmonic of H
    real(8), intent(in) :: alpha(3)  ! Invariants J_perp, p_phi, H
    real(8), intent(in) :: mth       ! Canonical poloidal harmonic number
    real(8), intent(in) :: mph       ! Canonical toroidal harmonic number
    real(8) :: HmnReIm(2)

    call bounce_average(2, Hm_integrand, HmnReIm)

    contains

    subroutine Hm_integrand(z, Hmdot)
      real(8), intent(in)  :: z(5)
      real(8), intent(out) :: Hmdot(2)

      Hmdot = z(1)
    end subroutine Hm_integrand

  end function Hm

!   function flux_integrand(alpha, mth, mph)
!     real(8) :: flux_integrand
!     real(8), intent(in) :: alpha(3)  ! invariants J_perp, p_phi, H
!     real(8), intent(in) :: mth       ! canonical poloidal harmonic number
!     real(8), intent(in) :: mph       ! canonical toroidal harmonic number

!     real(8) :: Omth, dOmthdv, dOmthdeta
!     real(8) :: dummy, dummy2
!     real(8) :: OmtB
!     real(8) :: Hmn2
!     real(8) :: dpp, dhh, fpeff, dres, dnorm, thatt,& ! for nonlin
!          Omph, dOmphdv, dOmphdeta, dOmdv, dOmdeta, Ompr, dOmphds, dOmthds,&
!          dOmdpph
!     real(8) :: ma, mb, mc, md, me, mf
!     real(8) :: dvdJ, detadJ

!     v = ux*vth
!     eta = etax(1)
!     call Om_th(Omth, dOmthdv, dOmthdeta)
!     call Om_tB(OmtB, dummy, dummy2)
!     call bounce(2d0*pi/abs(Omth))
!     Hmn2 = (bounceavg(4)**2 + bounceavg(5)**2)*(mi*(ux*vth)**2/2d0)**2

!     thatt = 1d0
!     if (intoutput .or. nonlin) then
!        call Om_ph(Omph, dOmphdv, dOmphdeta)
!        call d_Om_ds(dOmthds, dOmphds)
!        dOmdv = mth*dOmthdv + mph*dOmphdv
!        dOmdeta = mth*dOmthdeta + mph*dOmphdeta
!        dOmdpph = -(qi/c*iota*psi_pr)**(-1)*(mth*dOmthds+mph*dOmphds)

!        ma = mi*(ux*vth)*mi*c/qi*eta
!        mb = mi*(ux*vth)**2/2*mi*c/qi
!        mc = mi/(2d0*Omth)*(ux*vth)*(1d0-eta*bounceavg(7))
!        md = mi*(ux*vth)**2/2d0*Omth
!        me = -mth/mph
!        mf = 1d0/mph

!        dvdJ = mb*me/(ma*md*mf-mb*mc*mf)
!        detadJ = ma*me/(mb*mc*mf-ma*md*mf)

!        !Ompr=mth*(eta*dOmdeta-ux*vth/2*dOmdv)/(mi*(ux*vth)**2/(2d0*Omth))+dOmdpph

!        Ompr = dOmdv*dvdJ + dOmdeta*detadJ + mph*dOmdpph

!        if (intoutput) then
!           ! 0:n, 1:l, 2:Eth, 3:Jperp_tp, 4:drphi/dpphi, 5:E/Eth, 6:Jperp/Jperp_tp, 7:rphi,
!           ! 8:|Hmn|, 9:Omth, 10:Omph, 11:Ombarprime, 12:dOmdv, 13:dOmdeta, 14:dOmdpphi, 15:sigma
!           ! 16:iota=1/q, 17: Om_tE, 18: dOmthds, 19: dOmphds
!           write(11, *) mth, mph, mi*vth**2/2d0, mi*(ux*vth)**2/2d0*mi*c/qi*etatp,&
!                -(qi/c*iota*psi_pr)**(-1), ux**2, eta/etatp, s, sqrt(Hmn2), Omth, Omph,&
!                Ompr, dOmdv, dOmdeta, dOmdpph, sigv, iota, Om_tE, dOmthds, dOmphds
!        end if
!        if (nonlin) then
!           call coleff(ux,dpp,dhh,fpeff)
!           dhh = vth*dhh
!           dpp = vth**3*dpp
!           dres = dpp*(dOmdv/Ompr)**2 + dhh*eta*(bounceavg(6)-eta)*(dOmdeta/Ompr)**2
!           dnorm = dres*sqrt(abs(Ompr))/sqrt(abs(Hmn2))**(3d0/2d0)
!           call attenuation_factor(dnorm,thatt)
!        end if
!     end if

!     Tphi_int = -pi**(3d0/2d0)*mph*ni1*c*vth/qi*ux**3*exp(-ux**2)*taub*Hmn2*thatt*(A1+A2*ux**2)*1d0/abs(etax(2))
!   end function flux_integrand

!   function Tphi_int_u(ux)
!     real(8) :: ux
!     real(8) :: Tphi_int_u, dTphi_int_u
!     real(8) :: eta_res(2)
!     real(8) :: roots(nlev, 3)
!     integer :: nroots, kr

!     v = ux*vth
!     Tphi_int_u = 0d0

!     call driftorbit_coarse(etamin, etamax, roots, nroots)
!     if(nroots == 0) return

!     do kr = 1,nroots
!        eta_res = driftorbit_root(1d-8*abs(Om_tE), roots(kr,1), roots(kr,2))
!        eta = eta_res(1)
!        dTphi_int_u = Tphi_int(ux, eta_res)
!        Tphi_int_u = Tphi_int_u + dTphi_int_u

!        if(orbit_mode_transp>0) then
!           call taurel
!           torque_int_box = torque_int_box + dTphi_int_u*taubins
!        end if
!     end do

!   end function Tphi_int_u

!   subroutine comp_transp_harm(mth, sigv)
!     integer(4), intent(in) :: mth


!   end subroutine comp_transp_harm
end module transport
