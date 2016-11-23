!
  module collis_alp
    integer, parameter :: nsorts=3, ns=10000
    integer :: iswmod
    logical :: swcoll=.false.
    double precision, dimension(nsorts)    :: efcolf,velrat,enrat
!    double precision, dimension(nsorts,ns) :: efcolf_arr,velrat_arr,enrat_arr

  contains
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine coleff(p,dpp,dhh,fpeff)
!
!  Computes local values of dimensionless contravariant components
!  of collisional diffusion tensor and friction force for nonrelativistic
!  plasma. Backgound temperature is the same for all sorts.
!
!     Input variables:
!        formal: p      - dimensionless momentum module (p/(sqrt(2)*p_T)
!        common: efcolf - dmls collision frequencies
!                velrat - ratio of test species thermal velocity to
!                         background species thermal velocity
!     Output variables:
!        formal: dpp    - dimensionless momentum module diffusion
!                         coefficient
!                dhh    - dimensionless pitch angle diffusion coeff.
!                fpeff  - effective dimensionless drag force (prop. to linear
!                         deviation in Fokker-Planck eq.)
!
!
  implicit none
!
  integer i
  double precision :: p,dpp,dhh,fpeff,plim,xbeta,dp,dh,dpd
  !
  plim=max(p,1.d-8)
!
  dpp=0.0d0
  dhh=0.0d0
  fpeff=0.0d0
!
  do i=1,nsorts
    xbeta=p*velrat(i)
!
    call onseff(xbeta,dp,dh,dpd)
!
    dpp=dpp+dp*efcolf(i)
    dhh=dhh+dh*efcolf(i)
    fpeff=fpeff+(dpd/plim-2.0*dp*p*enrat(i))*efcolf(i)
  enddo
!
  dhh=dhh/plim**2
!
  return
  end subroutine coleff
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine onseff(v,dp,dh,dpd)
!
!  dp - dimensionless dpp
!  dh - dhh*p^2     (p - dmls)
!  dpd - (1/p)(d/dp)p^2*dp   (p - dmls)
!
  implicit none
!
! square root of pi
  double precision, parameter :: sqp=1.7724538d0
! cons=4./(3.*sqrt(pi))
  double precision, parameter :: cons=.75225278d0
  double precision :: v,dp,dh,dpd,v2,v3,ex,er
!
  v2=v**2
  v3=v2*v
  if(v.lt.0.01d0) then
    dp=cons*(1.d0-0.6d0*v2)
    dh=cons*(1.d0-0.2d0*v2)
    dpd=2.d0*cons*(1.d0-1.2d0*v2)
  elseif(v.gt.6.d0) then
    dp=1.d0/v3
    dh=(1.d0-0.5d0/v2)/v
    dpd=-1.d0/v3
  else
    ex=exp(-v2)/sqp
    er=erf(v)
    dp=er/v3-2.d0*ex/v2
    dh=er*(1.d0-0.5d0/v2)/v+ex/v2
    dpd=4.d0*ex-dp
  endif
!
  return
  end subroutine onseff
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      FUNCTION ERF(X)
!      PARAMETER  ( A1 = 0.07052 30784, A2 = 0.04228 20123,
!     ,             A3 = 0.00927 05272, A4 = 0.00015 10143,
!     ,             A5 = 0.00027 65672, A6 = 0.00004 30638 )
!      F(T) = 1./((1.+T*(A1+T*(A2+T*(A3+T*(A4+T*(A5+T*A6))))))**4)**4
!      W = 1. - F(ABS(X))
!      ERF = SIGN(W,X)
!      RETURN
!      END
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine loacol_nbi(amb,am1,am2,Zb,Z1,Z2,densi1,densi2,tempi1,tempi2,tempe,ebeam, &
                          v0,dchichi,slowrate,dchichi_norm,slowrate_norm)
!
!   Performs precomputation of the constants for Coulomb collision
!   operator for alpha-particles colliding with 2 sorts of ions and electrons
!
!
!   Normalisation: test particle (alpha) velocity is normalized by v0, alpha-particle
!   birth velocity, time is multiplied with v0 and has a meaning of free path of alpha
!   particle with v0. Respectively, collision frequencies have the meaning of inverse
!   mean free paths.
!
!   Input variables:
!        formal: am1,am2       - mass numbers of the first and second background ion species
!                Z1,Z2         - charge numbers of these species
!                densi1,densi2 - densities of these species, 1/cm**3
!                tempi1,tempi2,tempe - temperatures of two ion species and electrons, eV
!                ebeam         - initial beam particle energy, eV
!   Output variables:
!        formal: v0            - initial alpha particle velocity, cm/s
!                dchichi       - pitch angle scattering frequency, $D^{\chi\chi}$, of alpha
!                                particle with initial velocity, 1/s
!                slowrate      - slowing down rate, $F^v / v_0$, of alpha particle with 
!                                initial velocity, 1/s
!                dchichi_norim - normalized pitch angle scattering frequency, 1/cm
!                slowrate_norm - normalized slowing down rate, 1/cm
!        module collis_alp: 
!                efcolf - normalized collision frequencies
!                velrat - ratio of initial alpha particle velocity v0 to the
!                         specific background particle thermal velocity $v_{t}=\sqrt(2T/m)$
!                enrat  - ratio of initial alpha particle energy to the background species
!                         energy 
!
!
  implicit none
!
  double precision :: amb,am1,am2,Zb,Z1,Z2,densi1,densi2,tempi1,tempi2,tempe,ebeam,dense
  double precision :: v0,dchichi,slowrate,dchichi_norm,slowrate_norm,vti1,vti2,vte
  double precision :: pi,pmass,emass,e,ev,alame,frecol_base,alami1,alami2
!
  pi=3.14159265358979d0
  pmass=1.6726d-24
  emass=9.1094d-28
  e=4.8032d-10
  ev=1.6022d-12
!
  enrat(1)=ebeam/tempi1
  enrat(2)=ebeam/tempi2
  enrat(3)=ebeam/tempe
!
  v0=sqrt(2.d0*ebeam*ev/(amb*pmass))
  vti1=sqrt(2.d0*tempi1*ev/(pmass*am1))
  vti2=sqrt(2.d0*tempi2*ev/(pmass*am2))
  vte=sqrt(2.d0*tempe*ev/emass)
!
  velrat(1)=v0/vti1
  velrat(2)=v0/vti2
  velrat(3)=v0/vte
!
  dense=densi1*Z1+densi2*Z2
  alami1=23.d0-log(max(epsilon(1.d0), &
         sqrt(densi1*Z1**2/tempi1)*Zb*Z1*(amb+am1)/(amb*tempi1+am1*ebeam)))
  alami2=23.d0-log(max(epsilon(1.d0), &
         sqrt(densi2*Z2**2/tempi2)*Zb*Z2*(amb+am2)/(amb*tempi2+am2*ebeam)))
  alame=24.d0-log(sqrt(dense)/tempe)
  frecol_base=2.d0*pi*dense*e**4*Zb**2/((amb*pmass)**2*v0**3) !usual
  frecol_base=frecol_base/v0                                     !normalized
!
  efcolf(1)=frecol_base*Z1**2*alami1*densi1/dense
  efcolf(2)=frecol_base*Z2**2*alami2*densi2/dense
  efcolf(3)=frecol_base*alame
!
  efcolf=efcolf*velrat
!
!  s=0.d0
!  p=1.d0
!!
!  call coleff(s,p,dpp,dhh,fpeff)
!!
!  dchichi=dhh*v0
!  slowrate=abs(fpeff)*v0
!  dchichi_norm=dhh
!  slowrate_norm=abs(fpeff)
!
  return
  end subroutine loacol_nbi
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!


end module collis_alp
