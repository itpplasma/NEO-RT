!
  subroutine phielec_of_psi(psi,phi_elec,dPhi_dpsi)
!
! Profile of normalized electrostatic potential $\hat\Phi$ and its derivative
! as functions of (dimensional) poloidal flux $\psi$
!
  use phielec_of_psi_mod, only : npolyphi,polyphi
  use field_eq_mod,       only : psi_sep,psi_axis
!
  implicit none
!
  logical :: prop=.true.
  integer :: i
  double precision :: psi,phi_elec,dPhi_dpsi,psihat
  double precision :: denom,expon
!  double precision, save :: ampl=1.12d0,width_phi=1d2
  double precision, save :: ampl=-1.12d0,width_phi=0.15d0
!
!  if(.true.) then
  if(.false.) then
    denom=width_phi*(psi_sep-psi_axis)
    expon=exp((psi_axis-psi)/denom)
    phi_elec=ampl*width_phi*(1.d0-expon)
    dPhi_dpsi=ampl*expon/(psi_sep-psi_axis)
    return
  endif
!
  if(prop) then
    prop=.false.
!    ampl=0.d0 
!    ampl=1.12d0  !negative electric field
    ampl=-1.12d0 !positive electric field
    polyphi=0.d0
    polyphi(1)=ampl/(psi_sep-psi_axis)
polyphi(2)=-polyphi(1)/(psi_sep-psi_axis)/2.d0
polyphi(3)=-polyphi(2)/(psi_sep-psi_axis)/2.d0
  endif
!
  phi_elec=0.d0
  dPhi_dpsi=0.d0
!
  do i=npolyphi,0,-1
    phi_elec=polyphi(i)+phi_elec*(psi-psi_axis)
  enddo
!
  do i=npolyphi,1,-1
    dPhi_dpsi=polyphi(i)*dble(i)+dPhi_dpsi*(psi-psi_axis)
  enddo
!
  end subroutine phielec_of_psi
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine denstemp_of_psi(psi,dens,temp,ddens,dtemp)
!
! Profile of density $n_\alpha$ (dens) and normalized temperature $\hat T_\alpha$ (temp) 
! as functions of poloidal flux $\psi$ and their derivatives over this flux (ddens,dtemp)
!
!
  use field_eq_mod,       only : psi_sep,psi_axis
!
  implicit none
!
  double precision :: psi,dens,temp,ddens,dtemp
!
  dens=1.d0-(psi-psi_axis)/(psi_sep-psi_axis)
  temp=1.d0
!
  ddens=-1.d0/(psi_sep-psi_axis)
  dtemp=0.d0
!
! set realistic denisty value:
  dens=dens*5d13
  ddens=ddens*5d13
  end subroutine denstemp_of_psi
