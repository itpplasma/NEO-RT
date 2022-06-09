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
  double precision spol
  double precision :: psi,phi_elec,dPhi_dpsi

  spol = (psi-psi_axis)/(psi_sep-psi_axis)

  phi_elec=0.d0
  dPhi_dpsi=0.d0

  do i=npolyphi,0,-1
    phi_elec = phi_elec + polyphi(npolyphi-i)*spol**i
  enddo

  do i=npolyphi,1,-1
    dPhi_dpsi = dPhi_dpsi + dble(i)*polyphi(npolyphi-i)*spol**(i-1)
  enddo
  dPhi_dpsi = dPhi_dpsi/psi_sep

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
  use phielec_of_psi_mod, only : npolyphi,polydens,polytemp
  use field_eq_mod,       only : psi_sep,psi_axis
!
  implicit none
!
  integer :: i
  double precision :: spol
  double precision, intent(out) :: psi,dens,temp,ddens,dtemp
!
  spol = (psi-psi_axis)/(psi_sep-psi_axis)

  dens = 0d0
  temp = 0d0
  ddens = 0d0
  dtemp = 0d0

  do i=npolyphi,0,-1
    dens = dens + polydens(npolyphi-i)*spol**i
    temp = temp + polytemp(npolyphi-i)*spol**i
  enddo

  do i=npolyphi,1,-1
    ddens = ddens + dble(i)*polydens(npolyphi-i)*spol**(i-1)
    dtemp = dtemp + dble(i)*polytemp(npolyphi-i)*spol**(i-1)
  enddo
  ddens = ddens/psi_sep
  dtemp = dtemp/psi_sep

end subroutine denstemp_of_psi
