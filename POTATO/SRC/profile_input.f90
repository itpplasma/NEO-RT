  subroutine phielec_of_psi(psi,phi_elec,dPhi_dpsi)
!
! Profile of normalized electrostatic potential $\hat\Phi$ and its derivative
! as functions of (dimensional) poloidal flux $\psi$
!
  use phielec_of_psi_mod, only : npolyphi,polyphi,int_powers
  use potato_profile_data_mod, only : tabulated_profiles,evaluate_profile_table
  use field_eq_mod,       only : psi_sep,psi_axis
!
  implicit none
!
  logical :: prop=.true.
  integer :: i
  double precision spol,pw(0:npolyphi)
  double precision, intent(in) :: psi
  double precision, intent(out) :: phi_elec,dPhi_dpsi

  spol = (psi-psi_axis)/(psi_sep-psi_axis)

  if (tabulated_profiles) then
    call evaluate_profile_table(spol,3,phi_elec,dPhi_dpsi)
    dPhi_dpsi = dPhi_dpsi/(psi_sep-psi_axis)
    return
  endif

  call int_powers(spol,pw)

  phi_elec=0.d0
  dPhi_dpsi=0.d0

! novector: with the powers precomputed gfortran would vectorize the products
! and add them separately; scalar, it contracts each term to a fused
! multiply-add exactly as it did for the former spol**i, so results are unchanged.
!GCC$ novector
  do i=npolyphi,0,-1
    phi_elec = phi_elec + polyphi(npolyphi-i)*pw(i)
  enddo

!GCC$ novector
  do i=npolyphi,1,-1
    dPhi_dpsi = dPhi_dpsi + dble(i)*polyphi(npolyphi-i)*pw(i-1)
  enddo
  dPhi_dpsi = dPhi_dpsi/(psi_sep-psi_axis)

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
  use phielec_of_psi_mod, only : npolyphi,polydens,polytemp,int_powers
  use potato_profile_data_mod, only : tabulated_profiles,evaluate_profile_table
  use field_eq_mod,       only : psi_sep,psi_axis
!
  implicit none
!
  integer :: i
  double precision :: spol,pw(0:npolyphi)
  double precision, intent(in) :: psi
  double precision, intent(out) :: dens,temp,ddens,dtemp
!
  spol = (psi-psi_axis)/(psi_sep-psi_axis)

  if (tabulated_profiles) then
    call evaluate_profile_table(spol,1,dens,ddens)
    call evaluate_profile_table(spol,2,temp,dtemp)
    ddens = ddens/(psi_sep-psi_axis)
    dtemp = dtemp/(psi_sep-psi_axis)
    return
  endif

  call int_powers(spol,pw)

  dens = 0d0
  temp = 0d0
  ddens = 0d0
  dtemp = 0d0

!GCC$ novector
  do i=npolyphi,0,-1
    dens = dens + polydens(npolyphi-i)*pw(i)
    temp = temp + polytemp(npolyphi-i)*pw(i)
  enddo

!GCC$ novector
  do i=npolyphi,1,-1
    ddens = ddens + dble(i)*polydens(npolyphi-i)*pw(i-1)
    dtemp = dtemp + dble(i)*polytemp(npolyphi-i)*pw(i-1)
  enddo
  ddens = ddens/(psi_sep-psi_axis)
  dtemp = dtemp/(psi_sep-psi_axis)

end subroutine denstemp_of_psi
