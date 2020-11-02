!
module phielec_of_psi_mod
  integer, parameter :: npolyphi=1
  double precision, dimension(0:npolyphi) :: polyphi
end module phielec_of_psi_mod
!
!ccccccccc
!
subroutine phielec_of_psi(psi,phi_elec,dPhi_dpsi)
  !
    use phielec_of_psi_mod, only : npolyphi,polyphi
    use field_eq_mod,       only : psi_sep
  !
    implicit none
  !
    logical :: prop=.true.
    integer :: i
    double precision :: psi,phi_elec,dPhi_dpsi
    double precision :: ampl
  !
    ampl=1.12d0
  !
    if(prop) then
      prop=.false.
      polyphi(0)=0.d0
      polyphi(1)=-ampl/psi_sep
    endif
  !
    phi_elec=0.d0
    dPhi_dpsi=0.d0
  !
    do i=npolyphi,0,-1
      phi_elec=polyphi(i)+phi_elec*psi
    enddo
  !
    do i=npolyphi,1,-1
      dPhi_dpsi=polyphi(i)*dble(i)+dPhi_dpsi*psi
    enddo
  !
end subroutine phielec_of_psi
