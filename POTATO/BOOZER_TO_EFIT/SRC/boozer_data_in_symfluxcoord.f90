!
  subroutine boozer_data_in_symfluxcoord(s, theta, phi, R, Z, dR_ds, dR_dt, dZ_ds, dZ_dt)
!
  use mag_surf_module,  only : npoi,nsurf,s_b,m_pol_b,n_tor_b,rho_tor
  use read_file_module, only : R_mn_c,R_mn_s,Z_mn_c,Z_mn_s
!
  implicit none
!
    integer, parameter :: nder=1
    integer :: ibeg,iend
    double precision, intent (in) :: s, theta, phi
    double precision, dimension(0:nder,npoi) :: coef
    double precision, intent (out) :: R, Z, dR_ds, dR_dt, dZ_ds, dZ_dt
    double precision :: rho !<=NEW
!
    call intsrc(nsurf,npoi,s,0.d0,1.d0/dfloat(nsurf),ibeg,iend)
    rho=sqrt(s)
    call plag_coeff(npoi,nder,rho,rho_tor(ibeg:iend),coef)
!
    R = sum(cos(m_pol_b*theta + n_tor_b*phi)*matmul(coef(0,:), R_mn_c(ibeg:iend, :))) & 
          + sum(sin(m_pol_b*theta + n_tor_b*phi)*matmul(coef(0,:), R_mn_s(ibeg:iend, :)))
          
    Z = sum(cos(m_pol_b*theta + n_tor_b*phi)*matmul(coef(0,:), Z_mn_c(ibeg:iend, :))) & 
          + sum(sin(m_pol_b*theta + n_tor_b*phi)*matmul(coef(0,:), Z_mn_s(ibeg:iend, :)))
          
    dR_ds = sum(cos(m_pol_b*theta + n_tor_b*phi)*matmul(coef(1,:), R_mn_c(ibeg:iend, :))) & 
          + sum(sin(m_pol_b*theta + n_tor_b*phi)*matmul(coef(1,:), R_mn_s(ibeg:iend, :)))
          
    dZ_ds = sum(cos(m_pol_b*theta + n_tor_b*phi)*matmul(coef(1,:), Z_mn_c(ibeg:iend, :))) & 
              + sum(sin(m_pol_b*theta + n_tor_b*phi)*matmul(coef(1,:), Z_mn_s(ibeg:iend, :)))  
    dR_ds = dR_ds*0.5d0/rho !<=NEW
    dZ_ds = dZ_ds*0.5d0/rho !<=NEW
       
    dR_dt = -sum(m_pol_b*sin(m_pol_b*theta + n_tor_b*phi)*matmul(coef(0,:), R_mn_c(ibeg:iend, :))) & 
               + sum(m_pol_b*cos(m_pol_b*theta + n_tor_b*phi)*matmul(coef(0,:), R_mn_s(ibeg:iend, :)))
         
    dZ_dt = -sum(m_pol_b*sin(m_pol_b*theta + n_tor_b*phi)*matmul(coef(0,:), Z_mn_c(ibeg:iend, :))) & 
               + sum(m_pol_b*cos(m_pol_b*theta + n_tor_b*phi)*matmul(coef(0,:), Z_mn_s(ibeg:iend, :)))  

!
  end subroutine boozer_data_in_symfluxcoord
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine pert_modB(s,theta,amplbn)
!
  use mag_surf_module, only : ntor,npoi,nsurf,s_b,rho_tor,bmn_c,bmn_s,mvec
!
  implicit none
!
  integer, parameter :: nder=0
  double complex, parameter :: imun=(0.d0,1.d0)
  integer            :: ibeg,iend,n
  double precision   :: s,theta,rho
  double complex,   dimension(ntor) :: amplbn
  double precision, dimension(0:nder,npoi) :: coef
!
  rho=sqrt(s)
!
  call intsrc(nsurf,npoi,s,0.d0,1.d0/dfloat(nsurf),ibeg,iend)
!
  call plag_coeff(npoi,nder,rho,rho_tor(ibeg:iend),coef)
!
  do n=1,ntor
    amplbn(n)=sum(matmul(coef(0,:),bmn_c(ibeg:iend,:,n)-imun*bmn_s(ibeg:iend,:,n))*exp(imun*mvec*theta))   &
             +sum(matmul(coef(0,:),bmn_c(ibeg:iend,:,-n)+imun*bmn_s(ibeg:iend,:,-n))*exp(-imun*mvec*theta))
  enddo
!
  end subroutine pert_modB
