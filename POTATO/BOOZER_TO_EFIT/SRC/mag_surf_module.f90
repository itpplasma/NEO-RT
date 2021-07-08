!
  module mag_surf_module
!  
    implicit none
!    
    integer :: mpol, ntor, nsurf, nmodes, nper, inp_label, iunit
    double precision, dimension(:), allocatable :: rmnc, rmns, zmnc, zmns
    double precision, dimension(:), allocatable :: iota_arr,rho_tor
    integer, parameter :: npoi = 4
    integer, parameter :: nder = 0
    double precision, parameter :: mag_const = 1.25663706212d-6
    double precision :: pi, flux
    integer, dimension(:), allocatable :: mvec
    double precision, dimension(:), allocatable :: s_b, iota_b, Jpol_b, Itor_b, pprime_b, sqrtg_b, m_pol_b, n_tor_b, psi, B_tor!, rmnc_b, rmns_b, zmnc_b, zmns_b
    double precision, dimension(:,:),   allocatable :: R_mn_c, R_mn_s, Z_mn_c, Z_mn_s
    double precision, dimension(:,:,:), allocatable :: bmn_c, bmn_s
!    
  end module mag_surf_module
!  
!<= I removed file, where this module was described, and put it here (also removed unused variables)
!
  module pert_mag_field_mod
!
    implicit none
!    
    character(1) :: smb
    integer :: nr,nz,ir,iz,n_tor
    double precision :: rmin,zmin,hr,hz
    double precision, dimension(:), allocatable :: Rdata,Zdata,Rdata_src,Zdata_src
!    
  end module pert_mag_field_mod
!  
