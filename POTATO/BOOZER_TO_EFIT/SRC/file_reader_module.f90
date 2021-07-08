module read_file_module
!  
  use mag_surf_module, only : mpol, ntor, nsurf, nmodes, nper, rmnc, rmns, &
                              zmnc, zmns, npoi, nder, mag_const, pi, flux, &
                              rho_tor, iota_arr, R_mn_c, R_mn_s, Z_mn_c, Z_mn_s, &
                              bmn_c,bmn_s,mvec
                                
                                
!  
  implicit none
!
  contains
! 
  subroutine read_boozer_file (unitNum, fileName, s, iota, Jpol, Itor, pprime, sqrtg, m_pol, n_tor, &
                               R_mn_c, R_mn_s, Z_mn_c, Z_mn_s, psi, B_tor, Rzero)
!
    integer, intent (in) :: unitNum
    character (len=*), intent (in) :: fileName
    double precision :: s_0, iota_0, Jpol_0, Itor_0, pprime_0, sqrtg_0, sum_s, s_s, hs, s_start, s_end, integral
    double precision :: Rzero,dummy,rmnc_loc,rmns_loc,zmnc_loc,zmns_loc,bmnc_loc,bmns_loc
    double precision, dimension(:),   allocatable :: s, iota, Jpol, Itor, pprime, sqrtg, m_pol, n_tor, psi
    double precision, dimension(:),   allocatable, intent(out) :: B_tor
    double precision, dimension(:,:), allocatable, intent(out) :: R_mn_c, R_mn_s, Z_mn_c, Z_mn_s
    integer :: i, j, is, int_step
    double precision :: result, abserr, resabs, resasc, f, h_iota, start_iota, end_iota, sum_int, iota_test, iota_need, s_test
    double precision, dimension(:), allocatable :: iota_i
    integer :: ibeg,iend,ierr,nmodes_all,m,n,jn0
    double precision, dimension(0:nder,npoi) :: coef
!
    int_step = 100
!    
    open (unit=unitNum, file=fileName)
    read(unitNum,*)
    read(unitNum,*)
    read(unitNum,*)
    read(unitNum,*)
    read(unitNum,*)
    read(unitNum,*) mpol,ntor,nsurf,nper,flux,dummy,Rzero
    nmodes=mpol+1
!    nmodes_all=(mpol+1)*(2*ntor+1)
    nmodes_all=(mpol+1)*(2*ntor+1)-ntor
!    
    allocate(iota_i(0:int_step))
    allocate(s(0:nsurf), iota(0:nsurf), Jpol(nsurf), Itor(nsurf), pprime(nsurf), sqrtg(nsurf), psi(0:nsurf))
    allocate(B_tor(nsurf), iota_arr(nsurf))
    allocate(m_pol(nmodes), n_tor(nmodes))
    allocate(rmnc(nmodes), rmns(nmodes), zmnc(nmodes), zmns(nmodes))
    allocate(R_mn_c(nsurf, nmodes), R_mn_s(nsurf, nmodes), Z_mn_c(nsurf, nmodes), Z_mn_s(nsurf, nmodes))
!
    allocate(bmn_c(nsurf,0:mpol,-ntor:ntor),bmn_s(nsurf,0:mpol,-ntor:ntor),mvec(0:mpol))
!
    do i=0,mpol
      mvec(i)=i
    enddo
!
    pi = 4.d0 * atan(1.d0)
!
    s(0) = 0.d0
    psi(0) = 0.d0
    do i = 1, nsurf
      read(unitNum,*)
      read(unitNum,*)
      read(unitNum,*) s_0, iota_0, Jpol_0, Itor_0, pprime_0, sqrtg_0
      read(unitNum,*)
      jn0=0
      do j = 1, nmodes_all
!        read(unitNum,*) m_pol(j), n_tor(j), rmnc(j), rmns(j), zmnc(j), zmns(j)
!        R_mn_c(i, j) = rmnc(j)
!        R_mn_s(i, j) = rmns(j)
!        Z_mn_c(i, j) = zmnc(j)
!        Z_mn_s(i, j) = zmns(j)
        read(unitNum,*) m,n,rmnc_loc,rmns_loc,zmnc_loc,zmns_loc,dummy,dummy,bmnc_loc,bmns_loc
        if(n.eq.0) then
          jn0=jn0+1
          m_pol(jn0) = m
          n_tor(jn0) = n
          R_mn_c(i, jn0) = rmnc_loc
          R_mn_s(i, jn0) = rmns_loc
          Z_mn_c(i, jn0) = zmnc_loc
          Z_mn_s(i, jn0) = zmns_loc
        else
          bmn_c(i,m,n) = bmnc_loc
          bmn_s(i,m,n) = bmns_loc
        endif
      enddo
      s(i) = s_0
      iota(i) = iota_0
      Jpol(i) = Jpol_0
      Itor(i) = Itor_0
      pprime(i) = pprime_0
      sqrtg(i) = sqrtg_0
      B_tor(i) = - (mag_const / (2.d0 * pi)) * Jpol_0
    enddo
    close(unitNum)
    
    allocate(rho_tor(0:nsurf))   !<=NEW
    rho_tor=sqrt(s)              !<=NEW
    
    call intsrc(nsurf,npoi,0.d0,s(1),1.d0/dfloat(nsurf),ibeg,iend)
    
    call plag_coeff(npoi,nder,0.d0,s(ibeg:iend),coef)
    
    iota(0) = sum(iota(1:npoi)*coef(0,1:npoi))
!        
    call integrator_stw(nsurf-1,s(1),s(nsurf),iota(1:nsurf),psi(1:nsurf),ierr)
!
    call plag_coeff(npoi,nder,0.d0,s(1:npoi),coef)
!
    psi(0) = sum(psi(1:npoi)*coef(0,1:npoi))
    psi=psi-psi(0)
    flux=flux/(2.d0*pi)
    psi=psi*flux
!
!    start_iota = iota(0)
!    psi(0) = 0.d0 
!    do is = 1, nsurf
!      end_iota = iota(is)
!      h_iota = (end_iota - start_iota) / dfloat(int_step)
!      
!      do i = 1, int_step
!        iota_test = start_iota + h_iota * dfloat(i)
!        call intsrc(nsurf,npoi,iota_test,iota(0),1.d0/dfloat(nsurf),ibeg,iend)
!        call plag_coeff(npoi,nder,iota_test,iota(ibeg:iend),coef)
!        iota_i(i) = sum(iota(ibeg:iend)*coef(0,1:npoi))
!      enddo
!!      
!      call integrator(int_step,s(0),s(is),iota_i(0:int_step),integral)
!!      
!      psi(is) = integral 
!!
!do is = 0, nsurf
!write(1234,*) s(is),psi(is),iota(is)
!    enddo
        
    return
  end subroutine read_boozer_file
!
end module read_file_module
! 
