program boozer_to_efit
!
  use read_file_module,   only : read_boozer_file
!  
  use mag_surf_module,    only : mpol, ntor, nsurf, nmodes, nper, inp_label, iunit, &
                                 pi, s_b, iota_b, Jpol_b, Itor_b, pprime_b, sqrtg_b, &
                                 m_pol_b, n_tor_b, psi, B_tor, R_mn_c, R_mn_s, Z_mn_c, Z_mn_s, &
                                 bmn_c,bmn_s
!                              
  use pert_mag_field_mod, only : smb, nr, nz, rmin, zmin, hr, hz, &
                                 Rdata, Zdata, Rdata_src, Zdata_src
!                                 
  use signal_array_mod,   only : signal_arr, s_guess, theta_guess, s_grid, theta_grid, &
                                 R_mag_data, Z_mag_data
!                               
  use extender_mod,       only : s_bou !<=input for extender
!  
  implicit none
!
  logical :: err_newt
  logical :: err_grad
  logical :: filter,prop,outside
!
  integer, parameter :: nplag=4, nder=1
!  integer, parameter :: n = 10000
  integer, parameter :: ns = 10000
  integer, parameter :: grad_iter = 10000
  integer, parameter :: niter = 10
!
  integer :: i,j,nt,modfactor,t,m,n
  integer :: ir_beg, iz_beg, iter, ibeg,iend
  integer :: nwind=4
  integer :: nwEQD,nhEQD
!
  double precision, parameter :: eps_grad = 1d-4 !1d-9
  double precision, parameter :: eps_newt = 1d-9
  double precision, parameter :: sigma = 1.d0, rpowmin=1.d-6
!
  double precision :: s_test, theta_test, phi, htheta, R, Z, dR_ds, dR_dt, dZ_ds, dZ_dt
  double precision :: s_new, theta_new, R_src, Z_src
  double precision :: rr, zz, hs, ht, s, theta, s_0, theta_0, a_grad_s, a_grad_t, det_J
  double precision :: power,hmin,dldsmax,smax,del_s,del_t,dist_pol,dist_pol0,det_J0, & !<=NEW
                      s_out,theta_out,rmax,zmax,psiAxis,psiSep,bt0,rzero,psihat        !<=NEW
  double precision, dimension(0:nder,nplag)     :: coef
  integer,          dimension(:,:),   allocatable :: idummy
  double precision, dimension(:),     allocatable :: rpower, arr_in, arr_out, fpol, dummy1d
  double precision, dimension(:,:),   allocatable :: psi_grid
  double precision, dimension(:,:,:), allocatable :: bmod_n_re,bmod_n_im
  double complex,   dimension(:),     allocatable :: amplbn
!
  iunit=71
!  
  call read_boozer_file(iunit, 'boozfile_to_plot.bc', s_b, iota_b, Jpol_b, Itor_b, &
                        pprime_b, sqrtg_b, m_pol_b, n_tor_b,                       & 
                        R_mn_c, R_mn_s, Z_mn_c, Z_mn_s, psi, B_tor, Rzero)
!
!*********************
! Filter Fourier amplitudes:
!
!  filter=.false.
  filter=.true.
!
  if(filter) then
    allocate(rpower(nsurf),arr_in(nsurf),arr_out(nsurf))
!
    do m=1,nmodes
      power=0.5d0*dble(m_pol_b(m))
      rpower=s_b(1:nsurf)**power
      arr_in=R_mn_c(:,m)/max(rpower,rpowmin)
!
      call gauss_filter(nsurf,sigma,arr_in,arr_out)
!
      R_mn_c(:,m)=arr_out*rpower
      arr_in=R_mn_s(:,m)/max(rpower,rpowmin)
!
      call gauss_filter(nsurf,sigma,arr_in,arr_out)
!
      R_mn_s(:,m)=arr_out*rpower
      arr_in=Z_mn_c(:,m)/max(rpower,rpowmin)
!
      call gauss_filter(nsurf,sigma,arr_in,arr_out)
!
      Z_mn_c(:,m)=arr_out*rpower
      arr_in=Z_mn_s(:,m)/max(rpower,rpowmin)
!
      call gauss_filter(nsurf,sigma,arr_in,arr_out)
!
      Z_mn_s(:,m)=arr_out*rpower
    enddo
!
    do m=0,mpol
      power=0.5d0*dble(m)
      rpower=s_b(1:nsurf)**power
      do n=-ntor,ntor
        if(n.eq.0) cycle
        arr_in=bmn_c(:,m,n)/max(rpower,rpowmin)
!
        call gauss_filter(nsurf,sigma,arr_in,arr_out)
!
        bmn_c(:,m,n)=arr_out*rpower
        arr_in=bmn_s(:,m,n)/max(rpower,rpowmin)
!
        call gauss_filter(nsurf,sigma,arr_in,arr_out)
!
        bmn_s(:,m,n)=arr_out*rpower
      enddo
    enddo
!
    deallocate(rpower,arr_in,arr_out)
    print *,'Input filtered: sigma = ',sigma,' rpowmin = ',rpowmin
  endif
!
! End filter Fourier amplitudes
!*********************
!
! Determine box size
!
  phi=0.d0
!
  if(.false.) then
!
    open(8,file='out1I.dat')
      read (8,'(1A,2(1x,I4))')    smb, nr, nz
      read (8,'(1A,2(1x,e14.6))') smb, rmin,zmin
      read (8,'(1A,2(1x,e14.6))') smb, hr, hz
    close(8)
!
    rmin = rmin / 100.d0
    zmin = zmin / 100.d0
    hr = hr / 100.d0
    hz = hz / 100.d0  
!
  else
!
    nr=200
    nz=200
!
    nt=100
    htheta=2.d0*pi/dble(nt)
    s=1.d0
!
    do i=1,nt
      theta=htheta*dble(i)
!
      call boozer_data_in_symfluxcoord(s, theta, phi, R, Z, dR_ds, dR_dt, dZ_ds, dZ_dt)
!
      if(i.eq.1) then
        rmin=R
        rmax=R
        zmin=Z
        zmax=Z
      else
        rmin=min(R,rmin)
        rmax=max(R,rmax)
        zmin=min(Z,zmin)
        zmax=max(Z,zmax)
      endif
    enddo
!
    hr=0.1d0*(rmax-rmin)
    rmin=max(0.1d0*rmin,rmin-hr)
    rmax=rmax+hr
    hz=0.1d0*(zmax-zmin)
    zmin=zmin-hz
    zmax=zmax+hz
!
    hr=(rmax-rmin)/dble(nr)    
    hz=(zmax-zmin)/dble(nz)    
  endif
!
! End determine box size
!
  allocate(Rdata(0:nr))
  allocate(Zdata(0:nz))
  allocate(Rdata_src(0:nr))
  allocate(Zdata_src(0:nz))
!
  do t=0,nr
    rr = rmin + hr*dfloat(t)
    Rdata(t) = rr
  enddo  
!
  do m=0,nz
    zz = zmin + hz*dfloat(m)
    Zdata(m) = zz
  enddo  
!
!-----------Search grid----------------
!
  do i=0,nr
    R_src = Rdata(i)! - hr/2.d0
    Rdata_src(i) = R_src 
  enddo  
!  
  do j=0,nz
    Z_src = Zdata(j)! - hz/2.d0
    Zdata_src(j) = Z_src
  enddo  
!
!  hs = 0.97d0/dfloat(ns)
!  ht = 2.d0*pi/dfloat(n)
!
  allocate(signal_arr(0:nr,0:nz))
  allocate(s_guess(0:nr,0:nz),theta_guess(0:nr,0:nz))
  allocate(s_grid(0:nr,0:nz), theta_grid(0:nr,0:nz))
  allocate(R_mag_data(0:nr, 0:nz), Z_mag_data(0:nr, 0:nz))
  allocate(idummy(0:nr,0:nz),amplbn(ntor))
  allocate(psi_grid(0:nr,0:nz),bmod_n_re(0:nr,0:nz,ntor),bmod_n_im(0:nr,0:nz,ntor))
!
  do i=0,nr
    do j=0,nz
      signal_arr(i,j) = 0
      s_guess(i,j) = -1.d0
      theta_guess(i,j) = -1.d0
    enddo
  enddo  
!

  modfactor=30
  nt=mpol*modfactor
  htheta=(2.d0*pi)/dble(nt)
  
  hmin=min(hr,hz)/2.d0 
  smax=1.d0
  s=(hmin/rmin)**2
  prop=.true.
  det_J=1.d0
!
  do while(s.lt.smax.and.det_J.gt.0.d0)
    dldsmax=0.d0
    theta=0.d0
!
    do while(theta.lt.2.d0*pi)
!      
      call boozer_data_in_symfluxcoord(s, theta, phi, R, Z, dR_ds, dR_dt, dZ_ds, dZ_dt)
!      
      det_J = (dR_ds*dZ_dt - dR_dt*dZ_ds)
      if(prop) then
        prop=.false.
        det_J0=det_J
      endif
      det_J=det_J/det_J0
      if(det_J.lt.0.1d0) then
        dldsmax=1.d0
        exit
      endif
!
      ir_beg=nint((R-Rdata_src(0))/hr) 
      iz_beg=nint((Z-Zdata_src(0))/hz)
      if(ir_beg.ge.0.and.ir_beg.le.nr.and.iz_beg.ge.0.and.iz_beg.le.nz) then
        signal_arr(ir_beg,iz_beg) = 1
        theta_guess(ir_beg,iz_beg) = theta
        s_guess(ir_beg,iz_beg) = s
      else
        print *,'Warning: inner point outside computation box, R = ',R,'  Z = ',Z
      endif
      theta=theta+hmin/sqrt(dR_dt**2+dZ_dt**2)
      dldsmax=max(dldsmax,sqrt(dR_ds**2+dZ_ds**2))
    enddo
!
    if(theta.gt.2.d0*pi) s_bou=s
    s=s+hmin/dldsmax
  enddo
!
  print *,'Real last closed surface s = ',s_bou
!
  idummy=signal_arr
!
  do t=0,nr
    do m=0,nz
!
      call extender(Rdata(t),Zdata(m),outside,s_out,theta_out)
!
      if(outside) then
        signal_arr(t,m)=0
        s_grid(t,m) = s_out
        theta_grid(t,m) = theta_out
      endif
!
    enddo  
  enddo  
!
!-------------------------Gradient descent-------------------------
!
  do t=0,nr
    do m=0,nz        
      if(signal_arr(t,m) == 0) then 
        cycle
      endif    
      s_0 = s_guess(t,m)
      theta_0 = theta_guess(t,m)   
!
      err_grad = .true.
!
      a_grad_s = 1.d0
      a_grad_t = 1.d0
!      
      call boozer_data_in_symfluxcoord(s_0, theta_0, phi, R, Z, dR_ds, dR_dt, dZ_ds, dZ_dt)
!
      dist_pol0=abs(R - Rdata(t))+abs(Z - Zdata(m))
      det_J = (dR_ds*dZ_dt - dR_dt*dZ_ds)        
      del_s=(dR_dt*(Z-Zdata(m)) - dZ_dt*(R-Rdata(t)))/det_J
      del_t=(dZ_ds*(R-Rdata(t)) - dR_ds*(Z-Zdata(m)))/det_J  
!
      do iter=1,grad_iter
        s_new = abs(s_0 + a_grad_s*del_s)
        theta_new = theta_0 + a_grad_t*del_t
!      
        call boozer_data_in_symfluxcoord(s_new, theta_new, phi, R, Z, dR_ds, dR_dt, dZ_ds, dZ_dt)
!        
        if(abs(R - Rdata(t)) .lt. eps_grad .and. abs(Z - Zdata(m)) .lt. eps_grad) then
          err_grad = .false.
          exit
        endif
        dist_pol=abs(R - Rdata(t))+abs(Z - Zdata(m))
        if(dist_pol.lt.dist_pol0.and.s_new.lt.1.d0) then
          det_J = (dR_ds*dZ_dt - dR_dt*dZ_ds)        
          del_s=(dR_dt*(Z-Zdata(m)) - dZ_dt*(R-Rdata(t)))/det_J
          del_t=(dZ_ds*(R-Rdata(t)) - dR_ds*(Z-Zdata(m)))/det_J  
          dist_pol0=dist_pol
          s_0 = s_new 
          theta_0 = theta_new
          a_grad_s=min(1.d0,2.d0*a_grad_s)
          a_grad_t=a_grad_s
        else
          a_grad_s=0.5d0*a_grad_s
          a_grad_t=a_grad_s
        endif
      enddo  
      if(s_0 .gt. 1.d0 .or. s_0 .lt. 0.d0) then
        s_grid(t,m) = -1.d0
        theta_grid(t,m) = -1.d0
        cycle
      endif
      if(err_grad) print *,'error in Grad',t,m
         
!      
!------------------------Newton----------------------------
!        
      err_newt = .true.
!      
      do iter=1,niter
!      
        call boozer_data_in_symfluxcoord(s_0, theta_0, phi, R, Z, dR_ds, dR_dt, dZ_ds, dZ_dt)
!        
        det_J = (dR_ds*dZ_dt - dR_dt*dZ_ds)
        s_new = s_0 + (dR_dt*(Z-Zdata(m)) - dZ_dt*(R-Rdata(t)))/det_J
        theta_new = theta_0 + (dZ_ds*(R-Rdata(t)) - dR_ds*(Z-Zdata(m)))/det_J         
        if(abs(R - Rdata(t)) .lt. eps_newt .and. abs(Z - Zdata(m)) .lt. eps_newt) then
          err_newt = .false.
          exit
        endif 
        s_0 = abs(s_new)
        theta_0 = theta_new      
      enddo      
      if(s_0 .gt. 1.d0 .or. s_0 .lt. 0.d0) then
        s_grid(t,m) = -1.d0
        theta_grid(t,m) = -1.d0
        cycle
      endif
      if(err_newt) print *,'error in Newton',t,m
      s_grid(t,m) = s_new
      theta_grid(t,m) = theta_new     
    enddo
  enddo  
!  
!----------------------------------------------------------
!
! Compute psi-function and perturbation of mod-B
!
!
  do i=0,nr 
    do j=0,nz 
      s_0=s_grid(i,j)
!
      call binsrc(s_b(0:nsurf),0,nsurf,s_0,ibeg)
!
      ibeg=max(0,ibeg-nplag/2)
      iend=ibeg+nplag-1
      if(iend.gt.nsurf) then
        iend=nsurf
        ibeg=iend+1-nplag
      endif
!
      call plag_coeff(nplag,nder,s_0,s_b(ibeg:iend),coef)
!
      psi_grid(i,j)=sum(coef(0,:)*psi(ibeg:iend))
!
      if(s_0.gt.s_b(nsurf)) then
!
        call plag_coeff(nplag,nder,s_b(nsurf),s_b(ibeg:iend),coef)
!
        psi_grid(i,j)=sum(coef(0,:)*psi(ibeg:iend))                 &
                     +sum(coef(1,:)*psi(ibeg:iend))*(s_0-s_b(nsurf))
      endif
!
      theta_0=theta_grid(i,j)
!
      call pert_modB(s_0,theta_0,amplbn)
!
      bmod_n_re(i,j,:)=dble(amplbn)
      bmod_n_im(i,j,:)=dimag(amplbn)
    enddo
  enddo   
!  
! End compute psi-function and perturbation of mod-B
!
!----------------------------------------------------------
!
! Formation of data for EFIT file and for B-mod perturbation file
!
! Number of grid points per box width:
!
  nwEQD=nr+1
!
! Number of grid points per box height:
!
  nhEQD=nz+1
!
! poloidal flux at the axis:
!
  psiAxis=psi(0)
!
! poloidal flux at the "separatrix":
!
  call binsrc(s_b(0:nsurf),0,nsurf,s_bou,ibeg)
!
  ibeg=max(0,ibeg-nplag/2)
  iend=ibeg+nplag-1
  if(iend.gt.nsurf) then
    iend=nsurf
    ibeg=iend+1-nplag
  endif
!
  call plag_coeff(nplag,nder,s_bou,s_b(ibeg:iend),coef)
!
  psiSep=sum(coef(0,:)*psi(ibeg:iend))
!
! physical toroidal magnetic field component at the axis:
!
  call plag_coeff(nplag,nder,0.d0,s_b(1:nplag),coef)
!
  bt0=sum(coef(0,:)*B_tor(1:nplag))/rzero
!
! covariant toroidal componet of the magnetic field:
!
  allocate(fpol(nwEQD),dummy1d(0:nsurf))
  dummy1d=psi/(psiSep-psiAxis)
!
  do i=1,nwEQD
    psihat=dble(i-1)/dble(nr)
!
    call binsrc(dummy1d(0:nsurf),0,nsurf,psihat,ibeg)
!
    ibeg=max(1,ibeg-nplag/2)
    iend=ibeg+nplag-1
    if(iend.gt.nsurf) then
      iend=nsurf
      ibeg=iend+1-nplag
    endif
!
    call plag_coeff(nplag,nder,psihat,dummy1d(ibeg:iend),coef)
!
    fpol(i)=sum(coef(0,:)*B_tor(ibeg:iend))
  enddo
!
! write data to EQDSK file:
!
  call write_eqfile(nwEQD,nhEQD,psiAxis,psiSep,bt0,rzero,fpol,  &
                    Rdata(0:nr),Zdata(0:nz),psi_grid(0:nr,0:nz))
!
! write perturbed mod-B files:
!
  do n=1,ntor
    open(iunit,form='unformatted',file='bmod_n.'//char(48+n))
    write(iunit) nwEQD,nhEQD
    write(iunit) 1.d2*Rdata,1.d2*Zdata
    write(iunit) 1.d4*bmod_n_re(:,:,n),1.d4*bmod_n_im(:,:,n)
    close(iunit)
  enddo
!
! End formation of data for EFIT file and for B-mod perturbation files
!
!----------------------------------------------------------
!
! Convexwall:
!
  nt=100
  htheta=2.d0*pi/dble(nt)
  open(iunit,file='convexwall.boz')
!
  do i=0,nt
    theta=htheta*dble(i)
!
    call boozer_data_in_symfluxcoord(s_bou, theta, phi, R, Z, dR_ds, dR_dt, dZ_ds, dZ_dt)
!
    write(iunit,*) 1.d2*R,1.d2*Z
  enddo
!
  close(iunit)
!
! End convexwall
!
!----------------------------------------------------------
!
! Data for plotting:
!
  do i=0,nr 
    do j=0,nz 
      if(s_grid(i,j) .eq. -1.d0 .or. theta_grid(i,j) .eq. -1.d0) then
        R_mag_data(i,j) = -0.d0
        Z_mag_data(i,j) = -0.d0
        cycle
      endif  
      s = s_grid(i,j)
      theta = theta_grid(i,j)
!
      call boozer_data_in_symfluxcoord(s, theta, phi, R, Z, dR_ds, dR_dt, dZ_ds, dZ_dt)
!
      R_mag_data(i,j) = R
      Z_mag_data(i,j) = Z
    enddo
  enddo   
!  
!  do i=0,nr 
!    do j=0,nz 
!      if(s_grid(i,j) .eq. -1.d0 .or. theta_grid(i,j) .eq. -1.d0) cycle
!      s = s_grid(i,j)
!      theta = theta_grid(i,j)
!!
!      call boozer_data_in_symfluxcoord(s, theta, phi, R, Z, dR_ds, dR_dt, dZ_ds, dZ_dt)
!!
!      write(1488, *) s, theta, R, Z, abs(R - R_mag_data(i,j)), abs(Z - Z_mag_data(i,j))
!    enddo
!    write(1488, *) ' '
!  enddo 
!
  open(100,file='PLOT/psi_of_R_Z.dat')
  do i=0,nr
    write(100,*) psi_grid(i,:)
  enddo
  close(100)
!
  do n=1,ntor
    open(100,file='PLOT/bmod_n_re.'//char(48+n)//'.dat')
    open(101,file='PLOT/bmod_n_im.'//char(48+n)//'.dat')
    do i=0,nr
      write(100,*) bmod_n_re(i,:,n)
      write(101,*) bmod_n_im(i,:,n)
    enddo
    close(100)
    close(101)
  enddo
!
!
!
  open(100,file='PLOT/s_of_R_Z.dat')
  do i=0,nr
    write(100,*) s_grid(i,:)
  enddo
  close(100)
!  
  open(101,file='PLOT/theta_of_R_Z.dat')
  do i=0,nr
    write(101,*) theta_grid(i,:)
  enddo
  close(101)
!  
  open(102,file='PLOT/R_of_s_theta.dat')
  do i=0,nr
    write(102,*) R_mag_data(i,:)
  enddo
  close(102)
!  
  open(104,file='PLOT/Z_of_s_theta.dat')
  do i=0,nr
    write(104,*) Z_mag_data(i,:)
  enddo
  close(104)
!
end program boozer_to_efit
!
