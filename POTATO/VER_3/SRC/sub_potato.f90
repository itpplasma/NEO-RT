!
  module vparzero_line_mod
    logical :: load_psi_line=.true.
    integer :: nvpline,npsib
    double precision :: rinb,zinb,psi_bound,sigpsi
    double precision, dimension(:),   allocatable :: s_line
    double precision, dimension(:,:), allocatable :: vpar_line,psi_line
  end module vparzero_line_mod
!
!------------------------------------------------------
!
  module orbit_dim_mod
    integer, parameter  :: neqm=5, next=3, ndim=neqm+next
    logical             :: write_orb=.false.
    integer             :: iunit1
    double precision    :: Rorb_max
  end module orbit_dim_mod
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine find_vparzero_line(nline,toten,perpinv,rmin,rmax,zmin,zmax,icase)
!
! Looks for the boundary of the "forbidden" region for particles with normalized total
! energy "toten" and perpendicular invariant "perpinv", i.e. a line "v_par^2=0" within
! the volume limites by the flux surface psi=psi_bound which touches one of the boundaries
! of the box given by "rmin<R<rmax", "zmin<Z<zmax".
! Identifies 4 cases:
! icase=1 - line is outside the volume, all the volume is in the forbidden region (no orbits)
! icase=2 - line crosses (twice) volume boundary
! icase=3 - line is fully inside the volume
! icase=4 - line is absent in the volume, all the volume allows for the orbits
!
  use vparzero_line_mod, only : load_psi_line,npsib,rinb,zinb,psi_bound,sigpsi,psi_line, &
                                nvpline,s_line,vpar_line
  use field_eq_mod,      only : psif,psi_axis,psi_sep,dpsidr,dpsidz
!
  implicit none
!
  double precision, parameter :: relerr=1d-12 !8
!
  logical :: prop
  integer :: nline,icase,i,j,nr,nz,ir,iz
  double precision :: toten,perpinv,rmin,rmax,zmin,zmax
  double precision :: err_dist,vpar2,R,Z,dist
  double precision :: vpar2in,vpar2out,gradmod,w
  double precision :: Rin,Rout,delRZ,h_r,h_z,det,delR,delZ
  double precision :: rss,zss,rse,zse,h_in,h_out,delth_rec
  double precision, dimension(2)            :: gradvpar2,y,yst,dely
  integer, dimension(1) :: indrmin
!
! Find R,Z and psi where psi-contour touches the box:
!
  if(load_psi_line) then
    load_psi_line=.false.
    sigpsi=sign(1.d0,psi_sep-psi_axis)
    psi_bound=psi_sep
!
    nr=nline
    nz=nline
    h_r=(rmax-rmin)/dble(nr)
    h_z=(zmax-zmin)/dble(nz)
!
    do ir=0,nr
      R=rmin+h_r*dble(ir)
      Z=zmin
!
      call vparzero_vec(R,Z,vpar2in,gradvpar2)
!
      if(sigpsi.gt.0.d0) then
        psi_bound=min(psi_bound,psif)
        if(psi_bound.eq.psif) then
          rss=R
          zss=Z
        endif
      else
        psi_bound=max(psi_bound,psif)
        if(psi_bound.eq.psif) then
          rss=R
          zss=Z
        endif
      endif
      Z=zmax
!
      call vparzero_vec(R,Z,vpar2in,gradvpar2)
!
      if(sigpsi.gt.0.d0) then
        psi_bound=min(psi_bound,psif)
        if(psi_bound.eq.psif) then
          rss=R
          zss=Z
        endif
      else
        psi_bound=max(psi_bound,psif)
        if(psi_bound.eq.psif) then
          rss=R
          zss=Z
        endif
      endif
    enddo
!
    do iz=0,nz
      Z=zmin+h_z*dble(iz)
      R=rmin
!
      call vparzero_vec(R,Z,vpar2in,gradvpar2)
!
      if(sigpsi.gt.0.d0) then
        psi_bound=min(psi_bound,psif)
        if(psi_bound.eq.psif) then
          rss=R
          zss=Z
        endif
      else
        psi_bound=max(psi_bound,psif)
        if(psi_bound.eq.psif) then
          rss=R
          zss=Z
        endif
      endif
      R=rmax
!
      call vparzero_vec(R,Z,vpar2in,gradvpar2)
!
      if(sigpsi.gt.0.d0) then
        psi_bound=min(psi_bound,psif)
        if(psi_bound.eq.psif) then
          rss=R
          zss=Z
        endif
      else
        psi_bound=max(psi_bound,psif)
        if(psi_bound.eq.psif) then
          rss=R
          zss=Z
        endif
      endif
    enddo
!
! End find R,Z and psi where psi-contour touches the box:
!
!plot the forbidden region
do ir=0,nr
  R=rmin+h_r*dble(ir)
  do iz=0,nz
    Z=zmin+h_z*dble(iz)
!
    call vparzero_vec(R,Z,vpar2in,gradvpar2)
!
    if(vpar2in.gt.0.d0) then
      write(1001,*) R,Z
    else
      write(1002,*) R,Z
    endif
!
    if(sigpsi*(psif-psi_bound).lt.0.d0) write(1003,*) R,Z
  enddo
enddo
!
! Find the contour line psi=psi_bound
!
    y(1)=rss
    y(2)=zss
    h_in=rmax/dble(nline)
    delth_rec=6.28d-3
!
    call choose_step(psi_cont,y,delth_rec,h_in,h_out)
!
    err_dist=rmax*relerr
    yst=y
!
    call level_set_step_2D(psi_cont, h_out, y)
    call newt_adjust_root_2d(psi_cont,psi_bound,err_dist,y)
!
    npsib=1
!
    do
!
      call level_set_step_2D(psi_cont, h_out, y)
      call newt_adjust_root_2d(psi_cont,psi_bound,err_dist,y)
!
      npsib=npsib+1
      if(sum((y-yst)**2).lt.h_out**2*1.25d0) exit
    enddo
!
    npsib=npsib+1
    allocate(psi_line(2,0:npsib))
    y=yst
!
    do i=0,npsib
      psi_line(:,i)=y
!
      call level_set_step_2D(psi_cont, h_out, y)
      call newt_adjust_root_2d(psi_cont,psi_bound,err_dist,y)
!
!write(101,*) psi_line(:,i)
    enddo
!
    psi_line(:,npsib)=yst
!
    indrmin=minloc(psi_line(1,:))
    Rinb=psi_line(1,indrmin(1))
    Zinb=psi_line(2,indrmin(1))
  endif
!
! End find the countor line psi=psi_bound
!
! Find the intersections of v_par^2=0 line with psi=psi_bound line
!
  prop=.false.
!
  call vparzero_vec(psi_line(1,0),psi_line(2,0),vpar2in,gradvpar2)
!
  do i=1,npsib
!
    call vparzero_vec(psi_line(1,i),psi_line(2,i),vpar2out,gradvpar2)
!
    if(vpar2in*vpar2out.lt.0.d0) then
      if(prop) then
! second intersection found:
        rss=psi_line(1,i)
        zss=psi_line(2,i)
        exit
      endif
! first intersection found:
      rse=psi_line(1,i)
      zse=psi_line(2,i)
      vpar2in=-vpar2in
      prop=.true.
    endif
  enddo
!
  if(prop) then
! there are (two) intersections, case 2
    R=rss
    Z=zss
!
    call adjust_crossing
!
    rss=R
    zss=Z
    R=rse
    Z=zse
!
    call adjust_crossing
!
    rse=R
    zse=Z
!
  else
    if(vpar2in.lt.0.d0) then
! case 1, no orbits
      print *,'find_vparzero_line: case 1, whole volume is forbidden'
      icase=1
      return
    endif
  endif
!
! End find the intersections of v_par^2=0 line with psi=psi_bound line
!
! At this point, if prop=.true. this is case 2, if .false., these are cases 3 or 4
! Process cases 3 and 4
!
  if(prop) then
    icase=2
    print *,'find_vparzero_line: case 2, forbidden boundary crosses volume boundary'
  else
! Find minimum value of parallel energy using the steepest decent method
    y(1)=sum(psi_line(1,1:npsib))/dble(npsib)
    y(2)=sum(psi_line(2,1:npsib))/dble(npsib)
!
    call vpar_cont(y,vpar2,gradvpar2)
!
    gradmod=sqrt(sum(gradvpar2**2))
    delth_rec=0.1d0
!
    call choose_step(vpar_cont,y,delth_rec,h_in,h_out)
!
    do
!
      call vpar_cont(y,vpar2,gradvpar2)
!
      dely=h_out*gradvpar2/gradmod
      y=y-dely
      if(sum(dely**2).lt.err_dist**2) exit
      if(sigpsi*(psif-psi_bound).gt.0.d0) exit
    enddo
!
    vpar2out=vpar2
! End find minimum value of parallel energy using the steepest decent method
!
    if(vpar2out.ge.0.d0) then
! case 4
      icase=4
      print *,'find_vparzero_line: case 4, no forbidden region'
      return
    endif
!
! case 3
    w=vpar2in/(vpar2in-vpar2out)
    y=(1.d0-w)*psi_line(:,0)+w*y
!
    call newt_adjust_root_2d(vpar_cont,0.d0,err_dist,y)
!
    rss=y(1)
    zss=y(2)
    rse=y(1)
    zse=y(2)
    icase=3
    print *,'find_vparzero_line: case 3, whole forbidden boundary is within the volume'
!
  endif
!
! End process cases 3 and 4
!
! Find the countor line v_par^2=0
!
  y(1)=rss
  y(2)=zss
  h_in=rmax/dble(nline)
  delth_rec=6.28d-3
!
  call choose_step(vpar_cont,y,delth_rec,h_in,h_out)
!
  yst(1)=rse
  yst(2)=zse
!
  call level_set_step_2D(vpar_cont, h_out, y)
  call newt_adjust_root_2d(vpar_cont,0.d0,err_dist,y)
!
  nvpline=1
!
  do
!
    call level_set_step_2D(vpar_cont, h_out, y)
    call newt_adjust_root_2d(vpar_cont,0.d0,err_dist,y)
!
    nvpline=nvpline+1
    if(sum((y-yst)**2).lt.h_out**2*1.25d0) exit
  enddo
!
  nvpline=nvpline+1
  if(allocated(s_line)) then
    deallocate(s_line,vpar_line)
  else
    allocate(s_line(0:nvpline),vpar_line(2,0:nvpline))
  endif
  y(1)=rss
  y(2)=zss
!
  do i=0,nvpline
    s_line(i)=dble(i)
    vpar_line(:,i)=y
!
    call level_set_step_2D(vpar_cont, h_out, y)
    call newt_adjust_root_2d(vpar_cont,0.d0,err_dist,y)
!
  enddo
!
  s_line(nvpline)=s_line(nvpline-1)                                 &
                 +sqrt(sum((yst-vpar_line(:,nvpline-1))**2))/h_out
  vpar_line(:,nvpline)=yst
!
! End find the countor line v_par^2=0
!
!------------
!
  contains
!
!------------
!
  subroutine vparzero_vec(R,Z,vpar2,gradvpar2)
!
  implicit none
!
  double precision :: bmod,sqrtg,phi_elec,R,Z,vpar2
  double precision, dimension(3) :: x,bder,hcovar,hctrvr,hcurl,derphi
  double precision, dimension(2) :: gradvpar2
!
  x(1)=R
  x(2)=0.d0
  x(3)=Z
!
  call magfie(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
!
  call elefie(x,phi_elec,derphi)
!
  vpar2=toten-phi_elec-perpinv*bmod
  gradvpar2(1)=-derphi(1)-perpinv*bmod*bder(1)
  gradvpar2(2)=-derphi(3)-perpinv*bmod*bder(3)
!
  end subroutine vparzero_vec
!
!------------
!
  subroutine psi_cont(y,f,df)
!
  implicit none
!
  double precision :: f
  double precision, dimension(2) :: y,df
!
  R=y(1)
  Z=y(2)
!
  call vparzero_vec(R,Z,vpar2,gradvpar2)
!
  f=psif
  df(1)=dpsidr
  df(2)=dpsidz
!
  end subroutine psi_cont
!
!------------
!
  subroutine vpar_cont(y,vpar2,gradvpar2)
!
  implicit none
!
  double precision :: vpar2
  double precision, dimension(2) :: y,gradvpar2
!
  call vparzero_vec(y(1),y(2),vpar2,gradvpar2)
!
  end subroutine vpar_cont
!
!------------
!
  subroutine choose_step(fun,y0,delth_rec,h_in,h_out)
!
  implicit none
!
  external :: fun
  double precision :: delth_rec,h_in,h_out,delth
  double precision, dimension(2) :: y0,y1,y2
!
  h_out=h_in
  delth=3.d0*delth_rec
  do while(delth.ge.2.d0*delth_rec)
    y1=y0
!
    call level_set_step_2D(fun, h_out, y1)
!
    y2=y1
!
    call level_set_step_2D(fun, h_out, y2)
!
    delth=abs((y1(1)-y0(1))*(y2(2)-y1(2))-(y1(2)-y0(2))*(y2(1)-y1(1))) &
         /sum((y1-y0)**2)
    h_out=min(h_out*delth_rec/max(delth,epsilon(1.d0)),h_in)
  enddo
!
  end subroutine choose_step
!
!------------
!
  subroutine adjust_crossing
!
  implicit none
!
  do
!
    call vparzero_vec(R,Z,vpar2,gradvpar2)
!
    det=dpsidr*gradvpar2(2)-dpsidz*gradvpar2(1)
    delR=((psif-psi_bound)*gradvpar2(2)-dpsidz*vpar2)/det
    delZ=(dpsidr*vpar2-(psif-psi_bound)*gradvpar2(1))/det
    R=R-delR
    Z=Z-delZ
    if(delR**2+delZ**2.lt.err_dist) exit
  enddo
!
  end subroutine adjust_crossing
!
!------------
!
  end subroutine find_vparzero_line
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine vparzero_line(s,R,Z,dR_ds,dZ_ds)
!
! Interpolates v_par^2=0 line: computes R(Z) and dR(Z)/dZ using the data ponints
! pre-computed by "find_vparzero_line" and stored in the module "vparzero_line_mod"
!
  use vparzero_line_mod, only : nvpline,s_line,vpar_line
!
  implicit none
!
  integer, parameter :: nplag=4, nder=1
  integer :: ibeg,iend
  double precision :: s,Z,R,dR_ds,dZ_ds
  double precision, dimension(0:nder,nplag) :: coef
!
  ibeg=ceiling(s)
  ibeg=max(0,ibeg-nplag/2)
  iend=ibeg+nplag-1
  if(iend.gt.nvpline) then
    iend=nvpline
    ibeg=iend+1-nplag
  endif
!
  call plag_coeff(nplag,nder,s,s_line(ibeg:iend),coef)
!
  R=sum(coef(0,:)*vpar_line(1,ibeg:iend))
  dR_ds=sum(coef(1,:)*vpar_line(1,ibeg:iend))
  Z=sum(coef(0,:)*vpar_line(2,ibeg:iend))
  dZ_ds=sum(coef(1,:)*vpar_line(2,ibeg:iend))
!
  end subroutine vparzero_line
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine stagnation_point(icase,toten,perpinv,             &
                              sst,rst,zst,sst2,rst2,zst2,ierr)
!
! Finds the point where parallel acceleration $\dot \lambda$ is zero
! on the forbidden boundary $v_\parallel^2=0$. Straigh lines connecting
! this point (rst,zst) with each of two B^* axes are used as Poincare
! cuts for co- and counter-passing particles
!
  use vparzero_line_mod, only : nvpline,s_line
!
  implicit none
!
  double precision, parameter :: relerr=1d-10 !8
  logical :: prop
  integer :: icase,ierr,i
  double precision :: toten,perpinv,sst,rst,zst,alamdot,dalamdot_ds
  double precision :: sst2,rst2,zst2
  double precision :: Z_u,Z_d,ald_d,ald_u,dels,err_dist
  double precision :: vpar2,gradvpar2,delR
!
  ierr=0
!
  call par_accel(s_line(0),rst,zst,ald_d,dalamdot_ds)
!
  prop=.false.
!
  do i=1,nvpline
!
    call par_accel(s_line(i),rst,zst,alamdot,dalamdot_ds)
!
    if(alamdot*ald_d.lt.0.d0) then
      if(prop) then
        sst2=s_line(i)
        prop=.false.
        exit
      endif
      sst=s_line(i)
      ald_d=-ald_d
      prop=.true.
    endif
  enddo
!
  err_dist=rst*relerr
!
  do
!
    call par_accel(sst,rst,zst,alamdot,dalamdot_ds)
!
    dels=alamdot/dalamdot_ds
    sst=sst-dels
    if(abs(dels).lt.relerr) exit
  enddo
!
  do
!
    call vparzero(rst,zst,vpar2,gradvpar2)
!
    delR=vpar2/gradvpar2
    rst=rst-delR
    if(abs(delR).lt.err_dist) exit
  enddo
!
  if(icase.eq.3) then
!
    do
!
      call par_accel(sst2,rst2,zst2,alamdot,dalamdot_ds)
!
      dels=alamdot/dalamdot_ds
      sst2=sst2-dels
      if(abs(dels).lt.relerr) exit
    enddo
!
    do
!
      call vparzero(rst2,zst2,vpar2,gradvpar2)
!
      delR=vpar2/gradvpar2
      rst2=rst2-delR
      if(abs(delR).lt.err_dist) exit
    enddo
!
  endif
!------
!
  contains
!
!------
!
  subroutine par_accel(s,Rl,Zl,alamdot,dalamdot_ds)
!
  implicit none
!
  integer, parameter          :: ndim=5
  double precision, parameter :: eps_dif=1.d-6
  double precision :: s,Zl,Rl,alamdot,dalamdot_ds,dR_ds,dZ_ds,dtau,phi_elec,p2
  double precision :: hdif
  double precision, dimension(ndim) :: z,vz
  double precision, dimension(3) :: derphi
!
  call vparzero_line(s,Rl,Zl,dR_ds,dZ_ds)
!
! normal velocity:
!
  z(1)=Rl
  z(2)=0.d0
  z(3)=Zl
  z(5)=0.d0
!
  call elefie(z(1:3),phi_elec,derphi)
!
  p2=toten-phi_elec
  if(p2.le.0.d0) then
    print *,'tangential_velocity: negative kinetic energy '
    return
  endif
  z(4)=sqrt(p2)
!
  call velo(dtau,z,vz)
!
  alamdot=vz(5)
!
! end normal velocity
!
! derivative of normal velocity over Z:
!
  hdif=eps_dif
!
  call vparzero_line(s-hdif,Rl,Zl,dR_ds,dZ_ds)
!
  z(1)=Rl
  z(2)=0.d0
  z(3)=Zl
  z(5)=0.d0
!
  call elefie(z(1:3),phi_elec,derphi)
!
  p2=toten-phi_elec
  if(p2.le.0.d0) then
    print *,'tangential_velocity: negative kinetic energy '
    return
  endif
  z(4)=sqrt(p2)
!
  call velo(dtau,z,vz)
!
  dalamdot_ds=vz(5)
!
  call vparzero_line(s+hdif,Rl,Zl,dR_ds,dZ_ds)
!
  z(1)=Rl
  z(2)=0.d0
  z(3)=Zl
  z(5)=0.d0
!
  call elefie(z(1:3),phi_elec,derphi)
!
  p2=toten-phi_elec
  if(p2.le.0.d0) then
    print *,'tangential_velocity: negative kinetic energy '
    return
  endif
  z(4)=sqrt(p2)
!
  call velo(dtau,z,vz)
!
  dalamdot_ds=(vz(5)-dalamdot_ds)/(2.d0*hdif)
!
! end derivative of normal velocity over Z
!
  end subroutine par_accel
!
!------
!
  subroutine vparzero(R,Z,vpar2,gradvpar2)
!
  implicit none
!
  double precision :: bmod,sqrtg,phi_elec,R,Z,vpar2,gradvpar2
  double precision, dimension(3)    :: x,bder,hcovar,hctrvr,hcurl,derphi
!
  x(1)=R
  x(2)=0.d0
  x(3)=Z
!
  call magfie(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
!
  call elefie(x,phi_elec,derphi)
!
  vpar2=toten-phi_elec-perpinv*bmod
  gradvpar2=-derphi(1)-perpinv*bmod*bder(1)
!
  end subroutine vparzero
!
!------------
!
  end subroutine stagnation_point
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine find_bounce(next,velo_ext,dtau,z_eqm,taub,delphi,extraset)
!
! Integrates the orbit over one bounce time (finds this time). If needed
! (write_orb=.true.) writes it to the file with unit number "iunit1".
! Besides orbit equations integrates along the orbit and extra set of
! functions of phase space coordinates.
! Agruments:
! velo_ext       - external routine to integrate (input)
! dtau           - maximum value of time step (input)
! z_eqm(5)       - phase space variables (input/output)
! taub           - bounce time (output)
! delphi         - toroidal shift per bounce time (output)
! extraset(next) - extra integrals along the orbit
!
  use orbit_dim_mod, only : neqm,write_orb,iunit1,Rorb_max
!
  implicit none
!
  double precision, parameter :: relerr=1d-10 !8
!
  integer :: next, ndim
  double precision :: dtau,taub,delphi
  double precision :: dL2_pol,dL2_pol_start,dtau_newt,r_prev,z_prev
  double precision :: tau0,RNorm,ZNorm,vnorm,dnorm,vel_pol,dL2_pol_min
  double precision, dimension(neqm) :: z_eqm
  double precision, dimension(next) :: extraset
  double precision, dimension(neqm+next) :: z,z_start,vz
!
  external velo_ext
!
  ndim = neqm+next
!
  z(1:neqm)=z_eqm
  z(neqm+1:ndim)=extraset
  Rorb_max=z(1)
!
! Primary search:
!
  z_start=z
!
  call velo_ext(dtau,z,vz)
!
  vel_pol=sqrt(vz(1)**2+vz(3)**2)
  RNorm=vz(1)/vel_pol
  ZNorm=vz(3)/vel_pol
!
  tau0=0.d0
!
  if(write_orb) write (iunit1,*) z(1:neqm),vz(5)
!
  call odeint_allroutines(z,ndim,tau0,dtau,relerr,velo_ext)
!
  taub=dtau
  dL2_pol=2.d0*(z(1)-z_start(1))**2+(z(3)-z_start(3))**2
!
  if(write_orb) then
!
    call velo_ext(dtau,z,vz)
!
    write (iunit1,*) z(1:neqm),vz(5)
  endif
!
  do
    r_prev=z(1)
    z_prev=z(3)
!
    call odeint_allroutines(z,ndim,tau0,dtau,relerr,velo_ext)
!
    taub=taub+dtau
!    dL2_pol=max(dL2_pol,2.d0*((z(1)-r_prev)**2+(z(3)-z_prev)**2))
    dL2_pol=2.d0*((z(1)-r_prev)**2+(z(3)-z_prev)**2)
    dL2_pol_start=(z(1)-z_start(1))**2+(z(3)-z_start(3))**2
    if(dL2_pol_start.lt.dL2_pol) exit
    if(write_orb) then
!
      call velo_ext(dtau,z,vz)
!
      write (iunit1,*) z(1:neqm),vz(5)
    endif
    Rorb_max=max(Rorb_max,z(1))
  enddo
!
! End primary search
!
! Newton adjustment
!
  do
!
    call velo_ext(dtau,z,vz)
!
    vnorm=vz(1)*RNorm+vz(3)*ZNorm
    dnorm=(z_start(1)-z(1))*RNorm+(z_start(3)-z(3))*ZNorm
    if(dnorm**2.lt.dL2_pol*relerr) exit
    dtau_newt=dnorm/vnorm
!
    call odeint_allroutines(z,ndim,tau0,dtau_newt,relerr,velo_ext)
!
    taub=taub+dtau_newt
  enddo
!
  extraset=z(neqm+1:ndim)
  if(write_orb) then
!
    call velo_ext(dtau,z,vz)
!
    write (iunit1,*) z(1:neqm),vz(5)
  endif
!
! End Newton adjustment
!
  delphi=z(2)-z_start(2)
!
  end subroutine find_bounce

!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine find_bstar_axes(icase,toten,perpinv,                           &
                             sst,rst,zst,sst2,rst2,zst2,dtau,               &
                             trapped,copass,ctrpass,triplet,                &
                             raxis_co,zaxis_co,raxis_ctr,zaxis_ctr,         &
                             rxpoi,zxpoi,sigma_x,s_tp,                      &
                             rxpoi2,zxpoi2,sigma_x2,s_tp2,ierr)
!
! Finds stable points of the orbit map: magnetic axes of B^* field for co-
! and counter-passing B^* fields (if exist) and X-point if exists. In case
! X-point exists, finds the intersection of the separatrix orbit with v_par^2=0
! line "z_tp" - Poincare cut for trapped orbits is the segment of this line
! in the interval z_tp < Z < Z_max
!
  use field_eq_mod,      only : psif,dpsidr,dpsidz
  use vparzero_line_mod, only : nvpline,npsib,s_line,psi_line,rinb,zinb
  use orbit_dim_mod,     only : Rorb_max
!
  implicit none
!
  integer, parameter :: neqm=5,next=3,n_search=1000
  double precision, parameter :: relerr=1d-10, relerr_ax=1.d-4 !8
!
  logical :: trapped,copass,ctrpass,triplet,trigger,prop
  integer :: icase,ierr,iter,i
  double precision :: toten,perpinv,sst,rst,zst,sst2,rst2,zst2,dtau, &
                      raxis_co,zaxis_co,raxis_ctr,zaxis_ctr, &
                      rxpoi,zxpoi,sigma_x,s_tp,rxpoi2,zxpoi2,sigma_x2,s_tp2
  double precision :: bmod,phi_elec,p2,alam2,taub,delphi,sigma,Rorb_in
  double precision :: h_s,h_R,h_Z,vnorm,vnorm_prev,raxis,zaxis,p_phi
  double precision :: s,dR_ds,dZ_ds,delta_psif,delta_s,err_dist
  double precision, dimension(neqm)    :: z,vz
  double precision, dimension(next)    :: extraset
!
  ierr=0
!
  if(icase.eq.1) then
    print *,'find_bstar_axes: icase=1, no orbits'
    trapped=.false.
    copass=.false.
    ctrpass=.false.
    return
  elseif(icase.eq.4) then
    trapped=.false.
    z(1)=sum(psi_line(1,1:npsib))/dble(npsib)
    z(2)=0.d0
    z(3)=sum(psi_line(2,1:npsib))/dble(npsib)
    sigma=1.d0
!
    call find_axis(toten,perpinv,sigma,z,ierr)
!
    copass=.true.
    raxis_co=z(1)
    zaxis_co=z(3)
    sigma=-1.d0
!
    call find_axis(toten,perpinv,sigma,z,ierr)
!
    ctrpass=.true.
    raxis_ctr=z(1)
    zaxis_ctr=z(3)
    return
  endif
!
  trapped=.true.
  copass=.false.
  ctrpass=.false.
!
! Search for the main, non-vanishing axis:
!
  Rorb_in=rst
  prop=.false.
!
  z(1)=rst
  z(2)=0.d0
  z(3)=zst
!
  call get_bmod_and_Phi(z(1:3),bmod,phi_elec)
!
  p2=toten-phi_elec
  if(p2.le.0.d0) then
    print *,'find_bstar_axes: negative kinetic energy '
    ierr=1
    return
  endif
  z(4)=sqrt(p2)
  if(trapped) then
    z(5)=0.d0
  else
    alam2=1.d0-perpinv*bmod/p2
    alam2=max(0.d0,alam2)
    z(5)=sqrt(alam2)
  endif
  extraset=0.d0
!
  call find_bounce(next,velo_ext,dtau,z,taub,delphi,extraset)
!
  sigma=sign(1.d0,extraset(3))
  z(1)=extraset(1)/taub
  z(2)=0.d0
  z(3)=extraset(2)/taub
!
  do iter=1,30
!
    call get_bmod_and_Phi(z(1:3),bmod,phi_elec)
!
    p2=toten-phi_elec
    if(p2.le.0.d0) then
      print *,'find_bstar_axes: negative kinetic energy '
      ierr=1
      return
    endif
    z(4)=sqrt(p2)
    alam2=1.d0-perpinv*bmod/p2
    if(alam2.lt.0.d0) then
      Rorb_in=0.5d0*(Rorb_in+Rorb_max)
      z(1)=Rorb_in
      z(2)=0.d0
      z(3)=zst
      if(prop) then
        print *,'find_bstar_axes: negative lambda^2 '
        ierr=1
        return
      endif
      prop=.true.
      cycle
    endif
    prop=.false.
    z(5)=sign(sqrt(alam2),sigma)
    extraset=0.d0
!
    call find_bounce(next,velo_ext,dtau,z,taub,delphi,extraset)
!
    z(1)=extraset(1)/taub
    z(2)=0.d0
    z(3)=extraset(2)/taub
    if(abs(z(1)/Rorb_max-1.d0).lt.relerr_ax) exit
  enddo
!
  call find_axis(toten,perpinv,sigma,z,ierr)
!
  raxis=z(1)
  zaxis=z(3)
!
  if(sigma.gt.0.d0) then
    if(ierr.eq.0) then
      copass=.true.
      raxis_co=z(1)
      zaxis_co=z(3)
    else
      copass=.false.
    endif
  else
    if(ierr.eq.0) then
      ctrpass=.true.
      raxis_ctr=z(1)
      zaxis_ctr=z(3)
    else
      ctrpass=.false.
    endif
  endif
!
! End search for the main, non-vanishing axis
!
  sigma=-sigma
  sigma_x=sigma
!
  triplet=.false.
!
! Rough search for X-point and second axis:
!
  if(trapped) then
    trigger=.true.
    h_R=(z(1)-rst)/dble(n_search**2)
    h_Z=(z(3)-zst)/dble(n_search**2)
    do i=0,n_search
      z(1)=rst+h_R*dfloat(i**2)
      z(2)=0.d0
      z(3)=zst+h_Z*dfloat(i**2)
!
      call get_bmod_and_Phi(z(1:3),bmod,phi_elec)
!
      p2=toten-phi_elec
      if(p2.le.0.d0) then
        print *,'find_bstar_axes: negative kinetic energy '
        ierr=1
        return
      endif
      z(4)=sqrt(p2)
      alam2=1.d0-perpinv*bmod/p2
      if(alam2.lt.0.d0) then
        if(i.eq.0) then
          alam2=0.d0
        else
          print *,'find_bstar_axes: negative lambda^2 = ',alam2
          ierr=1
          return
        endif
      endif
      z(5)=sign(sqrt(alam2),sigma)
!
      call velo(dtau,z,vz)
!
      vnorm=vz(1)*h_Z-vz(3)*h_R
      if(i.eq.0) then
        vnorm_prev=vnorm
      else
        if(vnorm*vnorm_prev.lt.0.d0) then
          vnorm_prev=-vnorm_prev
          if(trigger) then
            rxpoi=z(1)
            zxpoi=z(3)
            trigger=.false.
            triplet=.true.
          else
            raxis=z(1)
            zaxis=z(3)
          endif
        endif
      endif
    enddo
  endif
!
! End rough search for X-point and second axis
!
! Refinement of the second axis:
!
  if(triplet.or..not.trapped) then
    z(1)=raxis
    z(2)=0.d0
    z(3)=zaxis
!
    call find_axis(toten,perpinv,sigma,z,ierr)
!
    if(sigma.gt.0.d0) then
      if(ierr.eq.0) then
        copass=.true.
        raxis_co=z(1)
        zaxis_co=z(3)
      else
        copass=.false.
      endif
    else
      if(ierr.eq.0) then
        ctrpass=.true.
        raxis_ctr=z(1)
        zaxis_ctr=z(3)
      else
        ctrpass=.false.
      endif
    endif
  endif
!
! End refinement of the second axis
!
! Refinement of the X-point:
!
  if(triplet) then
    z(1)=rxpoi
    z(2)=0.d0
    z(3)=zxpoi
!
    call find_axis(toten,perpinv,sigma,z,ierr)
!
    rxpoi=z(1)
    zxpoi=z(3)
!
! End refinement of the X-point
!
! Find trapped-passing boundary on the v_par=0 line:
!
    z(1)=rxpoi
    z(2)=0.d0
    z(3)=zxpoi
!
    call get_bmod_and_Phi(z(1:3),bmod,phi_elec)
!
    p2=toten-phi_elec
    z(4)=sqrt(p2)
    alam2=1.d0-perpinv*bmod/p2
    z(5)=sign(sqrt(alam2),sigma)
!
    call get_tormom(z,p_phi)
!
! rough search:
!
    h_s=(s_line(nvpline)-sst)/dble(n_search)**2
!
    do i=0,n_search
      s=sst+h_s*dble(i)**2
!
      call vparzero_line(s,z(1),z(3),dR_ds,dZ_ds)
!
      call get_bmod_and_Phi(z(1:3),bmod,phi_elec)
!
      if(i.eq.0) then
        delta_psif=psif-p_phi
      elseif((psif-p_phi)*delta_psif.lt.0.d0) then
        exit
      endif
    enddo
!
! refinement:
!
    err_dist=relerr
!
    do
!
      call vparzero_line(s,z(1),z(3),dR_ds,dZ_ds)
!
      call get_bmod_and_Phi(z(1:3),bmod,phi_elec)
!
      delta_s=(p_phi-psif)/(dpsidr*dR_ds+dpsidz*dZ_ds)
      s=s+delta_s
      if(abs(delta_s).lt.err_dist) exit
    enddo
!
    s_tp=s
!
! End find trapped-passing boundary on the v_par=0 line
!
  else
    sigma_x=0.d0
    s_tp=sst
  endif
!
  if(icase.eq.3) then
    sigma_x2=0.d0
    sigma=-sigma
    trigger=.true.
    h_R=(Rinb-rst2)/dble(n_search**2)
    h_Z=(Zinb-zst2)/dble(n_search**2)
    do i=0,n_search
      z(1)=rst2+h_R*dfloat(i**2)
      z(2)=0.d0
      z(3)=zst2+h_Z*dfloat(i**2)
!
      call get_bmod_and_Phi(z(1:3),bmod,phi_elec)
!
      p2=toten-phi_elec
      if(p2.le.0.d0) then
        print *,'find_bstar_axes: negative kinetic energy '
        ierr=1
        return
      endif
      z(4)=sqrt(p2)
      alam2=1.d0-perpinv*bmod/p2
      if(alam2.lt.0.d0) then
        if(i.eq.0) then
          alam2=0.d0
        else
          print *,'find_bstar_axes: negative lambda^2 = ',alam2
          ierr=1
          return
        endif
      endif
      z(5)=sign(sqrt(alam2),sigma)
!
      call velo(dtau,z,vz)
!
      vnorm=vz(1)*h_Z-vz(3)*h_R
      if(i.eq.0) then
        vnorm_prev=vnorm
      else
        if(vnorm*vnorm_prev.lt.0.d0) then
          vnorm_prev=-vnorm_prev
          rxpoi2=z(1)
          zxpoi2=z(3)
          sigma_x2=sigma
        endif
      endif
    enddo
!
    z(1)=rxpoi2
    z(2)=0.d0
    z(3)=zxpoi2
!
    call find_axis(toten,perpinv,sigma,z,ierr)
!
    rxpoi2=z(1)
    zxpoi2=z(3)
!
! Find trapped-passing boundary on the v_par=0 line:
!
    call get_bmod_and_Phi(z(1:3),bmod,phi_elec)
!
    p2=toten-phi_elec
    z(4)=sqrt(p2)
    alam2=1.d0-perpinv*bmod/p2
    z(5)=sign(sqrt(alam2),sigma)
!
    call get_tormom(z,p_phi)
!
! rough search:
!
    h_s=(sst-sst2)/dble(n_search)**2
!
    do i=0,n_search
      s=sst2+h_s*dble(i)**2
!
      call vparzero_line(s,z(1),z(3),dR_ds,dZ_ds)
!
      call get_bmod_and_Phi(z(1:3),bmod,phi_elec)
!
      if(i.eq.0) then
        delta_psif=psif-p_phi
      elseif((psif-p_phi)*delta_psif.lt.0.d0) then
        exit
      endif
    enddo
!
! refinement:
!
    err_dist=relerr
!
    do
!
      call vparzero_line(s,z(1),z(3),dR_ds,dZ_ds)
!
      call get_bmod_and_Phi(z(1:3),bmod,phi_elec)
!
      delta_s=(p_phi-psif)/(dpsidr*dR_ds+dpsidz*dZ_ds)
      s=s+delta_s
      if(abs(delta_s).lt.err_dist) exit
    enddo
!
    s_tp2=s
!
  endif
!
  contains
!
!---------------------
!
  subroutine velo_ext(dtau2,z2,vz2)
    !
      use orbit_dim_mod, only : neqm,next,ndim
    !
      double precision :: dtau2
      double precision, dimension(ndim) :: z2,vz2
    !
      call velo(dtau2,z2(1:neqm),vz2(1:neqm))
    !
      vz2(neqm+1)=z2(1)
      vz2(neqm+2)=z2(3)
      vz2(neqm+3)=z2(4)*z2(5)
    !
      end subroutine velo_ext

  end subroutine find_bstar_axes
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine find_axis(toten,perpinv,sigma,z,ierr)
!
! Finds invariant axes by Newton method solving equations $\dot R=0$ and $\dot Z = 0$
!
  implicit none
!
  integer, parameter          :: ndim=5
  double precision, parameter :: relerr=1d-10 !8
!
  integer :: ierr,i
  double precision :: toten,perpinv,sigma,toten_tmp,perpinv_tmp
  double precision :: err_dist,facnewt,bmod,phi_elec,p2,alam2
  double precision :: tau0,RNorm,ZNorm,vnorm,dnorm,vel_pol
  double precision, dimension(ndim) :: z,z_start,vz,z_restart
!
  ierr=0
!
  call get_bmod_and_Phi(z(1:3),bmod,phi_elec)
!
  p2=toten-phi_elec
  if(p2.le.0.d0) then
    print *,'find_axis: negative kinetic energy '
    return
  endif
  z(4)=sqrt(p2)
!
  alam2=1.d0-perpinv*bmod/p2
  if(alam2.lt.0.d0) then
    print *,'find_axis: negative lambda^2 '
    return
  endif
  z(5)=sign(sqrt(alam2),sigma)
!
  err_dist=z(1)*relerr
!
  facnewt=0.d0
!
  z_restart=z
!
  do i=1,1000
    z_start=z
!
    call next_iter_axis(facnewt,z_start,z,toten_tmp,perpinv_tmp,ierr)
!
    if(ierr.ne.0) exit
    if(max(abs(z(1)-z_start(1)),abs(z(3)-z_start(3))).lt.err_dist) exit
  enddo
!
  end subroutine find_axis
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine next_iter_axis(facnewt,z_start,z,toten,perpinv,ierr)
!
! Service routine for "find_axis"
!
  implicit none
!
  integer, parameter          :: ndim=5,nmat=4
  double precision, parameter :: eps_dif=1.d-6, relerr=1d-10 !8
!
  integer :: ierr,i
  double precision :: toten,perpinv,sigma,bmod,sqrtg,phi_elec
  double precision :: hdif,dtau,p2,alam2,facnewt,stepmod
  double precision, dimension(ndim) :: z,z_start,vz
  double precision, dimension(nmat)       :: bvec
  double precision, dimension(nmat,nmat)  :: amat_jac, amat
  integer,          dimension(nmat)       :: ipiv
  double precision, dimension(3)    :: x,bder,hcovar,hctrvr,hcurl,derphi
!
! y^i = (R,Z,p,lambda)
! F_1=V^R, F_2=V^Z, F^3=p^2+Phi-w, F^3=p^2(1-lambda^2)/B-J_perp (all normalized)
!
  ierr=0
!
  x=z(1:3)
!
  call magfie(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
!
  call elefie(x,phi_elec,derphi)
!
  toten=z(4)**2+phi_elec
!
  perpinv=z(4)**2*(1.d0-z(5)**2)/bmod
!
  sigma=sign(1.d0,z(5))
!
  hdif=x(1)*eps_dif
  z=z_start
!
  call velo(dtau,z,vz)
!
  bvec(1)=vz(1)
  bvec(2)=vz(3)
  bvec(3)=0.d0
  bvec(4)=0.d0
!
! derivative over R:
  z(1)=z_start(1)+hdif
!
  call velo(dtau,z,vz)
!
  amat_jac(1,1)=vz(1)
  amat_jac(2,1)=vz(3)
  z(1)=z_start(1)-hdif
!
  call velo(dtau,z,vz)
!
  amat_jac(1,1)=(amat_jac(1,1)-vz(1))/(2.d0*hdif)
  amat_jac(2,1)=(amat_jac(2,1)-vz(3))/(2.d0*hdif)
  z(1)=z_start(1)
  amat_jac(3,1)=derphi(1)
  amat_jac(4,1)=-perpinv*bder(1)
! end derivative over R
!
! derivative over Z:
  z(3)=z_start(3)+hdif
!
  call velo(dtau,z,vz)
!
  amat_jac(1,2)=vz(1)
  amat_jac(2,2)=vz(3)
  z(3)=z_start(3)-hdif
!
  call velo(dtau,z,vz)
!
  amat_jac(1,2)=(amat_jac(1,2)-vz(1))/(2.d0*hdif)
  amat_jac(2,2)=(amat_jac(2,2)-vz(3))/(2.d0*hdif)
  z(3)=z_start(3)
  amat_jac(3,2)=derphi(3)
  amat_jac(4,2)=-perpinv*bder(3)
! end derivative over Z
!
! derivative over p:
  z(4)=z_start(4)+eps_dif
!
  call velo(dtau,z,vz)
!
  amat_jac(1,3)=vz(1)
  amat_jac(2,3)=vz(3)
  z(4)=z_start(4)-eps_dif
!
  call velo(dtau,z,vz)
!
  amat_jac(1,3)=(amat_jac(1,3)-vz(1))/(2.d0*eps_dif)
  amat_jac(2,3)=(amat_jac(2,3)-vz(3))/(2.d0*eps_dif)
  z(4)=z_start(4)
  amat_jac(3,3)=2.d0*z(4)
  amat_jac(4,3)=2.d0*perpinv/z(4)
! end derivative over p
!
! derivative over lambda:
  z(5)=z_start(5)+eps_dif
!
  call velo(dtau,z,vz)
!
  amat_jac(1,4)=vz(1)
  amat_jac(2,4)=vz(3)
  z(5)=z_start(5)-eps_dif
!
  call velo(dtau,z,vz)
!
  amat_jac(1,4)=(amat_jac(1,4)-vz(1))/(2.d0*eps_dif)
  amat_jac(2,4)=(amat_jac(2,4)-vz(3))/(2.d0*eps_dif)
  z(5)=z_start(5)
  amat_jac(3,4)=0.d0
  amat_jac(4,4)=-2.d0*z(4)**2*z(5)/bmod
! end derivative over lambda
!
  call dgesv(nmat,1,amat_jac,nmat,ipiv,bvec,nmat,ierr)
!
  if(ierr.ne.0) then
    print *,'next_iter_axis: error in dgesv = ',ierr
    ierr=1
    return
  endif
!
  if(facnewt.eq.0.d0) then
    z(1)=z_start(1)-bvec(1)
    z(3)=z_start(3)-bvec(2)
  else
    stepmod=sqrt(bvec(1)**2+bvec(2)**2)
    z(1)=z_start(1)-bvec(1)*facnewt/stepmod
    z(3)=z_start(3)-bvec(2)*facnewt/stepmod
  endif
  x=z(1:3)
!
  call magfie(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
!
  call elefie(x,phi_elec,derphi)
!
  p2=toten-phi_elec
  if(p2.le.0.d0) then
    print *,'next_iter_axis: negative kinetic energy '
    ierr=2
    return
  endif
  z(4)=sqrt(p2)
!
  alam2=1.d0-perpinv*bmod/p2
  if(alam2.lt.0.d0) then
    print *,'next_iter_axis: negative lambda^2 '
    ierr=3
    return
  endif
  z(5)=sign(sqrt(alam2),sigma)
!
  end subroutine next_iter_axis
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine get_bmod_and_Phi(x,bmod,phi_elec)
!
  implicit none
!
  double precision :: bmod,sqrtg,phi_elec
  double precision, dimension(3)    :: x,bder,hcovar,hctrvr,hcurl,derphi
!
  call magfie(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
!
  call elefie(x,phi_elec,derphi)
!
  end subroutine get_bmod_and_Phi
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine get_tormom(z,p_phi)
!
  use field_eq_mod, only : psif
  use parmot_mod, only : ro0
!
  implicit none
!
  integer, parameter :: neqm=5
  double precision :: p_phi
  double precision :: bmod,sqrtg,phi_elec
  double precision, dimension(3)    :: bder,hcovar,hctrvr,hcurl,derphi
  double precision, dimension(neqm) :: z
!
  call magfie(z(1:3),bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
!
  p_phi=ro0*z(4)*z(5)*hcovar(2)+psif
!
  end subroutine get_tormom
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
