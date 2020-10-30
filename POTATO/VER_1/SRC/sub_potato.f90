!
  module vparzero_line_mod
    integer :: nzline
    double precision :: hz,Zmid,Z_up
    double precision, dimension(:), allocatable :: rline,zline
  end module vparzero_line_mod
!
!------------------------------------------------------
!
  module orbit_dim_mod
    integer, parameter  :: neqm=5, next=3, ndim=neqm+next
    logical             :: write_orb=.false.
    integer             :: iunit1
  end module orbit_dim_mod
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine find_vparzero_line(nline,toten,perpinv,rmin,rmax,zmin,zmax, &
                                R_line,Z_line,ierr)
!
! Looks for the forbidden line "v_par^2=0" in the box given by "rmin<R<rmax", "zmin<Z<zmax"
! for particles with normalized total energy "toten" and perpendicular invariant "perpinv"
! Points on the line, "Z_line", are equidistant in Z and are numbered as [-nline:nline].
! If this line is outside or crosses HFS boundary R=rmin, ierr=1 - most of the box is occupied
! by passing orbits, "R_line" is set to "rmin" if the actual line is fully outside the box
! (only passing orbits will be computed in the computation box), if actulal line crosses
! the inner boundary of the box, this boundary is replaced by the forbidden line.
! If this line is outside or crosses LFS boundary R=rmax, ierr=2 - most of the box is forbiden,
! no orbits will be computed, "R_line" is set to "rmax" for the points outside the box
!
  use vparzero_line_mod, only : nzline,hz,Zmid,rline,zline,Z_up
!
  implicit none
!
  integer, parameter :: nbsrc=10
  double precision, parameter :: relerr=1d-12 !8
!
  logical :: prop
  integer :: nline,ierr,i,j
  double precision :: toten,perpinv,rmin,rmax,zmin,zmax
  double precision :: err_dist,vpar2,gradvpar2,R,Z
  double precision :: vpar2in,vpar2out
  double precision :: Rin,Rout,delR
  double precision, dimension(-nline:nline) :: R_line,Z_line
!
  ierr=0
  err_dist=rmax*relerr
!
  Z_up=Zmax
  hz=(Zmax-Zmin)/dble(2*nline)
  Zmid=Zmin+hz*dble(nline)
  nzline=nline
  allocate(rline(-nzline:nzline),zline(-nzline:nzline))
  prop=.true.
!
  do i=-nline,nline
    Z_line(i)=Zmid+hz*dble(i)
!    Rin=Rmin
    Rin=0.8d0*Rmin
    Rout=Rmax
!
    call vparzero(Rmin,Z_line(i),vpar2in,gradvpar2)
!
    if(vpar2in.lt.0.d0) prop=.false.
!
    call vparzero(Rin,Z_line(i),vpar2in,gradvpar2)
    call vparzero(Rout,Z_line(i),vpar2out,gradvpar2)
!
    if(vpar2in*vpar2out.gt.0.d0) then
      if(vpar2out.lt.0.d0) then
        ierr=2
        R_line(i)=rmax
        cycle
      else
        ierr=1
        R_line(i)=rmin
        cycle
      endif
    endif
!
    do j=1,nbsrc
      R=0.5*(Rin+Rout)
!
      call vparzero(R,Z_line(i),vpar2,gradvpar2)
!
      if(vpar2*vpar2in.gt.0.d0) then
        Rin=R
      else
        Rout=R
      endif
    enddo
!
    do
!
      call vparzero(R,Z_line(i),vpar2,gradvpar2)
!
      delR=vpar2/gradvpar2
      R=R-delR
      if(abs(delR).lt.err_dist) exit
    enddo
!
    R_line(i)=R
  enddo
!
  if(prop) then
    R_line=rmin
    ierr=1
  endif
!
  rline=R_line
  zline=Z_line
!------------
!
  contains
!
!------------
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
  end subroutine find_vparzero_line
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine vparzero_line(Z,R,dR_dZ)
!
! Interpolates v_par^2=0 line: computes R(Z) and dR(Z)/dZ using the data ponints
! pre-computed by "find_vparzero_line" and stored in the module "vparzero_line_mod"
!
  use vparzero_line_mod, only : nzline,hz,Zmid,rline,zline
!
  implicit none
!
  integer, parameter :: nplag=4, nder=1
  integer :: ibeg,iend
  double precision :: Z,R,dR_dZ
  double precision, dimension(0:nder,nplag) :: coef
!
  ibeg=ceiling(Z/hz)
  ibeg=max(-nzline,ibeg-nplag/2)
  iend=ibeg+nplag-1
  if(iend.gt.nzline) then
    iend=nzline
    ibeg=iend+1-nplag
  endif
!
  call plag_coeff(nplag,nder,Z,zline(ibeg:iend),coef)
!
  R=sum(coef(0,:)*rline(ibeg:iend))
  dR_dZ=sum(coef(1,:)*rline(ibeg:iend))
!
  end subroutine vparzero_line
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine stagnation_point(toten,perpinv,rst,zst,ierr)
!
! Finds the point where parallel acceleration $\dot \lambda$ is zero
! on the forbidden boundary $v_\parallel^2=0$. Straigh lines connecting
! this point (rst,zst) with each of two B^* axes are used as Poincare
! cuts for co- and counter-passing particles
!
  use vparzero_line_mod, only : nzline,zline
!
  implicit none
!
  integer, parameter :: nbsrc=10
  double precision, parameter :: relerr=1d-10 !8
  integer :: ierr,i
  double precision :: toten,perpinv,rst,zst,alamdot,dalamdot_dZ
  double precision :: Z_u,Z_d,ald_d,ald_u,delZ,err_dist
  double precision :: vpar2,gradvpar2,delR
!
  ierr=0
!
  Z_d=zline(-nzline)
  Z_u=zline(nzline)
!
  call par_accel(rst,Z_d,ald_d,dalamdot_dZ)
  call par_accel(rst,Z_u,ald_u,dalamdot_dZ)
!
  if(ald_d*ald_u.gt.0.d0) then
    ierr=1
    return
  endif
!
  do i=1,nbsrc
    zst=0.5*(Z_d+Z_u)
!
    call par_accel(rst,zst,alamdot,dalamdot_dZ)
!
    if(alamdot*ald_d.gt.0.d0) then
      Z_d=zst
    else
      Z_u=zst
    endif
  enddo
!
  err_dist=rst*relerr
!
  do
!
    call par_accel(rst,zst,alamdot,dalamdot_dZ)
!
    delZ=alamdot/dalamdot_dZ
    zst=zst-delZ
    if(abs(delZ).lt.err_dist) exit
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
!------
!
  contains
!
!------
!
  subroutine par_accel(Rl,Zl,alamdot,dalamdot_dZ)
!
  implicit none
!
  integer, parameter          :: ndim=5
  double precision, parameter :: eps_dif=1.d-6
  double precision :: Zl,Rl,alamdot,dalamdot_dZ,dR_dZ,dtau,phi_elec,p2
  double precision :: hdif
  double precision, dimension(ndim) :: z,vz
  double precision, dimension(3) :: derphi
!
  call vparzero_line(Zl,Rl,dR_dZ)
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
  hdif=Rl*eps_dif
!
  call vparzero_line(Zl-hdif,Rl,dR_dZ)
!
  z(1)=Rl
  z(2)=0.d0
  z(3)=Zl-hdif
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
  dalamdot_dZ=vz(5)
!
  call vparzero_line(Zl+hdif,Rl,dR_dZ)
!
  z(1)=Rl
  z(2)=0.d0
  z(3)=Zl+hdif
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
  dalamdot_dZ=(vz(5)-dalamdot_dZ)/(2.d0*hdif)
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
  subroutine velo_ext(dtau,z,vz)
!
  use orbit_dim_mod, only : neqm,next,ndim
!
  double precision :: dtau
  double precision, dimension(ndim) :: z,vz
!
  call velo(dtau,z(1:neqm),vz(1:neqm))
!
  vz(neqm+1)=z(1)
  vz(neqm+2)=z(3)
  vz(neqm+3)=z(4)*z(5)
!
  end subroutine velo_ext
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine find_bounce(next, dtau,z_eqm,taub,delphi,extraset)
!
! Integrates the orbit over one bounce time (finds this time). If needed
! (write_orb=.true.) writes it to the file with unit number "iunit1".
! Besides orbit equations integrates along the orbit and extra set of
! functions of phase space coordinates.
! Agruments:
! naux           - auxiliary quantities to integrate
! dtau           - maximum value of time step (input)
! z_eqm(5)       - phase space variables (input/output)
! taub           - bounce time (output)
! delphi         - toroidal shift per bounce time (output)
! extraset(next) - extra integrals along the orbit
!
  use orbit_dim_mod, only : neqm,write_orb,iunit1
!
  implicit none
!
  double precision, parameter :: relerr=1d-10 !8
!
  integer :: next
  integer :: ndim
  double precision :: dtau,taub,delphi
  double precision :: dL2_pol,dL2_pol_start,dtau_newt,r_prev,z_prev
  double precision :: tau0,RNorm,ZNorm,vnorm,dnorm,vel_pol
  double precision, dimension(neqm) :: z_eqm
  double precision, dimension(next) :: extraset
  double precision, dimension(neqm+next) :: z,z_start,vz
!
  external velo_ext
!
  ndim = neqm + next
  z(1:neqm)=z_eqm
  z(neqm+1:ndim)=extraset
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
    dL2_pol=max(dL2_pol,2.d0*((z(1)-r_prev)**2+(z(3)-z_prev)**2))
    dL2_pol_start=(z(1)-z_start(1))**2+(z(3)-z_start(3))**2
    if(dL2_pol_start.lt.dL2_pol) exit
    if(write_orb) then
!
      call velo_ext(dtau,z,vz)
!
      write (iunit1,*) z(1:neqm),vz(5)
    endif
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
  subroutine find_bstar_axes(toten,perpinv,rst,zst,dtau,trapped,    &
                             copass,ctrpass,triplet,                &
                             raxis_co,zaxis_co,raxis_ctr,zaxis_ctr, &
                             rxpoi,zxpoi,sigma_x,z_tp,ierr)
!
! Finds stable points of the orbit map: magnetic axes of B^* field for co-
! and counter-passing B^* fields (if exist) and X-point if exists. In case
! X-point exists, finds the intersection of the separatrix orbit with v_par^2=0
! line "z_tp" - Poincare cut for trapped orbits is the segment of this line
! in the interval z_tp < Z < Z_max
!
  use field_eq_mod,      only : psif,dpsidr,dpsidz
  use vparzero_line_mod, only : nzline,hz,Zmid
!
  implicit none
!
  integer, parameter :: neqm=5,next=3,n_search=1000
  double precision, parameter :: relerr=1d-10 !8
!
  logical :: trapped,copass,ctrpass,triplet,trigger
  integer :: ierr,iter,i
  double precision :: toten,perpinv,rst,zst,dtau,            &
                      raxis_co,zaxis_co,raxis_ctr,zaxis_ctr, &
                      rxpoi,zxpoi,sigma_x,z_tp
  double precision :: bmod,phi_elec,p2,alam2,taub,delphi,sigma
  double precision :: h_R,h_Z,vnorm,vnorm_prev,raxis,zaxis,p_phi
  double precision :: dR_dZ,delta_psif,delta_Z,err_dist
  double precision, dimension(neqm)    :: z,vz
  double precision, dimension(next)    :: extraset
!
  ierr=0
  copass=.false.
  ctrpass=.false.
!
! Search for the main, non-vanishing axis:
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
  call find_bounce(dtau,z,taub,delphi,extraset)
!
  sigma=sign(1.d0,extraset(3))
  z(1)=extraset(1)/taub
  z(2)=0.d0
  z(3)=extraset(2)/taub
!
  do iter=1,3
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
      print *,'find_bstar_axes: negative lambda^2 '
      ierr=1
      return
    endif
    z(5)=sign(sqrt(alam2),sigma)
    extraset=0.d0
!
    call find_bounce(dtau,z,taub,delphi,extraset)
!
    z(1)=extraset(1)/taub
    z(2)=0.d0
    z(3)=extraset(2)/taub
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
        print *,'find_bstar_axes: negative lambda^2 = ',alam2
        if(i.eq.0) then
          alam2=0.d0
        else
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
    h_Z=(Zmid+hz*dble(nzline)-zst)/dble(n_search)**2
!
    do i=0,n_search
      z(3)=zst+h_Z*dble(i)**2
!
      call vparzero_line(z(3),z(1),dR_dZ)
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
    err_dist=z(1)*relerr
!
    do
!
      call vparzero_line(z(3),z(1),dR_dZ)
!
      call get_bmod_and_Phi(z(1:3),bmod,phi_elec)
!
      delta_Z=(p_phi-psif)/(dpsidr*dR_dZ+dpsidz)
      z(3)=z(3)+delta_Z
      if(abs(delta_Z).lt.err_dist) exit
    enddo
!
    z_tp=z(3)
!
! End find trapped-passing boundary on the v_par=0 line
!
  else
    sigma_x=0.d0
    z_tp=zst
  endif
!
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
