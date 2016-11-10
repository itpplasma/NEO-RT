MODULE magfield_mod
  INTEGER :: ierrfield=0
  INTEGER :: input_format,npmid,nr,np,nz,npoint
  INTEGER,          DIMENSION(:,:,:), ALLOCATABLE :: ipoint
  DOUBLE PRECISION, DIMENSION(:),     ALLOCATABLE :: rad,phi,zet,Brs,Bzs,Bps
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: Bx,By,Bz,Br,Bp
END MODULE
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
SUBROUTINE field(rrr,pp,zzz,Brad,Bphi,Bzet,dBrdR,dBrdp,dBrdZ  &             
     ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
!
  USE magfield_mod
!  use polylag_5, only : mp,indef, plag3d
  USE polylag_3, ONLY : mp,indef, plag3d
!
  IMPLICIT DOUBLE PRECISION (a-h,o-z)
!
! 22.11.2008  parameter (mp=4) ! power of Lagrange's polynomial =3
! 01.11.2015  parameter (mp=6) ! power of Lagrange's polynomial =5
  DIMENSION xp(mp),yp(mp),zp(mp),fp(mp,mp,mp)
  INTEGER indx(mp), indy(mp), indz(mp)
  CHARACTER *30 namedim(3), namevar(3) 
  CHARACTER *80 bez
  DATA icall/0/
  SAVE
!
!-------first call: read data from disk-------------------------------
  IF(icall .EQ. 0) THEN
    icall = 1
    iunit1=77
!
    OPEN(iunit1,file='MESH3D/field.format')
    READ(iunit1,*) input_format
    READ(iunit1,*) input_units
    IF(input_format.EQ.0) THEN
      READ(iunit1,*) nr
      READ(iunit1,*) np
      READ(iunit1,*) nz
    ENDIF
    CLOSE(iunit1)
!
    IF(input_format.EQ.0) THEN
!
       OPEN(iunit1,file='MESH3D/field.dat')
       READ(iunit1,*)   
       READ(iunit1,*)   
       READ(iunit1,*) bez 
       PRINT *,'description:',bez
       !
       ALLOCATE(Bx(nr,np,nz),By(nr,np,nz)) 
       ALLOCATE(Br(nr,np,nz),Bp(nr,np,nz),Bz(nr,np,nz))
       ALLOCATE(rad(nr),phi(np),zet(nz))
       !
       READ(iunit1,*)   
       READ(iunit1,*)        
       READ(iunit1,*) ndim 
       !      
       READ(iunit1,*)   
       READ(iunit1,*)        
       READ(iunit1,*) nvar          
       !      
       READ(iunit1,*)   
       READ(iunit1,*)        
       DO i=1,ndim
          READ(iunit1,*) namedim(i)            
       ENDDO
       PRINT *,'dimension:',ndim,namedim(1),namedim(2),namedim(3)
       DO i=1,nvar
          READ(iunit1,*) namevar(i)            
       ENDDO
       PRINT *,'function:',nvar,namevar(1),namevar(2),namevar(3)
       !      
       READ(iunit1,*)   
       READ(iunit1,*)        
       READ(iunit1,*) nwert 
       PRINT *,'number of values:',nwert
       READ(iunit1,*)   
       READ(iunit1,*)    
       !
       !---Input B      -->T = V*s/m/m
       DO i=1,nr
          DO j=1,np
             DO k=1,nz
                READ(iunit1,*) aaa,bbb,ccc,Bx(i,j,k),By(i,j,k),Bz(i,j,k)
                IF(i.EQ.1 .AND. j.EQ.1 .AND. k.EQ.1) THEN
                   rmin = aaa * 100.d0 !cm
                   pmin = bbb
                   zmin = ccc * 100.d0 !cm
                ELSEIF(i.EQ.1 .AND. j.EQ.2 .AND. k.EQ.1) THEN
                   delp = bbb
                ELSEIF(i.EQ.nr .AND. j.EQ.np .AND. k.EQ.nz) THEN
                   rmax = aaa * 100.d0 !cm
                   pmax = bbb
                   zmax = ccc * 100.d0 !cm
                ENDIF
             ENDDO
          ENDDO
       ENDDO
!
       CLOSE(iunit1)
!
       DO j=1,np
          ptmp = (j-1)*hphi
          sinp = SIN(ptmp)
          cosp = COS(ptmp)
          DO i=1,nr
             DO k=1,nz
                Br(i,j,k) = Bx(i,j,k)*cosp + By(i,j,k)*sinp
                Bp(i,j,k) = -Bx(i,j,k)*sinp + By(i,j,k)*cosp
             ENDDO
          ENDDO
       ENDDO
!
       DEALLOCATE(Bx,By)
!
       hrad = (rmax - rmin)/(nr-1)  
       hphi = (pmax - pmin)/(np-1)
       hzet = (zmax - zmin)/(nz-1)
       np = np + 1
       pmax = pmax + delp
       DO i=1,nr
          rad(i) = rmin + hrad*(i-1)
       ENDDO
       DO i=1,np
          phi(i) = pmin + hphi*(i-1)
       ENDDO
       DO i=1,nz
          zet(i) = zmin + hzet*(i-1)
       ENDDO
!
       DO i=1,nr
          DO k=1,nz
             Br(i,np,k) = Br(i,1,k)
             Bp(i,np,k) = Bp(i,1,k)
             Bz(i,np,k) = Bz(i,1,k)
          ENDDO
       ENDDO
!
       pmin_0=pmin
!
     ELSEIF(input_format.EQ.1) THEN
!
       OPEN(iunit1,file='MESH3D/field.dat')
       READ(iunit1,*) nr,np,nz
       READ(iunit1,*) rmin,rmax
       READ(iunit1,*) pmin,pmax
       READ(iunit1,*) zmin,zmax
       IF(input_units.EQ.0) THEN
         rmin = rmin * 100.d0 !cm
         rmax = rmax * 100.d0 !cm
         zmin = zmin * 100.d0 !cm
         zmax = zmax * 100.d0 !cm
       ENDIF
       np=np+4
       ALLOCATE(Br(nr,np,nz),Bp(nr,np,nz),Bz(nr,np,nz))
       ALLOCATE(rad(nr),phi(np),zet(nz))
       DO i=1,nr
          DO j=3,np-2
             DO k=1,nz
                READ(iunit1,*) Br(i,j,k),Bp(i,j,k),Bz(i,j,k)
             ENDDO
          ENDDO
       ENDDO
       CLOSE(iunit1)
!
       hrad = (rmax - rmin)/(nr-1)  
       hphi = (pmax - pmin)/(np-5)
       hzet = (zmax - zmin)/(nz-1)
!
       DO i=1,nr
          rad(i) = rmin + hrad*(i-1)
       ENDDO
       DO i=1,np
          phi(i) = pmin + hphi*(i-3)
       ENDDO
       DO i=1,nz
          zet(i) = zmin + hzet*(i-1)
       ENDDO
!
       Br(:,1,:)=Br(:,np-4,:)
       Br(:,2,:)=Br(:,np-3,:)
       Br(:,np-1,:)=Br(:,4,:)
       Br(:,np,:)=Br(:,5,:)
       Bp(:,1,:)=Bp(:,np-4,:)
       Bp(:,2,:)=Bp(:,np-3,:)
       Bp(:,np-1,:)=Bp(:,4,:)
       Bp(:,np,:)=Bp(:,5,:)
       Bz(:,1,:)=Bz(:,np-4,:)
       Bz(:,2,:)=Bz(:,np-3,:)
       Bz(:,np-1,:)=Bz(:,4,:)
       Bz(:,np,:)=Bz(:,5,:)
!
       pmin_0=pmin-2.d0*hphi
!
     ELSEIF(input_format.EQ.3) THEN
!
       PRINT *,'input_format = 3 : sparse field with stellarator symmetry'
       OPEN(iunit1,file='MESH3D/field_sparse.dat')
       READ(iunit1,*) 
       READ(iunit1,*)
       READ(iunit1,*)
       READ(iunit1,*)
       npoint=0
       DO
         READ(iunit1,*,END=11)
         npoint=npoint+1
       ENDDO
11     CLOSE(iunit1)
       OPEN(iunit1,file='MESH3D/field_sparse.dat')
       READ(iunit1,*) nr,np,nz
       READ(iunit1,*) rmin,rmax
       READ(iunit1,*) pmin,pmax
       READ(iunit1,*) zmin,zmax
!
       pmax=2*pmax-pmin
!
       npmid=np+2
       np=2*np+3
!
       ALLOCATE(ipoint(nr,np,nz))
       ipoint=0
!
       ALLOCATE(Brs(0:npoint),Bps(0:npoint),Bzs(0:npoint))
       Brs(0)=0.d0
       Bps(0)=0.d0
       Bzs(0)=0.d0
!
       DO i=1,npoint
         READ(iunit1,*) ir,ip,iz,Brs(i),Bps(i),Bzs(i)
         ipoint(ir,ip+2,iz)=i
       ENDDO
!
       CLOSE(iunit1)
!
       DO ir=1,nr
         DO iz=1,nz
           ipoint(ir,1,iz)=ipoint(ir,5,nz+1-iz)
           ipoint(ir,2,iz)=ipoint(ir,4,nz+1-iz)
         ENDDO
       ENDDO
!
       DO ir=1,nr
         DO ip=npmid+1,np
           DO iz=1,nz
             ipoint(ir,ip,iz)=ipoint(ir,np+1-ip,nz+1-iz)
           ENDDO
         ENDDO
       ENDDO
!
       ALLOCATE(rad(nr),phi(np),zet(nz))
!
       hrad = (rmax - rmin)/(nr-1)
       hphi = (pmax - pmin)/(np-5)
       hzet = (zmax - zmin)/(nz-1)
!
       DO i=1,nr
          rad(i) = rmin + hrad*(i-1)
       ENDDO
       DO i=1,np
          phi(i) = pmin + hphi*(i-3)
       ENDDO
       DO i=1,nz
          zet(i) = zmin + hzet*(i-1)
       ENDDO
!
       pmin_0=pmin-2.d0*hphi
!
     ELSE
!
       PRINT *,'unknown field format'
       STOP
!
     ENDIF
     !
     hrm1 = 1.d0/hrad
     hpm1 = 1.d0/hphi
     hzm1 = 1.d0/hzet
     phi_period=pmax - pmin
     !
!
  ENDIF
  !------- end first call ----------------------------------------------
!
  ierrfield=0
  !
  IF(pp.GT.pmax) THEN
    ppp=pp-phi_period*(INT((pp-pmax)/phi_period)+1)
  ELSEIF(pp.LT.pmin) THEN
    ppp=pp+phi_period*(INT((pmin-pp)/phi_period)+1)
  ELSE
    ppp=pp
  ENDIF
  !
  !
  CALL indef(rrr,rmin,hrm1,nr,indx)
  !      
  CALL indef(ppp,pmin_0,hpm1,np,indy)
  !      
  CALL indef(zzz,zmin,hzm1,nz,indz)
  !
  DO i=1,mp
     xp(i) = rad(indx(i))
     yp(i) = phi(indy(i))
     zp(i) = zet(indz(i))
  ENDDO
  !
  DO k=1,mp
    DO j=1,mp
      IF(input_format.EQ.3) THEN
        IF(indy(j).LE.2) THEN
          sigbr=-1.d0
        ELSEIF(indy(j).LE.npmid) THEN
          sigbr=1.d0
        ELSEIF(indy(j).LE.np-2) THEN
          sigbr=-1.d0
        ELSE
          sigbr=1.d0
        ENDIF
      ENDIF
      DO i=1,mp
        IF(input_format.EQ.3) THEN
          fp(i,j,k) = Brs(ipoint(indx(i),indy(j),indz(k)))*sigbr
        ELSE
          fp(i,j,k) = Br(indx(i),indy(j),indz(k))
        ENDIF
      ENDDO
    ENDDO
  ENDDO
  CALL plag3d(rrr,ppp,zzz,fp,hrm1,hpm1,hzm1,xp,yp,zp, &
       polylag,poly1x,poly1y,poly1z)
  Brad = polylag
  dBrdR = poly1x
  dBrdp = poly1y
  dBrdZ = poly1z
  !
  DO k=1,mp
    DO j=1,mp
      DO i=1,mp
        IF(input_format.EQ.3) THEN
          fp(i,j,k) = Bps(ipoint(indx(i),indy(j),indz(k)))
        ELSE
          fp(i,j,k) = Bp(indx(i),indy(j),indz(k))
        ENDIF
        IF(fp(i,j,k).EQ.0.d0) THEN
          ierrfield=1
PRINT *,'boundary touched'
!stop
          RETURN
        ENDIF
      ENDDO
    ENDDO
  ENDDO
  CALL plag3d(rrr,ppp,zzz,fp,hrm1,hpm1,hzm1,xp,yp,zp &
       ,polylag,poly1x,poly1y,poly1z)
  Bphi = polylag
  dBpdR = poly1x
  dBpdp = poly1y
  dBpdZ = poly1z
  !
  DO k=1,mp
    DO j=1,mp
      DO i=1,mp
        IF(input_format.EQ.3) THEN
          fp(i,j,k) = Bzs(ipoint(indx(i),indy(j),indz(k)))
        ELSE
          fp(i,j,k) = Bz(indx(i),indy(j),indz(k))
        ENDIF
      ENDDO
    ENDDO
  ENDDO
  CALL plag3d(rrr,ppp,zzz,fp,hrm1,hpm1,hzm1,xp,yp,zp &
       ,polylag,poly1x,poly1y,poly1z)
  Bzet = polylag
  dBzdR = poly1x
  dBzdp = poly1y
  dBzdZ = poly1z
  !
  RETURN
  END 
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
! dummies:
  MODULE magfie_mod
  END MODULE
!!$  MODULE neo_magfie_mod
!!$    INTEGER :: magfie_spline
!!$    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: magfie_sarray
!!$  END MODULE
  SUBROUTINE magfie_deallocate
  END
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
! 19.03.2010  SUBROUTINE stevvo(RT0,R0i,L1i,cbfi,BY0i,bf0)
!
! 19.03.2010  integer :: L1i,iunit1
! 19.03.2010  double precision :: RT0,R0i,cbfi,BY0i,bf0
!
! 19.03.2010  iunit1=77
! 19.03.2010  open(iunit1,file='MESH3D/stevvo_stuff.dat')
! 19.03.2010  read(iunit1,*) RT0
! 19.03.2010  read(iunit1,*) L1i
! 19.03.2010  read(iunit1,*) bf0
! 19.03.2010  read(iunit1,*) R0i
! 19.03.2010  read(iunit1,*) cbfi
! 19.03.2010  read(iunit1,*) BY0i
! 19.03.2010  close(iunit1)
!
! 19.03.2010  return
! 19.03.2010  end
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
! 19.03.2010  subroutine magfie(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
!
! Computes magnetic field module in units of the magnetic code  - bmod,
! square root of determinant of the metric tensor               - sqrtg,
! derivatives of the logarythm of the magnetic field module
! over coordinates                                              - bder,
! covariant componets of the unit vector of the magnetic
! field direction                                               - hcovar,
! contravariant components of this vector                       - hctrvr,
! contravariant component of the curl of this vector            - hcurl
! Order of coordinates is the following: x(1)=R (big radius), 
! x(2)=phi (toroidal angle), x(3)=Z (altitude).
!
!  Input parameters:
!            formal:  x                -    array of coordinates
!  Output parameters:
!            formal:  bmod
!                     sqrtg
!                     bder
!                     hcovar
!                     hctrvr
!                     hcurl
!
!  Called routines:  field
!
! 19.03.2010  double precision x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl
! 19.03.2010  double precision hr,hf,hz
!
! 19.03.2010  double precision ri,fii,zi,br,bf,bz,      &
! 19.03.2010      BRR,BRF,BRZ,BFR,BFF,BFZ,BZR,BZF,BZZ,  &
! 19.03.2010      BRK,BZK,BRRK,BRZK,BZRK,BZZK
!
! 19.03.2010  dimension x(3),bder(3),hcovar(3),hctrvr(3),hcurl(3)
!
! 19.03.2010  rbig=max(x(1),1d-12)
!
! 19.03.2010  ri=rbig
! 19.03.2010  fii=x(2)
! 19.03.2010  zi=x(3)
! 19.03.2010  call field(ri,fii,zi,br,bf,bz, &
! 19.03.2010             BRR,BRF,BRZ,BFR,BFF,BFZ,BZR,BZF,BZZ)
!
! 19.03.2010  bmod=dsqrt(br**2+bf**2+bz**2)
! 19.03.2010  sqrtg=rbig
! 19.03.2010  hr=br/bmod
! 19.03.2010  hf=bf/bmod
! 19.03.2010  hz=bz/bmod
!
! 19.03.2010  bder(1)=(brr*hr+bfr*hf+bzr*hz)/bmod
! 19.03.2010  bder(2)=(brf*hr+bff*hf+bzf*hz)/bmod
! 19.03.2010  bder(3)=(brz*hr+bfz*hf+bzz*hz)/bmod
!
! 19.03.2010  hcovar(1)=hr
! 19.03.2010  hcovar(2)=hf*rbig
! 19.03.2010  hcovar(3)=hz
!
! 19.03.2010  hctrvr(1)=hr
! 19.03.2010  hctrvr(2)=hf/rbig
! 19.03.2010  hctrvr(3)=hz
!
! 19.0 hcurl(1)=((bzf-rbig*bfz)/bmod+hcovar(2)*bder(3)-hcovar(3)*bder(2))/sqrtg
! 19.03  hcurl(2)=((brz-bzr)/bmod+hcovar(3)*bder(1)-hcovar(1)*bder(3))/sqrtg
!19 hcurl(3)=((bf+rbig*bfr-brf)/bmod+hcovar(1)*bder(2)-hcovar(2)*bder(1))/sqrtg
!
! 19.03.2010  return
! 19.03.2010  end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!

