  module chamb_new_mod
    logical :: prop=.true.
    integer :: iunit=71
    integer :: nr,np,nz
! 04.05.2014    double precision :: rmin,pmin,zmin,hr,hp,hz,pmax
    double precision :: rmin,pmin,zmin,hr,hp,hz,pmax,fac_ch
    logical, dimension(:,:,:), allocatable :: in_chamb
  end module chamb_new_mod
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine chamb(y,phi,ierr)
!
! checks whether the point is inside the vacuum chamber
!  Input parameters:
!            formal: y(i) - coordinates on the poloidal cross section
!                    phi  - toroidal angle
! Outout parameters:
!            formal: ierr -error code (0 if the point is inside 1 - othervice)
!
      use chamb_new_mod
!
      implicit none
!
!
      integer :: ierr,ir,ip,iz
      double precision :: phi
      double precision, dimension(2) :: y
! 23.11.2010
      logical first
      integer nier
! 04.05.2014
      double precision :: facch
! 05.05.2014      save nier
      save nier,facch
      data first/.true./
ierr=0
return
      if(first) then
      open(22,form='FORMATTED',file='chamb0.dat')
! 04.05.2014
      facch=fac_ch
      print*,facch,'=facch'
      write(22,*)facch,'=facch'
! 04.05.2014 end
      nier=0
      first=.false.
      endif
! 23.11.2010 end
!
      if(prop) then
        prop=.false.
        open(iunit,form='unformatted',file='chamb.dat')
        read (iunit) nr,np,nz
        read (iunit) rmin,pmin,zmin,hr,hp,hz,pmax
        allocate(in_chamb(nr,np,nz))
        read (iunit) in_chamb
        close(iunit)
      endif
!
      ip=nint(modulo(phi-pmin,pmax-pmin)/hp)
      if(ip.le.0) then 
        ip=np
      elseif(ip.gt.np) then
        ip=1
      endif
      ir=int((y(1)-rmin)/hr)+1
! 04.05.2014      iz=int((y(2)-zmin)/hz)+1
      iz=int((facch*y(2)-zmin)/hz)+1
!
      if(ir.lt.1.or.ir.gt.nr.or.iz.lt.1.or.iz.gt.nz) then
        ierr=1
        return
      endif
!
      if(in_chamb(ir,ip,iz)) then
        ierr=0
      else
        ierr=1
      endif
! 23.11.2010
      if(ierr.ne.0) then
      nier=nier+1
      print*,y(1),y(2),nier,'=nier, phi=',phi
      write(22,*)y(1),y(2),nier,'=nier, phi=',phi
      endif
! 23.11.2010 end
!
      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
