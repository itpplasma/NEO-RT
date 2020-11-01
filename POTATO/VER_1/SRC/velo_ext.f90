!------------------------------------------------------
!
  module velo_ext_mod
    integer, parameter  :: neqm=5, next=3, ndim=neqm+next
    logical             :: write_orb=.false.
    integer             :: iunit1

    contains

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !
    subroutine velo_ext(dtau,z,vz)
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

  end module velo_ext_mod
