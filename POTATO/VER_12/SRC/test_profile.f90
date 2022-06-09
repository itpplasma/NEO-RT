program test_profile
    use phielec_of_psi_mod, only : npolyphi, polyphi, polydens, polytemp
    use field_eq_mod, only : psif,dpsidr,dpsidz,psi_axis,psi_sep
    implicit none

    integer, parameter :: nsurf = 100
    integer :: iunit, i
    double precision :: bmod,sqrtg
    double precision, dimension(3) :: x,bder,hcovar,hctrvr,hcurl
    double precision :: psi, phi_elec, dens, temp, dPhi_dpsi, ddens, dtemp

    associate ( E_alpha => 5d3, A_alpha => 2d0, Z_alpha => 1d0, &
        e_charge => 4.8032d-10, ev => 1.6022d-12 )

    open(newunit=iunit,file='profile_poly.in')
        read(iunit,*)  ! dummy
        read(iunit,*)  ! dummy
        read(iunit,*) polydens
        read(iunit,*) polytemp  ! dummy
        read(iunit,*) polytemp
        read(iunit,*) polyphi
    close(iunit)
    polytemp = polytemp/E_alpha
    polyphi = polyphi*Z_alpha*e_charge/(E_alpha*ev)

    x(1) = 170.0d0
    x(2) = 0.0d0
    x(3) = 0.0d0
    call magfie(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
    close(iunit)

    open(newunit=iunit,file='profile.out')
        do i=1,nsurf
            psi = psi_axis + i*1d0/nsurf*(psi_sep-psi_axis)
            call phielec_of_psi(psi, phi_elec, dPhi_dpsi)
            call denstemp_of_psi(psi, dens, temp, ddens, dtemp)
            write(iunit, *) psi, phi_elec, dens, temp, dPhi_dpsi, ddens, dtemp
        end do
    end associate
    close(iunit)
end program test_profile
