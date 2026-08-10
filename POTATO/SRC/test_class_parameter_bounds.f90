program test_class_parameter_bounds
    implicit none

    double precision, parameter :: relmargin=1.d-7
    double precision :: xbeg,xend,xi,dxi_dx
    integer :: map_status
    external :: classbounds,xi_func

    call check_margin(12,1.d0)
    call check_margin(21,1.d0)
    call check_margin(22,1.d0)
    call check_margin(23,1.d0)
    call check_margin(24,1.d0)
    call check_margin(32,1.d0)
    call check_margin(42,1.d0)

contains

    subroutine check_margin(ifuntype,widthclass)
        integer, intent(in) :: ifuntype
        double precision, intent(in) :: widthclass
        double precision :: xleft,xright

        call classbounds(ifuntype,relmargin,widthclass,xleft,xright, &
                         map_status)
        if(map_status.ne.0) error stop 'generated class bounds rejected valid input'
        if(ifuntype.eq.12 .or. ifuntype.eq.32 .or. ifuntype.eq.42) then
            call xi_func(ifuntype,xright,xi,dxi_dx,relmargin,widthclass, &
                         map_status)
            if(abs(1.d0-xi-relmargin).gt.1.d-12) then
                error stop 'inner endpoint margin is not physical for chart'
            endif
        elseif(ifuntype.eq.22) then
            call xi_func(ifuntype,xleft,xi,dxi_dx,relmargin,widthclass, &
                         map_status)
            if(abs(xi-relmargin).gt.1.d-12) then
                error stop 'left inner endpoint margin is not physical for chart'
            endif
            call xi_func(ifuntype,xright,xi,dxi_dx,relmargin,widthclass, &
                         map_status)
            if(abs(1.d0-xi-relmargin).gt.1.d-12) then
                error stop 'right inner endpoint margin is not physical for chart'
            endif
        else
            call xi_func(ifuntype,xleft,xi,dxi_dx,relmargin,widthclass, &
                         map_status)
            if(abs(xi-relmargin).gt.1.d-12) then
                error stop 'inner endpoint margin is not physical for chart'
            endif
        endif
    end subroutine check_margin

end program test_class_parameter_bounds

subroutine phielec_of_psi(psi,phi_elec,dPhi_dpsi)
    double precision, intent(in) :: psi
    double precision, intent(out) :: phi_elec,dPhi_dpsi

    phi_elec=0.d0
    dPhi_dpsi=0.d0
end subroutine phielec_of_psi

subroutine denstemp_of_psi(psi,dens,temp,ddens,dtemp)
    double precision, intent(in) :: psi
    double precision, intent(out) :: dens,temp,ddens,dtemp

    dens=0.d0
    temp=1.d0
    ddens=0.d0
    dtemp=0.d0
end subroutine denstemp_of_psi
