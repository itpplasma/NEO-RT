program test_reslines
    use iso_fortran_env, only: dp => real64
    use util
    use do_magfie_mod, only: s, R0, inp_swi, do_magfie_init
    use do_magfie_pert_mod, only: mph
    use neort, only: init, set_to_trapped_region, set_to_passing_region
    use driftorbit, only: mth, M_t, vth, &
        etadt, etatp, Om_tE, sign_vpar, etamin, etamax, nlev
    use neort_resonance, only: driftorbit_coarse, driftorbit_root
    implicit none

    real(dp), parameter  :: vnorm = sqrt(0.15_dp)

    integer :: fid, ios, ns
    integer :: ks
    real(dp) :: row(3)
    real(dp), allocatable :: profile(:,:)

    integer, parameter :: nm2 = 11, nm3 = 2
    integer :: km2, km3
    integer :: m2(nm2) = [-5,-4,-3,-2,-1,0,1,2,3,4,5]
    integer :: m3(nm3) = [-1,1]

    inp_swi = 9 ! ASDEX Upgrade format
    call do_magfie_init("in_file")

    open(newunit=fid, file='profile.in', action='read', status='old', iostat=ios)
    if (ios /= 0) error stop 'could not open profile.in'
    ns = 0
    do
        read(fid,*,iostat=ios) row
        if (ios > 0) error stop 'invalid row in profile.in'
        if (ios /= 0) exit
        ns = ns + 1
    end do
    if (ns < 1) error stop 'profile.in has no rows'
    rewind(fid)
    allocate(profile(ns, 3))
    do ks=1,ns
        read(fid,*,iostat=ios) profile(ks,:)
        if (ios /= 0) error stop 'profile.in ended before all rows were read'
    end do
    close(fid)

    open(newunit=fid, file='reslines.out', action='write')
    write(fid,*)  '%# s, v/vth, eta [1/G], m2, m3, root number, class'
    do ks=1, ns
        s = profile(ks,1)
        M_t = profile(ks,2)
        vth = profile(ks,3)
        print *
        print *, ks, '/', ns
        print *, 's = ', s, 'M_t = ', M_t, 'vth = ', vth
        Om_tE = vth*M_t/R0 ! ExB precession frequency
        call init
        do km3=1,nm3
            mph = m3(km3)
            do km2=1,nm2
                mth = m2(km2)
                call resline
            end do
        end do
    end do
    close(fid)

contains

    subroutine resline

        real(dp) :: v, eta_res(2)
        real(dp) :: roots(nlev, 3)
        integer :: nroots, kr

        v = vnorm*vth

        ! Trapped
        sign_vpar = 1.0_dp
        call set_to_trapped_region(etamin, etamax)
        call driftorbit_coarse(v, etamin, etamax, roots, nroots)
        if(nroots == 0) return
        do kr = 1,nroots
            eta_res = driftorbit_root(v, 1.0e-8_dp*abs(Om_tE), roots(kr,1), roots(kr,2))
            write(fid,*) s, v/vth, eta_res(1), mth, mph, kr, 0
        end do
        ! Co-passing
        sign_vpar = 1.0_dp
        call set_to_passing_region(etamin, etamax)
        call driftorbit_coarse(v, etamin, etamax, roots, nroots)
        if(nroots == 0) return
        do kr = 1,nroots
            eta_res = driftorbit_root(v, 1.0e-8_dp*abs(Om_tE), roots(kr,1), roots(kr,2))
            write(fid,*) s, v/vth, eta_res(1), mth, mph, kr, 1
        end do
        ! Ctr-passing
        sign_vpar = -1.0_dp
        call set_to_passing_region(etamin, etamax)
        call driftorbit_coarse(v, etamin, etamax, roots, nroots)
        if(nroots == 0) return
        do kr = 1,nroots
            eta_res = driftorbit_root(v, 1.0e-8_dp*abs(Om_tE), roots(kr,1), roots(kr,2))
            write(fid,*) s, v/vth, eta_res(1), mth, mph, kr, -1
        end do

    end subroutine resline

end program test_reslines
