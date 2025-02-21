program test_reslines
  use util
  use do_magfie_mod, only: s, psi_pr, Bthcov, Bphcov, dBthcovds, dBphcovds,&
  q, dqds, iota, R0, a, eps, inp_swi, do_magfie_init, do_magfie
  use do_magfie_pert_mod, only: mph
  use driftorbit, only: init, epsmn, mth, M_t, vth, qi, mi, Jperp, &
    etadt, etatp, Om_th, Om_ph, Om_tE, sigv, &
    epst, epsp, etamin, etamax, nlev, driftorbit_coarse, driftorbit_root

  implicit none

  real(8), parameter  :: vnorm = sqrt(0.15)

  integer :: fid

  integer, parameter :: ns = 250
  integer :: ks
  real(8) :: profile(ns, 3)

  integer, parameter :: nm2 = 11, nm3 = 2
  integer :: km2, km3
  integer :: m2(nm2) = [-5,-4,-3,-2,-1,0,1,2,3,4,5]
  integer :: m3(nm3) = [-1,1]

  inp_swi = 9  ! ASDEX Upgrade format
  call do_magfie_init

  open(newunit=fid, file='profile.in', action='read')
    do ks=1,ns
      read(fid,*) profile(ks,:)
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
      Om_tE = vth*M_t/R0  ! ExB precession frequency
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

    real(8) :: v, eta, eta_res(2)
    real(8) :: roots(nlev, 3)
    integer :: nroots, kr

    v = vnorm*vth

    ! Trapped
    sigv = 1d0
    etamin = (1+epst)*etatp
    etamax = (1-epst)*etadt
    call driftorbit_coarse(v, etamin, etamax, roots, nroots)
    if(nroots == 0) return
    do kr = 1,nroots
      eta_res = driftorbit_root(v, 1d-8*abs(Om_tE), roots(kr,1), roots(kr,2))
      write(fid,*) s, v/vth, eta_res(1), mth, mph, kr, 0
      end do
    ! Co-passing
    sigv = 1d0
    etamin = epsp*etatp
    etamax = (1-epsp)*etatp
    call driftorbit_coarse(v, etamin, etamax, roots, nroots)
    if(nroots == 0) return
    do kr = 1,nroots
      eta_res = driftorbit_root(v, 1d-8*abs(Om_tE), roots(kr,1), roots(kr,2))
      write(fid,*) s, v/vth, eta_res(1), mth, mph, kr, 1
      end do
    ! Ctr-passing
    sigv = -1d0
    etamin = epsp*etatp
    etamax = (1-epsp)*etatp
    call driftorbit_coarse(v, etamin, etamax, roots, nroots)
    if(nroots == 0) return
    do kr = 1,nroots
      eta_res = driftorbit_root(v, 1d-8*abs(Om_tE), roots(kr,1), roots(kr,2))
      write(fid,*) s, v/vth, eta_res(1), mth, mph, kr, -1
    end do

  end subroutine resline

end program test_reslines
