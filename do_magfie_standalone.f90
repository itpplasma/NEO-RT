module do_magfie_mod
  
  use common
  use spline
  
  implicit none
  save

  real(8), public :: s, psi_pr, Bthcov, Bphcov, dBthcovds, dBphcovds,&
       q, dqds, iota, R0, a, eps, B0, B00
  
  real(8), protected :: config0(6), config(6)
  real(8), allocatable, protected :: params0(:,:), modes0(:,:,:)  
  real(8), allocatable, protected :: params(:,:),  modes(:,:,:)
  integer, protected :: m0b, n0b, nflux, nfp, nmode
  
  real(8), allocatable, protected :: spl_coeff1(:,:,:), spl_coeff2(:,:,:,:)
  real(8), allocatable, protected :: eps_spl(:,:)

  real(8), parameter :: ItoB = 2.0d-1 ! Factor for covar. field (cgs) from I(SI)
                                      ! Bcov=mu0/2pi*I,mu0->4pi/c,I->10^(-1)*c*I
contains
  
  subroutine do_magfie_init
    integer :: j,k

    call boozer_read('tok-synch2-n0.bc') ! TODO: general filename

    if (.not. allocated(spl_coeff1)) then
       allocate(spl_coeff1(nflux-1, 5, 5))
       allocate(spl_coeff2(nflux-1, 5, 4, nmode))
       allocate(eps_spl(nflux-1, 5))
    end if
    
    B00 = 1.0d4*modes0(1,1,6)
    
    ! calculate spline coefficients
    do k=1,5
       ! first column is s, so start with second column
       spl_coeff1(:, :, k) = spline_coeff(params0(:,1), params0(:,k+1))
    end do
    
    do j=1,nmode
       do k=1,4
          ! first two columns are mode numbers, so start with third column
          spl_coeff2(:, :, k, j) = spline_coeff(params0(:,1), modes0(:,j,k+2))
       end do
    end do

    eps_spl = spline_coeff(params0(:,1), abs(modes0(:,2,6)/modes0(:,1,6)))

  end subroutine do_magfie_init
  
  subroutine do_magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
    real(8), dimension(:),       intent(in)         :: x
    real(8),                     intent(out)        :: bmod
    real(8),                     intent(out)        :: sqrtg
    real(8), dimension(size(x)), intent(out)        :: bder
    real(8), dimension(size(x)), intent(out)        :: hcovar
    real(8), dimension(size(x)), intent(out)        :: hctrvr
    real(8), dimension(size(x)), intent(out)        :: hcurl

    integer :: j
    real(8) :: spl_val(3)
    real(8) :: B0mn(nmode), dB0dsmn(nmode)

    spl_val = spline_val_0(spl_coeff1(:,:,3), s)
    Bthcov = -ItoB*spl_val(1)
    dBthcovds = -ItoB*spl_val(2)
    spl_val = spline_val_0(spl_coeff1(:,:,2), s)
    Bphcov = -ItoB*spl_val(1)
    dBphcovds = -ItoB*spl_val(2)
    spl_val = spline_val_0(spl_coeff1(:,:,1), s)
    iota = spl_val(1)
    q = 1/iota
    dqds = -spl_val(2)/iota**2
    spl_val = spline_val_0(eps_spl, s)
    eps = spl_val(1)

    ! calculate B-field from modes
    do j=1,nmode
       spl_val = spline_val_0(spl_coeff2(:, :, 4, j), s)
       B0mn(j) = 1d4*spl_val(1)
       dB0dsmn(j) = 1d4*spl_val(2)
    end do
    B0 = B0mn(1)
    
    bmod = sum(B0mn*cos(modes0(1,:,1)*x(3)))
    sqrtg = psi_pr*(iota*Bthcov + Bphcov)/bmod**2
    bder(1) = sum(dB0dsmn*cos(modes0(1,:,1)*x(3)))/bmod
    bder(2) = 0d0 ! TODO 3: toroidal symmetry assumed 
    bder(3) = sum(-modes0(1,:,1)*B0mn*sin(modes0(1,:,1)*x(3)))/bmod 
    hcovar(1) = 0d0 ! TODO 2: get this from geometry
    hcovar(2) = Bphcov/bmod
    hcovar(3) = Bthcov/bmod
    hctrvr(1) = 0d0
    hctrvr(2) = psi_pr/(sqrtg*bmod)
    hctrvr(3) = iota*psi_pr/(sqrtg*bmod)
    hcurl(1) = 0d0 ! TODO
    hcurl(2) = 0d0 ! TODO
    hcurl(3) = 0d0 ! TODO
  end subroutine do_magfie  
  
  subroutine boozer_read(filename)
    integer :: ksurf, kmode
    real(8) :: flux
    character(len=*) :: filename
    open (unit=18, file=filename)
    read(18, '(////)')
    read(18, *) m0b, n0b, nflux, nfp, flux, a, R0
    a = 100*a   ! m -> cm
    R0 = 100*R0 ! m -> cm
    psi_pr = 1.0d8*abs(flux)/(2*pi) ! T -> Gauss, m -> cm
    nmode = (m0b+1)*(n0b+1)
    allocate(params0(nflux, 6))
    allocate(modes0(nflux, nmode, 6))
    do ksurf=1,nflux
       read(18, '(/)')
       read(18, *) params0(ksurf,:)
       read(18, *)
       do kmode=1,nmode
          read(18, *) modes0(ksurf,kmode,:)
       end do
    end do
    ! Set R0 to first harmonic instead
    R0 = modes0(1,1,3)*100
  end subroutine boozer_read
  
end module do_magfie_mod
