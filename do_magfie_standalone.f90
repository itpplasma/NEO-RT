module do_magfie_mod

  use common
  use spline

  implicit none
  save

  real(8), public :: s, psi_pr, Bthcov, Bphcov, dBthcovds, dBphcovds,&
       q, dqds, iota, R0, a, eps, B0h, B00
  real(8), public :: bfac = 1.0d0
  ! B0h is the 0th theta harmonic of bmod on current flux surface
  ! and B00 the 0th theta harmonic of bmod on the innermost flux surface

  real(8), allocatable, protected :: params0(:,:), modes0(:,:,:)
  integer, protected :: m0b, n0b, nflux, nfp, nmode

  real(8), allocatable, protected :: spl_coeff1(:,:,:), spl_coeff2(:,:,:,:)
  real(8), allocatable, protected :: eps_spl(:,:)

  real(8), parameter :: ItoB = 2.0d-1 ! Factor for covar. field (cgs) from I(SI)
                                      ! Bcov=mu0/2pi*I,mu0->4pi/c,I->10^(-1)*c*I
  integer :: ncol1, ncol2 ! number of columns in input file

  integer :: inp_swi ! type of input file
contains

  subroutine do_magfie_init
  ! Initializes spline and Fourier coefficients for later evaluation of
  ! unperturbed axisymmetric magnetic field in do_magfie
    integer :: j,k
    real(8)                      :: x(3)
    real(8)                      :: bmod
    real(8)                      :: sqrtg
    real(8), dimension(size(x))   :: bder
    real(8), dimension(size(x))   :: hcovar
    real(8), dimension(size(x))   :: hctrvr
    real(8), dimension(size(x))   :: hcurl

    ncol1 = 5
    if(inp_swi == 8) ncol2=4 ! tok_circ
    if(inp_swi == 9) ncol2=8 ! ASDEX
    call boozer_read('in_file') ! TODO: general filename

    if (.not. allocated(spl_coeff1)) then
       allocate(spl_coeff1(nflux-1, 5, ncol1))
       allocate(spl_coeff2(nflux-1, 5, ncol2, nmode))
       allocate(eps_spl(nflux-1, 5))
    end if

    B00 = 1.0d4*modes0(1,1,6)*bfac

    ! calculate spline coefficients
    do k=1,ncol1
       ! first column is s, so start with second column
       spl_coeff1(:, :, k) = spline_coeff(params0(:,1), params0(:,k+1))
    end do

    do j=1,nmode
       do k=1,ncol2
          ! first two columns are mode numbers, so start with third column
          spl_coeff2(:, :, k, j) = spline_coeff(params0(:,1), modes0(:,j,k+2))
       end do
    end do

    if (inp_swi == 8) eps_spl = spline_coeff(params0(:,1), abs(modes0(:,2,6)/modes0(:,1,6)))
    if (inp_swi == 9) eps_spl = spline_coeff(params0(:,1), abs(modes0(:,2,9)/modes0(:,1,9)))

    x(1) = s
    x(2) = 0.0
    x(3) = 0.0
    call do_magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)

  end subroutine do_magfie_init

  subroutine do_magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
  ! Evaluate unperturbed axisymmetric magnetic field in point x = (s, ph, th)
    real(8), dimension(:),       intent(in)         :: x
    real(8),                     intent(out)        :: bmod
    real(8),                     intent(out)        :: sqrtg
    real(8), dimension(size(x)), intent(out)        :: bder
    real(8), dimension(size(x)), intent(out)        :: hcovar
    real(8), dimension(size(x)), intent(out)        :: hctrvr
    real(8), dimension(size(x)), intent(out)        :: hcurl

    integer :: j
    real(8) :: spl_val(3), spl_val_c(3), spl_val_s(3)
    real(8) :: B0mnc(nmode), dB0dsmnc(nmode), B0mns(nmode), dB0dsmns(nmode)
    real(8) :: sqgbmod2 ! sqg*B^2

    real(8) :: x1

    ! safety measure in order not to extrapolate
    x1 = max(params0(1,1),x(1))
    x1 = min(params0(nflux,1),x1)

    spl_val = spline_val_0(abs(spl_coeff1(:,:,3)), x1)
    Bthcov = -ItoB*spl_val(1)*bfac
    dBthcovds = -ItoB*spl_val(2)*bfac
    spl_val = spline_val_0(abs(spl_coeff1(:,:,2)), x1)
    Bphcov = -ItoB*spl_val(1)*bfac
    dBphcovds = -ItoB*spl_val(2)*bfac
    spl_val = spline_val_0(abs(spl_coeff1(:,:,1)), x1)
    iota = spl_val(1)
    q = 1/iota
    dqds = -spl_val(2)/iota**2
    ! spl_val = spline_val_0(eps_spl, x1)
    ! eps = spl_val(1)

    ! calculate B-field from modes
    if (inp_swi == 8) then
       do j=1,nmode
          spl_val_c = spline_val_0(spl_coeff2(:, :, 4, j), x1)
          B0mnc(j) = 1d4*spl_val_c(1)*bfac
          dB0dsmnc(j) = 1d4*spl_val_c(2)*bfac
       end do
       B0h = B0mnc(1)

       bmod = sum(B0mnc*cos(modes0(1,:,1)*x(3)))
       sqrtg = psi_pr*(iota*Bthcov + Bphcov)/bmod**2
       bder(1) = sum(dB0dsmnc*cos(modes0(1,:,1)*x(3)))/bmod
       bder(2) = 0d0 ! TODO 3: toroidal symmetry assumed
       bder(3) = sum(-modes0(1,:,1)*B0mnc*sin(modes0(1,:,1)*x(3)))/bmod
    else if (inp_swi == 9) then
       Bthcov = -Bthcov
       dBthcovds = -dBthcovds
       Bphcov = -Bphcov
       dBphcovds = -dBphcovds
       do j=1,nmode
          spl_val_c = spline_val_0(spl_coeff2(:, :, 7, j), x1)
          B0mnc(j) = 1d4*spl_val_c(1)*bfac
          dB0dsmnc(j) = 1d4*spl_val_c(2)*bfac
          spl_val_s = spline_val_0(spl_coeff2(:, :, 8, j), x1)
          B0mns(j) = 1d4*spl_val_s(1)*bfac
          dB0dsmns(j) = 1d4*spl_val_s(2)*bfac
       end do
       B0h = B0mnc(1)

       bmod = sum(B0mnc*cos(modes0(1,:,1)*x(3))+B0mns*sin(modes0(1,:,1)*x(3)))
       sqrtg = psi_pr*(iota*Bthcov + Bphcov)/bmod**2
       bder(1) = sum(dB0dsmnc*cos(modes0(1,:,1)*x(3))+dB0dsmns*sin(modes0(1,:,1)*x(3)))/bmod
       bder(2) = 0d0 ! TODO 3: toroidal symmetry assumed
       bder(3) = sum(-modes0(1,:,1)*B0mnc*sin(modes0(1,:,1)*x(3))&
            +modes0(1,:,1)*B0mns*cos(modes0(1,:,1)*x(3)))/bmod
    end if

    sqgbmod2 = psi_pr*(Bphcov+iota*Bthcov)

    hcovar(1) = 0d0 ! TODO 2: get this from geometry
    hcovar(2) = Bphcov/bmod
    hcovar(3) = Bthcov/bmod

    hctrvr(1) = 0d0
    hctrvr(2) = psi_pr/(sqrtg*bmod)
    hctrvr(3) = iota*psi_pr/(sqrtg*bmod)

    ! Common factor bmod since bder(k) = (dB/dx(k))/bmod
    hcurl(1) = (bder(2)*Bthcov - bder(3)*Bphcov)*bmod/sqgbmod2
    hcurl(3) = (-dBphcovds+bder(1)*Bphcov)*bmod/sqgbmod2 ! TODO: B_s
    hcurl(2) = (dBthcovds-bder(1)*Bthcov)*bmod/sqgbmod2 ! TODO: B_s

  end subroutine do_magfie

  subroutine boozer_read(filename)
    ! Reads Boozer in_file and converts SI to CGS

    integer :: ksurf, kmode
    real(8) :: flux
    character(len=*) :: filename
    open (unit=18, file=filename)
    read(18, '(////)')
    read(18, *) m0b, n0b, nflux, nfp, flux, a, R0
    a = 100*a   ! m -> cm
    R0 = 100*R0 ! m -> cm

    ! psi_pr = dpsi_tor/ds = psi_tor/s is constantly spaced
    ! psi_tor = psi_pr*s
    ! A_theta = -psi_tor
    psi_pr = 1.0d8*abs(flux)/(2*pi)*bfac ! T -> Gauss, m -> cm

    nmode = (m0b+1)*(n0b+1)
    allocate(params0(nflux, ncol1+1))
    allocate(modes0(nflux, nmode, ncol2+2))
    do ksurf=1,nflux
       read(18, '(/)')
       read(18, *) params0(ksurf,:)
       read(18, *)
       do kmode=1,nmode
          read(18, *) modes0(ksurf,kmode,:)
       end do
    end do
    close(unit=18)
    ! Set R0 to first harmonic
    R0 = modes0(1,1,3)*100
  end subroutine boozer_read

  subroutine booz_to_cyl(x, r)

   real(8), intent(in)  :: x(3)  ! Boozer coordinates (s, ph, th)
   real(8), intent(out) :: r(3)  ! Cylindrical coordinates (R, phi, Z)

   real(8) :: spl_val(3), x1
   real(8) :: rmnc(nmode), rmns(nmode), zmnc(nmode), zmns(nmode)

   integer :: j

   if(inp_swi /= 9) error stop  ! Only implemented for ASDEX-U type of data

   x1 = max(params0(1,1),x(1))
   x1 = min(params0(nflux,1),x1)

   do j=1,nmode
      spl_val = spline_val_0(spl_coeff2(:, :, 1, j), x1)
      rmnc(j) = 1.0d2*spl_val(1)
      spl_val = spline_val_0(spl_coeff2(:, :, 2, j), x1)
      rmns(j) = 1.0d2*spl_val(1)
      spl_val = spline_val_0(spl_coeff2(:, :, 3, j), x1)
      zmnc(j) = 1.0d2*spl_val(1)
      spl_val = spline_val_0(spl_coeff2(:, :, 4, j), x1)
      zmns(j) = 1.0d2*spl_val(1)
   end do

   r(1) = sum(rmnc*cos(modes0(1,:,1)*x(3))+rmns*sin(modes0(1,:,1)*x(3)))
   r(2) = 0.0d0  ! TODO: phi
   r(3) = sum(zmnc*cos(modes0(1,:,1)*x(3))+zmns*sin(modes0(1,:,1)*x(3)))

  end subroutine booz_to_cyl

end module do_magfie_mod

module do_magfie_pert_mod

  use common
  use spline
  use do_magfie_mod, only: s, bfac, inp_swi

  implicit none
  save

  real(8), allocatable, protected :: params(:,:),  modes(:,:,:)
  integer, protected :: mb, nb, nflux, nfp, nmode

  real(8), allocatable, protected :: spl_coeff1(:,:,:), spl_coeff2(:,:,:,:)

  integer :: ncol1, ncol2 ! number of columns in input file
  real(8) :: mph ! toroidal perturbation mode
contains

  subroutine do_magfie_pert_init
    integer :: j,k

    ncol1 = 5
    if(inp_swi == 8) ncol2=4 ! tok_circ
    if(inp_swi == 9) ncol2=8 ! ASDEX
    call boozer_read_pert('in_file_pert') ! TODO: general filename

    mph = nfp*modes(1,1,2)

    if (.not. allocated(spl_coeff1)) then
       allocate(spl_coeff1(nflux-1, 5, ncol1))
       allocate(spl_coeff2(nflux-1, 5, ncol2, nmode))
    end if

    ! calculate spline coefficients
    do k=1,ncol1
       ! first column is s, so start with second column
       spl_coeff1(:, :, k) = spline_coeff(params(:,1), params(:,k+1))
    end do

    do j=1,nmode
       do k=1,ncol2
          ! first two columns are mode numbers, so start with third column
          spl_coeff2(:, :, k, j) = spline_coeff(params(:,1), modes(:,j,k+2))
       end do
    end do

  end subroutine do_magfie_pert_init

  subroutine do_magfie_pert_amp(x, bamp)
    real(8), dimension(:),       intent(in)         :: x
    complex(8),                  intent(out)        :: bamp

    integer :: j
    real(8) :: spl_val_c(3), spl_val_s(3)
    real(8) :: Bmnc(nmode), Bmns(nmode)
    real(8) :: x1

    ! safety measure in order not to extrapolate
    x1   = max(params(1,1),x(1))
    x1   = min(params(nflux,1),x1)

    ! calculate B-field from modes
    if (inp_swi == 8) then
        do j=1,nmode
          spl_val_c = spline_val_0(spl_coeff2(:, :, 4, j), x1)
          Bmnc(j) = 1d4*spl_val_c(1)*bfac
       end do
       bamp = sum(Bmnc*cos(modes(1,:,1)*x(3)))
    else if (inp_swi == 9) then
       do j=1,nmode
          spl_val_c = spline_val_0(spl_coeff2(:, :, 7, j), x1)
          Bmnc(j) = 1d4*spl_val_c(1)*bfac
          spl_val_s = spline_val_0(spl_coeff2(:, :, 8, j), x1)
          Bmns(j) = 1d4*spl_val_s(1)*bfac
       end do

       bamp = sum((Bmnc-imun*Bmns)*exp(imun*modes(1,:,1)*x(3)))
    end if
  end subroutine do_magfie_pert_amp

  subroutine do_magfie_pert(x, bmod)
   real(8), dimension(:),       intent(in)         :: x
   real(8),                     intent(out)        :: bmod
   complex(8)                                      :: bamp

   call do_magfie_pert_amp(x, bamp)
   bmod = real(sum(bamp*exp(imun*nfp*modes(1,:,2)*x(2))))

  end subroutine do_magfie_pert

  subroutine boozer_read_pert(filename)
    integer :: ksurf, kmode
    real(8) :: flux, dummy
    character(len=*) :: filename
    character(len=128) :: buffer
    open (unit=18, file=filename)
    read(18, '(////)')
    read(18, *) mb, nb, nflux, nfp, flux, dummy, dummy
    nmode = (mb+1)*(nb+1)
    allocate(params(nflux, ncol1+1))
    allocate(modes(nflux, nmode, ncol2+2))
    do ksurf=1,nflux
       read(18, '(/)')
       read(18, *) params(ksurf,:)
       read(18, *)
       do kmode=1,nmode
          read(18, *) modes(ksurf,kmode,:)
       end do
    end do
    close(unit=18)
  end subroutine boozer_read_pert
end module do_magfie_pert_mod
