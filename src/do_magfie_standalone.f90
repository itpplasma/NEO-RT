module do_magfie_mod

    use util
    use spline

    implicit none

    real(8), private :: s_prev = -1.0d0
    real(8), private, allocatable :: spl_val_c(:,:), spl_val_s(:,:)
    ! Work arrays previously automatic on stack (size=nmode)
    real(8), private, allocatable :: B0mnc(:), dB0dsmnc(:), B0mns(:), dB0dsmns(:)
    real(8), private, allocatable :: costerm(:), sinterm(:)

    real(8), parameter :: sign_theta = -1.0d0  ! negative for left-handed

    real(8) :: s = 0d0, psi_pr = 0d0, Bthcov = 0d0, Bphcov = 0d0, &
                       dBthcovds = 0d0, dBphcovds = 0d0, &
                       q = 0d0, dqds = 0d0, iota = 0d0, R0 = 0d0, a = 0d0, &
                       eps = 0d0, B0h = 0d0, B00 = 0d0
    real(8) :: bfac = 1.0d0
    ! B0h is the 0th theta harmonic of bmod on current flux surface
    ! and B00 the 0th theta harmonic of bmod on the innermost flux surface

    real(8), allocatable, protected :: params0(:, :), modes0(:, :, :)
    integer, protected :: m0b, n0b, nflux, nfp, nmode

    real(8), allocatable, protected :: spl_coeff1(:, :, :), spl_coeff2(:, :, :, :)
    ! Work arrays for booz_to_cyl (size=nmode)
    real(8), private, allocatable :: rmnc(:), rmns(:), zmnc(:), zmns(:)

    real(8), parameter :: ItoB = 2.0d-1*sign_theta ! Covarient B (cgs) from I (SI)
    ! Bcov=mu0/2pi*I,mu0->4pi/c,I->10^(-1)*c*I

    integer :: ncol1 = 0, ncol2 = 0 ! number of columns in input file

    integer :: inp_swi = 0 ! type of input file
    contains

    subroutine do_magfie_init(path)
        character(len=*), intent(in) :: path

        ! Initializes spline and Fourier coefficients for later evaluation of
        ! unperturbed axisymmetric magnetic field in do_magfie
        integer :: j, k
        real(8) :: x(3)
        real(8) :: bmod
        real(8) :: sqrtg
        real(8), dimension(size(x)) :: bder
        real(8), dimension(size(x)) :: hcovar
        real(8), dimension(size(x)) :: hctrvr
        real(8), dimension(size(x)) :: hcurl

        ncol1 = 5
        if (inp_swi == 8) ncol2 = 4 ! tok_circ
        if (inp_swi == 9) ncol2 = 8 ! ASDEX
        call boozer_read(path)

        ! Allocate spline coefficient arrays (deallocate first if size changed)
        if (allocated(spl_coeff1)) then
            if (size(spl_coeff1, 1) /= nflux - 1) deallocate(spl_coeff1, spl_coeff2)
        end if
        if (.not. allocated(spl_coeff1)) then
            allocate (spl_coeff1(nflux - 1, 5, ncol1))
            allocate (spl_coeff2(nflux - 1, 5, ncol2, nmode))
        end if

        ! Allocate work arrays (deallocate first if size changed)
        if (allocated(B0mnc)) then
            if (size(B0mnc) /= nmode) then
                deallocate(B0mnc, dB0dsmnc, costerm, sinterm)
                if (allocated(B0mns)) deallocate(B0mns, dB0dsmns)
            end if
        end if
        if (.not. allocated(B0mnc)) then
            allocate(B0mnc(nmode), dB0dsmnc(nmode))
            if (ncol2 >= 8) allocate(B0mns(nmode), dB0dsmns(nmode))
            allocate(costerm(nmode), sinterm(nmode))
        end if

        if (allocated(rmnc)) then
            if (size(rmnc) /= nmode) deallocate(rmnc, rmns, zmnc, zmns)
        end if
        if (.not. allocated(rmnc)) then
            allocate(rmnc(nmode), rmns(nmode), zmnc(nmode), zmns(nmode))
        end if

        B00 = 1.0d4*modes0(1, 1, 6)*bfac

        ! calculate spline coefficients
        do k = 1, ncol1
            ! first column is s, so start with second column
            spl_coeff1(:, :, k) = spline_coeff(params0(:, 1), params0(:, k + 1))
        end do

        do j = 1, nmode
            do k = 1, ncol2
                ! first two columns are mode numbers, so start with third column
                spl_coeff2(:, :, k, j) = spline_coeff(params0(:, 1), modes0(:, j, k + 2))
            end do
        end do

        ! Allocate cached spline values (deallocate first if size changed)
        if (allocated(spl_val_c)) then
            if (size(spl_val_c, 2) /= nmode) deallocate(spl_val_c, spl_val_s)
        end if
        if (.not. allocated(spl_val_c)) then
            allocate (spl_val_c(3, nmode))
            allocate (spl_val_s(3, nmode))
        end if

        x(1) = s
        x(2) = 0.0
        x(3) = 0.0
        call do_magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
    end subroutine do_magfie_init

    subroutine do_magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
        ! Evaluate unperturbed axisymmetric magnetic field in point x = (s, ph, th)
        real(8), dimension(:), intent(in) :: x
        real(8), intent(out) :: bmod
        real(8), intent(out) :: sqrtg
        real(8), dimension(size(x)), intent(out) :: bder
        real(8), dimension(size(x)), intent(out) :: hcovar
        real(8), dimension(size(x)), intent(out) :: hctrvr
        real(8), dimension(size(x)), intent(out) :: hcurl

        real(8) :: spl_val(3)
        real(8) :: sqgbmod, sqgbmod2  ! sqg*B, sqg*B^2

        real(8) :: x1

        ! safety measure in order not to extrapolate
        ! note: this is s
        x1 = max(params0(1, 1), x(1))
        x1 = min(params0(nflux, 1), x1)

        spl_val = spline_val_0(spl_coeff1(:, :, 3), x1)
        Bthcov = ItoB*spl_val(1)*bfac
        dBthcovds = ItoB*spl_val(2)*bfac
        spl_val = spline_val_0(spl_coeff1(:, :, 2), x1)
        Bphcov = ItoB*spl_val(1)*bfac
        dBphcovds = ItoB*spl_val(2)*bfac
        spl_val = spline_val_0(spl_coeff1(:, :, 1), x1)
        iota = spl_val(1)
        q = 1/iota
        dqds = -spl_val(2)/iota**2

        call fast_sin_cos(modes0(1, :, 1), x(3), sinterm, costerm)

        ! calculate B-field from modes
        if (inp_swi == 8) then
            call cached_spline(x1, s_prev, spl_coeff2(:, :, 4, :), spl_val_c)
            B0mnc(:) = 1d4*spl_val_c(1, :)*bfac
            dB0dsmnc(:) = 1d4*spl_val_c(2, :)*bfac
            B0h = B0mnc(1)

            bmod = sum(B0mnc*costerm)
            bder(1) = sum(dB0dsmnc*costerm)/bmod
            bder(2) = 0d0
            bder(3) = sum(-modes0(1, :, 1)*B0mnc*sinterm)/bmod
        else if (inp_swi == 9) then
            call cached_spline(x1, s_prev, spl_coeff2(:, :, 7, :), spl_val_c)
            B0mnc(:) = 1d4*spl_val_c(1, :)*bfac
            dB0dsmnc(:) = 1d4*spl_val_c(2, :)*bfac
            call cached_spline(x1, s_prev, spl_coeff2(:, :, 8, :), spl_val_s)
            B0mns(:) = 1d4*spl_val_s(1, :)*bfac
            dB0dsmns(:) = 1d4*spl_val_s(2, :)*bfac
            B0h = B0mnc(1)

            bmod = sum(B0mnc*costerm + B0mns*sinterm)
            bder(1) = sum(dB0dsmnc*costerm + dB0dsmns*sinterm)/bmod
            bder(2) = 0d0
            bder(3) = sum(-modes0(1, :, 1)*B0mnc*sinterm &
                          + modes0(1, :, 1)*B0mns*costerm)/bmod
        end if

        sqgbmod2 = sign_theta*psi_pr*(Bphcov + iota*Bthcov)
        sqgbmod = sqgbmod2/bmod
        sqrtg = sqgbmod/bmod

        hcovar(1) = 0d0  ! TODO
        hcovar(2) = Bphcov/bmod
        hcovar(3) = Bthcov/bmod

        hctrvr(1) = 0d0
        hctrvr(2) = sign_theta*psi_pr/sqgbmod
        hctrvr(3) = sign_theta*iota*psi_pr/sqgbmod

        hcurl(1) = 0d0  ! TODO
        hcurl(3) = 0d0  ! TODO
        hcurl(2) = 0d0  ! TODO

        s_prev = x1

    end subroutine do_magfie

    subroutine boozer_read(filename)
        ! Reads Boozer in_file and converts SI to CGS

        integer :: ksurf, kmode
        real(8) :: flux
        character(len=*) :: filename
        open (unit=18, file=filename, action='read', status='old')
        read (18, '(////)')
        read (18, *) m0b, n0b, nflux, nfp, flux, a, R0
        a = 100*a   ! m -> cm
        R0 = 100*R0 ! m -> cm

        psi_pr = 1.0d8*flux/(2*pi)*bfac ! T -> Gauss, m -> cm

        nmode = (m0b + 1)*(n0b + 1)

        ! Allocate params and modes (deallocate first if size changed)
        if (allocated(params0)) then
            if (size(params0, 1) /= nflux) deallocate(params0, modes0)
        end if
        if (.not. allocated(params0)) then
            allocate (params0(nflux, ncol1 + 1))
            allocate (modes0(nflux, nmode, ncol2 + 2))
        end if
        do ksurf = 1, nflux
            read (18, '(/)')
            read (18, *) params0(ksurf, :)
            read (18, *)
            do kmode = 1, nmode
                read (18, *) modes0(ksurf, kmode, :)
            end do
        end do
        close (unit=18)
        ! Set R0 to first harmonic
        R0 = modes0(1, 1, 3)*100
    end subroutine boozer_read

    subroutine booz_to_cyl(x, r)

        real(8), intent(in) :: x(3)  ! Boozer coordinates (s, ph, th)
        real(8), intent(out) :: r(3)  ! Cylindrical coordinates (R, phi, Z)

        real(8) :: spl_val(3), x1

        integer :: j

        if (inp_swi /= 9) error stop  ! Only implemented for ASDEX-U type of data

        x1 = max(params0(1, 1), x(1))
        x1 = min(params0(nflux, 1), x1)

        do j = 1, nmode
            spl_val = spline_val_0(spl_coeff2(:, :, 1, j), x1)
            rmnc(j) = 1.0d2*spl_val(1)
            spl_val = spline_val_0(spl_coeff2(:, :, 2, j), x1)
            rmns(j) = 1.0d2*spl_val(1)
            spl_val = spline_val_0(spl_coeff2(:, :, 3, j), x1)
            zmnc(j) = 1.0d2*spl_val(1)
            spl_val = spline_val_0(spl_coeff2(:, :, 4, j), x1)
            zmns(j) = 1.0d2*spl_val(1)
        end do

        r(1) = sum(rmnc*cos(modes0(1, :, 1)*x(3)) + rmns*sin(modes0(1, :, 1)*x(3)))
        r(2) = 0.0d0  ! TODO: phi
        r(3) = sum(zmnc*cos(modes0(1, :, 1)*x(3)) + zmns*sin(modes0(1, :, 1)*x(3)))

    end subroutine booz_to_cyl

    subroutine fast_sin_cos(m, x, sinterm_, costerm_)
        ! Fast sine and cosine that assumes equally spaced ascending mode numbers
        real(8), intent(in) :: m(:), x
        real(8), intent(out) :: sinterm_(:), costerm_(:)

        real(8) :: dm
        complex(8) :: fourier_factor, rotation
        integer :: j

        dm = m(2) - m(1)
        fourier_factor  = exp(imun*m(1)*x)
        rotation = exp(imun*dm*x)

        costerm_ = (0.0d0, 0.0d0)
        sinterm_ = (0.0d0, 0.0d0)
        do j = 1, size(m)
            costerm_(j) = real(fourier_factor)
            sinterm_(j) = imag(fourier_factor)
            fourier_factor = fourier_factor*rotation
        end do
    end subroutine fast_sin_cos

end module do_magfie_mod

module do_magfie_pert_mod

    use util
    use spline
    use do_magfie_mod, only: s, bfac, inp_swi

    implicit none

    real(8), private :: s_prev = -1.0d0
    real(8), private, allocatable :: spl_val_c(:,:), spl_val_s(:,:)

    real(8), allocatable, protected :: params(:, :), modes(:, :, :)
    integer, protected :: mb, nb, nflux, nfp, nmode

    real(8), allocatable, protected :: spl_coeff1(:, :, :), spl_coeff2(:, :, :, :)

    ! Work arrays (size=nmode)
    real(8), private, allocatable :: Bmnc(:), Bmns(:)

    integer :: ncol1, ncol2 ! number of columns in input file
    real(8) :: mph ! toroidal perturbation mode
    contains

    subroutine do_magfie_pert_init(path)
        character(len=*), intent(in) :: path

        integer :: j, k
        real(8) :: x(3)
        complex(8) :: dummy

        ncol1 = 5
        if (inp_swi == 8) ncol2 = 4 ! tok_circ
        if (inp_swi == 9) ncol2 = 8 ! ASDEX
        call boozer_read_pert(path)

        mph = nfp*modes(1, 1, 2)

        ! Allocate spline coefficient arrays (deallocate first if size changed)
        if (allocated(spl_coeff1)) then
            if (size(spl_coeff1, 1) /= nflux - 1) deallocate(spl_coeff1, spl_coeff2)
        end if
        if (.not. allocated(spl_coeff1)) then
            allocate (spl_coeff1(nflux - 1, 5, ncol1))
            allocate (spl_coeff2(nflux - 1, 5, ncol2, nmode))
        end if

        ! Allocate work arrays (deallocate first if size changed)
        if (allocated(Bmnc)) then
            if (size(Bmnc) /= nmode) then
                deallocate(Bmnc)
                if (allocated(Bmns)) deallocate(Bmns)
            end if
        end if
        if (.not. allocated(Bmnc)) then
            allocate(Bmnc(nmode))
            if (ncol2 >= 8) allocate(Bmns(nmode))
        end if

        ! calculate spline coefficients
        do k = 1, ncol1
            ! first column is s, so start with second column
            spl_coeff1(:, :, k) = spline_coeff(params(:, 1), params(:, k + 1))
        end do

        do j = 1, nmode
            do k = 1, ncol2
                ! first two columns are mode numbers, so start with third column
                spl_coeff2(:, :, k, j) = spline_coeff(params(:, 1), modes(:, j, k + 2))
            end do
        end do

        ! Allocate cached spline values (deallocate first if size changed)
        if (allocated(spl_val_c)) then
            if (size(spl_val_c, 2) /= nmode) deallocate(spl_val_c, spl_val_s)
        end if
        if (.not. allocated(spl_val_c)) then
            allocate (spl_val_c(3, nmode))
            allocate (spl_val_s(3, nmode))
        end if

        x(1) = s
        x(2) = 0.0
        x(3) = 0.0
        call do_magfie_pert_amp(x, dummy)
    end subroutine do_magfie_pert_init

    subroutine do_magfie_pert_amp(x, bamp)
        real(8), dimension(:), intent(in) :: x
        complex(8), intent(out) :: bamp

        real(8) :: x1

        ! safety measure in order not to extrapolate
        x1 = max(params(1, 1), x(1))
        x1 = min(params(nflux, 1), x1)

        ! calculate B-field from modes
        if (inp_swi == 8) then
            call cached_spline(x1, s_prev, spl_coeff2(:, :, 4, :), spl_val_c)
            Bmnc(:) = 1d4*spl_val_c(1, :)*bfac
            bamp = sum(Bmnc*cos(modes(1, :, 1)*x(3)))
        else if (inp_swi == 9) then
            call cached_spline(x1, s_prev, spl_coeff2(:, :, 7, :), spl_val_c)
            call cached_spline(x1, s_prev, spl_coeff2(:, :, 8, :), spl_val_s)
            Bmnc(:) = 1d4*spl_val_c(1, :)*bfac
            Bmns(:) = 1d4*spl_val_s(1, :)*bfac
            bamp = fast_fourier_sum(Bmnc, Bmns, modes(1, :, 1), x(3))
        end if

        s_prev = x1
    end subroutine do_magfie_pert_amp

    subroutine do_magfie_pert(x, bmod)
        real(8), dimension(:), intent(in) :: x
        real(8), intent(out) :: bmod
        complex(8) :: bamp

        call do_magfie_pert_amp(x, bamp)
        bmod = real(sum(bamp*exp(imun*nfp*modes(1, :, 2)*x(2))))

    end subroutine do_magfie_pert

    subroutine boozer_read_pert(filename)
        integer :: ksurf, kmode
        real(8) :: flux, dummy
        character(len=*) :: filename
        open (unit=18, file=filename)
        read (18, '(////)')
        read (18, *) mb, nb, nflux, nfp, flux, dummy, dummy
        nmode = (mb + 1)*(nb + 1)

        ! Allocate params and modes (deallocate first if size changed)
        if (allocated(params)) then
            if (size(params, 1) /= nflux) deallocate(params, modes)
        end if
        if (.not. allocated(params)) then
            allocate (params(nflux, ncol1 + 1))
            allocate (modes(nflux, nmode, ncol2 + 2))
        end if
        do ksurf = 1, nflux
            read (18, '(/)')
            read (18, *) params(ksurf, :)
            read (18, *)
            do kmode = 1, nmode
                read (18, *) modes(ksurf, kmode, :)
                call check_equal_space_ascending(modes(ksurf, :, 1), kmode)
                call check_equal(modes(ksurf, :, 2), kmode)
            end do
        end do
        close (unit=18)
    end subroutine boozer_read_pert

    subroutine check_equal_space_ascending(m, kmode)
        real(8), dimension(:), intent(in) :: m
        integer, intent(in) :: kmode

        real(8), parameter :: tol = 1.0d-10
        real(8) :: diff, diff2

        if (kmode < 3) return

        diff = m(kmode - 1) - m(kmode - 2)
        diff2 = m(kmode) - m(kmode - 1)
        if (abs(diff - diff2) > tol) then
            error stop "harmonics not equally spaced and ascending"
        end if
    end subroutine check_equal_space_ascending

    subroutine check_equal(m, kmode)
        real(8), dimension(:), intent(in) :: m
        integer, intent(in) :: kmode

        real(8), parameter :: tol = 1.0d-10

        if (kmode < 2) return

        if (abs(m(kmode) - m(kmode - 1)) > tol) then
            error stop "harmonics not equal"
        end if
    end subroutine check_equal

    function fast_fourier_sum(fmnc, fmns, m, x)
        ! Fast Fourier sum that assumes equally spaced ascending mode numbers
        real(8), dimension(:), intent(in) :: fmnc, fmns, m
        real(8), intent(in) :: x
        complex(8) :: fast_fourier_sum

        real(8) :: dm
        complex(8) :: fourier_factor, rotation
        integer :: j

        dm = m(2) - m(1)
        fourier_factor  = exp(imun*m(1)*x)
        rotation = exp(imun*dm*x)

        fast_fourier_sum = (0.0d0, 0.0d0)
        do j = 1, size(m)
            fast_fourier_sum = fast_fourier_sum + &
                (fmnc(j) - imun*fmns(j))*fourier_factor
            fourier_factor = fourier_factor*rotation
        end do
    end function fast_fourier_sum
end module do_magfie_pert_mod
