module do_magfie_mod
    use iso_fortran_env, only: dp => real64
    use util
    use spline

    implicit none

    real(dp), private :: s_prev = -1.0_dp
    real(dp), private, allocatable :: spl_val_c(:, :), spl_val_s(:, :)
    ! Work arrays previously automatic on stack (size=nmode)
    real(dp), private, allocatable :: B0mnc(:), dB0dsmnc(:), B0mns(:), dB0dsmns(:)
    real(dp), private, allocatable :: costerm(:), sinterm(:)

    real(dp), parameter :: sign_theta = -1.0_dp ! negative for left-handed

    real(dp) :: s = 0.0_dp, psi_pr = 0.0_dp, Bthcov = 0.0_dp, Bphcov = 0.0_dp, &
        dBthcovds = 0.0_dp, dBphcovds = 0.0_dp, &
        q = 0.0_dp, dqds = 0.0_dp, iota = 0.0_dp, R0 = 0.0_dp, a = 0.0_dp, &
        eps = 0.0_dp, B0h = 0.0_dp, B00 = 0.0_dp
    real(dp) :: bfac = 1.0_dp
    ! B0h is the 0th theta harmonic of bmod on current flux surface
    ! and B00 the 0th theta harmonic of bmod on the innermost flux surface

    real(dp), allocatable, protected :: params0(:, :), modes0(:, :, :)
    integer, protected :: m0b, n0b, nflux, nfp, nmode

    real(dp), allocatable, protected :: spl_coeff1(:, :, :), spl_coeff2(:, :, :, :)
    ! Work arrays for booz_to_cyl (size=nmode)
    real(dp), private, allocatable :: rmnc(:), rmns(:), zmnc(:), zmns(:)

    real(dp), parameter :: ItoB = 2.0e-1_dp * sign_theta ! Covarient B (cgs) from I (SI)
    ! Bcov=mu0/2pi*I,mu0->4pi/c,I->10^(-1)*c*I

    integer :: ncol1 = 0, ncol2 = 0 ! number of columns in input file

    integer :: inp_swi = 0 ! type of input file

    ! Initialization flag for threadprivate allocatable arrays
    logical, save :: magfie_arrays_initialized = .false.

    ! Magic sentinel for auto-initializing threadprivate state in worker threads
    ! DATA/initializers don't run in worker threads, so we use a magic value
    integer, parameter :: MAGFIE_INIT_SENTINEL = 314159265
    integer, save :: magfie_thread_init_state = 0

    ! Working buffers (per-thread cache and temporary arrays)
    !$omp threadprivate (magfie_arrays_initialized, magfie_thread_init_state)
    !$omp threadprivate (s_prev, spl_val_c, spl_val_s)
    !$omp threadprivate (B0mnc, dB0dsmnc, B0mns, dB0dsmns, costerm, sinterm)
    !$omp threadprivate (rmnc, rmns, zmnc, zmns)

    ! Flux-surface dependent quantities (computed per-thread for each s)
    !$omp threadprivate (s, Bthcov, Bphcov, dBthcovds, dBphcovds)
    !$omp threadprivate (q, dqds, iota, eps, B0h)

    ! Shared data (NOT threadprivate): bfac, params0, modes0, m0b, n0b, nflux,
    ! nfp, nmode, spl_coeff1, spl_coeff2, ncol1, ncol2, inp_swi, a, B00, psi_pr, R0

    ! inp_swi == 10: Boozer chartmap (NetCDF) input via libneo reader.
    ! Shared read-only arrays; populated once by read_boozer_chartmap_file.
    integer :: cm_n_rho = 0, cm_n_theta = 0
    real(dp) :: cm_torflux = 0.0_dp, cm_h_theta = 0.0_dp
    real(dp) :: cm_rho_min = 0.0_dp, cm_rho_max = 0.0_dp
    real(dp), allocatable, protected :: cm_rho(:)
    real(dp), allocatable, protected :: cm_s_aphi(:) ! s abscissa for A_phi
    real(dp), allocatable, protected :: cm_spl_aphi(:, :) ! spline of A_phi(s)
    real(dp), allocatable, protected :: cm_spl_Bth(:, :) ! spline of B_theta(rho)
    real(dp), allocatable, protected :: cm_spl_Bph(:, :) ! spline of B_phi(rho)
    ! (n_rho-1, 5, n_theta): for each theta grid index, spline of Bmod vs rho
    real(dp), allocatable, protected :: cm_spl_bmod(:, :, :)

contains

    subroutine set_s(s_)
        real(dp), intent(in) :: s_

        s = s_
    end subroutine set_s

    subroutine magfie_thread_init()
        ! Initialize threadprivate variables for this thread
        ! Must be called once per thread before using magfie routines
        magfie_arrays_initialized = .false.
        s_prev = -1.0_dp
    end subroutine magfie_thread_init

    subroutine read_boozer_file(path)
        ! Main thread only: Read Boozer file and prepare shared spline data
        ! This should be called ONCE before parallel region
        character(len=*), intent(in) :: path
        integer :: j, k

        if (inp_swi == 10) then
            call read_boozer_chartmap_file(path)
            return
        end if

        ncol1 = 5
        if (inp_swi == 8) ncol2 = 4 ! tok_circ
        if (inp_swi == 9) ncol2 = 8 ! ASDEX

        ! allocates params0, modes0 - SHARED arrays
        call boozer_read(path)

        ! Allocate SHARED spline coefficient arrays
        if (allocated(spl_coeff1)) then
            if (size(spl_coeff1, 1) /= nflux - 1) deallocate(spl_coeff1, spl_coeff2)
        end if
        if (.not. allocated(spl_coeff1)) then
            allocate (spl_coeff1(nflux - 1, 5, ncol1))
            allocate (spl_coeff2(nflux - 1, 5, ncol2, nmode))
        end if

        do k = 1, ncol1
            spl_coeff1(:, :, k) = spline_coeff(params0(:, 1), params0(:, k + 1))
        end do

        do j = 1, nmode
            do k = 1, ncol2
                spl_coeff2(:, :, k, j) = spline_coeff(params0(:, 1), modes0(:, j, k + 2))
            end do
        end do

        ! Set B00 from first mode
        B00 = 1.0e4_dp * modes0(1, 1, 6) * bfac
    end subroutine read_boozer_file

    subroutine init_magfie_at_s()
        ! Per-thread: Initialize magnetic field at current s value
        ! Assumes s is already set via set_s()
        ! Allocates threadprivate working buffers and computes s-dependent values
        real(dp) :: x(3), bmod, sqrtg
        real(dp), dimension(3) :: bder, hcovar, hctrvr, hcurl

        ! Auto-initialize threadprivate state if not yet done for this thread
        ! (DATA/initializers don't run in worker threads)
        if (magfie_thread_init_state /= MAGFIE_INIT_SENTINEL) then
            call magfie_thread_init()
            magfie_thread_init_state = MAGFIE_INIT_SENTINEL
        end if

        ! Allocate threadprivate working buffers (safe for undefined allocation status)
        if (.not. magfie_arrays_initialized) then
            if (inp_swi /= 10) then
                if (allocated(B0mnc)) deallocate(B0mnc)
                if (allocated(dB0dsmnc)) deallocate(dB0dsmnc)
                allocate(B0mnc(nmode), dB0dsmnc(nmode))
                if (ncol2 >= 8) then
                    if (allocated(B0mns)) deallocate(B0mns)
                    if (allocated(dB0dsmns)) deallocate(dB0dsmns)
                    allocate(B0mns(nmode), dB0dsmns(nmode))
                end if
                if (allocated(costerm)) deallocate(costerm)
                if (allocated(sinterm)) deallocate(sinterm)
                allocate(costerm(nmode), sinterm(nmode))

                if (allocated(rmnc)) deallocate(rmnc)
                if (allocated(rmns)) deallocate(rmns)
                if (allocated(zmnc)) deallocate(zmnc)
                if (allocated(zmns)) deallocate(zmns)
                allocate(rmnc(nmode), rmns(nmode), zmnc(nmode), zmns(nmode))

                if (allocated(spl_val_c)) deallocate(spl_val_c)
                if (allocated(spl_val_s)) deallocate(spl_val_s)
                allocate(spl_val_c(3, nmode), spl_val_s(3, nmode))
            end if

            magfie_arrays_initialized = .true.
        end if

        ! Initialize cache
        s_prev = -1.0_dp

        x(1) = s
        x(2) = 0.0
        x(3) = 0.0
        call do_magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
    end subroutine init_magfie_at_s

    subroutine do_magfie_init(path)
        ! init axisymmetric part of field from infile
        ! Backward compatibility wrapper: calls read_boozer_file + init_magfie_at_s
        character(len=*), intent(in) :: path

        if (inp_swi == 10) then
            call read_boozer_chartmap_file(path)
        else
            call read_boozer_file(path)
        end if
        call init_magfie_at_s()
    end subroutine do_magfie_init

    subroutine do_magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
        ! Evaluate unperturbed axisymmetric magnetic field in point x = (s, ph, th)
        real(dp), dimension(:), intent(in) :: x
        real(dp), intent(out) :: bmod
        real(dp), intent(out) :: sqrtg
        real(dp), dimension(size(x)), intent(out) :: bder
        real(dp), dimension(size(x)), intent(out) :: hcovar
        real(dp), dimension(size(x)), intent(out) :: hctrvr
        real(dp), dimension(size(x)), intent(out) :: hcurl

        real(dp) :: spl_val(3)
        real(dp) :: sqgbmod, sqgbmod2 ! sqg*B, sqg*B^2
        real(dp) :: x1

        if (inp_swi == 10) then
            call do_magfie_chartmap(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
            return
        end if

        ! safety measure in order not to extrapolate
        ! note: this is s
        x1 = max(params0(1, 1), x(1))
        x1 = min(params0(nflux, 1), x1)

        spl_val = spline_val_0(spl_coeff1(:, :, 3), x1)
        Bthcov = ItoB * spl_val(1) * bfac
        dBthcovds = ItoB * spl_val(2) * bfac
        spl_val = spline_val_0(spl_coeff1(:, :, 2), x1)
        Bphcov = ItoB * spl_val(1) * bfac
        dBphcovds = ItoB * spl_val(2) * bfac
        spl_val = spline_val_0(spl_coeff1(:, :, 1), x1)
        iota = spl_val(1)
        q = 1 / iota
        dqds = -spl_val(2) / iota**2

        call fast_sin_cos(modes0(1, :, 1), x(3), sinterm, costerm)

        ! calculate B-field from modes
        if (inp_swi == 8) then
            call cached_spline(x1, s_prev, spl_coeff2(:, :, 4, :), spl_val_c)
            B0mnc(:) = 1.0e4_dp * spl_val_c(1, :) * bfac
            dB0dsmnc(:) = 1.0e4_dp * spl_val_c(2, :) * bfac
            B0h = B0mnc(1)

            bmod = sum(B0mnc * costerm)
            bder(1) = sum(dB0dsmnc * costerm) / bmod
            bder(2) = 0.0_dp
            bder(3) = sum(-modes0(1, :, 1) * B0mnc * sinterm) / bmod
        else if (inp_swi == 9) then
            call cached_spline(x1, s_prev, spl_coeff2(:, :, 7, :), spl_val_c)
            B0mnc(:) = 1.0e4_dp * spl_val_c(1, :) * bfac
            dB0dsmnc(:) = 1.0e4_dp * spl_val_c(2, :) * bfac
            call cached_spline(x1, s_prev, spl_coeff2(:, :, 8, :), spl_val_s)
            B0mns(:) = 1.0e4_dp * spl_val_s(1, :) * bfac
            dB0dsmns(:) = 1.0e4_dp * spl_val_s(2, :) * bfac
            B0h = B0mnc(1)

            bmod = sum(B0mnc * costerm + B0mns * sinterm)
            bder(1) = sum(dB0dsmnc * costerm + dB0dsmns * sinterm) / bmod
            bder(2) = 0.0_dp
            bder(3) = sum(-modes0(1, :, 1) * B0mnc * sinterm &
                + modes0(1, :, 1) * B0mns * costerm) / bmod
        end if

        sqgbmod2 = sign_theta * psi_pr * (Bphcov + iota * Bthcov)
        sqgbmod = sqgbmod2 / bmod
        sqrtg = sqgbmod / bmod

        hcovar(1) = 0.0_dp ! TODO
        hcovar(2) = Bphcov / bmod
        hcovar(3) = Bthcov / bmod

        hctrvr(1) = 0.0_dp
        hctrvr(2) = sign_theta * psi_pr / sqgbmod
        hctrvr(3) = sign_theta * iota * psi_pr / sqgbmod

        hcurl(1) = 0.0_dp ! TODO
        hcurl(3) = 0.0_dp ! TODO
        hcurl(2) = 0.0_dp ! TODO

        s_prev = x1

    end subroutine do_magfie

    subroutine do_magfie_chartmap(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
        ! Evaluate field from Boozer chartmap data (inp_swi == 10).
        ! x = (s, phi, theta); Boozer convention theta is the poloidal angle.
        real(dp), dimension(:), intent(in) :: x
        real(dp), intent(out) :: bmod
        real(dp), intent(out) :: sqrtg
        real(dp), dimension(size(x)), intent(out) :: bder
        real(dp), dimension(size(x)), intent(out) :: hcovar
        real(dp), dimension(size(x)), intent(out) :: hctrvr
        real(dp), dimension(size(x)), intent(out) :: hcurl

        real(dp) :: spl_val(3), rho1, s1, dAphi_ds
        real(dp) :: bmod0, bmod1, dBmod_drho, dBmod_dth
        real(dp) :: sqgbmod, sqgbmod2
        real(dp) :: alpha, th_shift
        integer :: k_th, k_th1

        s1 = max(cm_s_aphi(1), x(1))
        s1 = min(cm_s_aphi(size(cm_s_aphi)), s1)
        rho1 = max(cm_rho_min, sqrt(s1))
        rho1 = min(cm_rho_max, rho1)

        ! Radial profiles from chartmap splines.
        spl_val = spline_val_0(cm_spl_aphi, s1)
        dAphi_ds = spl_val(2)
        ! iota = -dA_phi/ds / torflux  (A_phi = -torflux * integral iota ds)
        iota = -dAphi_ds / cm_torflux
        q = 1.0_dp / iota
        dqds = spl_val(3) / (cm_torflux * iota**2)

        ! Chartmap stores B_theta/B_phi in G*cm with the sign from buco/bvco (no
        ! sign_theta flip).  The .bc path applies ItoB = sign_theta * mu0/(2*pi)*1e6
        ! which introduces an extra sign_theta factor.  Apply it here so that Bthcov
        ! and Bphcov carry the same sign convention as the .bc path.
        spl_val = spline_val_0(cm_spl_Bth, rho1)
        Bthcov = sign_theta * spl_val(1) * bfac
        ! chain rule: dBthcov/ds = dBthcov/drho * drho/ds, drho/ds = 1/(2*rho)
        dBthcovds = sign_theta * spl_val(2) * bfac / (2.0_dp * rho1)

        spl_val = spline_val_0(cm_spl_Bph, rho1)
        Bphcov = sign_theta * spl_val(1) * bfac
        dBphcovds = sign_theta * spl_val(2) * bfac / (2.0_dp * rho1)

        ! psi_pr matches the .bc path: .bc stores flux < 0 for this orientation so
        ! psi_pr = flux*1e8/(2*pi) < 0, while torflux in the chartmap = |psi_pr|.
        ! Set psi_pr = -cm_torflux so sign_theta*psi_pr has the same sign in both paths.
        psi_pr = -cm_torflux * bfac

        ! Bmod on (rho, theta) grid: wrap theta to [0, 2*pi) then interpolate.
        th_shift = x(3) - floor(x(3) / (2.0_dp * pi)) * 2.0_dp * pi
        k_th = max(0, min(cm_n_theta - 1, int(th_shift / cm_h_theta)))
        k_th1 = mod(k_th + 1, cm_n_theta)
        alpha = (th_shift - real(k_th, dp) * cm_h_theta) / cm_h_theta

        ! Spline in rho at k_th grid point (1-based index k_th+1).
        spl_val = spline_val_0(cm_spl_bmod(:, :, k_th + 1), rho1)
        bmod0 = spl_val(1) * bfac
        dBmod_drho = spl_val(2) * bfac

        ! Spline in rho at next theta grid point (periodic wrap).
        spl_val = spline_val_0(cm_spl_bmod(:, :, k_th1 + 1), rho1)
        bmod1 = spl_val(1) * bfac

        bmod = (1.0_dp - alpha) * bmod0 + alpha * bmod1
        B0h = bmod

        ! dBmod/ds via chain rule through rho.
        spl_val = spline_val_0(cm_spl_bmod(:, :, k_th + 1), rho1)
        dBmod_drho = (1.0_dp - alpha) * spl_val(2) * bfac
        spl_val = spline_val_0(cm_spl_bmod(:, :, k_th1 + 1), rho1)
        dBmod_drho = dBmod_drho + alpha * spl_val(2) * bfac
        bder(1) = dBmod_drho / (2.0_dp * rho1 * bmod)
        ! x(2) is the toroidal angle phi; Bmod has no phi dependence in axisymmetric
        ! (nfp=1, phi=0 slice) Boozer coordinates, so bder(2)=0 is exact.
        bder(2) = 0.0_dp
        ! dBmod/dtheta by finite difference across the theta grid cell.
        dBmod_dth = (bmod1 - bmod0) / (cm_h_theta * bmod)
        bder(3) = dBmod_dth

        sqgbmod2 = sign_theta * psi_pr * (Bphcov + iota * Bthcov)
        sqgbmod = sqgbmod2 / bmod
        sqrtg = sqgbmod / bmod

        hcovar(1) = 0.0_dp
        hcovar(2) = Bphcov / bmod
        hcovar(3) = Bthcov / bmod

        hctrvr(1) = 0.0_dp
        hctrvr(2) = sign_theta * psi_pr / sqgbmod
        hctrvr(3) = sign_theta * iota * psi_pr / sqgbmod

        hcurl(1) = 0.0_dp
        hcurl(3) = 0.0_dp
        hcurl(2) = 0.0_dp

    end subroutine do_magfie_chartmap

    subroutine read_boozer_chartmap_file(path)
        ! Read a Boozer chartmap NetCDF file via libneo's boozer_chartmap_io reader
        ! and build the shared spline arrays for do_magfie_chartmap.
        use boozer_chartmap_io, only: read_boozer_chartmap
        use boozer_chartmap_types, only: boozer_chartmap_data_t

        character(len=*), intent(in) :: path

        type(boozer_chartmap_data_t) :: d
        integer :: j, n_s_aphi

        call read_boozer_chartmap(path, d)

        cm_n_rho = d%n_rho
        cm_n_theta = d%n_theta
        cm_torflux = d%torflux
        cm_h_theta = d%h_theta
        cm_rho_min = d%rho(1)
        cm_rho_max = d%rho(d%n_rho)
        n_s_aphi = d%n_s
        R0 = d%rmajor * 100.0_dp ! m -> cm
        a = d%aminor * 100.0_dp ! m -> cm
        if (a <= 0.0_dp) then
            error stop "Boozer chartmap has no usable minor radius"
        end if
        ! psi_pr stored globally; negate torflux to match the .bc sign convention
        ! (chartmap torflux > 0, .bc psi_pr = flux*1e8/(2*pi) < 0 for typical orientation).
        psi_pr = -d%torflux * bfac

        if (allocated(cm_rho)) deallocate(cm_rho)
        if (allocated(cm_s_aphi)) deallocate(cm_s_aphi)
        if (allocated(cm_spl_aphi)) deallocate(cm_spl_aphi)
        if (allocated(cm_spl_Bth)) deallocate(cm_spl_Bth)
        if (allocated(cm_spl_Bph)) deallocate(cm_spl_Bph)
        if (allocated(cm_spl_bmod)) deallocate(cm_spl_bmod)

        allocate(cm_rho(cm_n_rho))
        cm_rho = d%rho

        allocate(cm_s_aphi(n_s_aphi))
        cm_s_aphi = d%s

        allocate(cm_spl_aphi(n_s_aphi - 1, 5))
        cm_spl_aphi = spline_coeff(d%s, d%A_phi)

        allocate(cm_spl_Bth(cm_n_rho - 1, 5))
        cm_spl_Bth = spline_coeff(d%rho, d%B_theta)

        allocate(cm_spl_Bph(cm_n_rho - 1, 5))
        cm_spl_Bph = spline_coeff(d%rho, d%B_phi)

        ! Bmod grid: (n_rho, n_theta, n_phi). For the do_magfie evaluation we
        ! use the phi=0 slice (index 1), which equals the phi=2pi/nfp plane
        ! for a field with nfp-fold symmetry; works for axisymmetric too.
        allocate(cm_spl_bmod(cm_n_rho - 1, 5, cm_n_theta))
        do j = 1, cm_n_theta
            cm_spl_bmod(:, :, j) = spline_coeff(d%rho, d%Bmod(:, j, 1))
        end do

        ! B00: (0,0) mode of Bmod = rho=0 surface average; approximate with rho=rho_min.
        B00 = d%Bmod(1, 1, 1) * bfac

    end subroutine read_boozer_chartmap_file

    subroutine boozer_read(filename)
        ! Reads Boozer in_file and converts SI to CGS

        integer :: ksurf, kmode
        real(dp) :: flux
        character(len=*) :: filename
        open (unit=18, file=filename, action='read', status='old')
        read (18, '(////)')
        read (18, *) m0b, n0b, nflux, nfp, flux, a, R0
        a = 100 * a ! m -> cm
        R0 = 100 * R0 ! m -> cm

        psi_pr = 1.0e8_dp * flux / (2 * pi) * bfac ! T -> Gauss, m -> cm

        nmode = (m0b + 1) * (n0b + 1)

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

        real(dp), intent(in) :: x(3) ! Boozer coordinates (s, ph, th)
        real(dp), intent(out) :: r(3) ! Cylindrical coordinates (R, phi, Z)

        real(dp) :: spl_val(3), x1

        integer :: j

        if (inp_swi /= 9) error stop ! Only implemented for ASDEX-U type of data

        x1 = max(params0(1, 1), x(1))
        x1 = min(params0(nflux, 1), x1)

        do j = 1, nmode
            spl_val = spline_val_0(spl_coeff2(:, :, 1, j), x1)
            rmnc(j) = 1.0e2_dp * spl_val(1)
            spl_val = spline_val_0(spl_coeff2(:, :, 2, j), x1)
            rmns(j) = 1.0e2_dp * spl_val(1)
            spl_val = spline_val_0(spl_coeff2(:, :, 3, j), x1)
            zmnc(j) = 1.0e2_dp * spl_val(1)
            spl_val = spline_val_0(spl_coeff2(:, :, 4, j), x1)
            zmns(j) = 1.0e2_dp * spl_val(1)
        end do

        r(1) = sum(rmnc * cos(modes0(1, :, 1) * x(3)) + rmns * sin(modes0(1, :, 1) * x(3)))
        r(2) = 0.0_dp ! TODO: phi
        r(3) = sum(zmnc * cos(modes0(1, :, 1) * x(3)) + zmns * sin(modes0(1, :, 1) * x(3)))

    end subroutine booz_to_cyl

    subroutine fast_sin_cos(m, x, sinterm_, costerm_)
        ! Fast sine and cosine that assumes equally spaced ascending mode numbers
        real(dp), intent(in) :: m(:), x
        real(dp), intent(out) :: sinterm_(:), costerm_(:)

        real(dp) :: dm
        complex(dp) :: fourier_factor, rotation
        integer :: j

        dm = m(2) - m(1)
        fourier_factor  = exp(imun*m(1)*x)
        rotation = exp(imun*dm*x)

        costerm_ = (0.0_dp, 0.0_dp)
        sinterm_ = (0.0_dp, 0.0_dp)
        do j = 1, size(m)
            costerm_(j) = real(fourier_factor)
            sinterm_(j) = imag(fourier_factor)
            fourier_factor = fourier_factor*rotation
        end do
    end subroutine fast_sin_cos

end module do_magfie_mod

module do_magfie_pert_mod
    use iso_fortran_env, only: dp => real64
    use util
    use spline
    use do_magfie_mod, only: s, bfac

    implicit none

    real(dp), private :: s_prev = -1.0_dp
    real(dp), private, allocatable :: spl_val_c(:, :), spl_val_s(:, :)

    real(dp), allocatable, protected :: params(:, :), modes(:, :, :)
    integer, protected :: mb, nb, nflux, nfp, nmode

    real(dp), allocatable, protected :: spl_coeff1(:, :, :), spl_coeff2(:, :, :, :)

    ! Work arrays (size=nmode)
    real(dp), private, allocatable :: Bmnc(:), Bmns(:)

    integer :: ncol1, ncol2 ! number of columns in input file
    integer :: inp_swi_pert = -1 ! perturbation .bc format, independent of axisymmetric input
    integer :: mph ! toroidal perturbation mode (threadprivate)
    integer :: mph_shared = 0 ! shared copy for namelist input (when pertfile=.false.)

    ! Initialization flag for threadprivate allocatable arrays
    logical, save :: magfie_pert_arrays_initialized = .false.

    ! Magic sentinel for auto-initializing threadprivate state in worker threads
    integer, parameter :: MAGFIE_PERT_INIT_SENTINEL = 271828182
    integer, save :: magfie_pert_thread_init_state = 0

    ! Working buffers (per-thread cache and temporary arrays)
    !$omp threadprivate (magfie_pert_arrays_initialized, magfie_pert_thread_init_state)
    !$omp threadprivate (s_prev, spl_val_c, spl_val_s, Bmnc, Bmns)

    ! Flux-surface dependent quantities
    !$omp threadprivate (mph)

    ! Shared data (NOT threadprivate): params, modes, mb, nb, nflux, nfp, nmode,
    ! spl_coeff1, spl_coeff2, ncol1, ncol2, mph_shared

contains

    subroutine magfie_pert_thread_init()
        ! Initialize threadprivate variables for this thread
        ! Must be called once per thread before using magfie_pert routines
        magfie_pert_arrays_initialized = .false.
        s_prev = -1.0_dp
    end subroutine magfie_pert_thread_init

    subroutine read_boozer_pert_file(path)
        ! Main thread only: Read perturbation Boozer file and prepare shared spline data
        ! This should be called ONCE before parallel region
        character(len=*), intent(in) :: path
        integer :: j, k

        ncol1 = 5
        if (inp_swi_pert == 8) ncol2 = 4 ! tok_circ
        if (inp_swi_pert == 9) ncol2 = 8 ! ASDEX

        ! allocates params, modes - SHARED arrays
        call boozer_read_pert(path)

        ! Update mph_shared from perturbation file (for use in parallel regions)
        mph_shared = nint(nfp * modes(1, 1, 2))

        ! Allocate shared spline coefficient arrays
        if (allocated(spl_coeff1)) then
            if (size(spl_coeff1, 1) /= nflux - 1) deallocate(spl_coeff1, spl_coeff2)
        end if
        if (.not. allocated(spl_coeff1)) then
            allocate (spl_coeff1(nflux - 1, 5, ncol1))
            allocate (spl_coeff2(nflux - 1, 5, ncol2, nmode))
        end if

        do k = 1, ncol1
            spl_coeff1(:, :, k) = spline_coeff(params(:, 1), params(:, k + 1))
        end do

        do j = 1, nmode
            do k = 1, ncol2
                spl_coeff2(:, :, k, j) = spline_coeff(params(:, 1), modes(:, j, k + 2))
            end do
        end do
    end subroutine read_boozer_pert_file

    subroutine init_magfie_pert_at_s()
        ! Per-thread: Initialize perturbation field at current s value
        ! Assumes s is already set via set_s()
        ! Allocates threadprivate working buffers and computes s-dependent values
        real(dp) :: x(3)
        complex(dp) :: dummy

        ! Auto-initialize threadprivate state if not yet done for this thread
        if (magfie_pert_thread_init_state /= MAGFIE_PERT_INIT_SENTINEL) then
            call magfie_pert_thread_init()
            magfie_pert_thread_init_state = MAGFIE_PERT_INIT_SENTINEL
        end if

        ! Allocate threadprivate working buffers (safe for undefined allocation status)
        if (.not. magfie_pert_arrays_initialized) then
            if (allocated(Bmnc)) deallocate(Bmnc)
            allocate(Bmnc(nmode))
            if (ncol2 >= 8) then
                if (allocated(Bmns)) deallocate(Bmns)
                allocate(Bmns(nmode))
            end if

            if (allocated(spl_val_c)) deallocate(spl_val_c)
            if (allocated(spl_val_s)) deallocate(spl_val_s)
            allocate(spl_val_c(3, nmode), spl_val_s(3, nmode))

            magfie_pert_arrays_initialized = .true.
        end if

        ! Initialize cache
        s_prev = -1.0_dp

        ! Compute mph at current s
        mph = nint(nfp * modes(1, 1, 2))

        x(1) = s
        x(2) = 0.0
        x(3) = 0.0
        call do_magfie_pert_amp(x, dummy)
    end subroutine init_magfie_pert_at_s

    subroutine set_mph(mph_value)
        integer, intent(in) :: mph_value
        mph = mph_value
        mph_shared = mph_value
    end subroutine set_mph

    subroutine init_mph_from_shared()
        mph = mph_shared
    end subroutine init_mph_from_shared

    subroutine do_magfie_pert_init(path)
        ! Backward compatibility wrapper: calls read_boozer_pert_file + init_magfie_pert_at_s
        character(len=*), intent(in) :: path

        call read_boozer_pert_file(path)
        call init_magfie_pert_at_s()
    end subroutine do_magfie_pert_init

    subroutine do_magfie_pert_amp(x, bamp)
        real(dp), dimension(:), intent(in) :: x
        complex(dp), intent(out) :: bamp

        real(dp) :: x1

        ! safety measure in order not to extrapolate
        x1 = max(params(1, 1), x(1))
        x1 = min(params(nflux, 1), x1)

        ! calculate B-field from modes
        if (inp_swi_pert == 8) then
            call cached_spline(x1, s_prev, spl_coeff2(:, :, 4, :), spl_val_c)
            Bmnc(:) = 1.0e4_dp * spl_val_c(1, :) * bfac
            bamp = sum(Bmnc * cos(modes(1, :, 1) * x(3)))
        else if (inp_swi_pert == 9) then
            call cached_spline(x1, s_prev, spl_coeff2(:, :, 7, :), spl_val_c)
            call cached_spline(x1, s_prev, spl_coeff2(:, :, 8, :), spl_val_s)
            Bmnc(:) = 1.0e4_dp * spl_val_c(1, :) * bfac
            Bmns(:) = 1.0e4_dp * spl_val_s(1, :) * bfac
            bamp = fast_fourier_sum(Bmnc, Bmns, modes(1, :, 1), x(3))
        end if

        s_prev = x1
    end subroutine do_magfie_pert_amp

    subroutine do_magfie_pert(x, bmod)
        real(dp), dimension(:), intent(in) :: x
        real(dp), intent(out) :: bmod
        complex(dp) :: bamp

        call do_magfie_pert_amp(x, bamp)
        bmod = real(sum(bamp*exp(imun*nfp*modes(1, :, 2)*x(2))))

    end subroutine do_magfie_pert

    subroutine boozer_read_pert(filename)
        integer :: ksurf, kmode
        real(dp) :: flux, dummy
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
        real(dp), dimension(:), intent(in) :: m
        integer, intent(in) :: kmode

        real(dp), parameter :: tol = 1.0e-10_dp
        real(dp) :: diff, diff2

        if (kmode < 3) return

        diff = m(kmode - 1) - m(kmode - 2)
        diff2 = m(kmode) - m(kmode - 1)
        if (abs(diff - diff2) > tol) then
            error stop "harmonics not equally spaced and ascending"
        end if
    end subroutine check_equal_space_ascending

    subroutine check_equal(m, kmode)
        real(dp), dimension(:), intent(in) :: m
        integer, intent(in) :: kmode

        real(dp), parameter :: tol = 1.0e-10_dp

        if (kmode < 2) return

        if (abs(m(kmode) - m(kmode - 1)) > tol) then
            error stop "harmonics not equal"
        end if
    end subroutine check_equal

    function fast_fourier_sum(fmnc, fmns, m, x)
        ! Fast Fourier sum that assumes equally spaced ascending mode numbers
        real(dp), dimension(:), intent(in) :: fmnc, fmns, m
        real(dp), intent(in) :: x
        complex(dp) :: fast_fourier_sum

        real(dp) :: dm
        complex(dp) :: fourier_factor, rotation
        integer :: j

        dm = m(2) - m(1)
        fourier_factor  = exp(imun*m(1)*x)
        rotation = exp(imun*dm*x)

        fast_fourier_sum = (0.0_dp, 0.0_dp)
        do j = 1, size(m)
            fast_fourier_sum = fast_fourier_sum + &
                (fmnc(j) - imun*fmns(j))*fourier_factor
            fourier_factor = fourier_factor*rotation
        end do
    end function fast_fourier_sum
end module do_magfie_pert_mod
