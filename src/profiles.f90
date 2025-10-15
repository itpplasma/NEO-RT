module neort_profiles
    use iso_fortran_env, only: dp => real64
    use spline, only: spline_coeff, spline_val_0
    use util, only: readdata, ev, qi, mi, mu, qe, c
    use collis_alp, only: loacol_nbi
    use do_magfie_mod, only: sign_theta

    implicit none

    ! Thermal velocity, Mach number, and their derivatives over s_tor
    real(dp) :: vth = 0d0, dvthds = 0d0, M_t = 0d0, dM_tds = 0d0

    ! Electric precession frequency and its derivative over s_tor
    real(8) :: Om_tE = 0d0, dOm_tEds = 0d0

    ! density and temperature profiles and their derivatives over s_tor
    real(dp) :: ni1 = 0d0, ni2 = 0d0, Ti1 = 0d0, Ti2 = 0d0, Te = 0d0, dni1ds = 0d0, &
                dni2ds = 0d0, dTi1ds = 0d0, dTi2ds = 0d0, dTeds = 0d0

    ! Thermodynamic forces in radial variable s_tor
    real(dp) :: A1 = 0d0, A2 = 0d0

contains

    subroutine init_profiles(R0)
        real(dp), intent(in) :: R0
        Om_tE = vth*M_t/R0                   ! toroidal ExB drift frequency
        dOm_tEds = 0d0
    end subroutine init_profiles

    subroutine read_plasma_input(path, nplasma, am1, am2, Z1, Z2, plasma)
        character(len=*), intent(in) :: path
        integer, intent(out) :: nplasma
        real(dp), intent(out) :: am1, am2, Z1, Z2
        real(dp), allocatable, intent(out) :: plasma(:, :)

        integer :: k
        integer, parameter :: NCOL = 6
        integer, parameter :: fd = 1

        open (fd, file=path, status="old")
        read (fd, *)
        read (fd, *) nplasma, am1, am2, Z1, Z2
        read (fd, *)
        allocate (plasma(nplasma, NCOL))
        do k = 1, nplasma
            read (fd, *) plasma(k, :)
        end do
        close (fd)

    end subroutine read_plasma_input

    subroutine init_plasma_input(s, nplasma, am1, am2, Z1, Z2, plasma)
        use spline, only: spline_coeff, spline_val_0
        real(dp), intent(in) :: s
        integer, intent(in) :: nplasma
        real(dp), intent(in) :: am1, am2, Z1, Z2
        real(dp), allocatable, intent(in) :: plasma(:, :)

        real(dp), parameter :: pmass = 1.6726d-24
        integer, parameter :: NCOL = 6

        real(dp) :: amb, Zb, dchichi, slowrate, dchichi_norm, slowrate_norm
        real(dp) :: v0, ebeam
        real(dp), allocatable :: spl_coeff(:, :, :)
        real(dp) :: spl_val(3)
        integer :: k

        allocate (spl_coeff(nplasma - 1, 5, NCOL))

        do k = 1, 5
            spl_coeff(:, :, k) = spline_coeff(plasma(:, 1), plasma(:, k + 1))
        end do

        spl_val = spline_val_0(spl_coeff(:, :, 1), s)
        ni1 = spl_val(1)
        dni1ds = spl_val(2)
        spl_val = spline_val_0(spl_coeff(:, :, 2), s)
        ni2 = spl_val(1)
        dni2ds = spl_val(2)
        spl_val = spline_val_0(spl_coeff(:, :, 3), s)
        Ti1 = spl_val(1)
        dTi1ds = spl_val(2)
        spl_val = spline_val_0(spl_coeff(:, :, 4), s)
        Ti2 = spl_val(1)
        dTi2ds = spl_val(2)
        spl_val = spline_val_0(spl_coeff(:, :, 5), s)
        Te = spl_val(1)
        dTeds = spl_val(2)

        qi = Z1*qe
        mi = am1*mu
        vth = sqrt(2d0*Ti1*ev/mi)
        dvthds = 0.5d0*sqrt(2d0*ev/(mi*Ti1))*dTi1ds
        v0 = vth
        amb = 2d0
        Zb = 1d0
        ebeam = amb*pmass*v0**2/(2d0*ev)

        call loacol_nbi(amb, am1, am2, Zb, Z1, Z2, ni1, ni2, Ti1, Ti2, Te, &
                        ebeam, v0, dchichi, slowrate, dchichi_norm, slowrate_norm)

    end subroutine init_plasma_input

    subroutine read_and_init_plasma_input(path, s)
        character(len=*), intent(in) :: path
        real(dp), intent(in) :: s

        integer :: nplasma
        real(dp) :: am1, am2, Z1, Z2
        real(dp), allocatable :: plasma(:, :)

        call read_plasma_input(path, nplasma, am1, am2, Z1, Z2, plasma)  ! allocates plasma
        call init_plasma_input(s, nplasma, am1, am2, Z1, Z2, plasma)

        deallocate (plasma)
    end subroutine read_and_init_plasma_input

    subroutine init_profile_input(s, R0, efac, bfac, data)
        ! Init s profile for finite orbit width boxes in radial s
        real(8), intent(in) :: s, R0, efac, bfac
        real(8), allocatable, intent(in) :: data(:, :)

        ! For splining electric precession frequency
        real(8), allocatable :: Mt_spl_coeff(:, :)

        real(8) :: splineval(3)

        allocate (Mt_spl_coeff(size(data(:, 1)), 5))

        Mt_spl_coeff = spline_coeff(data(:, 1), data(:, 2))
        splineval = spline_val_0(Mt_spl_coeff, s)

        M_t = splineval(1)*efac/bfac
        dM_tds = splineval(2)*efac/bfac

        Om_tE = vth*M_t/R0                   ! toroidal ExB drift frequency
        dOm_tEds = vth*dM_tds/R0 + M_t*dvthds/R0
    end subroutine init_profile_input

    subroutine read_and_init_profile_input(path, s, R0, efac, bfac)
        character(len=*), intent(in) :: path
        real(8), intent(in) :: s, R0, efac, bfac

        real(8), allocatable :: data(:, :)

        call readdata(path, 3, data)  ! allocates data
        call init_profile_input(s, R0, efac, bfac, data)

        deallocate (data)
    end subroutine read_and_init_profile_input

    subroutine init_thermodynamic_forces(psi_pr, q)
        real(dp), intent(in) :: psi_pr, q

        A1 = dni1ds/ni1 - qi/(Ti1*ev)*sign_theta*psi_pr/(q*c)*Om_tE - 3d0/2d0*dTi1ds/Ti1
        A2 = dTi1ds/Ti1
    end subroutine init_thermodynamic_forces

end module neort_profiles
