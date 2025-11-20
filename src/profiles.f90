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

    ! Plasma profile spline data
    integer :: nplasma_global = 0
    real(dp) :: am1_global = 0d0, am2_global = 0d0, Z1_global = 0d0, Z2_global = 0d0
    real(dp), allocatable :: plasma_spl_coeff(:,:,:)

    ! Rotation profile spline data
    real(dp), allocatable :: Mt_spl_coeff(:,:)
    ! Flux-surface dependent quantities (interpolated per-thread at each s)
    !$omp threadprivate (vth, dvthds, M_t, dM_tds, Om_tE, dOm_tEds)
    !$omp threadprivate (ni1, ni2, Ti1, Ti2, Te, dni1ds, dni2ds, dTi1ds, dTi2ds, dTeds)
    !$omp threadprivate (A1, A2)

contains

    subroutine init_profiles(R0)
        real(dp), intent(in) :: R0

        ! vth: set in read_and_set_control, overwritten in init_plasma_input
        ! M_t: set in read_and_set_control, overwritten in init_profile_input
        Om_tE = vth * M_t / R0  ! toroidal ExB drift frequency
        dOm_tEds = 0d0
    end subroutine init_profiles

    subroutine prepare_plasma_splines(nplasma, am1, am2, Z1, Z2, plasma)
        ! Main thread only: Prepare shared spline coefficients from plasma data
        ! This should be called ONCE before parallel region
        integer, intent(in) :: nplasma
        real(dp), intent(in) :: am1, am2, Z1, Z2
        real(dp), intent(in) :: plasma(:, :)

        integer, parameter :: NCOL = 6
        integer :: k

        nplasma_global = nplasma
        am1_global = am1
        am2_global = am2
        Z1_global = Z1
        Z2_global = Z2

        ! Allocate and compute spline coefficients (SHARED)
        if (allocated(plasma_spl_coeff)) deallocate(plasma_spl_coeff)
        allocate (plasma_spl_coeff(nplasma - 1, 5, NCOL))

        do k = 1, 5
            plasma_spl_coeff(:, :, k) = spline_coeff(plasma(:, 1), plasma(:, k + 1))
        end do
    end subroutine prepare_plasma_splines

    subroutine init_plasma_at_s()
        ! Per-thread: Interpolate plasma profiles at current s value
        ! Assumes s is already set via set_s() and prepare_plasma_splines() was called
        use do_magfie_mod, only: s
        real(dp) :: spl_val(3)
        real(dp), parameter :: pmass = 1.6726d-24
        real(dp) :: amb, Zb, dchichi, slowrate, dchichi_norm, slowrate_norm
        real(dp) :: v0, ebeam

        spl_val = spline_val_0(plasma_spl_coeff(:, :, 1), s)
        ni1 = spl_val(1)
        dni1ds = spl_val(2)

        spl_val = spline_val_0(plasma_spl_coeff(:, :, 2), s)
        ni2 = spl_val(1)
        dni2ds = spl_val(2)

        spl_val = spline_val_0(plasma_spl_coeff(:, :, 3), s)
        Ti1 = spl_val(1)
        dTi1ds = spl_val(2)

        spl_val = spline_val_0(plasma_spl_coeff(:, :, 4), s)
        Ti2 = spl_val(1)
        dTi2ds = spl_val(2)

        spl_val = spline_val_0(plasma_spl_coeff(:, :, 5), s)
        Te = spl_val(1)
        dTeds = spl_val(2)

        ! Compute derived quantities
        qi = Z1_global*qe
        mi = am1_global*mu
        vth = sqrt(2d0*Ti1*ev/mi)
        dvthds = 0.5d0*sqrt(2d0*ev/(mi*Ti1))*dTi1ds

        ! Call collision routine
        v0 = vth
        amb = 2d0
        Zb = 1d0
        ebeam = amb*pmass*v0**2/(2d0*ev)
        call loacol_nbi(amb, am1_global, am2_global, Zb, Z1_global, Z2_global, &
                        ni1, ni2, Ti1, Ti2, Te, ebeam, v0, &
                        dchichi, slowrate, dchichi_norm, slowrate_norm)
    end subroutine init_plasma_at_s

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

        deallocate(spl_coeff)
    end subroutine init_plasma_input

    subroutine read_and_init_plasma_input(path, s_in)
        ! Backward compatibility wrapper: Read plasma file and initialize at given s
        ! This combines read_plasma_input, prepare_plasma_splines, and init_plasma_at_s
        use do_magfie_mod, only: s
        character(len=*), intent(in) :: path
        real(dp), intent(in) :: s_in

        integer :: nplasma
        real(dp) :: am1, am2, Z1, Z2
        real(dp), allocatable :: plasma(:, :)

        call read_plasma_input(path, nplasma, am1, am2, Z1, Z2, plasma)

        ! Prepare shared spline coefficients (main thread)
        call prepare_plasma_splines(nplasma, am1, am2, Z1, Z2, plasma)

        s = s_in
        call init_plasma_at_s()

        deallocate (plasma)
    end subroutine read_and_init_plasma_input

    subroutine prepare_profile_splines(data)
        ! Main thread only: Prepare shared spline coefficients from rotation profile data
        ! This should be called ONCE before parallel region
        real(8), intent(in) :: data(:, :)

        ! Deallocate if already allocated
        if (allocated(Mt_spl_coeff)) deallocate(Mt_spl_coeff)

        ! Allocate and compute spline coefficients (SHARED)
        allocate (Mt_spl_coeff(size(data, 1) - 1, 5))
        Mt_spl_coeff = spline_coeff(data(:, 1), data(:, 2))
    end subroutine prepare_profile_splines

    subroutine init_profile_at_s(R0, efac, bfac)
        ! Per-thread: Interpolate rotation profile at current s value
        ! Assumes s is already set and prepare_profile_splines() was called
        use do_magfie_mod, only: s
        real(8), intent(in) :: R0, efac, bfac
        real(8) :: splineval(3)

        ! Interpolate M_t at s (writes threadprivate variables)
        splineval = spline_val_0(Mt_spl_coeff, s)

        M_t = splineval(1)*efac/bfac
        dM_tds = splineval(2)*efac/bfac

        ! Compute ExB drift frequency (same functionality as init_profiles)
        Om_tE = vth*M_t/R0
        dOm_tEds = vth*dM_tds/R0 + M_t*dvthds/R0
    end subroutine init_profile_at_s

    subroutine init_profile_input(s, R0, efac, bfac, data)
        ! Init s profile for finite orbit width boxes in radial s
        real(8), intent(in) :: s, R0, efac, bfac
        real(8), allocatable, intent(in) :: data(:, :)

        ! For splining electric precession frequency
        real(8), allocatable :: local_Mt_spl_coeff(:, :)

        real(8) :: splineval(3)

        allocate (local_Mt_spl_coeff(size(data, 1) - 1, 5))

        local_Mt_spl_coeff = spline_coeff(data(:, 1), data(:, 2))
        splineval = spline_val_0(local_Mt_spl_coeff, s)

        M_t = splineval(1)*efac/bfac
        dM_tds = splineval(2)*efac/bfac

        ! note: same functionality as init_profiles
        Om_tE = vth*M_t/R0                   ! toroidal ExB drift frequency
        dOm_tEds = vth*dM_tds/R0 + M_t*dvthds/R0

        deallocate(local_Mt_spl_coeff)
    end subroutine init_profile_input

    subroutine read_and_init_profile_input(path, s_in, R0, efac, bfac)
        ! Backward compatibility wrapper: Read profile file and initialize at given s
        ! This combines readdata, prepare_profile_splines, and init_profile_at_s
        use do_magfie_mod, only: s
        character(len=*), intent(in) :: path
        real(8), intent(in) :: s_in, R0, efac, bfac

        real(8), allocatable :: data(:, :)

        call readdata(path, 2, data)  ! allocates data

        ! Prepare shared spline coefficients (main thread)
        call prepare_profile_splines(data)

        s = s_in
        call init_profile_at_s(R0, efac, bfac)

        deallocate (data)
    end subroutine read_and_init_profile_input

    subroutine init_thermodynamic_forces(psi_pr, q)
        real(dp), intent(in) :: psi_pr  ! toroidal flux at plasma boundary == dpsi_tor/ds
        real(dp), intent(in) :: q  ! safety factor

        A1 = dni1ds/ni1 - qi/(Ti1*ev)*sign_theta*psi_pr/(q*c)*Om_tE - 3d0/2d0*dTi1ds/Ti1
        A2 = dTi1ds/Ti1
    end subroutine init_thermodynamic_forces

end module neort_profiles
