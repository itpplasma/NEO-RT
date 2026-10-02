module test_raw_fixture
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_raw_orbit, only: raw_background_t, raw_species_t
    use neort_raw_drive, only: RAW_OK, RAW_BAD_INPUT
    use neo_perturbation_field, only: UNITS_SI
    implicit none
    private

    real(dp), parameter, public :: speed_unit = 48242.666078326634_dp
    real(dp), parameter, public :: time_unit = 1.0_dp/speed_unit
    real(dp), parameter, public :: charge_si = 1.602176634e-19_dp
    real(dp), parameter, public :: mass_si = 2.0_dp*1.66053906660e-27_dp
    real(dp), parameter, public :: energy_unit = charge_si*speed_unit
    real(dp), parameter, public :: moment = 0.0005151488525985234_dp*energy_unit
    real(dp), parameter, public :: slope = 0.002942126607849932_dp
    real(dp), parameter, public :: epsilon = 1.0e-6_dp, baseline = 1.0e-5_dp
    real(dp), parameter, public :: chi_scale = 0.01_dp
    real(dp), public :: gauge_omega = 0.0_dp

    public :: circular_background, circular_perturbation, linear_gauge
    public :: read_reference, species_si

contains

    function species_si() result(species)
        type(raw_species_t) :: species

        species%mass = mass_si
        species%charge = charge_si
        species%units = UNITS_SI
    end function species_si

    subroutine circular_background(x, background, ierr, data)
        real(dp), intent(in) :: x(3)
        type(raw_background_t), intent(out) :: background
        integer, intent(out) :: ierr
        class(*), intent(in), optional :: data
        real(dp) :: r, cosine, sine, p, pr, length, flux

        background = raw_background_t()
        ierr = RAW_BAD_INPUT
        r = hypot(x(1) - 1.0_dp, x(3))
        if (r <= 0.001_dp .or. r >= 0.75_dp) return
        cosine = (x(1) - 1.0_dp)/r
        sine = x(3)/r
        p = r/(1.5_dp + 8.0_dp*r*r)
        pr = (1.5_dp - 8.0_dp*r*r)/(1.5_dp + 8.0_dp*r*r)**2
        length = sqrt(1.0_dp + p*p)
        flux = log(1.0_dp + 8.0_dp*r*r/1.5_dp)/16.0_dp
        call circular_geometry(x(1), r, cosine, sine, p, pr, length, background)
        background%grad_phi(1) = slope*p*cosine*speed_unit
        background%grad_phi(3) = slope*p*sine*speed_unit
        background%potential = slope*flux*speed_unit
        background%a(2) = flux/x(1)
        background%a(3) = -log(x(1))
        background%has_vector_potential = .true.
        if (present(data)) continue
        ierr = RAW_OK
    end subroutine circular_background

    subroutine circular_geometry(radius, r, cosine, sine, p, pr, length, background)
        real(dp), intent(in) :: radius, r, cosine, sine, p, pr, length
        type(raw_background_t), intent(inout) :: background
        real(dp) :: br, bt, fp, gp

        background%b(1) = -p*sine/radius
        background%b(2) = 1.0_dp/radius
        background%b(3) = p*cosine/radius
        br = p*pr/(length*radius) - length*cosine/radius**2
        bt = length*sine/radius**2
        background%grad_b(1) = br*cosine - bt*sine
        background%grad_b(3) = br*sine + bt*cosine
        fp = pr/length**3
        gp = -p*pr/length**3
        background%curl_bhat(1) = -gp*sine
        background%curl_bhat(2) = -(fp + p/(length*r))
        background%curl_bhat(3) = gp*cosine + 1.0_dp/(length*radius)
    end subroutine circular_geometry

    subroutine circular_perturbation(n, radius, z, a, potential)
        integer, intent(in) :: n
        real(dp), intent(in) :: radius, z
        complex(dp), intent(out) :: a(3), potential
        complex(dp) :: ztheta
        real(dp) :: r, center

        if (n /= 1) error stop "circular physical fixture requires n one"
        r = hypot(radius - 1.0_dp, z)
        if (r <= 0.0_dp .or. r >= 1.0_dp) error stop "invalid circular field grid"
        ztheta = cmplx((radius - 1.0_dp)/r + r, &
            sqrt(1.0_dp - r*r)*z/r, dp)/radius
        center = sqrt((sqrt(77.0_dp) - 7.0_dp)/32.0_dp)
        a = cmplx(0.0_dp, 0.0_dp, dp)
        a(2) = epsilon*exp(-((r - center)/0.12_dp)**2)*ztheta**(-2)/radius
        potential = baseline*speed_unit &
            *exp(-((r - center)/0.18_dp)**2)/ztheta
    end subroutine circular_perturbation

    subroutine linear_gauge(n, radius, z, a, potential)
        integer, intent(in) :: n
        real(dp), intent(in) :: radius, z
        complex(dp), intent(out) :: a(3), potential

        if (n /= 1) error stop "linear gauge fixture requires n one"
        if (abs(z) > 1.0_dp) error stop "unexpected gauge grid"
        a(1) = cmplx(chi_scale, 0.0_dp, dp)
        a(2) = cmplx(0.0_dp, chi_scale, dp)
        a(3) = cmplx(0.0_dp, 0.0_dp, dp)
        potential = cmplx(0.0_dp, gauge_omega*chi_scale*radius, dp)
    end subroutine linear_gauge

    subroutine read_reference(path, reference)
        character(len=*), intent(in) :: path
        real(dp), allocatable, intent(out) :: reference(:, :)
        integer :: unit, n, k, ios

        open (newunit=unit, file=path, status='old', action='read', iostat=ios)
        if (ios /= 0) error stop "independent raw GC reference missing"
        read (unit, *, iostat=ios) n
        if (ios /= 0) error stop "raw GC reference count missing"
        if (n < 3) error stop "raw GC reference needs full orbit"
        allocate (reference(21, n))
        do k = 1, n
            read (unit, *, iostat=ios) reference(:, k)
            if (ios /= 0) error stop "raw GC reference row truncated"
        end do
        close (unit)
    end subroutine read_reference

end module test_raw_fixture
