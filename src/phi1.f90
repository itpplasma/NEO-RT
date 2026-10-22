module neort_phi1
    use iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use util, only: imun, spline_cache_tol
    use spline, only: spline_coeff, spline_val_0
    use do_magfie_pert_mod, only: horner_fourier_sum, fourier_chains, mph_shared
    implicit none
    private
    public :: phi1, phi1file, phi1_scale, phi_mn, m0_phi, reset_phi1
    public :: read_phi1_file, do_elec_pert_amp, phi1_amplitude_bound

    logical :: phi1 = .false.
    character(1024) :: phi1file = ''
    real(dp) :: phi1_scale = 1.0_dp
    complex(dp) :: phi_mn = (0.0_dp, 0.0_dp)
    integer :: m0_phi = 0
    real(dp), allocatable :: surfaces(:), coefficients(:, :, :, :), values(:, :, :)
    real(dp) :: first_m, delta_m
    integer :: mode_count, generation = 0
    real(dp), allocatable :: cached_re(:), cached_im(:)
    real(dp) :: cached_s = -huge(1.0_dp)
    integer :: cached_generation = -1
    !$omp threadprivate(cached_re, cached_im, cached_s, cached_generation)
contains
    subroutine reset_phi1()
        phi1 = .false.
        phi1file = ''
        phi1_scale = 1.0_dp
        phi_mn = (0.0_dp, 0.0_dp)
        m0_phi = 0
        generation = generation + 1
        if (allocated(surfaces)) deallocate (surfaces, values)
        if (allocated(coefficients)) deallocate (coefficients)
    end subroutine reset_phi1

    subroutine read_phi1_file()
        integer :: unit, ns, nper, k, j
        real(dp) :: params(6), row(10)
        if (.not. phi1) return
        if (.not. all(ieee_is_finite([phi1_scale, real(phi_mn, dp), aimag(phi_mn)]))) &
            error stop 'nonfinite Phi_1 configuration'
        if (len_trim(phi1file) == 0) return
        if (allocated(surfaces)) deallocate (surfaces, values)
        if (allocated(coefficients)) deallocate (coefficients)
        open (newunit=unit, file=trim(phi1file), status='old', action='read')
        call read_phi1_header(unit, ns, nper)
        allocate (surfaces(ns), values(ns, 2, mode_count))
        do k = 1, ns
            read (unit, '(/)')
            read (unit, *) params
            surfaces(k) = params(1)
            read (unit, *)
            do j = 1, mode_count
                read (unit, *) row
                call check_phi1_row(row, k, j, nper)
                values(k, :, j) = row(9:10)
            end do
        end do
        close (unit)
        call prepare_phi1_splines()
        generation = generation + 1
    end subroutine read_phi1_file

    subroutine read_phi1_header(unit, ns, nper)
        integer, intent(in) :: unit
        integer, intent(out) :: ns, nper
        integer :: mb, nb, status
        real(dp) :: flux, a, radius
        character(1024) :: line
        logical :: quantity
        quantity = .false.
        do
            read (unit, '(a)', iostat=status) line
            if (status /= 0) error stop 'missing Phi_1 header'
            if (trim(line) == 'CC quantity: Phi_1 [statvolt] in bmnc/bmns columns') &
                quantity = .true.
            if (index(adjustl(line), 'CC') == 1) cycle
            if (index(line, 'm0b') > 0) cycle
            read (line, *, iostat=status) mb, nb, ns, nper, flux, a, radius
            if (status /= 0) error stop 'invalid Phi_1 header'
            exit
        end do
        if (.not. quantity) error stop 'missing mandatory Phi_1 quantity comment'
        if (mb < 0 .or. nb < 0 .or. ns < 1 .or. nper < 1) &
            error stop 'invalid Phi_1 dimensions'
        mode_count = (mb + 1)*(nb + 1)
    end subroutine read_phi1_header

    subroutine check_phi1_row(row, k, j, nper)
        real(dp), intent(in) :: row(10)
        integer, intent(in) :: k, j, nper
        if (.not. all(ieee_is_finite(row))) error stop 'nonfinite Phi_1 row'
        if (any(row(1:2) /= real(nint(row(1:2)), dp))) &
            error stop 'noninteger Phi_1 harmonic'
        if (nper*nint(row(2)) /= mph_shared) error stop 'Phi_1 n must equal mph'
        if (any(row(3:8) /= 0.0_dp)) error stop 'nonzero Phi_1 geometry columns'
        if (k == 1 .and. j == 1) then
            first_m = row(1)
            delta_m = 1.0_dp
        else if (k == 1 .and. j == 2) then
            delta_m = row(1) - first_m
            if (delta_m <= 0.0_dp) error stop 'Phi_1 m must ascend'
        end if
        if (row(1) /= first_m + (j - 1)*delta_m) &
            error stop 'Phi_1 harmonics must be equally spaced on every surface'
    end subroutine check_phi1_row

    subroutine prepare_phi1_splines()
        integer :: j, k, ns
        ns = size(surfaces)
        if (.not. all(ieee_is_finite(surfaces))) error stop 'nonfinite Phi_1 s'
        if (any(surfaces(2:) <= surfaces(:ns - 1))) error stop 'Phi_1 s must ascend'
        if (ns == 1) then
            print *, 'Phi_1: one surface, constant in s (clamped to supplied surface)'
            return
        end if
        allocate (coefficients(ns - 1, 5, 2, mode_count))
        do j = 1, mode_count
            do k = 1, 2
                if (ns == 2) then
                    coefficients(1, :, k, j) = [surfaces(1), values(1, k, j), &
                      (values(2, k, j) - values(1, k, j))/(surfaces(2) - surfaces(1)), &
                                                0.0_dp, 0.0_dp]
                else
                    coefficients(:, :, k, j) = spline_coeff(surfaces, values(:, k, j))
                end if
            end do
        end do
    end subroutine prepare_phi1_splines

    subroutine refresh_phi1(x)
        real(dp), intent(in) :: x(3)
        real(dp) :: surface, val(3)
        integer :: j, padded
        if (.not. allocated(surfaces)) error stop 'Phi_1 file was not initialized'
        surface = min(max(x(1), surfaces(1)), surfaces(size(surfaces)))
        if (cached_generation == generation) then
            if (abs(surface - cached_s) < spline_cache_tol) return
        else
            if (allocated(cached_re)) deallocate (cached_re, cached_im)
            padded = fourier_chains*((mode_count - 1)/fourier_chains + 1)
            allocate (cached_re(padded), cached_im(padded))
        end if
        cached_re = 0.0_dp
        cached_im = 0.0_dp
        do j = 1, mode_count
            cached_re(j) = values(1, 1, j)
            cached_im(j) = -values(1, 2, j)
            if (size(surfaces) == 1) cycle
            val = spline_val_0(coefficients(:, :, 1, j), surface)
            cached_re(j) = val(1)
            val = spline_val_0(coefficients(:, :, 2, j), surface)
            cached_im(j) = -val(1)
        end do
        cached_s = surface
        cached_generation = generation
    end subroutine refresh_phi1

    subroutine do_elec_pert_amp(x, amplitude)
        real(dp), intent(in) :: x(3)
        complex(dp), intent(out) :: amplitude
        amplitude = (0.0_dp, 0.0_dp)
        if (.not. phi1) return
        if (len_trim(phi1file) == 0) then
            amplitude = phi_mn*exp(imun*m0_phi*x(3))
        else
            call refresh_phi1(x)
            amplitude = horner_fourier_sum(cached_re, cached_im, first_m, delta_m, x(3))
        end if
        amplitude = phi1_scale*amplitude
    end subroutine do_elec_pert_amp

    subroutine phi1_amplitude_bound(x, amplitude)
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: amplitude
        amplitude = 0.0_dp
        if (.not. phi1) return
        if (len_trim(phi1file) == 0) then
            amplitude = abs(phi_mn)
        else
            call refresh_phi1(x)
            amplitude = sum(sqrt(cached_re**2 + cached_im**2))
        end if
        amplitude = abs(phi1_scale)*amplitude
    end subroutine phi1_amplitude_bound
end module neort_phi1
