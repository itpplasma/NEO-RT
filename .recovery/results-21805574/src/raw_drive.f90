module neort_raw_drive
    !! Retained-order Hamiltonian in physical cylindrical components:
    !! H = mu*b0.curl(A) + q*Phi - (q/c)*A.Xdot in Gaussian units.
    !! SI replaces q/c by q. Sources include exp(+i*n*phi-i*omega*t).
    !! Xdot must advance the same orbit used in harmonic projection.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use neo_perturbation_field, only: perturbation_field_t, PERTFIELD_OK, &
        UNITS_GAUSSIAN, UNITS_SI
    use neo_perturbation_field_netcdf, only: read_perturbation_field_netcdf
    implicit none
    private

    integer, parameter, public :: RAW_OK = 0, RAW_BAD_INPUT = 1, RAW_FIELD_ERROR = 2
    real(dp), parameter, public :: C_GAUSSIAN = 29979245800.0_dp

    type, public :: raw_source_t
        type(perturbation_field_t) :: field
        real(dp), allocatable :: omega(:)
    contains
        procedure :: load => load_source
        procedure :: set_frequencies
        procedure :: evaluate
        procedure :: ready
    end type raw_source_t

    type :: raw_sample_t
        real(dp) :: time, x(3), velocity(3), b0(3)
    end type raw_sample_t

    type :: raw_projection_t
        real(dp) :: moment, charge, omega_b, omega_phi
    end type raw_projection_t

    public :: raw_harmonics

contains

    subroutine load_source(self, path, frequencies, ierr)
        class(raw_source_t), intent(inout) :: self
        character(len=*), intent(in) :: path
        real(dp), intent(in) :: frequencies(:)
        integer, intent(out) :: ierr
        integer :: field_status

        if (allocated(self%omega)) deallocate (self%omega)
        call read_perturbation_field_netcdf(path, self%field, field_status)
        ierr = RAW_FIELD_ERROR
        if (field_status /= PERTFIELD_OK) return
        call self%set_frequencies(frequencies, ierr)
    end subroutine load_source

    subroutine set_frequencies(self, frequencies, ierr)
        class(raw_source_t), intent(inout) :: self
        real(dp), intent(in) :: frequencies(:)
        integer, intent(out) :: ierr

        ierr = RAW_BAD_INPUT
        if (allocated(self%omega)) deallocate (self%omega)
        if (self%field%n_modes < 1) return
        if (size(frequencies) /= self%field%n_modes) return
        if (.not. all(ieee_is_finite(frequencies))) return
        if (self%field%units /= UNITS_GAUSSIAN) then
            if (self%field%units /= UNITS_SI) return
        end if
        allocate (self%omega, source=frequencies)
        ierr = RAW_OK
    end subroutine set_frequencies

    logical function ready(self)
        class(raw_source_t), intent(in) :: self

        ready = .false.
        if (.not. allocated(self%omega)) return
        if (.not. allocated(self%field%ntor)) return
        if (self%field%n_modes < 1) return
        if (size(self%field%ntor) /= self%field%n_modes) return
        if (size(self%omega) /= self%field%n_modes) return
        if (self%field%units /= UNITS_GAUSSIAN) then
            if (self%field%units /= UNITS_SI) return
        end if
        ready = .true.
    end function ready

    subroutine evaluate(self, x, t, b0, velocity, moment, charge, h, ierr)
        class(raw_source_t), intent(in) :: self
        real(dp), intent(in) :: x(3), t, b0(3), velocity(3), moment, charge
        complex(dp), intent(out) :: h(:)
        integer, intent(out) :: ierr
        complex(dp) :: a(3, self%field%n_modes), db(3, self%field%n_modes)
        complex(dp) :: potential(self%field%n_modes), phase
        real(dp) :: coupling
        integer :: k, field_status
        h = cmplx(0.0_dp, 0.0_dp, dp)
        ierr = RAW_BAD_INPUT
        if (.not. self%ready()) return
        if (size(h) /= self%field%n_modes) return
        if (.not. valid_point(x, t, b0, velocity, moment, charge)) return
        call evaluate_potentials(self%field, x, a, db, potential, field_status)
        ierr = RAW_FIELD_ERROR
        if (field_status /= PERTFIELD_OK) return
        if (.not. finite_potentials(a, db, potential)) return
        coupling = charge
        if (self%field%units == UNITS_GAUSSIAN) coupling = charge/C_GAUSSIAN
        do k = 1, size(h)
            phase = exp(cmplx(0.0_dp, real(self%field%ntor(k), dp)*x(2) &
                - self%omega(k)*t, dp))
            h(k) = (moment*sum(b0*db(:, k))/norm2(b0) + charge*potential(k) &
                - coupling*sum(a(:, k)*velocity))*phase
        end do
        if (all(ieee_is_finite(real(h, dp)))) then
            if (all(ieee_is_finite(aimag(h)))) ierr = RAW_OK
        end if
    end subroutine evaluate

    logical function valid_point(x, t, b0, velocity, moment, charge)
        real(dp), intent(in) :: x(3), t, b0(3), velocity(3), moment, charge

        valid_point = .false.
        if (.not. all(ieee_is_finite(x))) return
        if (.not. all(ieee_is_finite(b0))) return
        if (.not. all(ieee_is_finite(velocity))) return
        if (.not. ieee_is_finite(t)) return
        if (.not. ieee_is_finite(moment)) return
        if (.not. ieee_is_finite(charge)) return
        if (x(1) <= 0.0_dp) return
        if (norm2(b0) <= 0.0_dp) return
        if (moment < 0.0_dp) return
        valid_point = .true.
    end function valid_point

    subroutine evaluate_potentials(field, x, a, db, potential, ierr)
        type(perturbation_field_t), intent(in) :: field
        real(dp), intent(in) :: x(3)
        complex(dp), intent(out) :: a(:, :), db(:, :), potential(:)
        integer, intent(out) :: ierr

        potential = cmplx(0.0_dp, 0.0_dp, dp)
        if (field%has_potential) then
            call field%eval_modes(x(1), x(3), a, db, ierr, dPhi=potential)
        else
            call field%eval_modes(x(1), x(3), a, db, ierr)
        end if
    end subroutine evaluate_potentials

    logical function finite_potentials(a, db, potential)
        complex(dp), intent(in) :: a(:, :), db(:, :), potential(:)

        finite_potentials = .false.
        if (.not. all(ieee_is_finite(real(a, dp)))) return
        if (.not. all(ieee_is_finite(aimag(a)))) return
        if (.not. all(ieee_is_finite(real(db, dp)))) return
        if (.not. all(ieee_is_finite(aimag(db)))) return
        if (.not. all(ieee_is_finite(real(potential, dp)))) return
        if (.not. all(ieee_is_finite(aimag(potential)))) return
        finite_potentials = .true.
    end function finite_potentials

    subroutine raw_harmonics(source, time, x, velocity, b0, moment, charge, &
            mb, omega_b, omega_phi, hmn, ierr)
        type(raw_source_t), intent(in) :: source
        real(dp), intent(in) :: time(:), x(:, :), velocity(:, :), b0(:, :)
        real(dp), intent(in) :: moment, charge, omega_b, omega_phi
        integer, intent(in) :: mb(:)
        complex(dp), intent(out) :: hmn(:, :)
        integer, intent(out) :: ierr
        type(raw_projection_t) :: projection

        ierr = RAW_BAD_INPUT
        hmn = cmplx(0.0_dp, 0.0_dp, dp)
        if (.not. source%ready()) return
        if (size(hmn, 2) /= source%field%n_modes) return
        if (.not. valid_samples(time, x, velocity, b0, mb, hmn)) return
        if (.not. ieee_is_finite(moment)) return
        if (.not. ieee_is_finite(charge)) return
        if (.not. ieee_is_finite(omega_b)) return
        if (.not. ieee_is_finite(omega_phi)) return
        projection%moment = moment
        projection%charge = charge
        projection%omega_b = omega_b
        projection%omega_phi = omega_phi
        call accumulate_harmonics(source, time, x, velocity, b0, mb, &
            projection, hmn, ierr)
    end subroutine raw_harmonics

    subroutine accumulate_harmonics(source, time, x, velocity, b0, mb, &
            projection, hmn, ierr)
        type(raw_source_t), intent(in) :: source
        real(dp), intent(in) :: time(:), x(:, :), velocity(:, :), b0(:, :)
        integer, intent(in) :: mb(:)
        type(raw_projection_t), intent(in) :: projection
        complex(dp), intent(out) :: hmn(:, :)
        integer, intent(out) :: ierr
        complex(dp) :: previous(size(mb), source%field%n_modes)
        complex(dp) :: current(size(mb), source%field%n_modes)
        integer :: k

        hmn = cmplx(0.0_dp, 0.0_dp, dp)
        call projected_point(source, raw_sample_t(time(1), x(:, 1), velocity(:, 1), &
            b0(:, 1)), mb, projection, previous, ierr)
        if (ierr /= RAW_OK) return
        do k = 2, size(time)
            call projected_point(source, raw_sample_t(time(k), x(:, k), &
                velocity(:, k), b0(:, k)), mb, projection, current, ierr)
            if (ierr /= RAW_OK) return
            hmn = hmn + 0.5_dp*(time(k) - time(k - 1))*(previous + current)
            previous = current
        end do
        hmn = hmn/(time(size(time)) - time(1))
    end subroutine accumulate_harmonics

    subroutine projected_point(source, sample, mb, projection, projected, ierr)
        type(raw_source_t), intent(in) :: source
        type(raw_sample_t), intent(in) :: sample
        integer, intent(in) :: mb(:)
        type(raw_projection_t), intent(in) :: projection
        complex(dp), intent(out) :: projected(:, :)
        integer, intent(out) :: ierr
        complex(dp) :: physical_h(source%field%n_modes)
        real(dp) :: detuning
        integer :: ih, mode

        call source%evaluate(sample%x, sample%time, sample%b0, sample%velocity, &
            projection%moment, projection%charge, physical_h, ierr)
        if (ierr /= RAW_OK) return
        do mode = 1, source%field%n_modes
            do ih = 1, size(mb)
                detuning = mb(ih)*projection%omega_b &
                    + source%field%ntor(mode)*projection%omega_phi - source%omega(mode)
                projected(ih, mode) = physical_h(mode) &
                    *exp(cmplx(0.0_dp, -detuning*sample%time, dp))
            end do
        end do
    end subroutine projected_point

    logical function valid_samples(time, x, velocity, b0, mb, hmn)
        real(dp), intent(in) :: time(:), x(:, :), velocity(:, :), b0(:, :)
        integer, intent(in) :: mb(:)
        complex(dp), intent(in) :: hmn(:, :)

        valid_samples = .false.
        if (size(time) < 2) return
        if (size(x, 1) /= 3 .or. size(x, 2) /= size(time)) return
        if (size(velocity, 1) /= 3 .or. size(velocity, 2) /= size(time)) return
        if (size(b0, 1) /= 3 .or. size(b0, 2) /= size(time)) return
        if (size(hmn, 1) /= size(mb)) return
        if (.not. all(ieee_is_finite(time))) return
        if (any(time(2:) <= time(:size(time) - 1))) return
        valid_samples = .true.
    end function valid_samples

end module neort_raw_drive
