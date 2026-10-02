module neort_raw_orbit
    !! Axisymmetric Littlejohn guiding-center orbit in physical (R,phi,Z).
    !! The identical physical velocity feeds the raw perturbation Hamiltonian.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
    use fortnum_ode_vode, only: vode_state_t, vode_init, vode_integrate_to
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    use neort_raw_drive, only: RAW_OK, RAW_BAD_INPUT, C_GAUSSIAN
    use neo_perturbation_field, only: UNITS_GAUSSIAN, UNITS_SI
    implicit none
    private

    integer, parameter, public :: RAW_ORBIT_ERROR = 3

    type, public :: raw_background_t
        real(dp) :: b(3) = 0.0_dp, grad_b(3) = 0.0_dp
        real(dp) :: curl_bhat(3) = 0.0_dp, grad_phi(3) = 0.0_dp
        real(dp) :: potential = 0.0_dp, a(3) = 0.0_dp
        logical :: has_vector_potential = .false.
    end type raw_background_t

    type, public :: raw_species_t
        real(dp) :: mass = 0.0_dp, charge = 0.0_dp
        character(len=16) :: units = ''
    end type raw_species_t

    type, public :: raw_trajectory_t
        real(dp), allocatable :: time(:), state(:, :), velocity(:, :), b(:, :)
        real(dp), allocatable :: energy(:), pphi(:), poloidal_oneform(:)
        real(dp), allocatable :: liouville_density(:)
        logical :: has_actions = .false.
    end type raw_trajectory_t

    abstract interface
        subroutine raw_background_i(x, background, ierr, data)
            import :: dp, raw_background_t
            real(dp), intent(in) :: x(3)
            type(raw_background_t), intent(out) :: background
            integer, intent(out) :: ierr
            class(*), intent(in), optional :: data
        end subroutine raw_background_i
    end interface

    type :: raw_gc_context_t
        procedure(raw_background_i), pointer, nopass :: background => null()
        class(*), pointer :: data => null()
        type(raw_species_t) :: species
        real(dp) :: moment
        integer, pointer :: error => null()
    end type raw_gc_context_t

    public :: raw_gc_point, raw_gc_trajectory, raw_gc_actions, raw_background_i

contains

    subroutine raw_gc_point(background, species, moment, vpar, velocity, &
            acceleration, bstar_parallel, ierr)
        type(raw_background_t), intent(in) :: background
        type(raw_species_t), intent(in) :: species
        real(dp), intent(in) :: moment, vpar
        real(dp), intent(out) :: velocity(3), acceleration, bstar_parallel
        integer, intent(out) :: ierr
        real(dp) :: bhat(3), bstar(3), force(3), drift(3), factor, bmag
        ierr = RAW_BAD_INPUT
        velocity = 0.0_dp
        acceleration = 0.0_dp
        bstar_parallel = 0.0_dp
        if (.not. valid_gc_input(background, species, moment, vpar)) return
        factor = 1.0_dp/gc_coupling(species)
        bmag = norm2(background%b)
        bhat = background%b/bmag
        bstar = background%b + species%mass*factor*vpar*background%curl_bhat
        bstar_parallel = sum(bhat*bstar)
        if (abs(bstar_parallel) < 0.5_dp*bmag) return
        force = moment*background%grad_b + species%charge*background%grad_phi
        drift(1) = bhat(2)*force(3) - bhat(3)*force(2)
        drift(2) = bhat(3)*force(1) - bhat(1)*force(3)
        drift(3) = bhat(1)*force(2) - bhat(2)*force(1)
        velocity = (vpar*bstar + factor*drift)/bstar_parallel
        acceleration = -sum(bstar*force)/(species%mass*bstar_parallel)
        if (.not. all(ieee_is_finite(velocity))) return
        if (.not. ieee_is_finite(acceleration)) return
        ierr = RAW_OK
    end subroutine raw_gc_point

    logical function valid_gc_input(background, species, moment, vpar)
        type(raw_background_t), intent(in) :: background
        type(raw_species_t), intent(in) :: species
        real(dp), intent(in) :: moment, vpar

        valid_gc_input = .false.
        if (.not. all(ieee_is_finite(background%b))) return
        if (.not. all(ieee_is_finite(background%grad_b))) return
        if (.not. all(ieee_is_finite(background%curl_bhat))) return
        if (.not. all(ieee_is_finite(background%grad_phi))) return
        if (.not. ieee_is_finite(background%potential)) return
        if (.not. ieee_is_finite(species%mass)) return
        if (.not. ieee_is_finite(species%charge)) return
        if (.not. ieee_is_finite(moment)) return
        if (.not. ieee_is_finite(vpar)) return
        if (species%mass <= 0.0_dp .or. species%charge == 0.0_dp) return
        if (moment < 0.0_dp .or. norm2(background%b) <= 0.0_dp) return
        if (species%units /= UNITS_GAUSSIAN) then
            if (species%units /= UNITS_SI) return
        end if
        if (background%has_vector_potential) then
            if (.not. all(ieee_is_finite(background%a))) return
        end if
        valid_gc_input = .true.
    end function valid_gc_input

    pure real(dp) function gc_coupling(species) result(coupling)
        type(raw_species_t), intent(in) :: species

        coupling = species%charge
        if (species%units == UNITS_GAUSSIAN) coupling = coupling/C_GAUSSIAN
    end function gc_coupling

    subroutine raw_gc_trajectory(background, species, moment, time, initial, &
            relative_tolerance, absolute_tolerance, orbit, ierr, data)
        procedure(raw_background_i) :: background
        type(raw_species_t), intent(in) :: species
        real(dp), intent(in) :: moment, time(:), initial(4)
        real(dp), intent(in) :: relative_tolerance, absolute_tolerance(4)
        type(raw_trajectory_t), intent(out) :: orbit
        integer, intent(out), target :: ierr
        class(*), intent(in), optional, target :: data
        type(raw_gc_context_t) :: context
        integer, target :: callback_error

        ierr = RAW_BAD_INPUT
        if (.not. valid_trajectory_input(time, initial, relative_tolerance, &
            absolute_tolerance)) return
        context%background => background
        context%species = species
        context%moment = moment
        callback_error = RAW_OK
        context%error => callback_error
        if (present(data)) context%data => data
        call initialize_trajectory(time, initial, orbit)
        ierr = RAW_OK
        call integrate_trajectory(context, relative_tolerance, absolute_tolerance, &
            orbit, ierr)
        if (ierr /= RAW_OK) orbit%has_actions = .false.
    end subroutine raw_gc_trajectory

    logical function valid_trajectory_input(time, initial, relative_tolerance, &
            absolute_tolerance)
        real(dp), intent(in) :: time(:), initial(4)
        real(dp), intent(in) :: relative_tolerance, absolute_tolerance(4)

        valid_trajectory_input = .false.
        if (size(time) < 2) return
        if (.not. all(ieee_is_finite(time))) return
        if (.not. all(ieee_is_finite(initial))) return
        if (.not. ieee_is_finite(relative_tolerance)) return
        if (.not. all(ieee_is_finite(absolute_tolerance))) return
        if (any(time(2:) <= time(:size(time) - 1))) return
        if (initial(1) <= 0.0_dp .or. relative_tolerance <= 0.0_dp) return
        if (any(absolute_tolerance <= 0.0_dp)) return
        valid_trajectory_input = .true.
    end function valid_trajectory_input

    subroutine initialize_trajectory(time, initial, orbit)
        real(dp), intent(in) :: time(:), initial(4)
        type(raw_trajectory_t), intent(out) :: orbit
        integer :: n

        n = size(time)
        allocate (orbit%time(n), source=time)
        allocate (orbit%state(4, n), orbit%velocity(3, n), orbit%b(3, n))
        allocate (orbit%energy(n), orbit%pphi(n), orbit%poloidal_oneform(n))
        allocate (orbit%liouville_density(n))
        orbit%state = 0.0_dp
        orbit%state(:, 1) = initial
        orbit%has_actions = .true.
    end subroutine initialize_trajectory

    subroutine integrate_trajectory(context, rtol, atol, orbit, ierr)
        type(raw_gc_context_t), intent(in) :: context
        real(dp), intent(in) :: rtol, atol(4)
        type(raw_trajectory_t), intent(inout) :: orbit
        integer, intent(out) :: ierr
        type(vode_state_t) :: integrator
        type(fortnum_status_t) :: status
        real(dp), allocatable :: yend(:)
        integer :: k

        ierr = RAW_OK
        call vode_init(integrator, 4, orbit%time(1), orbit%state(:, 1))
        call sample_trajectory(context, orbit, 1, ierr)
        if (ierr /= RAW_OK) return
        do k = 2, size(orbit%time)
            call vode_integrate_to(gc_rhs, integrator, orbit%time(k), rtol, atol, &
                yend, status, ctx=context)
            if (context%error /= RAW_OK) then
                ierr = context%error
                return
            end if
            if (status%code /= FORTNUM_OK) then
                ierr = RAW_ORBIT_ERROR
                return
            end if
            orbit%state(:, k) = yend
            call sample_trajectory(context, orbit, k, ierr)
            if (ierr /= RAW_OK) return
        end do
    end subroutine integrate_trajectory

    subroutine gc_rhs(t, state, derivative, data)
        real(dp), intent(in) :: t, state(:)
        real(dp), intent(out) :: derivative(:)
        class(*), intent(in), optional :: data
        type(raw_background_t) :: background
        real(dp) :: velocity(3), acceleration, bstar_parallel, x(3)
        integer :: ierr

        derivative = 0.0_dp
        if (.not. present(data)) return
        select type (data)
        type is (raw_gc_context_t)
            if (data%error /= RAW_OK) return
            x = state(1:3)
            call evaluate_background(data, x, background, ierr)
            if (ierr == RAW_OK) call raw_gc_point(background, data%species, &
                data%moment, state(4), velocity, acceleration, bstar_parallel, ierr)
            if (ierr /= RAW_OK) then
                data%error = ierr
                return
            end if
            derivative(1) = velocity(1)
            derivative(2) = velocity(2)/state(1)
            derivative(3) = velocity(3)
            derivative(4) = acceleration
        end select
        if (.not. ieee_is_finite(t)) derivative = 0.0_dp
    end subroutine gc_rhs

    subroutine evaluate_background(context, x, background, ierr)
        type(raw_gc_context_t), intent(in) :: context
        real(dp), intent(in) :: x(3)
        type(raw_background_t), intent(out) :: background
        integer, intent(out) :: ierr

        ierr = RAW_BAD_INPUT
        if (.not. all(ieee_is_finite(x))) return
        if (x(1) <= 0.0_dp) return
        if (associated(context%data)) then
            call context%background(x, background, ierr, context%data)
        else
            call context%background(x, background, ierr)
        end if
    end subroutine evaluate_background

    subroutine sample_trajectory(context, orbit, k, ierr)
        type(raw_gc_context_t), intent(in) :: context
        type(raw_trajectory_t), intent(inout) :: orbit
        integer, intent(in) :: k
        integer, intent(out) :: ierr
        type(raw_background_t) :: background
        real(dp) :: acceleration, bstar_parallel

        call evaluate_background(context, orbit%state(1:3, k), background, ierr)
        if (ierr /= RAW_OK) return
        call raw_gc_point(background, context%species, context%moment, &
            orbit%state(4, k), orbit%velocity(:, k), acceleration, bstar_parallel, ierr)
        if (ierr /= RAW_OK) return
        orbit%b(:, k) = background%b
        call sample_invariants(context, background, orbit, k, bstar_parallel)
    end subroutine sample_trajectory

    subroutine sample_invariants(context, background, orbit, k, bstar_parallel)
        type(raw_gc_context_t), intent(in) :: context
        type(raw_background_t), intent(in) :: background
        type(raw_trajectory_t), intent(inout) :: orbit
        integer, intent(in) :: k
        real(dp), intent(in) :: bstar_parallel
        real(dp) :: canonical(3), coupling, bmag, vp, nan

        bmag = norm2(background%b)
        vp = orbit%state(4, k)
        orbit%energy(k) = 0.5_dp*context%species%mass*vp**2 &
            + context%moment*bmag + context%species%charge*background%potential
        orbit%liouville_density(k) = context%species%mass**2 &
            *abs(bstar_parallel)*orbit%state(1, k)
        nan = ieee_value(0.0_dp, ieee_quiet_nan)
        orbit%pphi(k) = nan
        orbit%poloidal_oneform(k) = nan
        if (.not. background%has_vector_potential) then
            orbit%has_actions = .false.
            return
        end if
        coupling = gc_coupling(context%species)
        canonical = coupling*background%a + context%species%mass*vp*background%b/bmag
        orbit%pphi(k) = orbit%state(1, k)*canonical(2)
        orbit%poloidal_oneform(k) = canonical(1)*orbit%velocity(1, k) &
            + canonical(3)*orbit%velocity(3, k)
    end subroutine sample_invariants

    subroutine raw_gc_actions(orbit, closure_atol, cycle_orientation, action_b, &
            pphi, omega_b, &
            omega_phi, closure_error, ierr)
        !! Caller supplies one primitive meridional period, not multiple turns.
        !! Trapped cycles use +1; passing cycles use their poloidal winding sign.
        type(raw_trajectory_t), intent(in) :: orbit
        real(dp), intent(in) :: closure_atol(3)
        integer, intent(in) :: cycle_orientation
        real(dp), intent(out) :: action_b, pphi, omega_b, omega_phi, closure_error(3)
        integer, intent(out) :: ierr
        real(dp), parameter :: two_pi = 2.0_dp*acos(-1.0_dp)
        real(dp) :: period
        integer :: n

        ierr = RAW_BAD_INPUT
        if (cycle_orientation /= 1 .and. cycle_orientation /= -1) return
        if (.not. valid_action_samples(orbit, closure_atol)) return
        n = size(orbit%time)
        closure_error = abs(orbit%state([1, 3, 4], n) - orbit%state([1, 3, 4], 1))
        if (any(closure_error > closure_atol)) return
        period = orbit%time(n) - orbit%time(1)
        action_b = bounce_action(orbit, cycle_orientation)
        pphi = orbit%pphi(1)
        omega_b = cycle_orientation*two_pi/period
        omega_phi = (orbit%state(2, n) - orbit%state(2, 1))/period
        ierr = RAW_OK
    end subroutine raw_gc_actions

    real(dp) function bounce_action(orbit, cycle_orientation) result(action)
        type(raw_trajectory_t), intent(in) :: orbit
        integer, intent(in) :: cycle_orientation
        integer :: k

        action = 0.0_dp
        do k = 2, size(orbit%time)
            action = action + 0.5_dp*(orbit%time(k) - orbit%time(k - 1)) &
                *(orbit%poloidal_oneform(k) + orbit%poloidal_oneform(k - 1))
        end do
        action = cycle_orientation*action/(2.0_dp*acos(-1.0_dp))
    end function bounce_action

    logical function valid_action_samples(orbit, closure_atol)
        type(raw_trajectory_t), intent(in) :: orbit
        real(dp), intent(in) :: closure_atol(3)

        valid_action_samples = .false.
        if (.not. orbit%has_actions) return
        if (.not. allocated(orbit%time)) return
        if (.not. allocated(orbit%state)) return
        if (.not. allocated(orbit%pphi)) return
        if (.not. allocated(orbit%poloidal_oneform)) return
        if (size(orbit%time) < 2) return
        if (size(orbit%state, 1) /= 4) return
        if (size(orbit%state, 2) /= size(orbit%time)) return
        if (size(orbit%pphi) /= size(orbit%time)) return
        if (size(orbit%poloidal_oneform) /= size(orbit%time)) return
        if (.not. all(ieee_is_finite(orbit%time))) return
        if (.not. all(ieee_is_finite(orbit%state))) return
        if (.not. all(ieee_is_finite(orbit%pphi))) return
        if (.not. all(ieee_is_finite(orbit%poloidal_oneform))) return
        if (.not. all(ieee_is_finite(closure_atol))) return
        if (any(closure_atol <= 0.0_dp)) return
        if (any(orbit%time(2:) <= orbit%time(:size(orbit%time) - 1))) return
        valid_action_samples = .true.
    end function valid_action_samples

end module neort_raw_orbit
