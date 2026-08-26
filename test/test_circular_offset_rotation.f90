program test_circular_offset_rotation
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_datatypes, only: transport_data_t
    use neort_lib, only: config_t, neort_compute_no_splines, neort_init
    use neort_magfie, only: B0
    use neort_profiles, only: M_t, Om_tE, Ti1, dTi1ds, dni1ds, init_plasma_at_s, &
        ni1, prepare_plasma_splines, read_plasma_input, vth
    use do_magfie_mod, only: Bphcov, R0, psi_pr, q, set_s, sign_theta
    use omp_lib, only: omp_get_num_threads, omp_set_dynamic, omp_set_num_threads
    use util, only: c, ev, qi

    implicit none

    integer, parameter :: N_SURFACES = 4
    integer, parameter :: N_CASES = 3
    integer, parameter :: N_BISECT = 20
    integer, parameter :: N_THREADS = 4
    real(dp), parameter :: K_NC = 1.17_dp
    real(dp), parameter :: ROOT_BRACKET = 1.0_dp
    real(dp), parameter :: DELTA_OMEGA = 2.0e3_dp
    real(dp), parameter :: SURFACES(N_SURFACES) = &
        [0.045_dp, 0.075_dp, 0.105_dp, 0.135_dp]
    real(dp), parameter :: CASE_SIGN(N_CASES) = [-1.0_dp, 0.0_dp, 1.0_dp]

    type :: offset_state_t
        real(dp) :: D11
        real(dp) :: D12
        real(dp) :: geometry_factor
        real(dp) :: k_na
        real(dp) :: Om_tE
        real(dp) :: Vphi
        real(dp) :: Vphi_in
        real(dp) :: residual
        real(dp) :: torque
    end type offset_state_t

    real(dp) :: root_mt(N_SURFACES)
    type(offset_state_t) :: states(N_SURFACES, N_CASES)
    integer :: surface_index, case_index

    call check_electric_sign_oracle()
    call initialize_library()

    do surface_index = 1, N_SURFACES
        call find_offset_root(SURFACES(surface_index), root_mt(surface_index))
    end do

    call omp_set_dynamic(.false.)
    call omp_set_num_threads(N_THREADS)
    call verify_thread_count()
    !$omp parallel do collapse(2) default(shared) private(surface_index, case_index) &
    !$omp schedule(dynamic, 1)
    do case_index = 1, N_CASES
        do surface_index = 1, N_SURFACES
            call evaluate_delta(SURFACES(surface_index), root_mt(surface_index), &
                CASE_SIGN(case_index)*DELTA_OMEGA, states(surface_index, case_index))
        end do
    end do
    !$omp end parallel do

    do surface_index = 1, N_SURFACES
        call check_restoring_response(SURFACES(surface_index), states(surface_index, :))
    end do

    call write_showcase(states(2, :))
    print *, "Circular Kasilov offset response passes with k =", K_NC

contains

    subroutine initialize_library()
        type(config_t) :: config
        real(dp), allocatable :: plasma(:, :)
        integer :: nplasma
        real(dp) :: am1, am2, Z1, Z2

        config%s = SURFACES(1)
        config%M_t = 0.0_dp
        config%qs = 1.0_dp
        config%ms = 2.014_dp
        config%vth = 4.0e7_dp
        config%epsmn = 0.001_dp
        config%m0 = 0
        config%mph = 3
        config%comptorque = .true.
        config%supban = .false.
        config%magdrift = .true.
        config%magdrift_passing = 1
        config%nopassing = .false.
        config%noshear = .true.
        config%pertfile = .false.
        config%nonlin = .false.
        config%bfac = 1.0_dp
        config%efac = 1.0_dp
        config%inp_swi = 8
        config%inp_swi_pert = -1
        config%vsteps = 32
        config%mth_max_abs = -1
        config%vmax_over_vth = 4.0_dp
        config%log_level = -1
        config%output_format = "text"

        call neort_init(config, "in_file")
        call read_plasma_input("plasma.in", nplasma, am1, am2, Z1, Z2, plasma)
        call prepare_plasma_splines(nplasma, am1, am2, Z1, Z2, plasma)
        deallocate (plasma)
    end subroutine initialize_library

    subroutine check_electric_sign_oracle()
        real(dp), parameter :: phi_prime = 3.0_dp
        real(dp), parameter :: Bphi = 4.0_dp
        real(dp), parameter :: Btheta = 2.0_dp
        real(dp), parameter :: charge = 5.0_dp
        real(dp), parameter :: temperature = 7.0_dp
        real(dp) :: B2, electric(3), magnetic(3), exb(3)
        real(dp) :: omega_from_vector, omega_from_potential
        real(dp) :: force_from_omega, force_from_radial_field

        ! Use a right-handed Cartesian realization of (s, theta, phi), then
        ! form the field-line-label velocity V^phi - q V^theta used by the
        ! code, whose stored coordinate ordering is (s, phi, theta).
        electric = [-phi_prime, 0.0_dp, 0.0_dp]
        magnetic = [0.0_dp, Btheta, Bphi]
        B2 = dot_product(magnetic, magnetic)
        exb = c*cross_product(electric, magnetic)/B2

        omega_from_vector = exb(3) - (Bphi/Btheta)*exb(2)
        omega_from_potential = -c*phi_prime/Btheta
        call assert_close("E cross B frequency sign", omega_from_potential, &
            omega_from_vector, 1.0e-14_dp)

        force_from_omega = -charge*Btheta*omega_from_vector/(c*temperature)
        force_from_radial_field = charge*phi_prime/temperature
        call assert_close("A1 electric-force sign", force_from_radial_field, &
            force_from_omega, 1.0e-14_dp)
    end subroutine check_electric_sign_oracle

    subroutine verify_thread_count()
        integer :: observed_threads

        observed_threads = 0
        !$omp parallel default(shared)
        !$omp single
        observed_threads = omp_get_num_threads()
        !$omp end single
        !$omp end parallel
        if (observed_threads /= N_THREADS) then
            error stop "circular offset test did not create four OpenMP threads"
        end if
    end subroutine verify_thread_count

    pure function cross_product(a, b) result(cross)
        real(dp), intent(in) :: a(3), b(3)
        real(dp) :: cross(3)

        cross = [a(2)*b(3) - a(3)*b(2), &
            a(3)*b(1) - a(1)*b(3), &
            a(1)*b(2) - a(2)*b(1)]
    end function cross_product

    subroutine find_offset_root(surface, root)
        real(dp), intent(in) :: surface
        real(dp), intent(out) :: root
        type(offset_state_t) :: lower_state, middle_state, upper_state
        real(dp) :: lower, middle, upper
        integer :: iteration

        lower = -ROOT_BRACKET
        upper = ROOT_BRACKET
        call evaluate(surface, lower, lower_state)
        call evaluate(surface, upper, upper_state)
        if (lower_state%residual*upper_state%residual >= 0.0_dp) then
            print *, "surface, lower residual, upper residual:", surface, &
                lower_state%residual, upper_state%residual
            error stop "circular offset root is not bracketed"
        end if

        do iteration = 1, N_BISECT
            middle = 0.5_dp*(lower + upper)
            call evaluate(surface, middle, middle_state)
            if (lower_state%residual*middle_state%residual <= 0.0_dp) then
                upper = middle
                upper_state = middle_state
            else
                lower = middle
                lower_state = middle_state
            end if
        end do
        root = 0.5_dp*(lower + upper)
    end subroutine find_offset_root

    subroutine evaluate_delta(surface, root, delta_omega, state)
        real(dp), intent(in) :: surface, root, delta_omega
        type(offset_state_t), intent(out) :: state
        real(dp) :: trial_mt

        call set_s(surface)
        call init_plasma_at_s()
        trial_mt = root + delta_omega*R0/vth
        call evaluate(surface, trial_mt, state)
    end subroutine evaluate_delta

    subroutine evaluate(surface, trial_mt, state)
        real(dp), intent(in) :: surface, trial_mt
        type(offset_state_t), intent(out) :: state
        type(transport_data_t) :: result

        call set_s(surface)
        call init_plasma_at_s()
        M_t = trial_mt
        call neort_compute_no_splines(result)
        call derive_offset_state(result, state)
    end subroutine evaluate

    subroutine derive_offset_state(result, state)
        type(transport_data_t), intent(in) :: result
        type(offset_state_t), intent(out) :: state
        real(dp) :: chi_prime, gradient_scale, pressure_gradient_over_density

        state%D11 = sum([result%summary%Dco(1), result%summary%Dctr(1), &
            result%summary%Dt(1)])
        state%D12 = sum([result%summary%Dco(2), result%summary%Dctr(2), &
            result%summary%Dt(2)])
        if (state%D11 <= 0.0_dp) error stop "D11 must be positive at the offset"
        if (q == 0.0_dp .or. psi_pr == 0.0_dp) then
            error stop "invalid circular magnetic coordinates"
        end if
        if (B0 == 0.0_dp .or. R0 == 0.0_dp .or. ni1 == 0.0_dp) then
            error stop "invalid circular profile or geometry"
        end if

        chi_prime = sign_theta*psi_pr/q
        gradient_scale = c*ev/(qi*chi_prime)
        pressure_gradient_over_density = Ti1*dni1ds/ni1 + dTi1ds
        state%geometry_factor = (Bphcov/(B0*R0))**2
        state%k_na = state%D12/state%D11 - 2.5_dp + &
            state%geometry_factor*K_NC
        state%Om_tE = Om_tE
        state%Vphi_in = gradient_scale*state%k_na*dTi1ds
        state%Vphi = Om_tE - gradient_scale*pressure_gradient_over_density + &
            gradient_scale*state%geometry_factor*K_NC*dTi1ds
        state%residual = state%Vphi - state%Vphi_in
        state%torque = result%torque%Tco + result%torque%Tctr + result%torque%Tt
    end subroutine derive_offset_state

    subroutine check_restoring_response(surface, surface_states)
        real(dp), intent(in) :: surface
        type(offset_state_t), intent(in) :: surface_states(:)
        real(dp) :: side_torque_scale

        if (surface_states(1)%residual >= 0.0_dp) then
            error stop "below-offset circular case is not below the offset"
        end if
        if (surface_states(3)%residual <= 0.0_dp) then
            error stop "above-offset circular case is not above the offset"
        end if
        if (surface_states(1)%torque*surface_states(1)%residual >= 0.0_dp) then
            error stop "below-offset torque is not restoring"
        end if
        if (surface_states(3)%torque*surface_states(3)%residual >= 0.0_dp) then
            error stop "above-offset torque is not restoring"
        end if

        if (abs(surface_states(2)%residual) > 5.0e-4_dp*DELTA_OMEGA) then
            print *, "surface, root residual:", surface, surface_states(2)%residual
            error stop "circular offset fixed point did not converge"
        end if
        side_torque_scale = min(abs(surface_states(1)%torque), &
            abs(surface_states(3)%torque))
        if (abs(surface_states(2)%torque) > 1.0e-3_dp*side_torque_scale) then
            print *, "surface, root torque:", surface, surface_states(2)%torque
            error stop "torque does not vanish at the circular offset"
        end if
    end subroutine check_restoring_response

    subroutine write_showcase(surface_states)
        type(offset_state_t), intent(in) :: surface_states(:)
        character(len=1024) :: output_path
        integer :: k, length, status, unit

        call get_environment_variable("NEORT_OFFSET_SHOWCASE", output_path, &
            length=length, status=status)
        if (status /= 0 .or. length == 0) return

        open (newunit=unit, file=output_path(:length), status="replace", action="write")
        write (unit, '(A)') &
            "s,delta_omega_s-1,vphi_s-1,vphi_in_s-1,om_te_s-1,"// &
            "torque_native,residual_s-1,k_na,D11,D12"
        do k = 1, size(surface_states)
            write (unit, '(*(ES24.16,:,","))') SURFACES(2), &
                CASE_SIGN(k)*DELTA_OMEGA, surface_states(k)%Vphi, &
                surface_states(k)%Vphi_in, surface_states(k)%Om_tE, &
                surface_states(k)%torque, surface_states(k)%residual, &
                surface_states(k)%k_na, surface_states(k)%D11, surface_states(k)%D12
        end do
        close (unit)
    end subroutine write_showcase

    subroutine assert_close(label, expected, actual, tolerance)
        character(len=*), intent(in) :: label
        real(dp), intent(in) :: expected, actual, tolerance

        if (abs(actual - expected) > tolerance*max(1.0_dp, abs(expected))) then
            print *, trim(label), "expected", expected, "actual", actual
            error stop "manufactured electric-sign oracle failed"
        end if
    end subroutine assert_close

end program test_circular_offset_rotation
