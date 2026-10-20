program test_mixture_thread_invariance
    ! The velocity nodes of all resonance components of all species and
    ! harmonics are one OpenMP work list. A thread computing a node must hold
    ! the state the serial loops would hold, whatever it computed before (e.g.
    ! a counter-passing node of another species), so the result must not
    ! depend on the thread count at all. The argument selects the finite-width
    ! Krook ('krook') or the collisionless ('collisionless') velocity integral.
    use iso_fortran_env, only: dp => real64, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use neort_lib, only: config_t, mixture_t, species_state_t, mixture_result_t, &
        neort_init, neort_prepare_mixture, neort_compute_mixture_at_s
    use omp_lib, only: omp_set_num_threads, omp_set_dynamic
    implicit none
    type(config_t) :: config
    type(mixture_t) :: mixture
    type(mixture_result_t) :: serial, threaded
    integer, parameter :: thread_counts(4) = [2, 4, 8, 16]
    integer :: n, nthreads, failures
    character(len=32) :: mode

    call get_command_argument(1, mode)
    if (mode /= 'krook' .and. mode /= 'collisionless') &
        error stop 'usage: test_mixture_thread_invariance.x krook|collisionless'
    call setup_config(config, mode == 'krook')
    call neort_init(config, 'in_file')
    call setup_mixture(mixture)
    call neort_prepare_mixture(mixture)
    call omp_set_dynamic(.false.)
    call omp_set_num_threads(1)
    call neort_compute_mixture_at_s(config%s, serial)
    if (size(serial%species(1)%transport%harmonics) < 3) &
        error stop 'test needs several harmonics per species'
    if (all(serial%torque_total == 0.0_dp)) error stop 'test needs resonant torque'
    print *, trim(mode), ' torque_total', serial%torque_total
    failures = 0
    do n = 1, size(thread_counts)
        nthreads = thread_counts(n)
        call omp_set_num_threads(nthreads)
        call neort_compute_mixture_at_s(config%s, threaded)
        call compare(serial, threaded, nthreads, failures)
    end do
    if (failures > 0) error stop 'mixture transport depends on the thread count'
    print *, 'Mixture transport (', trim(mode), ') is independent of the thread count'
contains
    subroutine setup_config(value, krook)
        type(config_t), intent(out) :: value
        logical, intent(in) :: krook
        value%s = 0.09_dp
        value%epsmn = 1.0e-3_dp
        value%m0 = 0
        value%mph = 3
        value%comptorque = .true.
        value%nopassing = .false.
        value%magdrift = .true.
        value%inp_swi = 8
        value%vsteps = 4
        value%mth_max_abs = 1
        value%vmax_over_vth = 2.0_dp
        value%log_level = -1
        value%krook_collisions = krook
        value%krook_nu = 1.0e3_dp
        value%krook_eta_steps = 4
        value%krook_rtol = 1.0e-4_dp
    end subroutine setup_config

    subroutine setup_mixture(value)
        type(mixture_t), intent(out) :: value
        value%reference_s = 0.09_dp
        value%omega_e = -700.0_dp
        value%coulomb_log = 20.0_dp
        allocate(value%species(3))
        value%species(1) = species_state_t('D', 3.3435837724e-27_dp, 1.0_dp, &
            5.0e15_dp, 2000.0_dp, 0.0_dp, 0.0_dp)
        value%species(2) = species_state_t('T', 5.0073567446e-27_dp, 1.0_dp, &
            5.0e15_dp, 2000.0_dp, 0.0_dp, 0.0_dp)
        value%species(3) = species_state_t('e', 9.1093837139e-31_dp, -1.0_dp, &
            1.0e16_dp, 2000.0_dp, 0.0_dp, 0.0_dp)
    end subroutine setup_mixture

    subroutine compare(expected, actual, nthreads, failures)
        type(mixture_result_t), intent(in) :: expected, actual
        integer, intent(in) :: nthreads
        integer, intent(inout) :: failures
        integer :: k, j
        do k = 1, size(expected%species)
            associate (e => expected%species(k)%transport, &
                       a => actual%species(k)%transport)
                call check('torque co', k, nthreads, &
                    e%torque%Tco, a%torque%Tco, failures)
                call check('torque ctr', k, nthreads, &
                    e%torque%Tctr, a%torque%Tctr, failures)
                call check('torque t', k, nthreads, &
                    e%torque%Tt, a%torque%Tt, failures)
                call check_values('D co', k, nthreads, &
                    e%summary%Dco, a%summary%Dco, failures)
                call check_values('D ctr', k, nthreads, &
                    e%summary%Dctr, a%summary%Dctr, failures)
                call check_values('D t', k, nthreads, &
                    e%summary%Dt, a%summary%Dt, failures)
                do j = 1, size(e%harmonics)
                    call check_values('harmonic D co', k, nthreads, &
                        e%harmonics(j)%Dresco, a%harmonics(j)%Dresco, failures)
                    call check_values('harmonic D ctr', k, nthreads, &
                        e%harmonics(j)%Dresctr, a%harmonics(j)%Dresctr, failures)
                    call check_values('harmonic D t', k, nthreads, &
                        e%harmonics(j)%Drest, a%harmonics(j)%Drest, failures)
                    call check('harmonic co', k, nthreads, &
                        e%harmonics(j)%Tresco, a%harmonics(j)%Tresco, failures)
                    call check('harmonic ctr', k, nthreads, &
                        e%harmonics(j)%Tresctr, a%harmonics(j)%Tresctr, failures)
                    call check('harmonic t', k, nthreads, &
                        e%harmonics(j)%Trest, a%harmonics(j)%Trest, failures)
                end do
            end associate
        end do
        call check_values('total', 0, nthreads, &
            expected%torque_total, actual%torque_total, failures)
    end subroutine compare

    subroutine check_values(label, species, nthreads, reference, value, failures)
        character(len=*), intent(in) :: label
        integer, intent(in) :: species, nthreads
        real(dp), intent(in) :: reference(:), value(:)
        integer, intent(inout) :: failures
        integer :: i
        if (size(reference) /= size(value)) error stop 'comparison shape mismatch'
        do i = 1, size(reference)
            call check(label, species, nthreads, reference(i), value(i), failures)
        end do
    end subroutine check_values

    subroutine check(label, species, nthreads, reference, value, failures)
        character(len=*), intent(in) :: label
        integer, intent(in) :: species, nthreads
        real(dp), intent(in) :: reference, value
        integer, intent(inout) :: failures
        if (.not. ieee_is_finite(reference) .or. &
                .not. ieee_is_finite(value) .or. &
                transfer(reference, 0_int64) /= transfer(value, 0_int64)) then
            failures = failures + 1
            print *, label, ' species', species, ' threads', nthreads, reference, value
        end if
    end subroutine check
end program test_mixture_thread_invariance
