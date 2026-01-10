program test_neort_lib
    use iso_fortran_env, only: dp => real64
    use logger, only: set_log_level
    use neort_lib
    use omp_lib, only: omp_get_thread_num, omp_set_num_threads

    implicit none

    call set_log_level(-1)

    call test_single_s()
    call test_parallel_loop()
    call test_spline_update()

    print *, "============================================"
    print *, "ALL NEORT_LIB TESTS PASSED"
    print *, "============================================"

contains

    subroutine test_single_s()
        type(transport_data_t) :: result
        real(dp), parameter :: TEST_S = 0.5_dp

        print *, "Test: single s computation..."

        call neort_init("driftorbit.in", "in_file")
        call neort_prepare_splines("plasma.in", "profile.in")
        call neort_compute_at_s(TEST_S, result)

        if (abs(result%torque%s - TEST_S) > 1.0e-12_dp) then
            print *, "ERROR: s value mismatch, expected", TEST_S, "got", result%torque%s
            error stop "test_single_s failed"
        end if

        if (.not. result%torque%has_torque) then
            error stop "test_single_s: expected torque computation"
        end if

        if (result%torque%dVds <= 0.0_dp) then
            error stop "test_single_s: dVds should be positive"
        end if

        print *, "  OK"
    end subroutine test_single_s

    subroutine test_parallel_loop()
        integer, parameter :: N_S = 5
        integer, parameter :: N_THREADS = 4
        type(transport_data_t) :: results(N_S)
        real(dp) :: s_array(N_S)
        real(dp), parameter :: TOL = 1.0e-10_dp
        integer :: i
        logical :: all_positive

        print *, "Test: parallel loop over s values..."

        call omp_set_num_threads(N_THREADS)

        do i = 1, N_S
            s_array(i) = 0.2_dp + 0.1_dp * real(i, dp)
        end do

        call neort_init("driftorbit.in", "in_file")
        call neort_prepare_splines("plasma.in", "profile.in")

        !$omp parallel do schedule(dynamic)
        do i = 1, N_S
            call neort_compute_at_s(s_array(i), results(i))
        end do
        !$omp end parallel do

        all_positive = .true.
        do i = 1, N_S
            if (abs(results(i)%torque%s - s_array(i)) > TOL) then
                print *, "ERROR: s mismatch at index", i
                error stop "test_parallel_loop: s value mismatch"
            end if
            if (results(i)%torque%dVds <= 0.0_dp) then
                all_positive = .false.
            end if
        end do

        if (.not. all_positive) then
            error stop "test_parallel_loop: not all dVds positive"
        end if

        print *, "  OK"
    end subroutine test_parallel_loop

    subroutine test_spline_update()
        use neort_profiles, only: read_plasma_input
        use util, only: readdata

        type(transport_data_t) :: result1, result2
        real(dp), parameter :: TEST_S = 0.5_dp
        integer :: nplasma
        real(dp) :: am1, am2, Z1, Z2
        real(dp), allocatable :: plasma_data(:, :), profile_data(:, :)

        print *, "Test: spline update (KAMEL use case)..."

        call neort_init("driftorbit.in", "in_file")

        call read_plasma_input("plasma.in", nplasma, am1, am2, Z1, Z2, plasma_data)
        call readdata("profile.in", 2, profile_data)

        call neort_prepare_splines(nplasma, am1, am2, Z1, Z2, plasma_data, profile_data)
        call neort_compute_at_s(TEST_S, result1)

        call neort_prepare_splines(nplasma, am1, am2, Z1, Z2, plasma_data, profile_data)
        call neort_compute_at_s(TEST_S, result2)

        if (abs(result1%torque%s - result2%torque%s) > 1.0e-12_dp) then
            error stop "test_spline_update: results should match"
        end if

        if (abs(result1%torque%Tt - result2%torque%Tt) > 1.0e-10_dp) then
            error stop "test_spline_update: torque should match"
        end if

        print *, "  OK"
    end subroutine test_spline_update

end program test_neort_lib
