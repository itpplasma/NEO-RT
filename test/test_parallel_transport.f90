program test_parallel_transport
    use iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use neort_datatypes, only: transport_harmonic_t
    use neort_lib, only: config_t, transport_data_t, neort_compute_at_s, &
        neort_init, neort_prepare_splines
    use neort_profiles, only: A1, A2, Om_tE, Ti1, qi
    use do_magfie_mod, only: psi_pr, q, sign_theta
    use util, only: c, ev
    use omp_lib, only: omp_get_num_threads, omp_set_dynamic, &
        omp_set_num_threads

    implicit none

    integer, parameter :: N_SURFACES = 8
    integer, parameter :: N_THREADS = 4
    integer, parameter :: N_REPEATS = 4
    integer, parameter :: MAX_REPORTED_FAILURES = 20
    real(dp), parameter :: REL_TOL = 1.0e-10_dp
    real(dp), parameter :: ABS_TOL = 1.0e-12_dp
    real(dp), parameter :: SURFACES(N_SURFACES) = [ &
        0.030_dp, 0.045_dp, 0.060_dp, 0.075_dp, &
        0.090_dp, 0.105_dp, 0.120_dp, 0.135_dp]

    type(transport_data_t) :: serial_results(N_SURFACES)
    type(transport_data_t) :: parallel_results(N_SURFACES)
    integer :: surface_index, repeat_index, failures

    failures = 0

    call initialize_library()

    do surface_index = 1, N_SURFACES
        call neort_compute_at_s( &
            SURFACES(surface_index), serial_results(surface_index))
        call check_offset_response(serial_results(surface_index), surface_index)
    end do

    call omp_set_dynamic(.false.)
    call omp_set_num_threads(N_THREADS)
    call verify_thread_count()

    do repeat_index = 1, N_REPEATS
        !$omp parallel do default(shared) private(surface_index) &
        !$omp schedule(dynamic, 1)
        do surface_index = 1, N_SURFACES
            call neort_compute_at_s( &
                SURFACES(surface_index), parallel_results(surface_index))
        end do
        !$omp end parallel do

        do surface_index = 1, N_SURFACES
            call compare_result( &
                serial_results(surface_index), &
                parallel_results(surface_index), &
                surface_index, repeat_index, failures)
        end do
    end do

    if (failures > 0) then
        print *, "Parallel transport mismatches:", failures
        error stop "parallel transport differs from serial reference"
    end if

    print *, "Parallel transport matches serial reference"

contains

    subroutine initialize_library()
        type(config_t) :: config

        config%s = 0.5_dp
        config%M_t = 0.036_dp
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
        call neort_prepare_splines("plasma.in", "profile.in")
    end subroutine initialize_library

    subroutine verify_thread_count()
        integer :: observed_threads

        observed_threads = 0

        !$omp parallel default(shared)
        !$omp single
        observed_threads = omp_get_num_threads()
        !$omp end single
        !$omp end parallel

        if (observed_threads /= N_THREADS) then
            print *, "Expected threads:", N_THREADS
            print *, "Observed threads:", observed_threads
            error stop "OpenMP did not create the requested team"
        end if
    end subroutine verify_thread_count

    subroutine compare_result(expected, actual, surface, round, failures_)
        type(transport_data_t), intent(in) :: expected, actual
        integer, intent(in) :: surface, round
        integer, intent(inout) :: failures_

        call compare_summary(expected, actual, surface, round, failures_)
        call compare_torque(expected, actual, surface, round, failures_)
        call compare_harmonics(expected, actual, surface, round, failures_)
    end subroutine compare_result

    subroutine compare_summary(expected, actual, surface, round, failures_)
        type(transport_data_t), intent(in) :: expected, actual
        integer, intent(in) :: surface, round
        integer, intent(inout) :: failures_

        call check_real("summary M_t", expected%summary%M_t, &
            actual%summary%M_t, surface, round, failures_)
        call check_array("summary Dco", expected%summary%Dco, &
            actual%summary%Dco, surface, round, failures_)
        call check_array("summary Dctr", expected%summary%Dctr, &
            actual%summary%Dctr, surface, round, failures_)
        call check_array("summary Dt", expected%summary%Dt, &
            actual%summary%Dt, surface, round, failures_)
    end subroutine compare_summary

    subroutine check_offset_response(result, surface)
        type(transport_data_t), intent(in) :: result
        integer, intent(in) :: surface

        real(dp) :: D11, D12, k_na, delta_omega, total_torque
        real(dp) :: dA1dOm, expected_offset

        D11 = sum([result%summary%Dco(1), result%summary%Dctr(1), &
            result%summary%Dt(1)])
        D12 = sum([result%summary%Dco(2), result%summary%Dctr(2), &
            result%summary%Dt(2)])
        if (D11 == 0.0_dp) error stop "D11 vanished on circular torque test"
        if (Ti1 == 0.0_dp .or. q == 0.0_dp) then
            error stop "invalid circular profile for offset test"
        end if
        k_na = D12 / D11 - 2.5_dp
        dA1dOm = -qi / (Ti1 * ev) * sign_theta * psi_pr / (q * c)
        expected_offset = Om_tE + (-(k_na + 2.5_dp) * A2 - A1) / dA1dOm
        delta_omega = Om_tE - expected_offset
        total_torque = torque_total(result)
        if (abs(delta_omega) > 1.0e-10_dp * &
            max(1.0_dp, abs(Om_tE))) then
            if (abs(total_torque) > 1.0e-12_dp) then
                if (total_torque * delta_omega >= 0.0_dp) then
                    print *, "torque", total_torque, "delta omega", delta_omega
                    error stop "torque does not point toward the reported offset"
                end if
            end if
        end if
    end subroutine check_offset_response

    subroutine compare_torque(expected, actual, surface, round, failures_)
        type(transport_data_t), intent(in) :: expected, actual
        integer, intent(in) :: surface, round
        integer, intent(inout) :: failures_

        call check_logical("torque enabled", expected%torque%has_torque, &
            actual%torque%has_torque, surface, round, failures_)
        call check_real("torque s", expected%torque%s, actual%torque%s, &
            surface, round, failures_)
        call check_real("torque dVds", expected%torque%dVds, &
            actual%torque%dVds, surface, round, failures_)
        call check_real("torque M_t", expected%torque%M_t, &
            actual%torque%M_t, surface, round, failures_)
        call check_real("torque Tco", expected%torque%Tco, &
            actual%torque%Tco, surface, round, failures_)
        call check_real("torque Tctr", expected%torque%Tctr, &
            actual%torque%Tctr, surface, round, failures_)
        call check_real("torque Tt", expected%torque%Tt, &
            actual%torque%Tt, surface, round, failures_)
        call check_real("torque total", torque_total(expected), &
            torque_total(actual), surface, round, failures_)
    end subroutine compare_torque

    subroutine compare_harmonics(expected, actual, surface, round, failures_)
        type(transport_data_t), intent(in) :: expected, actual
        integer, intent(in) :: surface, round
        integer, intent(inout) :: failures_
        integer :: harmonic_index

        call check_integer("harmonic count", size(expected%harmonics), &
            size(actual%harmonics), surface, round, failures_)

        if (size(expected%harmonics) /= size(actual%harmonics)) return

        do harmonic_index = 1, size(expected%harmonics)
            call compare_harmonic( &
                expected%harmonics(harmonic_index), &
                actual%harmonics(harmonic_index), harmonic_index, &
                surface, round, failures_)
        end do
    end subroutine compare_harmonics

    subroutine compare_harmonic(expected, actual, index, surface, round, failures_)
        type(transport_harmonic_t), intent(in) :: expected, actual
        integer, intent(in) :: index, surface, round
        integer, intent(inout) :: failures_

        call check_integer(harmonic_label(index, "mth"), expected%mth, &
            actual%mth, surface, round, failures_)
        call check_array(harmonic_label(index, "Dresco"), expected%Dresco, &
            actual%Dresco, surface, round, failures_)
        call check_array(harmonic_label(index, "Dresctr"), expected%Dresctr, &
            actual%Dresctr, surface, round, failures_)
        call check_array(harmonic_label(index, "Drest"), expected%Drest, &
            actual%Drest, surface, round, failures_)
        call check_real(harmonic_label(index, "Tresco"), expected%Tresco, &
            actual%Tresco, surface, round, failures_)
        call check_real(harmonic_label(index, "Tresctr"), expected%Tresctr, &
            actual%Tresctr, surface, round, failures_)
        call check_real(harmonic_label(index, "Trest"), expected%Trest, &
            actual%Trest, surface, round, failures_)
        call compare_velocity_bounds(expected, actual, index, surface, &
            round, failures_)
    end subroutine compare_harmonic

    subroutine compare_velocity_bounds(expected, actual, index, surface, &
            round, failures_)
        type(transport_harmonic_t), intent(in) :: expected, actual
        integer, intent(in) :: index, surface, round
        integer, intent(inout) :: failures_

        call check_real(harmonic_label(index, "vminp"), &
            expected%vminp_over_vth, actual%vminp_over_vth, &
            surface, round, failures_)
        call check_real(harmonic_label(index, "vmaxp"), &
            expected%vmaxp_over_vth, actual%vmaxp_over_vth, &
            surface, round, failures_)
        call check_real(harmonic_label(index, "vmint"), &
            expected%vmint_over_vth, actual%vmint_over_vth, &
            surface, round, failures_)
        call check_real(harmonic_label(index, "vmaxt"), &
            expected%vmaxt_over_vth, actual%vmaxt_over_vth, &
            surface, round, failures_)
    end subroutine compare_velocity_bounds

    subroutine check_array(name, expected, actual, surface, round, failures_)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: expected(:), actual(:)
        integer, intent(in) :: surface, round
        integer, intent(inout) :: failures_
        integer :: component
        character(len=128) :: component_name

        do component = 1, size(expected)
            write (component_name, '(A,"(",I0,")")') trim(name), component
            call check_real(component_name, expected(component), &
                actual(component), surface, round, failures_)
        end do
    end subroutine check_array

    subroutine check_real(name, expected, actual, surface, round, failures_)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: expected, actual
        integer, intent(in) :: surface, round
        integer, intent(inout) :: failures_
        real(dp) :: difference, limit

        difference = abs(actual - expected)
        limit = ABS_TOL + REL_TOL*max(abs(expected), abs(actual))

        if (.not. ieee_is_finite(expected) .or. &
            .not. ieee_is_finite(actual) .or. difference > limit) then
            if (failures_ < MAX_REPORTED_FAILURES) then
                print *, trim(name), "surface", surface, "round", round
                print *, "expected", expected, "actual", actual
            end if
            failures_ = failures_ + 1
        end if
    end subroutine check_real

    subroutine check_integer(name, expected, actual, surface, round, failures_)
        character(len=*), intent(in) :: name
        integer, intent(in) :: expected, actual, surface, round
        integer, intent(inout) :: failures_

        if (expected /= actual) then
            if (failures_ < MAX_REPORTED_FAILURES) then
                print *, trim(name), "surface", surface, "round", round
                print *, "expected", expected, "actual", actual
            end if
            failures_ = failures_ + 1
        end if
    end subroutine check_integer

    subroutine check_logical(name, expected, actual, surface, round, failures_)
        character(len=*), intent(in) :: name
        logical, intent(in) :: expected, actual
        integer, intent(in) :: surface, round
        integer, intent(inout) :: failures_

        if (expected .neqv. actual) then
            if (failures_ < MAX_REPORTED_FAILURES) then
                print *, trim(name), "surface", surface, "round", round
                print *, "expected", expected, "actual", actual
            end if
            failures_ = failures_ + 1
        end if
    end subroutine check_logical

    pure function torque_total(result) result(total)
        type(transport_data_t), intent(in) :: result
        real(dp) :: total

        total = result%torque%Tco + result%torque%Tctr + result%torque%Tt
    end function torque_total

    function harmonic_label(index, field) result(label)
        integer, intent(in) :: index
        character(len=*), intent(in) :: field
        character(len=96) :: label

        write (label, '("harmonic(",I0,") ",A)') index, trim(field)
    end function harmonic_label

end program test_parallel_transport
