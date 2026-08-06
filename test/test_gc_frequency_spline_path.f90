program test_gc_frequency_spline_path
    !! End-to-end gate for the opt-in production frequency seam.  Independent
    !! checks use the historical zero-width bounce solver only for the lambda=0
    !! period and exact thin-limit velocity scalings for the drift terms.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use do_magfie_mod, only: inp_swi, read_boozer_file, set_s, &
        init_magfie_at_s
    use driftorbit, only: FREQUENCY_MODEL_GC_THIN, frequency_model, &
        magdrift, magdrift_passing, sign_vpar, etatp, etadt
    use neort_freq, only: init_canon_freq_trapped_spline, &
        init_canon_freq_passing_spline, Om_th, Om_ph, Om_tB
    use neort_gc_frequency_splines, only: gc_spline_q
    use neort_magfie, only: init_flux_surface_average
    use neort_orbit, only: bounce_time
    use neort_profiles, only: vth, Om_tE
    use util, only: pi, qe, mu, qi, mi
    use util_for_test, only: pass_test

    implicit none

    real(dp), parameter :: surface = 0.25_dp
    real(dp), parameter :: reference_velocity = 1.0e8_dp
    real(dp), parameter :: electric_frequency = 4.0e3_dp
    real(dp), parameter :: speed_ratio = 1.37_dp
    character(len=1024) :: eqdsk_file
    real(dp) :: eta_trapped, eta_passing, legacy_omega
    real(dp) :: omega_b_1, omega_b_2, omega_phi_1, omega_phi_2
    real(dp) :: omega_b_dv, omega_b_deta, omega_phi_dv, omega_phi_deta
    real(dp) :: omega_b_dv_2, omega_b_deta_2
    real(dp) :: omega_phi_dv_2, omega_phi_deta_2
    real(dp) :: omega_bmag_1, omega_bmag_2, omega_bmag_dv
    real(dp) :: omega_bmag_deta, electric_1, electric_2, q_gc

    call get_environment_variable('EQDSK_FILE', eqdsk_file)
    if (len_trim(eqdsk_file) == 0) then
        error stop 'GEQDSK test fixture is not configured'
    end if

    inp_swi = 11
    call read_boozer_file(trim(eqdsk_file))
    call set_s(surface)
    call init_magfie_at_s()
    call init_flux_surface_average(surface)
    qi = qe
    mi = 2.014_dp*mu
    vth = reference_velocity
    Om_tE = electric_frequency
    sign_vpar = 1.0_dp
    magdrift = .true.
    magdrift_passing = 1
    frequency_model = FREQUENCY_MODEL_GC_THIN

    call init_canon_freq_trapped_spline()
    call init_canon_freq_passing_spline()

    eta_trapped = etatp + 0.4_dp*(etadt - etatp)
    call evaluate_all(reference_velocity, eta_trapped, omega_b_1, &
        omega_phi_1, omega_bmag_1, omega_b_dv, omega_b_deta, &
        omega_phi_dv, omega_phi_deta, omega_bmag_dv, omega_bmag_deta)
    legacy_omega = 2.0_dp*pi/bounce_time(reference_velocity, eta_trapped)
    call require_relative('trapped lambda=0 bounce', omega_b_1, &
        legacy_omega, 4.0e-4_dp)
    electric_1 = omega_phi_1 - omega_bmag_1
    call require_relative('trapped electric frequency', electric_1, &
        electric_frequency, 6.0e-3_dp)

    call evaluate_all(speed_ratio*reference_velocity, eta_trapped, omega_b_2, &
        omega_phi_2, omega_bmag_2, omega_b_dv_2, omega_b_deta_2, &
        omega_phi_dv_2, omega_phi_deta_2, omega_bmag_dv, omega_bmag_deta)
    electric_2 = omega_phi_2 - omega_bmag_2
    call require_relative('bounce v scaling', omega_b_2, &
        speed_ratio*omega_b_1, 2.0e-12_dp)
    call require_relative('magnetic v-squared scaling', omega_bmag_2, &
        speed_ratio**2*omega_bmag_1, 2.0e-12_dp)
    call require_relative('electric v-zero scaling', electric_2, electric_1, &
        2.0e-12_dp)
    call require_relative('bounce velocity derivative', omega_b_dv, &
        omega_b_1/reference_velocity, 2.0e-12_dp)

    eta_passing = 0.4_dp*etatp
    call evaluate_all(reference_velocity, eta_passing, omega_b_1, &
        omega_phi_1, omega_bmag_1, omega_b_dv, omega_b_deta, &
        omega_phi_dv, omega_phi_deta, omega_bmag_dv, omega_bmag_deta)
    q_gc = gc_spline_q()
    electric_1 = omega_phi_1 - omega_bmag_1 - q_gc*omega_b_1
    call require_relative('passing electric frequency', electric_1, &
        electric_frequency, 6.0e-3_dp)
    legacy_omega = 2.0_dp*pi/bounce_time(reference_velocity, eta_passing)
    call require_relative('passing lambda=0 transit', omega_b_1, &
        legacy_omega, 4.0e-4_dp)

    call pass_test

contains

    subroutine evaluate_all(velocity, eta, omega_b, omega_phi, omega_bmag, &
            domega_b_dv, domega_b_deta, domega_phi_dv, domega_phi_deta, &
            domega_bmag_dv, domega_bmag_deta)
        real(dp), intent(in) :: velocity, eta
        real(dp), intent(out) :: omega_b, omega_phi, omega_bmag
        real(dp), intent(out) :: domega_b_dv, domega_b_deta
        real(dp), intent(out) :: domega_phi_dv, domega_phi_deta
        real(dp), intent(out) :: domega_bmag_dv, domega_bmag_deta

        call Om_th(velocity, eta, omega_b, domega_b_dv, domega_b_deta)
        call Om_ph(velocity, eta, omega_phi, domega_phi_dv, domega_phi_deta)
        call Om_tB(velocity, eta, omega_bmag, domega_bmag_dv, &
            domega_bmag_deta)
    end subroutine evaluate_all

    subroutine require_relative(label, value, reference, tolerance)
        character(*), intent(in) :: label
        real(dp), intent(in) :: value, reference, tolerance

        if (abs(value - reference) > tolerance*max(abs(reference), 1.0_dp)) then
            write(*, *) trim(label)//' value/reference:', value, reference
            error stop 'GC production frequency scaling test failed'
        end if
    end subroutine require_relative

end program test_gc_frequency_spline_path
