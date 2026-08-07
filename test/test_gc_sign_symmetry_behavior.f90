program test_gc_sign_symmetry_behavior
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_full_fow_sign_symmetry_symbolic, only: &
        evaluate_neort_full_fow_sign_symmetry_ledger
    use util_for_test, only: pass_test
    implicit none

    real(dp), parameter :: a_real = 1.25_dp, a_imag = -0.75_dp
    real(dp), parameter :: m_mode = 2.0_dp, n_mode = 3.0_dp
    real(dp), parameter :: theta = 0.37_dp, phi = -0.22_dp
    real(dp), parameter :: global_phase = 0.41_dp
    real(dp), parameter :: dpsi_dx = -1.4_dp, f_prime = 0.8_dp
    real(dp), parameter :: cdot = -2.3_dp, signed_measure = -1.1_dp
    real(dp), parameter :: h = 4.2_dp, charge = -2.2_dp
    real(dp), parameter :: electrostatic_potential = 0.6_dp
    real(dp), parameter :: temperature = 1.7_dp
    real(dp), parameter :: a1 = 0.4_dp, a2 = -0.25_dp, gauge_c = 0.9_dp
    real(dp), parameter :: omega_b = 1.7_dp, omega_phi = -0.6_dp
    real(dp), parameter :: torque_phi = -2.4_dp, rotation_phi = 0.35_dp
    real(dp) :: ledger(14)

    call evaluate_neort_full_fow_sign_symmetry_ledger(a_real, a_imag, &
        m_mode, n_mode, theta, phi, global_phase, dpsi_dx, f_prime, cdot, &
        signed_measure, h, charge, electrostatic_potential, temperature, &
        a1, a2, gauge_c, omega_b, omega_phi, torque_phi, rotation_phi, &
        ledger(1), ledger(2), ledger(3), ledger(4), ledger(5), ledger(6), &
        ledger(7), ledger(8), ledger(9), ledger(10), ledger(11), ledger(12), &
        ledger(13), ledger(14))

    call check_field_fixture(ledger)
    call check_section_fixture(ledger)
    call check_resonance_fixture(ledger)
    call check_potential_fixture(ledger)

    write (*, '(A)') 'test_gc_sign_symmetry_behavior OK'
    call pass_test

contains

    subroutine check_field_fixture(values)
        real(dp), intent(in) :: values(:)
        complex(dp) :: amplitude, rotated_amplitude
        real(dp) :: original, paired, fixed_conjugate
        real(dp) :: amplitude_squared, signed_squared, phase_squared

        amplitude = cmplx(a_real, a_imag, dp)
        rotated_amplitude = amplitude*cmplx(cos(global_phase), &
            sin(global_phase), dp)
        original = real_field(m_mode, n_mode, a_real, a_imag, theta, phi)
        paired = real_field(-m_mode, -n_mode, real(conjg(amplitude), dp), &
            aimag(conjg(amplitude)), theta, phi)
        fixed_conjugate = real_field(m_mode, n_mode, real(conjg(amplitude), dp), &
            aimag(conjg(amplitude)), theta, phi)
        amplitude_squared = amplitude_square(a_real, a_imag)
        signed_squared = amplitude_square(-a_real, -a_imag)
        phase_squared = amplitude_square(real(rotated_amplitude, dp), &
            aimag(rotated_amplitude))

        call require_close('paired real field', paired, original)
        call require_close('A to minus A modulus', signed_squared, &
            amplitude_squared)
        call require_close('global phase modulus', phase_squared, &
            amplitude_squared)
        call require_nonzero('fixed-mode conjugation field change', &
            fixed_conjugate - original)
        call require_close('ledger A sign residual', values(1), 0.0_dp)
        call require_close('ledger phase residual', values(2), 0.0_dp)
        call require_close('ledger paired-field residual', values(3), 0.0_dp)
        call require_close('ledger fixed-conjugation difference', values(4), &
            fixed_conjugate - original)
        call require_close('ledger conjugate modulus residual', values(5), &
            0.0_dp)
        call require_close('ledger gauge-force residual', values(8), 0.0_dp)
        call require_close('ledger gauge-outer residual', values(9), 0.0_dp)
    end subroutine check_field_fixture

    subroutine check_section_fixture(values)
        real(dp), intent(in) :: values(:)
        real(dp) :: root_weight, reversed_root_weight
        real(dp) :: section_signed, reversed_section_signed
        real(dp) :: signed_crossing, reversed_signed_crossing
        real(dp) :: positive_measure, reversed_positive_measure

        root_weight = abs(dpsi_dx/f_prime)
        reversed_root_weight = abs((-dpsi_dx)/(-f_prime))
        section_signed = dpsi_dx
        reversed_section_signed = -dpsi_dx
        signed_crossing = signed_measure*cdot
        reversed_signed_crossing = signed_measure*(-cdot)
        positive_measure = abs(signed_measure)*abs(cdot)
        reversed_positive_measure = abs(signed_measure)*abs(-cdot)

        call require_close('section signed orientation', &
            reversed_section_signed, -section_signed)
        call require_close('section positive root weight', &
            reversed_root_weight, root_weight)
        call require_close('cut positive measure', reversed_positive_measure, &
            positive_measure)
        call require_close('cut signed orientation', reversed_signed_crossing, &
            -signed_crossing)
        call require_nonzero('cut orientation changes signed quantity', &
            reversed_signed_crossing - signed_crossing)
        call require_close('ledger section residual', values(6), 0.0_dp)
        call require_close('ledger cut measure residual', values(7), 0.0_dp)
    end subroutine check_section_fixture

    subroutine check_resonance_fixture(values)
        real(dp), intent(in) :: values(:)
        real(dp) :: resonance, reversed_resonance
        real(dp) :: power, reversed_power, reversed_torque

        resonance = resonance_value(m_mode, n_mode, omega_b, omega_phi)
        reversed_resonance = resonance_value(m_mode, -n_mode, omega_b, &
            -omega_phi)
        power = torque_power(torque_phi, rotation_phi)
        reversed_torque = -torque_phi
        reversed_power = torque_power(reversed_torque, -rotation_phi)

        call require_close('toroidal resonance relabelling', &
            reversed_resonance, resonance)
        call require_close('toroidal torque reversal', reversed_torque, &
            -torque_phi)
        call require_close('toroidal torque power', reversed_power, power)
        call require_nonzero('toroidal torque component changes sign', &
            reversed_torque - torque_phi)
        call require_close('ledger resonance residual', values(10), 0.0_dp)
        call require_close('ledger power residual', values(11), 0.0_dp)
        call require_close('ledger torque difference', values(12), &
            reversed_torque - torque_phi)
    end subroutine check_resonance_fixture

    subroutine check_potential_fixture(values)
        real(dp), intent(in) :: values(:)
        real(dp) :: base_energy, charge_reversed, potential_reversed

        base_energy = charge*electrostatic_potential
        charge_reversed = (-charge)*electrostatic_potential
        potential_reversed = charge*(-electrostatic_potential)
        call require_nonzero('q reversal is not a symmetry', &
            charge_reversed - base_energy)
        call require_nonzero('Phi reversal is not a symmetry', &
            potential_reversed - base_energy)
        call require_close('ledger q reversal difference', values(13), &
            charge_reversed - base_energy)
        call require_close('ledger Phi reversal difference', values(14), &
            potential_reversed - base_energy)
    end subroutine check_potential_fixture

    pure real(dp) function real_field(m_value, n_value, real_part, imag_part, &
            theta_value, phi_value)
        real(dp), intent(in) :: m_value, n_value, real_part, imag_part
        real(dp), intent(in) :: theta_value, phi_value
        real(dp) :: phase

        phase = m_value*theta_value + n_value*phi_value
        real_field = real_part*cos(phase) - imag_part*sin(phase)
    end function real_field

    pure real(dp) function amplitude_square(real_part, imag_part)
        real(dp), intent(in) :: real_part, imag_part

        amplitude_square = real_part**2 + imag_part**2
    end function amplitude_square

    pure real(dp) function resonance_value(m_value, n_value, omega_b_value, &
            omega_phi_value)
        real(dp), intent(in) :: m_value, n_value, omega_b_value
        real(dp), intent(in) :: omega_phi_value

        resonance_value = m_value*omega_b_value + n_value*omega_phi_value
    end function resonance_value

    pure real(dp) function torque_power(torque_value, rotation_value)
        real(dp), intent(in) :: torque_value, rotation_value

        torque_power = torque_value*rotation_value
    end function torque_power

    subroutine require_close(label, actual, expected)
        character(*), intent(in) :: label
        real(dp), intent(in) :: actual, expected
        real(dp) :: scale

        scale = max(1.0_dp, max(abs(actual), abs(expected)))
        if (abs(actual - expected) > 2.0e-12_dp*scale) then
            write (*, '(A,2(1X,ES24.16))') trim(label), actual, expected
            error stop 'sign-symmetry behavioral oracle failed'
        end if
    end subroutine require_close

    subroutine require_nonzero(label, value)
        character(*), intent(in) :: label
        real(dp), intent(in) :: value

        if (abs(value) <= 2.0e-12_dp) then
            write (*, '(A,1X,ES24.16)') trim(label), value
            error stop 'sign-symmetry distinction was lost'
        end if
    end subroutine require_nonzero

end program test_gc_sign_symmetry_behavior
