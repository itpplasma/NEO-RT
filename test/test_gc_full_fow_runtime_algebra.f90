program test_gc_full_fow_runtime_algebra
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_full_fow_cycle_average_symbolic, only: &
        evaluate_neort_full_fow_cycle_average
    use neort_full_fow_cycle_frequency_symbolic, only: &
        evaluate_neort_full_fow_cycle_frequency
    use neort_full_fow_mode_mapping_symbolic, only: &
        evaluate_neort_full_fow_mode_mapping
    use neort_full_fow_resonance_scalar_symbolic, only: &
        evaluate_neort_full_fow_resonance_scalar
    use neort_full_fow_torque_assembly_symbolic, only: &
        evaluate_neort_full_fow_torque_assembly
    use util_for_test, only: pass_test
    implicit none

    call test_cycle_frequency_against_numerical_quotients
    call test_winding_and_eqdsk_conventions
    call test_signed_resonance_zero_set
    call test_cycle_average_against_midpoint_quadrature
    call test_signed_torque_assembly
    write (*, '(A)') 'test_gc_full_fow_runtime_algebra OK'
    call pass_test

contains

    subroutine test_cycle_frequency_against_numerical_quotients
        real(dp), parameter :: period = 2.4_dp, delta_phi = -0.73_dp
        real(dp), parameter :: q_fieldline = 1.7_dp, dperiod = 0.31_dp
        real(dp), parameter :: ddelta_phi = -0.22_dp, dq_fieldline = 0.08_dp
        real(dp), parameter :: hs(4) = [1.0e-2_dp, 3.0e-3_dp, 1.0e-3_dp, &
            3.0e-4_dp]
        real(dp) :: values(6), plus(3), minus(3), fd(3), h
        real(dp) :: error_history(3, 4)
        integer :: i

        call evaluate_neort_full_fow_cycle_frequency(period, -1.0_dp, delta_phi, &
            q_fieldline, dperiod, ddelta_phi, dq_fieldline, values(1), values(2), &
            values(3), values(4), values(5), values(6))
        call require_close('signed winding frequency', values(1), &
            -2.0_dp*acos(-1.0_dp)/period, 1.0e-13_dp)
        call require_close('signed toroidal frequency', values(2), &
            delta_phi/period, 1.0e-13_dp)
        call require_close('signed precession', values(3), &
            delta_phi/period-q_fieldline*values(1), 1.0e-13_dp)

        do i = 1, size(hs)
            h = hs(i)
            plus = direct_frequencies(period+h*dperiod, delta_phi+h*ddelta_phi, &
                q_fieldline+h*dq_fieldline, -1.0_dp)
            minus = direct_frequencies(period-h*dperiod, delta_phi-h*ddelta_phi, &
                q_fieldline-h*dq_fieldline, -1.0_dp)
            fd = (plus-minus)/(2.0_dp*h)
            error_history(:, i) = abs(fd-values(4:6))
        end do
        call require_less('omega_b derivative finest error', &
            error_history(1, 4), 2.0e-9_dp)
        call require_less('omega_phi derivative finest error', &
            error_history(2, 4), 2.0e-9_dp)
        call require_less('omega_prec derivative finest error', &
            error_history(3, 4), 2.0e-9_dp)
        do i = 1, 2
            call require_rate('omega_b centered-FD convergence', &
                error_history(1, i), error_history(1, i+1))
            call require_rate('omega_phi centered-FD convergence', &
                error_history(2, i), error_history(2, i+1))
            call require_rate('omega_prec centered-FD convergence', &
                error_history(3, i), error_history(3, i+1))
        end do
    end subroutine test_cycle_frequency_against_numerical_quotients

    subroutine test_winding_and_eqdsk_conventions
        real(dp) :: positive(6), negative(6)

        call evaluate_neort_full_fow_cycle_frequency(2.5_dp, 1.0_dp, -0.7_dp, &
            1.4_dp, 0.0_dp, 0.0_dp, 0.0_dp, positive(1), positive(2), &
            positive(3), positive(4), positive(5), positive(6))
        call evaluate_neort_full_fow_cycle_frequency(2.5_dp, -1.0_dp, -0.7_dp, &
            1.4_dp, 0.0_dp, 0.0_dp, 0.0_dp, negative(1), negative(2), &
            negative(3), negative(4), negative(5), negative(6))
        call require_nonzero('EQDSK positive winding is positive', positive(1))
        call require_close('negative winding reverses omega_b', negative(1), &
            -positive(1), 1.0e-13_dp)
        call require_close('delta phi remains signed', positive(2), -0.28_dp, &
            1.0e-13_dp)
        call require_close('winding does not alter omega_phi', negative(2), &
            positive(2), 1.0e-13_dp)
    end subroutine test_winding_and_eqdsk_conventions

    subroutine test_signed_resonance_zero_set
        real(dp) :: residual, derivative, reversed, unused_derivative
        real(dp) :: m_positive, n_positive, signed_map, positive_map

        call evaluate_neort_full_fow_resonance_scalar(2.0_dp, -3.0_dp, 3.0_dp, &
            2.0_dp, 0.4_dp, -0.1_dp, residual, derivative)
        call require_close('resonance zero', residual, 0.0_dp, 1.0e-13_dp)
        call require_close('resonance derivative', derivative, 0.5_dp, 1.0e-13_dp)
        call evaluate_neort_full_fow_resonance_scalar(2.0_dp, -3.0_dp, -3.0_dp, &
            -2.0_dp, -0.4_dp, 0.1_dp, reversed, unused_derivative)
        call require_close('frequency reversal preserves zero set', reversed, &
            0.0_dp, 1.0e-13_dp)
        call evaluate_neort_full_fow_resonance_scalar(-2.0_dp, 3.0_dp, 3.0_dp, &
            2.0_dp, 0.4_dp, -0.1_dp, reversed, unused_derivative)
        call require_close('Fourier pair reversal preserves zero set', reversed, &
            0.0_dp, 1.0e-13_dp)
        call evaluate_neort_full_fow_mode_mapping(2.0_dp, -3.0_dp, -1.0_dp, &
            3.0_dp, -2.0_dp, m_positive, n_positive, signed_map, positive_map)
        call require_close('counterpassing mode m relabelling', m_positive, &
            -2.0_dp, 1.0e-13_dp)
        call require_close('mode mapping signed residual', signed_map, &
            positive_map, 1.0e-13_dp)
        call require_close('mapped resonance zero set', positive_map, 0.0_dp, &
            1.0e-13_dp)
        call evaluate_neort_full_fow_resonance_scalar(1.0_dp, 1.0_dp, -1.0_dp, &
            0.2_dp, 0.0_dp, 0.0_dp, reversed, unused_derivative)
        call require_nonzero('resonance scalar remains signed', reversed)
        call evaluate_neort_full_fow_resonance_scalar(2.0_dp, -3.0_dp, -3.0_dp, &
            2.0_dp, 0.0_dp, 0.0_dp, reversed, unused_derivative)
        call require_nonzero('unmapped omega_b reversal moves off-resonance root', &
            reversed)
    end subroutine test_signed_resonance_zero_set

    subroutine test_cycle_average_against_midpoint_quadrature
        integer, parameter :: panels = 16384
        real(dp), parameter :: period = 2.7_dp
        real(dp) :: accum_real, accum_imag, accum_residence, accum_field
        real(dp) :: average_real, average_imag, average_residence, average_field
        real(dp) :: dt, time, expected_real, expected_imag
        real(dp) :: expected_residence, expected_field
        integer :: i

        dt = period/real(panels, dp)
        accum_real = 0.0_dp
        accum_imag = 0.0_dp
        accum_residence = 0.0_dp
        accum_field = 0.0_dp
        do i = 1, panels
            time = (real(i, dp)-0.5_dp)*dt
            accum_real = accum_real + exp(-time)*dt
            accum_imag = accum_imag + cos(1.3_dp*time)*dt
            accum_residence = accum_residence + (1.0_dp+sin(time))*dt
            accum_field = accum_field + (2.0_dp+time*time)*dt
        end do
        call evaluate_neort_full_fow_cycle_average(accum_real, accum_imag, &
            accum_residence, accum_field, period, average_real, average_imag, &
            average_residence, average_field)
        expected_real = (1.0_dp-exp(-period))/period
        expected_imag = sin(1.3_dp*period)/(1.3_dp*period)
        expected_residence = 1.0_dp+(1.0_dp-cos(period))/period
        expected_field = 2.0_dp+period**2/3.0_dp
        call require_close('real cycle average midpoint oracle', average_real, &
            expected_real, 2.0e-7_dp)
        call require_close('imaginary cycle average midpoint oracle', average_imag, &
            expected_imag, 2.0e-8_dp)
        call require_close('residence cycle average midpoint oracle', &
            average_residence, expected_residence, 2.0e-8_dp)
        call require_close('field cycle average midpoint oracle', average_field, &
            expected_field, 2.0e-8_dp)
    end subroutine test_cycle_average_against_midpoint_quadrature

    subroutine test_signed_torque_assembly
        real(dp) :: torque, reversed_force, reversed_n, reversed_outer
        real(dp) :: weight_ratios(3)

        ! Independent manufactured-root oracle: g_phase=x and
        ! F_freq=(n/tau)g_phase.  The phase weight is |n| tau^2 delta(g),
        ! while the frequency weight is tau delta(F); their ratio tends to
        ! n^2 before the emitted torque assembly is exercised.
        weight_ratios(1) = delta_weight_ratio(4096, 0.4_dp)
        weight_ratios(2) = delta_weight_ratio(8192, 0.2_dp)
        weight_ratios(3) = delta_weight_ratio(16384, 0.1_dp)
        call require_less('phase/frequency weight ratio converges', &
            abs(weight_ratios(3)-9.0_dp), 5.0e-2_dp)
        call require_less('phase/frequency ratio improves', &
            abs(weight_ratios(3)-9.0_dp), abs(weight_ratios(1)-9.0_dp))

        call evaluate_neort_full_fow_torque_assembly(3.0_dp, 2.5_dp, -0.8_dp, &
            torque)
        call evaluate_neort_full_fow_torque_assembly(3.0_dp, 2.5_dp, 0.8_dp, &
            reversed_force)
        call evaluate_neort_full_fow_torque_assembly(-3.0_dp, 2.5_dp, -0.8_dp, &
            reversed_n)
        call evaluate_neort_full_fow_torque_assembly(3.0_dp, -2.5_dp, -0.8_dp, &
            reversed_outer)
        call require_close('signed force controls torque sign', torque, -18.0_dp, &
            1.0e-13_dp)
        call require_close('force reversal flips torque', reversed_force, -torque, &
            1.0e-13_dp)
        call require_close('n reversal preserves n squared torque', reversed_n, &
            torque, 1.0e-13_dp)
        call require_close('signed outer factor controls torque sign', &
            reversed_outer, -torque, 1.0e-13_dp)
    end subroutine test_signed_torque_assembly

    real(dp) function delta_weight_ratio(panels, epsilon_value)
        integer, intent(in) :: panels
        real(dp), intent(in) :: epsilon_value
        real(dp), parameter :: n_mode = 3.0_dp, tau = 2.5_dp
        real(dp), parameter :: lower = -4.0_dp, upper = 4.0_dp
        real(dp) :: dx, x, g_phase, f_frequency, phase_delta, frequency_delta
        real(dp) :: phase_integral, frequency_integral
        integer :: i

        dx = (upper-lower)/real(panels, dp)
        phase_integral = 0.0_dp
        frequency_integral = 0.0_dp
        do i = 1, panels
            x = lower+(real(i, dp)-0.5_dp)*dx
            g_phase = x
            f_frequency = (n_mode/tau)*g_phase
            phase_delta = epsilon_value/(acos(-1.0_dp)* &
                (g_phase**2+epsilon_value**2))
            frequency_delta = epsilon_value/(acos(-1.0_dp)* &
                (f_frequency**2+epsilon_value**2))
            phase_integral = phase_integral+phase_delta*dx
            frequency_integral = frequency_integral+frequency_delta*dx
        end do
        delta_weight_ratio = abs(n_mode)*tau**2*phase_integral/ &
            (tau*frequency_integral)
    end function delta_weight_ratio

    pure function direct_frequencies(period, delta_phi, q_fieldline, winding) &
            result(values)
        real(dp), intent(in) :: period, delta_phi, q_fieldline, winding
        real(dp) :: values(3)

        values(1) = 2.0_dp*acos(-1.0_dp)*winding/period
        values(2) = delta_phi/period
        values(3) = values(2)-q_fieldline*values(1)
    end function direct_frequencies

    subroutine require_close(label, value, expected, tolerance)
        character(*), intent(in) :: label
        real(dp), intent(in) :: value, expected, tolerance

        if (abs(value-expected) > tolerance) then
            write (*, '(A,2ES24.16)') trim(label)//' failed: ', value, expected
            error stop 1
        end if
    end subroutine require_close

    subroutine require_nonzero(label, value)
        character(*), intent(in) :: label
        real(dp), intent(in) :: value

        if (abs(value) <= 1.0e-13_dp) then
            write (*, '(A)') trim(label)//' failed: unexpected zero'
            error stop 1
        end if
    end subroutine require_nonzero

    subroutine require_less(label, value, bound)
        character(*), intent(in) :: label
        real(dp), intent(in) :: value, bound

        if (value >= bound) then
            write (*, '(A,2ES24.16)') trim(label)//' failed: ', value, bound
            error stop 1
        end if
    end subroutine require_less

    subroutine require_rate(label, coarse_error, fine_error)
        character(*), intent(in) :: label
        real(dp), intent(in) :: coarse_error, fine_error
        real(dp) :: observed_rate

        if (coarse_error <= 0.0_dp .or. fine_error <= 0.0_dp) then
            write (*, '(A)') trim(label)//' failed: zero error hides rate'
            error stop 1
        end if
        observed_rate = coarse_error/fine_error
        if (observed_rate < 7.0_dp .or. observed_rate > 16.0_dp) then
            write (*, '(A,ES24.16)') trim(label)//' failed: ', observed_rate
            error stop 1
        end if
    end subroutine require_rate

end program test_gc_full_fow_runtime_algebra
