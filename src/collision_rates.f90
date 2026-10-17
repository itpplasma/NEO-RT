module neort_collision_rates
    use iso_fortran_env, only: dp => real64
    use ieee_arithmetic, only: ieee_is_finite
    implicit none
    private

    public :: collision_rates_t, initialize_collision_rates
    public :: thermal_deflection_rates, pair_deflection_rates
    public :: deflection_frequency, reference_frequency
    public :: radial_diffusion, pair_radial_diffusion
    public :: collision_elementary_charge, collision_epsilon0

    ! SI defining charge and 2022 CODATA recommended epsilon0 (not exact).
    ! https://physics.nist.gov/cuu/pdf/wall_2022.pdf
    real(dp), parameter :: collision_elementary_charge = 1.602176634e-19_dp
    real(dp), parameter :: collision_epsilon0 = 8.8541878188e-12_dp
    real(dp), parameter :: pi = acos(-1.0_dp), sqrt_pi = sqrt(pi)

    type collision_rates_t
        private
        real(dp), allocatable, public :: vth_m_s(:)
        real(dp), allocatable :: gamma_ab(:, :)
    end type collision_rates_t

contains

    subroutine initialize_collision_rates(mass_kg, charge_e, density_m3, &
                                          temperature_ev, coulomb_log, state)
        real(dp), intent(in) :: mass_kg(:), charge_e(:), density_m3(:)
        real(dp), intent(in) :: temperature_ev(:), coulomb_log
        type(collision_rates_t), intent(out) :: state
        integer :: n, a, b

        n = size(mass_kg)
        if (n == 0) error stop 'collision species list must not be empty'
        if (size(charge_e) /= n .or. size(density_m3) /= n .or. &
            size(temperature_ev) /= n) error stop 'collision species shapes differ'
        if (.not. all(ieee_is_finite(mass_kg))) &
            error stop 'collision masses must be finite'
        if (.not. all(ieee_is_finite(charge_e))) &
            error stop 'collision charges must be finite'
        if (.not. all(ieee_is_finite(density_m3))) &
            error stop 'collision densities must be finite'
        if (.not. all(ieee_is_finite(temperature_ev))) &
            error stop 'collision temperatures must be finite'
        if (.not. ieee_is_finite(coulomb_log)) &
            error stop 'Coulomb logarithm must be finite'
        if (any(mass_kg <= 0.0_dp)) error stop 'collision masses must be positive'
        if (any(abs(charge_e) <= 0.0_dp)) error stop 'collision charges must be nonzero'
        if (any(density_m3 < 0.0_dp)) &
            error stop 'collision densities must be nonnegative'
        if (any(temperature_ev <= 0.0_dp)) &
            error stop 'collision temperatures must be positive'
        if (coulomb_log <= 0.0_dp) error stop 'Coulomb logarithm must be positive'

        allocate (state%vth_m_s(n), state%gamma_ab(n, n))
        state%vth_m_s = sqrt(2.0_dp*collision_elementary_charge* &
                            temperature_ev/mass_kg)
        if (.not. all(ieee_is_finite(state%vth_m_s))) &
            error stop 'collision thermal speeds exceed numerical range'
        if (any(state%vth_m_s <= 0.0_dp)) &
            error stop 'collision thermal speeds are below numerical range'
        do b = 1, n
            do a = 1, n
                if (density_m3(b) <= 0.0_dp) then
                    state%gamma_ab(a, b) = 0.0_dp
                    cycle
                end if
                call collision_prefactor(pi, density_m3(b), coulomb_log, &
                    charge_e(a), charge_e(b), collision_elementary_charge, &
                    collision_epsilon0, mass_kg(a), state%gamma_ab(a, b))
            end do
        end do
        if (.not. all(ieee_is_finite(state%gamma_ab))) &
            error stop 'collision prefactors exceed numerical range'
    end subroutine initialize_collision_rates

    function reference_frequency(state, test_index) result(frequency)
        type(collision_rates_t), intent(in) :: state
        integer, intent(in) :: test_index
        real(dp) :: frequency

        call validate_index(state, test_index)
        ! Reference Coulomb amplitude Gamma_aa/v_ta^3, without shape B(1).
        frequency = state%gamma_ab(test_index, test_index)/state%vth_m_s(test_index)
        frequency = frequency/state%vth_m_s(test_index)/state%vth_m_s(test_index)
    end function reference_frequency

    subroutine thermal_deflection_rates(state, pair_rates, total_rates)
        type(collision_rates_t), intent(in) :: state
        real(dp), intent(out) :: pair_rates(:, :), total_rates(:)
        integer :: n, a

        call validate_state(state)
        n = size(state%vth_m_s)
        if (size(pair_rates, 1) /= n .or. size(pair_rates, 2) /= n .or. &
            size(total_rates) /= n) error stop 'thermal collision output shapes differ'
        do a = 1, n
            call pair_deflection_rates(state, a, state%vth_m_s(a), pair_rates(a, :))
            total_rates(a) = sum(pair_rates(a, :))
        end do
    end subroutine thermal_deflection_rates

    subroutine pair_deflection_rates(state, test_index, speed_m_s, pair_rates)
        type(collision_rates_t), intent(in) :: state
        integer, intent(in) :: test_index
        real(dp), intent(in) :: speed_m_s
        real(dp), intent(out) :: pair_rates(:)
        integer :: b

        call validate_index(state, test_index)
        call validate_speed(speed_m_s)
        if (size(pair_rates) /= size(state%vth_m_s)) &
            error stop 'pair collision output shape differs'
        do b = 1, size(state%vth_m_s)
            pair_rates(b) = pair_frequency(state%gamma_ab(test_index, b), &
                                          state%vth_m_s(b), speed_m_s)
        end do
    end subroutine pair_deflection_rates

    function deflection_frequency(state, test_index, speed_m_s) result(frequency)
        type(collision_rates_t), intent(in) :: state
        integer, intent(in) :: test_index
        real(dp), intent(in) :: speed_m_s
        real(dp) :: frequency
        integer :: b

        call validate_index(state, test_index)
        call validate_speed(speed_m_s)
        frequency = 0.0_dp
        do b = 1, size(state%vth_m_s)
            frequency = frequency + pair_frequency(state%gamma_ab(test_index, b), &
                                                   state%vth_m_s(b), speed_m_s)
        end do
        if (.not. ieee_is_finite(frequency)) &
            error stop 'total collision frequency exceeds numerical range'
    end function deflection_frequency

    subroutine pair_radial_diffusion(state, test_index, speed_m_s, pair_coefficients)
        type(collision_rates_t), intent(in) :: state
        integer, intent(in) :: test_index
        real(dp), intent(in) :: speed_m_s
        real(dp), intent(out) :: pair_coefficients(:)
        integer :: b

        call validate_index(state, test_index)
        call validate_radial_speed(speed_m_s)
        if (size(pair_coefficients) /= size(state%vth_m_s)) &
            error stop 'pair radial diffusion output shape differs'
        do b = 1, size(state%vth_m_s)
            pair_coefficients(b) = pair_radial_coefficient( &
                state%gamma_ab(test_index, b), state%vth_m_s(b), speed_m_s)
        end do
    end subroutine pair_radial_diffusion

    function radial_diffusion(state, test_index, speed_m_s) result(coefficient)
        type(collision_rates_t), intent(in) :: state
        integer, intent(in) :: test_index
        real(dp), intent(in) :: speed_m_s
        real(dp) :: coefficient
        integer :: b

        call validate_index(state, test_index)
        call validate_radial_speed(speed_m_s)
        coefficient = 0.0_dp
        do b = 1, size(state%vth_m_s)
            coefficient = coefficient + pair_radial_coefficient( &
                state%gamma_ab(test_index, b), state%vth_m_s(b), speed_m_s)
        end do
        if (.not. ieee_is_finite(coefficient)) &
            error stop 'total radial diffusion exceeds numerical range'
    end function radial_diffusion

    function pair_radial_coefficient(gamma_ab, background_vth, speed_m_s) &
        result(coefficient)
        real(dp), intent(in) :: gamma_ab, background_vth, speed_m_s
        real(dp) :: coefficient, x, inverse_x, g_over_x

        if (gamma_ab <= 0.0_dp) then
            coefficient = 0.0_dp
            return
        end if
        ! Dvv = Gamma*G(v/vtb)/v [m^2/s^3], with a finite isotropic origin.
        ! C_rad[f] = v^-2 d_v[v^2 Dvv (d_v f + ma*v*f/Tb)].
        ! This radial Fick friction differs from Cartesian slowing-down drag.
        if (speed_m_s <= 0.5_dp*background_vth) then
            x = speed_m_s/background_vth
            if (x < sqrt(epsilon(1.0_dp))) x = 0.0_dp
            call collision_small_parallel_shape(x*x, sqrt_pi, g_over_x)
        else
            inverse_x = background_vth/speed_m_s
            if (inverse_x < 0.125_dp) then
                call collision_large_parallel_shape(inverse_x, g_over_x)
            else
                x = speed_m_s/background_vth
                call collision_middle_parallel_shape(x, inverse_x, sqrt_pi, g_over_x)
            end if
        end if
        coefficient = gamma_ab/background_vth*g_over_x
        if (.not. ieee_is_finite(coefficient)) &
            error stop 'pair radial diffusion exceeds numerical range'
        if (coefficient < 0.0_dp) error stop 'negative pair radial diffusion'
    end function pair_radial_coefficient

    function pair_frequency(gamma_ab, background_vth, speed_m_s) result(frequency)
        real(dp), intent(in) :: gamma_ab, background_vth, speed_m_s
        real(dp) :: frequency, x, inverse_x, b_over_x

        if (gamma_ab <= 0.0_dp) then
            frequency = 0.0_dp
            return
        end if
        ! Small x uses B(x)/x directly, preserving nu_D ~ v^-2 as v -> 0.
        if (speed_m_s <= 0.5_dp*background_vth) then
            x = speed_m_s/background_vth
            if (x < sqrt(epsilon(1.0_dp))) x = 0.0_dp
            call collision_small_shape(x*x, sqrt_pi, b_over_x)
        else
            inverse_x = background_vth/speed_m_s
            if (inverse_x < 0.125_dp) then
                call collision_large_shape(inverse_x, b_over_x)
            else
                x = speed_m_s/background_vth
                call collision_middle_shape(x, inverse_x, sqrt_pi, b_over_x)
            end if
        end if
        ! Sequential division avoids forming v^3 or v^2 at extreme speeds.
        frequency = gamma_ab/background_vth*b_over_x
        frequency = frequency/speed_m_s/speed_m_s
        if (.not. ieee_is_finite(frequency)) &
            error stop 'pair collision frequency exceeds numerical range'
        if (frequency < 0.0_dp) error stop 'negative pair collision frequency'
    end function pair_frequency

    subroutine validate_state(state)
        type(collision_rates_t), intent(in) :: state

        if (.not. allocated(state%vth_m_s)) &
            error stop 'collision rate state is not initialized'
        if (.not. allocated(state%gamma_ab)) &
            error stop 'collision rate state is not initialized'
    end subroutine validate_state

    subroutine validate_index(state, test_index)
        type(collision_rates_t), intent(in) :: state
        integer, intent(in) :: test_index

        call validate_state(state)
        if (test_index < 1 .or. test_index > size(state%vth_m_s)) &
            error stop 'collision test species index is out of range'
    end subroutine validate_index

    subroutine validate_speed(speed_m_s)
        real(dp), intent(in) :: speed_m_s

        if (.not. ieee_is_finite(speed_m_s)) &
            error stop 'collision speed must be finite'
        ! For a Maxwellian bath nu_D diverges as v^-2 at zero speed.
        if (speed_m_s <= 0.0_dp) error stop 'collision speed must be positive'
    end subroutine validate_speed

    subroutine validate_radial_speed(speed_m_s)
        real(dp), intent(in) :: speed_m_s

        if (.not. ieee_is_finite(speed_m_s)) &
            error stop 'radial diffusion speed must be finite'
        if (speed_m_s < 0.0_dp) &
            error stop 'radial diffusion speed must be nonnegative'
    end subroutine validate_radial_speed

    include 'generated/collision_rate_leaf.inc'

end module neort_collision_rates
