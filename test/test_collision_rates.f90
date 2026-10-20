program test_collision_rates
    use iso_fortran_env, only: dp => real64
    use ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf
    use neort_collision_rates, only: collision_rates_t, initialize_collision_rates, &
        thermal_deflection_rates, pair_deflection_rates, deflection_frequency, &
        reference_frequency, collision_elementary_charge
    implicit none

    type(collision_rates_t) :: state, altered, single, empty
    real(dp) :: masses(3), charges(3), densities(3), temperatures(3)
    real(dp) :: pairs(3, 3), totals(3), pair_altered(3, 3), total_altered(3)
    real(dp) :: pair_vector(3), parallel_values(3), split_mass(4), split_density(4)
    real(dp) :: low, high, gamma, x, speed, loglam, bad, single_pair(1)
    integer :: a, b, order(3), k
    character(len=64) :: mode
    real(dp), parameter :: pi = acos(-1.0_dp), sqrt_pi = sqrt(pi)

    call get_command_argument(1, mode)
    masses = [3.345228e-27_dp, 5.017842e-27_dp, 9.1091e-31_dp]
    charges = [1.0_dp, 1.0_dp, -1.0_dp]
    densities = [5.0e15_dp, 5.0e15_dp, 1.0e16_dp]
    temperatures = 2000.0_dp
    loglam = 20.0_dp
    if (len_trim(mode) > 0) then
        call rejected_input(trim(mode))
        print *, 'Unexpected acceptance of invalid collision input: ', trim(mode)
        stop 0
    end if

    call initialize_collision_rates(masses, charges, densities, temperatures, &
                                    loglam, state)
    call thermal_deflection_rates(state, pairs, totals)
    call check_golden_mixture()
    do a = 1, 3
        call pair_deflection_rates(state, a, state%vth_m_s(a), pair_vector)
        do b = 1, 3
            call check(pair_vector(b), pairs(a, b), 2.0e-15_dp, 'pair API')
        end do
        call check(deflection_frequency(state, a, state%vth_m_s(a)), &
                   totals(a), 2.0e-15_dp, 'total API')
    end do
    !$omp parallel do default(none) shared(state, parallel_values)
    do a = 1, 3
        parallel_values(a) = deflection_frequency(state, a, state%vth_m_s(a))
    end do
    !$omp end parallel do
    do a = 1, 3
        call check(parallel_values(a), totals(a), 2.0e-15_dp, 'parallel state reads')
    end do

    ! Reference amplitudes have simple mass scaling; shaped thermal sums do not.
    call check(reference_frequency(state, 1)/reference_frequency(state, 2), &
               sqrt(masses(2)/masses(1)), 3.0e-15_dp, 'equal-T reference scaling')
    order = [3, 1, 2]
    call initialize_collision_rates(masses(order), charges(order), &
        densities(order), temperatures(order), loglam, altered)
    call thermal_deflection_rates(altered, pair_altered, total_altered)
    do a = 1, 3
        call check(total_altered(a), totals(order(a)), 3.0e-15_dp, 'permutation sum')
        do b = 1, 3
            call check(pair_altered(a, b), pairs(order(a), order(b)), &
                       3.0e-15_dp, 'permutation pairs')
        end do
    end do
    call initialize_collision_rates(masses, charges, 2.0_dp*densities, &
                                    temperatures, loglam, altered)
    call thermal_deflection_rates(altered, pair_altered, total_altered)
    do a = 1, 3
        call check(total_altered(a), 2.0_dp*totals(a), &
                   3.0e-15_dp, 'background density scaling')
    end do
    call initialize_collision_rates(masses, -charges, densities, &
                                    temperatures, loglam, altered)
    call thermal_deflection_rates(altered, pair_altered, total_altered)
    do a = 1, 3
        call check(total_altered(a), totals(a), 3.0e-15_dp, 'charge-sign invariance')
    end do

    split_mass = [masses(1), masses(2), masses(2), masses(3)]
    split_density = [densities(1), densities(2)/2, densities(2)/2, densities(3)]
    call initialize_collision_rates(split_mass, [1.0_dp, 1.0_dp, 1.0_dp, -1.0_dp], &
        split_density, [2000.0_dp, 2000.0_dp, 2000.0_dp, 2000.0_dp], loglam, altered)
    do a = 1, 3
        b = a
        if (a == 3) b = 4
        call check(deflection_frequency(altered, b, state%vth_m_s(a)), &
                   totals(a), 3.0e-15_dp, 'identical background splitting')
    end do
    densities(2) = 0.0_dp
    call initialize_collision_rates(masses, charges, densities, temperatures, &
                                    loglam, altered)
    do a = 1, 3
        call pair_deflection_rates(altered, a, state%vth_m_s(a), pair_vector)
        if (abs(pair_vector(2)) > 0.0_dp) error stop 'zero background contributed rate'
        call check(sum(pair_vector), totals(a) - pairs(a, 2), &
                   3.0e-15_dp, 'trace-background limit')
    end do
    if (abs(reference_frequency(altered, 2)) > 0.0_dp) &
        error stop 'trace reference amplitude is not zero'
    if (deflection_frequency(altered, 2, state%vth_m_s(2)) <= 0.0_dp) &
        error stop 'trace test species lost background collisions'
    call check_unequal_mixture()

    ! vth=1 exactly in the physical definition; these checks isolate B(x).
    call initialize_collision_rates([2.0_dp*collision_elementary_charge], &
        [2.0_dp], [1.0e10_dp], [1.0_dp], 10.0_dp, single)
    gamma = reference_frequency(single, 1)
    call check(single%vth_m_s(1), 1.0_dp, 2.0e-15_dp, 'single thermal speed')
    call check(deflection_frequency(single, 1, 1.0_dp)/gamma, &
        0.6289041451851547863407444422619733123325_dp, &
        3.0e-15_dp, 'unit-speed Gaussian integral')
    low = deflection_frequency(single, 1, 1.0e-8_dp)
    call check(low*1.0e-16_dp/gamma, 4.0_dp/(3.0_dp*sqrt_pi), &
               3.0e-15_dp, 'Maxwellian low-speed inverse-square limit')
    high = deflection_frequency(single, 1, 1.0e8_dp)
    call check(high*1.0e24_dp/gamma, 1.0_dp, &
               3.0e-15_dp, 'high-speed inverse-cube limit')
    do k = 1, 3
        if (k == 1) x = 0.5_dp
        if (k == 2) x = 8.0_dp
        if (k == 3) x = 1.0e-8_dp
        low = deflection_frequency(single, 1, x*(1.0_dp - 1.0e-10_dp))
        high = deflection_frequency(single, 1, x*(1.0_dp + 1.0e-10_dp))
        if (abs(high/low - 1.0_dp) > 1.0e-8_dp) &
            error stop 'collision shape branch is discontinuous'
    end do
    call pair_deflection_rates(single, 1, 1.0_dp, single_pair)
    call check_shape_golden()
    print *, 'All multispecies scalar collision-rate oracles passed.'

contains

    subroutine check(actual, expected, tolerance, label)
        real(dp), intent(in) :: actual, expected, tolerance
        character(len=*), intent(in) :: label

        if (abs(actual - expected) <= tolerance*abs(expected)) return
        print *, trim(label), actual, expected, abs(actual/expected - 1.0_dp)
        error stop 'collision rate oracle failed'
    end subroutine check

    subroutine check_golden_mixture()
        ! Filled from an independent 100-digit positive Gaussian integral.
        ! SI constants, explicit legacy-comparator masses, 50/50 D/T/e.
        real(dp), parameter :: expected_pair(3, 3) = reshape([ &
            0.0448280934468414679464717610884294_dp, &
            0.0315856317587746833054069889753953_dp, &
            4.31898608228269459330976295309951_dp, &
            0.0508896991830302328916311166162081_dp, &
            0.0366019850288546596688706059297059_dp, &
            4.31918211990853620985990477030000_dp, &
            0.00176954181439292302818499151192114_dp, &
            0.00117971595750433711519916182298357_dp, &
            5.43319623354218156795027424862159_dp], [3, 3])
        real(dp), parameter :: expected_reference(3) = [ &
            0.0712796914920058157636013749566165_dp, &
            0.0581996243928058713133605699697751_dp, &
            8.63914839032043888592037680940196_dp]
        integer :: i, j

        do i = 1, 3
            call check(reference_frequency(state, i), expected_reference(i), &
                       8.0e-15_dp, 'D/T/e independent reference frequency')
            do j = 1, 3
                call check(pairs(i, j), expected_pair(i, j), &
                           8.0e-15_dp, 'D/T/e independent thermal pair')
            end do
        end do
    end subroutine check_golden_mixture

    subroutine check_unequal_mixture()
        ! Unequal-mass/charge/temperature independent Gaussian oracle.
        real(dp), parameter :: expected_pair(3, 3) = reshape([ &
            0.392945039427472378434952041324323_dp, &
            1.61943296239464178216908541935421_dp, &
            7.84848932681211255256447942813956_dp, &
            1.33173784571972988822361398176941_dp, &
            7.15738360236755855525562982322375_dp, &
            18.8370514332625441668916089134490_dp, &
            0.00633718775708738991731301306993546_dp, &
            0.0228142039144091090620507237419567_dp, &
            6.91064695327476169478966381862463_dp], [3, 3])
        integer :: i, j

        call initialize_collision_rates([3.0e-27_dp, 2.0e-26_dp, 8.0e-31_dp], &
            [1.0_dp, 6.0_dp, -1.0_dp], [1.2e16_dp, 8.0e14_dp, 1.68e16_dp], &
            [800.0_dp, 1200.0_dp, 2300.0_dp], 17.5_dp, altered)
        call thermal_deflection_rates(altered, pair_altered, total_altered)
        do i = 1, 3
            do j = 1, 3
                call check(pair_altered(i, j), expected_pair(i, j), &
                           8.0e-15_dp, 'unequal mixture independent thermal pair')
            end do
        end do
    end subroutine check_unequal_mixture

    subroutine check_shape_golden()
        real(dp), parameter :: speeds(8) = [1.0e-8_dp, 0.01_dp, 0.5_dp, &
            0.5000000000001_dp, 1.0_dp, 8.0_dp, 8.00000000001_dp, 1.0e8_dp]
        real(dp), parameter :: shapes(8) = [ &
            7.5225277806367503421905037414e-9_dp, &
            0.00752237733330501853339532408169_dp, &
            0.358282701122398256410977300932_dp, &
            0.358282701122463143281653558355_dp, &
            0.628904145185154786340744442262_dp, &
            0.992187500000000000000000000000_dp, &
            0.992187500000019531249999963379_dp, &
            0.99999999999999995_dp]
        integer :: i
        real(dp) :: value

        do i = 1, size(speeds)
            value = deflection_frequency(single, 1, speeds(i))
            value = value*speeds(i)*speeds(i)*speeds(i)/gamma
            call check(value, shapes(i), 8.0e-15_dp, &
                       'independent Gaussian shape across small/large branches')
        end do
    end subroutine check_shape_golden

    subroutine rejected_input(which)
        character(len=*), intent(in) :: which

        select case (which)
        case ('zero-mass')
            masses(1) = 0.0_dp
        case ('negative-mass')
            masses(1) = -masses(1)
        case ('nan-mass')
            masses(1) = ieee_value(1.0_dp, ieee_quiet_nan)
        case ('zero-charge')
            charges(1) = 0.0_dp
        case ('infinite-charge')
            charges(1) = ieee_value(1.0_dp, ieee_positive_inf)
        case ('negative-density')
            densities(1) = -densities(1)
        case ('nan-density')
            densities(1) = ieee_value(1.0_dp, ieee_quiet_nan)
        case ('zero-temperature')
            temperatures(1) = 0.0_dp
        case ('nan-temperature')
            temperatures(1) = ieee_value(1.0_dp, ieee_quiet_nan)
        case ('zero-log')
            loglam = 0.0_dp
        case ('nan-log')
            loglam = ieee_value(1.0_dp, ieee_quiet_nan)
        case ('zero-speed', 'negative-speed', 'nan-speed', 'infinite-speed')
            call initialize_collision_rates(masses, charges, densities, temperatures, &
                                            loglam, state)
            speed = 0.0_dp
            if (which == 'negative-speed') speed = -1.0_dp
            if (which == 'nan-speed') speed = ieee_value(1.0_dp, ieee_quiet_nan)
            if (which == 'infinite-speed') &
                speed = ieee_value(1.0_dp, ieee_positive_inf)
            bad = deflection_frequency(state, 1, speed)
            return
        case ('invalid-index')
            call initialize_collision_rates(masses, charges, densities, temperatures, &
                                            loglam, state)
            bad = reference_frequency(state, 0)
            return
        case ('uninitialized')
            bad = reference_frequency(empty, 1)
            return
        case ('shape-mismatch')
            call initialize_collision_rates(masses, charges(1:2), densities, &
                                            temperatures, loglam, state)
            return
        case default
            error stop 'unknown collision rejection case'
        end select
        call initialize_collision_rates(masses, charges, densities, temperatures, &
                                        loglam, state)
    end subroutine rejected_input
end program test_collision_rates
