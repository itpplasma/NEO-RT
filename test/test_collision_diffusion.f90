program test_collision_diffusion
    use iso_fortran_env, only: dp => real64
    use ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf
    use neort_collision_rates, only: collision_rates_t, initialize_collision_rates, &
        radial_diffusion, pair_radial_diffusion, deflection_frequency, &
        reference_frequency, collision_elementary_charge
    implicit none

    type(collision_rates_t) :: state, altered, single
    real(dp) :: masses(3), charges(3), densities(3), temperatures(3)
    real(dp) :: pairs(3, 3), zero_pairs(3, 3), values(3), thermal(3), origin(3)
    real(dp) :: parallel_values(3), speed, gamma, value
    integer :: a, b, order(3)
    character(len=64) :: mode

    masses = [3.345228e-27_dp, 5.017842e-27_dp, 9.1091e-31_dp]
    charges = [1.0_dp, 1.0_dp, -1.0_dp]
    densities = [5.0e15_dp, 5.0e15_dp, 1.0e16_dp]
    temperatures = 2000.0_dp
    call initialize_collision_rates(masses, charges, densities, temperatures, &
                                    20.0_dp, state)
    call get_command_argument(1, mode)
    if (len_trim(mode) > 0) then
        speed = 0.0_dp
        select case (trim(mode))
        case ('negative-speed')
            speed = -1.0_dp
        case ('nan-speed')
            speed = ieee_value(1.0_dp, ieee_quiet_nan)
        case ('infinite-speed')
            speed = ieee_value(1.0_dp, ieee_positive_inf)
        case ('pair-shape')
            call pair_radial_diffusion(state, 1, speed, values(1:2))
            stop 0
        case default
            error stop 'unknown diffusion rejection case'
        end select
        value = radial_diffusion(state, 1, speed)
        print *, 'Unexpected acceptance of invalid diffusion input: ', value
        stop 0
    end if

    do a = 1, 3
        call pair_radial_diffusion(state, a, state%vth_m_s(a), pairs(a, :))
        call pair_radial_diffusion(state, a, 0.0_dp, zero_pairs(a, :))
        thermal(a) = radial_diffusion(state, a, state%vth_m_s(a))
        origin(a) = radial_diffusion(state, a, 0.0_dp)
        call check(thermal(a), sum(pairs(a, :)), 'thermal pair sum')
        call check(origin(a), sum(zero_pairs(a, :)), 'origin pair sum')
        call check(radial_diffusion(state, a, 1.0e-200_dp), origin(a), &
                   'finite diffusion at arbitrarily small speed')
        speed = 1.0e-8_dp*minval(state%vth_m_s)
        call check(deflection_frequency(state, a, speed)*speed*speed/2.0_dp, &
                   origin(a), 'isotropic Cartesian diffusion at origin')
    end do
    call check_golden_mixtures()
    !$omp parallel do default(none) shared(state, parallel_values)
    do a = 1, 3
        parallel_values(a) = radial_diffusion(state, a, state%vth_m_s(a))
    end do
    !$omp end parallel do
    do a = 1, 3
        call check(parallel_values(a), thermal(a), 'parallel immutable state reads')
    end do

    order = [3, 1, 2]
    call initialize_collision_rates(masses(order), charges(order), &
        densities(order), temperatures(order), 20.0_dp, altered)
    do a = 1, 3
        call pair_radial_diffusion(altered, a, state%vth_m_s(order(a)), values)
        call check(radial_diffusion(altered, a, state%vth_m_s(order(a))), &
                   thermal(order(a)), 'species permutation sum')
        do b = 1, 3
            call check(values(b), pairs(order(a), order(b)), 'species permutation pair')
        end do
    end do
    call initialize_collision_rates(masses, charges, 3.0_dp*densities, &
                                    temperatures, 20.0_dp, altered)
    do a = 1, 3
        call check(radial_diffusion(altered, a, state%vth_m_s(a)), &
                   3.0_dp*thermal(a), 'linear background density scaling')
    end do
    call initialize_collision_rates(masses, charges, densities, &
                                    4.0_dp*temperatures, 20.0_dp, altered)
    do a = 1, 3
        call check(radial_diffusion(altered, a, 2.0_dp*state%vth_m_s(a)), &
                   thermal(a)/2.0_dp, 'temperature and speed scaling')
    end do
    call initialize_collision_rates([masses(1), masses(2), masses(2), masses(3)], &
        [1.0_dp, 1.0_dp, 1.0_dp, -1.0_dp], &
        [densities(1), densities(2)/2, densities(2)/2, densities(3)], &
        [2000.0_dp, 2000.0_dp, 2000.0_dp, 2000.0_dp], 20.0_dp, altered)
    do a = 1, 3
        b = a
        if (a == 3) b = 4
        call check(radial_diffusion(altered, b, state%vth_m_s(a)), thermal(a), &
                   'identical background splitting')
        call check(radial_diffusion(altered, b, 0.0_dp), origin(a), &
                   'identical background splitting at origin')
    end do
    densities(2) = 0.0_dp
    call initialize_collision_rates(masses, charges, densities, temperatures, &
                                    20.0_dp, altered)
    do a = 1, 3
        call pair_radial_diffusion(altered, a, state%vth_m_s(a), values)
        if (abs(values(2)) > 0.0_dp) error stop 'zero-density bath diffuses particles'
        call check(sum(values), thermal(a) - pairs(a, 2), 'trace background limit')
    end do
    if (radial_diffusion(altered, 2, state%vth_m_s(2)) <= 0.0_dp) &
        error stop 'trace test species lost speed diffusion'

    call initialize_collision_rates([2.0_dp*collision_elementary_charge], &
        [2.0_dp], [1.0e10_dp], [1.0_dp], 10.0_dp, single)
    gamma = reference_frequency(single, 1)
    call check_shapes()
    call check(radial_diffusion(single, 1, 1.0e8_dp)*1.0e24_dp/gamma, &
               0.5_dp, 'high-speed radial diffusion inverse-cube coefficient')
    print *, 'All multispecies radial-diffusion oracles passed.'

contains

    subroutine check(actual, expected, label)
        real(dp), intent(in) :: actual, expected
        character(len=*), intent(in) :: label

        if (abs(actual - expected) <= 1.2e-14_dp*abs(expected)) return
        print *, trim(label), actual, expected
        error stop 'radial diffusion oracle failed'
    end subroutine check

    subroutine check_golden_mixtures()
        ! 100-digit direct Landau-tensor Maxwellian shell integrals, not erf/G.
        ! Dvv=(Gamma/2)<u_perp^2/|v-u|^3>; positive radial quadrature.
        real(dp), parameter :: expected(3, 3) = reshape([ &
            2919519394.03069686891599692546380_dp, &
            1554088066.91795148982174052392841_dp, &
            413766765333.487551689948168211533_dp, &
            2769239080.35659918748808579883047_dp, &
            1589185068.78550036893073367517257_dp, &
            275844510222.325034459965445474355_dp, &
            169483799.879495783384216724850716_dp, &
            75330235.5661944285073141436349958_dp, &
            1299471444273419.44307722344698424_dp], [3, 3])
        real(dp), parameter :: expected_zero(3, 3) = reshape([ &
            5136227807.43717501012384238024813_dp, &
            2282767914.41652222672170772455472_dp, &
            69269935721693760.6451729833973468_dp, &
            6290568665.45754631144687549471675_dp, &
            2795808295.75890947175416688654078_dp, &
            84837998516769465.0432593783108410_dp, &
            169511492.108172375337523030646186_dp, &
            75338440.9369655001500102358427493_dp, &
            2286123319027867.18342435600293392_dp], [3, 3])
        real(dp), parameter :: unequal(3, 3) = reshape([ &
            11414512457.1460402727045845577073_dp, &
            14236530689.5324615728789399662174_dp, &
            335339984309.052475613511155899630_dp, &
            13971122461.3819161599079963636467_dp, &
            46780231637.4507916488218639330112_dp, &
            181083591526.888336831296024185800_dp, &
            270744465.486258827861247646146643_dp, &
            219312475.820908923940501591376102_dp, &
            2164279666045873.34258587596869718_dp], [3, 3])
        integer :: i, j

        do i = 1, 3
            do j = 1, 3
                call check(pairs(i, j), expected(i, j), 'D/T/e thermal Landau tensor')
                call check(zero_pairs(i, j), expected_zero(i, j), &
                           'D/T/e zero-speed Landau tensor')
            end do
        end do
        call initialize_collision_rates([3.0e-27_dp, 2.0e-26_dp, 8.0e-31_dp], &
            [1.0_dp, 6.0_dp, -1.0_dp], [1.2e16_dp, 8.0e14_dp, 1.68e16_dp], &
            [800.0_dp, 1200.0_dp, 2300.0_dp], 17.5_dp, altered)
        do i = 1, 3
            call pair_radial_diffusion(altered, i, altered%vth_m_s(i), values)
            do j = 1, 3
                call check(values(j), unequal(i, j), 'unequal mixture Landau tensor')
            end do
        end do
    end subroutine check_golden_mixtures

    subroutine check_shapes()
        real(dp), parameter :: speeds(9) = [0.0_dp, 1.0e-8_dp, 0.01_dp, 0.5_dp, &
            0.5000000000001_dp, 1.0_dp, 8.0_dp, 8.00000000001_dp, 1.0e8_dp]
        real(dp), parameter :: shapes(9) = [ &
            0.376126389031837524632052967707_dp, &
            0.376126389031837502064469625797_dp, &
            0.376103822254459838251419946638_dp, &
            0.324434353381296562543538705919_dp, &
            0.324434353381277658447297015714_dp, &
            0.213796647764560083000476192821_dp, &
            0.000976562499999999999999999999_dp, &
            0.000976562499996337890625009154_dp, &
            5.0e-25_dp]
        integer :: i

        do i = 1, size(speeds)
            call check(radial_diffusion(single, 1, speeds(i))/gamma, shapes(i), &
                       'Landau tensor across origin and numerical branches')
        end do
    end subroutine check_shapes
end program test_collision_diffusion
