program test_landau_projection
    use iso_fortran_env, only: dp => real64
    use ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use omp_lib, only: omp_set_num_threads
    use neort_landau_pair_quadrature, only: landau_quadrature_t, &
        initialize_landau_quadrature
    use neort_landau_projection, only: assemble_landau_projection
    use neort_collision_rates, only: collision_elementary_charge, collision_epsilon0
    implicit none

    type(landau_quadrature_t) :: quadrature
    real(dp) :: masses(3), charges(3), densities(3), vth(3), temperature
    complex(dp) :: matrix(9, 9), null_vector(9)
    integer :: a
    character(len=64) :: mode

    masses = [3.345228e-27_dp, 5.017842e-27_dp, 9.1091e-31_dp]
    charges = [1.0_dp, 1.0_dp, -1.0_dp]
    densities = [5.0e15_dp, 5.0e15_dp, 1.0e16_dp]
    temperature = 2000.0_dp
    vth = sqrt(2.0_dp*collision_elementary_charge*temperature/masses)
    call initialize_landau_quadrature(2, 2, 2, 4, quadrature)
    call get_command_argument(1, mode)
    if (len_trim(mode) > 0) then
        call rejected_input(trim(mode))
        stop 0
    end if
    call assemble_landau_projection(masses, charges, densities, temperature, &
        20.0_dp, [3, 3, 3], quadrature, invariant_gradients, matrix)
    if (maxval(abs(matrix - conjg(transpose(matrix)))) > &
        1.0e-14_dp*maxval(abs(matrix))) error stop 'Landau matrix is not Hermitian'
    do a = 1, 3
        null_vector = cmplx(0.0_dp, 0.0_dp, dp)
        null_vector(3*a - 2) = 1.0_dp
        call check_null(null_vector, 'individual density')
    end do
    null_vector = cmplx(0.0_dp, 0.0_dp, dp)
    null_vector(2:9:3) = sqrt(masses/masses(1))
    call check_null(null_vector, 'common parallel flow')
    null_vector = cmplx(0.0_dp, 0.0_dp, dp)
    null_vector(3:9:3) = 1.0_dp
    call check_null(null_vector, 'common temperature')
    call check_absolute_flow()
    call check_temperature_fixture()
    call check_permutation_and_scaling()
    call check_thread_counts()
    call check_complex_basis()
    call check_arbitrary_charge_fixture()
    call check_stress_and_splitting()
    call check_nonpolynomial_refinement()
    print *, 'All local coupled Landau projection oracles passed.'

contains

    pure subroutine invariant_gradients(species_index, velocity, gradients)
        integer, intent(in) :: species_index
        real(dp), intent(in) :: velocity(3)
        complex(dp), intent(out) :: gradients(:, :)

        gradients = cmplx(0.0_dp, 0.0_dp, dp)
        gradients(3, 2) = 1.0_dp/vth(species_index)
        gradients(:, 3) = 2.0_dp*velocity/vth(species_index)**2
    end subroutine invariant_gradients

    subroutine check_null(vector, label)
        complex(dp), intent(in) :: vector(:)
        character(len=*), intent(in) :: label

        if (maxval(abs(matmul(matrix, vector))) <= &
            5.0e-14_dp*maxval(abs(matrix))*maxval(abs(vector))) return
        print *, label
        error stop 'Landau collision invariant is not a null'
    end subroutine check_null

    subroutine check_absolute_flow()
        real(dp) :: expected(3, 3), pi, kab, reaction, aa, bb
        integer :: i, j

        pi = acos(-1.0_dp)
        expected = 0.0_dp
        do j = 2, 3
            do i = 1, j - 1
                kab = 2.0_dp*pi*20.0_dp*(charges(i)*charges(j)* &
                    collision_elementary_charge**2/(4*pi*collision_epsilon0))**2
                reaction = 4*kab*densities(i)*densities(j)/ &
                    (3*sqrt(pi)*sqrt(vth(i)**2 + vth(j)**2))
                aa = 1/(masses(i)*vth(i)); bb = 1/(masses(j)*vth(j))
                expected(i, i) = expected(i, i) + reaction*aa*aa
                expected(j, j) = expected(j, j) + reaction*bb*bb
                expected(i, j) = -reaction*aa*bb
                expected(j, i) = expected(i, j)
            end do
        end do
        if (maxval(abs(matrix(2:9:3, 2:9:3) - expected)) > &
            5.0e-14_dp*maxval(abs(expected))) error stop 'absolute Landau flow rate'
    end subroutine check_absolute_flow

    subroutine check_temperature_fixture()
        ! Independent 100-digit relative-Maxwellian heat-exchange derivative.
        real(dp), parameter :: expected(3, 3) = reshape([ &
            275738582494511.443540219276285602_dp, &
            -249204848232571.280161913646903611_dp, &
            -26533734261940.1633783056293819909_dp, &
            -249204848232571.280161913646903611_dp, &
            266896412414989.916126782442715896_dp, &
            -17691564182418.6359648687958122841_dp, &
            -26533734261940.1633783056293819909_dp, &
            -17691564182418.6359648687958122841_dp, &
            44225298444358.7993431744251942750_dp], [3, 3])

        if (maxval(abs(matrix(3:9:3, 3:9:3) - expected)) > &
            5.0e-14_dp*maxval(abs(expected))) error stop 'absolute Landau heat exchange'
        if (maxval(abs(matrix(2:9:3, 3:9:3))) > &
            5.0e-14_dp*maxval(abs(expected))) error stop 'odd/even moment coupling'
    end subroutine check_temperature_fixture

    subroutine check_permutation_and_scaling()
        complex(dp) :: changed(9, 9)
        real(dp) :: original_vth(3)
        integer :: order(3), i, j, global_order(9)

        original_vth = vth; order = [3, 1, 2]
        vth = original_vth(order)
        call assemble_landau_projection(masses(order), charges(order), &
            densities(order), temperature, 20.0_dp, [3, 3, 3], quadrature, &
            invariant_gradients, changed)
        do i = 1, 3
            global_order(3*i - 2:3*i) = [3*order(i) - 2, 3*order(i) - 1, 3*order(i)]
        end do
        do j = 1, 9
            do i = 1, 9
                if (abs(changed(i, j) - matrix(global_order(i), global_order(j))) > &
                    8.0e-14_dp*maxval(abs(matrix))) error stop 'species permutation'
            end do
        end do
        vth = original_vth
        call assemble_landau_projection(masses, -charges, 2*densities, temperature, &
            20.0_dp, [3, 3, 3], quadrature, invariant_gradients, changed)
        if (maxval(abs(changed - 4*matrix)) > 1.0e-13_dp*maxval(abs(matrix))) &
            error stop 'pair density scaling or charge-sign invariance'
    end subroutine check_permutation_and_scaling

    subroutine check_thread_counts()
        complex(dp) :: changed(9, 9)

        call omp_set_num_threads(1)
        call assemble_landau_projection(masses, charges, densities, temperature, &
            20.0_dp, [3, 3, 3], quadrature, invariant_gradients, changed)
        call omp_set_num_threads(4)
        if (maxval(abs(changed - matrix)) > 1.0e-14_dp*maxval(abs(matrix))) &
            error stop 'Landau pair assembly differs across thread counts'
    end subroutine check_thread_counts

    subroutine check_complex_basis()
        complex(dp) :: changed(9, 9), phases(9), vector(9), entropy
        integer :: i, j

        phases = cmplx(0.0_dp, 0.0_dp, dp)
        do i = 1, 9
            phases(i) = exp(cmplx(0.0_dp, real(modulo(i - 1, 3) + 1, dp)/7, dp))
            vector(i) = cmplx(sin(real(i, dp)), cos(real(2*i, dp)), dp)
        end do
        call assemble_landau_projection(masses, charges, densities, temperature, &
            20.0_dp, [3, 3, 3], quadrature, phased_gradients, changed)
        do j = 1, 9
            do i = 1, 9
                if (abs(changed(i, j) - conjg(phases(i))*matrix(i, j)*phases(j)) > &
                    2.0e-13_dp*maxval(abs(matrix))) error stop 'complex bilinear form'
            end do
        end do
        entropy = dot_product(vector, matmul(changed, vector))
        if (real(entropy, dp) < 0.0_dp .or. &
            abs(aimag(entropy)) > 1.0e-13_dp*abs(entropy)) &
            error stop 'complex entropy production is not real nonnegative'
    end subroutine check_complex_basis

    pure subroutine phased_gradients(species_index, velocity, gradients)
        integer, intent(in) :: species_index
        real(dp), intent(in) :: velocity(3)
        complex(dp), intent(out) :: gradients(:, :)
        integer :: j

        call invariant_gradients(species_index, velocity, gradients)
        do j = 1, size(gradients, 2)
            gradients(:, j) = gradients(:, j)*exp(cmplx(0.0_dp, real(j, dp)/7, dp))
        end do
    end subroutine phased_gradients

    subroutine check_arbitrary_charge_fixture()
        complex(dp) :: changed(6, 6)
        real(dp), parameter :: flow(2, 2) = reshape([ &
            6311385758592285.53045982427503904_dp, &
            -2444389193451858.62387596620802663_dp, &
            -2444389193451858.62387596620802663_dp, &
            946707863788842.829568973641255856_dp], [2, 2])
        real(dp), parameter :: heat = 4939345376289614.76296855812829142_dp

        masses(1:2) = [3.0e-27_dp, 2.0e-26_dp]
        charges(1:2) = [1.0_dp, 6.0_dp]
        densities(1:2) = [1.2e16_dp, 8.0e14_dp]
        temperature = 800.0_dp
        vth(1:2) = sqrt(2*collision_elementary_charge*temperature/masses(1:2))
        call assemble_landau_projection(masses(1:2), charges(1:2), densities(1:2), &
            temperature, 17.5_dp, [3, 3], quadrature, invariant_gradients, changed)
        if (maxval(abs(changed(2:6:3, 2:6:3) - flow)) > &
            8.0e-14_dp*maxval(flow)) error stop 'arbitrary-Z absolute flow fixture'
        if (abs(changed(3, 3)/heat - 1) > 8.0e-14_dp .or. &
            abs(changed(3, 6)/heat + 1) > 8.0e-14_dp) &
            error stop 'arbitrary-Z absolute heat fixture'
    end subroutine check_arbitrary_charge_fixture

    subroutine check_stress_and_splitting()
        complex(dp) :: single(1, 1), split(2, 2)

        masses(1:2) = 3.3435837768e-27_dp
        densities(1) = 5.0e19_dp
        temperature = 1000.0_dp
        vth(1:2) = sqrt(2*collision_elementary_charge*temperature/masses(1:2))
        call initialize_landau_quadrature(2, 2, 3, 6, quadrature)
        call assemble_landau_projection(masses(1:1), [1.0_dp], densities(1:1), &
            temperature, 17.0_dp, [1], quadrature, stress_gradients, single)
        if (abs(single(1, 1)/1.8235399640301296e22_dp - 1) > 8.0e-14_dp) &
            error stop 'self-collision traceless-stress absolute rate or self half'
        call assemble_landau_projection(masses(1:2), [1.0_dp, 1.0_dp], &
            [1.5e19_dp, 3.5e19_dp], temperature, 17.0_dp, [1, 1], quadrature, &
            stress_gradients, split)
        if (abs(sum(split)/single(1, 1) - 1) > 8.0e-14_dp) &
            error stop 'identical-species stress splitting'
    end subroutine check_stress_and_splitting

    pure subroutine stress_gradients(species_index, velocity, gradients)
        integer, intent(in) :: species_index
        real(dp), intent(in) :: velocity(3)
        complex(dp), intent(out) :: gradients(:, :)

        gradients(:, 1) = -2*velocity/(3*vth(species_index)**2)
        gradients(3, 1) = gradients(3, 1) + 2*velocity(3)/vth(species_index)**2
    end subroutine stress_gradients

    subroutine check_nonpolynomial_refinement()
        complex(dp) :: single(1, 1), mixture(2, 2)
        real(dp) :: previous(2), errors(2), cross
        integer :: k, radial(3), polar(3), azimuth(3)

        masses(1:2) = [3.3435837768e-27_dp, 9.1093837139e-31_dp]
        vth(1:2) = sqrt(2*collision_elementary_charge*temperature/masses(1:2))
        radial = [2, 8, 16]; polar = [3, 8, 12]; azimuth = [6, 16, 24]
        previous = huge(1.0_dp)
        do k = 1, 3
            call initialize_landau_quadrature(radial(k), radial(k), &
                polar(k), azimuth(k), quadrature)
            call assemble_landau_projection(masses(1:1), [1.0_dp], [5.0e19_dp], &
                temperature, 17.0_dp, [1], quadrature, exponential_gradients, single)
            call assemble_landau_projection(masses(1:2), [1.0_dp, -1.0_dp], &
                [5.0e19_dp, 1.0e20_dp], temperature, 17.0_dp, [1, 1], quadrature, &
                exponential_gradients, mixture)
            cross = real(mixture(1, 1) - single(1, 1), dp)
            errors = [abs(real(single(1, 1), dp)/1.0679630988348584e21_dp - 1), &
                      abs(cross/3.505532237843788e20_dp - 1)]
            print '(a,3i4,a,2es12.4)', 'Landau orders ', radial(k), polar(k), &
                azimuth(k), ' self/cross relative errors ', errors
            if (any(errors >= previous/4)) error stop 'Landau quadrature refinement'
            previous = errors
        end do
        if (maxval(errors) > 2.0e-8_dp) error stop 'Landau nonpolynomial reference'
    end subroutine check_nonpolynomial_refinement

    pure subroutine exponential_gradients(species_index, velocity, gradients)
        integer, intent(in) :: species_index
        real(dp), intent(in) :: velocity(3)
        complex(dp), intent(out) :: gradients(:, :)
        real(dp) :: u2

        gradients = cmplx(0.0_dp, 0.0_dp, dp)
        if (species_index /= 1) return
        u2 = sum(velocity*velocity)/vth(species_index)**2
        gradients(:, 1) = -1.4_dp*exp(-0.7_dp*u2)*velocity/vth(species_index)**2
    end subroutine exponential_gradients

    subroutine rejected_input(which)
        character(len=*), intent(in) :: which

        select case (which)
        case ('small-order')
            call initialize_landau_quadrature(1, 2, 2, 4, quadrature)
        case ('nan-node')
            quadrature%t(1) = ieee_value(1.0_dp, ieee_quiet_nan)
        case ('zero-relative-node')
            quadrature%t(1) = 0.0_dp
        case ('negative-weight')
            quadrature%wy(1) = -quadrature%wy(1)
        case ('nan-mass')
            masses(1) = ieee_value(1.0_dp, ieee_quiet_nan)
        case ('zero-temperature')
            temperature = 0.0_dp
        case ('negative-density')
            densities(1) = -densities(1)
        case ('zero-charge')
            charges(1) = 0.0_dp
        case ('nan-temperature')
            temperature = ieee_value(1.0_dp, ieee_quiet_nan)
        case ('nan-gradient')
            call assemble_landau_projection(masses, charges, densities, temperature, &
                20.0_dp, [3, 3, 3], quadrature, invalid_gradients, matrix)
            return
        case default
            error stop 'unknown Landau rejection case'
        end select
        call assemble_landau_projection(masses, charges, densities, temperature, &
            20.0_dp, [3, 3, 3], quadrature, invariant_gradients, matrix)
    end subroutine rejected_input

    pure subroutine invalid_gradients(species_index, velocity, gradients)
        integer, intent(in) :: species_index
        real(dp), intent(in) :: velocity(3)
        complex(dp), intent(out) :: gradients(:, :)

        call invariant_gradients(species_index, velocity, gradients)
        gradients(1, 1) = ieee_value(1.0_dp, ieee_quiet_nan)
    end subroutine invalid_gradients
end program test_landau_projection
