module neort_landau_projection
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use neort_collision_rates, only: collision_elementary_charge, collision_epsilon0
    use neort_landau_pair_quadrature, only: landau_quadrature_t, &
        validate_landau_quadrature
    implicit none
    private

    real(dp), parameter :: pi = acos(-1.0_dp)
    type :: landau_pair_t
        integer :: a, b, first_a, last_a, first_b, last_b
        real(dp) :: ma, mb, vcm, vrel, fa, fb, prefactor
    end type landau_pair_t
    type :: landau_workspace_t
        complex(dp), allocatable :: ga(:, :), gb(:, :), delta(:, :), transverse(:, :)
    end type landau_workspace_t

    abstract interface
        pure subroutine landau_gradient(species_index, velocity_m_s, gradient_columns)
            import dp
            integer, intent(in) :: species_index
            real(dp), intent(in) :: velocity_m_s(3)
            complex(dp), intent(out) :: gradient_columns(:, :)
        end subroutine landau_gradient
    end interface
    public :: assemble_landau_projection, landau_gradient

contains

    ! Local velocity-pair Dirichlet matrix in SI; no spatial Jacobian is included.
    ! Dimensionless basis functions give entries in m^-3 s^-1, not generator rates.
    ! Species-contiguous basis functions must be gyrotropic about the local z axis.
    ! A common spatial harmonic phase is permitted and cancels from the pair form.
    subroutine assemble_landau_projection(mass_kg, charge_e, density_m3, &
        temperature_ev, coulomb_log, ndof_species, quadrature, gradient, matrix, &
        pair_blocks, pair_species)
        real(dp), intent(in) :: mass_kg(:), charge_e(:), density_m3(:)
        real(dp), intent(in) :: temperature_ev, coulomb_log
        integer, intent(in) :: ndof_species(:)
        type(landau_quadrature_t), intent(in) :: quadrature
        procedure(landau_gradient) :: gradient
        complex(dp), intent(out) :: matrix(:, :)
        complex(dp), allocatable, intent(out), optional :: pair_blocks(:, :, :)
        integer, allocatable, intent(out), optional :: pair_species(:, :)
        type(landau_pair_t), allocatable :: pairs(:)
        complex(dp), allocatable :: blocks(:, :, :)
        integer :: k

        call validate_inputs(mass_kg, charge_e, density_m3, temperature_ev, &
                             coulomb_log, ndof_species, matrix)
        call validate_landau_quadrature(quadrature)
        call initialize_pairs(mass_kg, charge_e, density_m3, temperature_ev, &
                              coulomb_log, ndof_species, pairs)
        allocate (blocks(size(matrix, 1), size(matrix, 2), size(pairs)))
        !$omp parallel do default(none) shared(pairs, quadrature, blocks)
        do k = 1, size(pairs)
            call integrate_pair(pairs(k), quadrature, gradient, blocks(:, :, k))
        end do
        !$omp end parallel do
        call combine_blocks(blocks, matrix)
        if (present(pair_species)) then
            allocate(pair_species(2, size(pairs)))
            do k = 1, size(pairs)
                pair_species(:, k) = [pairs(k)%a, pairs(k)%b]
            end do
        end if
        if (present(pair_blocks)) call export_pair_blocks(blocks, pair_blocks)
    end subroutine assemble_landau_projection

    subroutine export_pair_blocks(blocks, pair_blocks)
        complex(dp), allocatable, intent(inout) :: blocks(:, :, :)
        complex(dp), allocatable, intent(out) :: pair_blocks(:, :, :)
        integer :: i, j, k
        ! Pair accumulation stores only the upper triangle. Complete each
        ! retained block exactly as combine_blocks completes the total.
        do k = 1, size(blocks, 3)
            do j = 1, size(blocks, 2)
                blocks(j, j, k) = cmplx(real(blocks(j, j, k), dp), 0.0_dp, dp)
                do i = j + 1, size(blocks, 1)
                    blocks(i, j, k) = conjg(blocks(j, i, k))
                end do
            end do
        end do
        call move_alloc(blocks, pair_blocks)
    end subroutine export_pair_blocks

    subroutine validate_inputs(mass, charge, density, temperature, loglam, ndof, matrix)
        real(dp), intent(in) :: mass(:), charge(:), density(:), temperature, loglam
        integer, intent(in) :: ndof(:)
        complex(dp), intent(in) :: matrix(:, :)
        integer :: n

        n = size(mass)
        if (n == 0) error stop 'Landau species list must not be empty'
        if (size(charge) /= n .or. size(density) /= n .or. size(ndof) /= n) &
            error stop 'Landau species input shapes differ'
        if (.not. all(ieee_is_finite(mass))) error stop 'Landau mass must be finite'
        if (.not. all(ieee_is_finite(charge))) error stop 'Landau charge must be finite'
        if (.not. all(ieee_is_finite(density))) &
            error stop 'Landau density must be finite'
        if (.not. ieee_is_finite(temperature)) error stop 'Landau T must be finite'
        if (.not. ieee_is_finite(loglam)) error stop 'Landau log must be finite'
        if (any(mass <= 0.0_dp)) error stop 'Landau mass must be positive'
        if (any(abs(charge) <= 0.0_dp)) error stop 'Landau charge must be nonzero'
        if (any(density < 0.0_dp)) error stop 'Landau density must be nonnegative'
        if (temperature <= 0.0_dp) error stop 'Landau common T must be positive'
        if (loglam <= 0.0_dp) error stop 'Landau Coulomb log must be positive'
        if (any(ndof < 1)) error stop 'Landau species basis must not be empty'
        if (sum(real(ndof, dp)) > real(huge(n), dp)) &
            error stop 'Landau basis size exceeds integer range'
        if (size(matrix, 1) /= sum(ndof) .or. size(matrix, 2) /= sum(ndof)) &
            error stop 'Landau matrix output shape differs'
    end subroutine validate_inputs

    subroutine initialize_pairs(mass, charge, density, temperature, loglam, ndof, pairs)
        real(dp), intent(in) :: mass(:), charge(:), density(:), temperature, loglam
        integer, intent(in) :: ndof(:)
        type(landau_pair_t), allocatable, intent(out) :: pairs(:)
        integer :: a, b, k, offsets(size(ndof) + 1)

        offsets(1) = 1
        do a = 1, size(ndof)
            offsets(a + 1) = offsets(a) + ndof(a)
        end do
        allocate (pairs(size(mass)*(size(mass) + 1)/2))
        k = 0
        do b = 1, size(mass)
            do a = 1, b
                k = k + 1
                pairs(k)%a = a; pairs(k)%b = b
                pairs(k)%first_a = offsets(a); pairs(k)%last_a = offsets(a + 1) - 1
                pairs(k)%first_b = offsets(b); pairs(k)%last_b = offsets(b + 1) - 1
                call initialize_pair_scales(mass(a), mass(b), charge(a), charge(b), &
                    density(a), density(b), temperature, loglam, pairs(k))
            end do
        end do
    end subroutine initialize_pairs

    subroutine initialize_pair_scales(ma, mb, za, zb, na, nb, temperature, loglam, pair)
        real(dp), intent(in) :: ma, mb, za, zb, na, nb, temperature, loglam
        type(landau_pair_t), intent(inout) :: pair

        pair%ma = ma; pair%mb = mb
        call landau_pair_scales(ma, mb, collision_elementary_charge*temperature, &
            za, zb, collision_elementary_charge, collision_epsilon0, loglam, &
            na, nb, pi, pair%vcm, pair%vrel, pair%fa, pair%fb, pair%prefactor)
        if (.not. all(ieee_is_finite([pair%vcm, pair%vrel, pair%prefactor]))) &
            error stop 'Landau pair scales exceed numerical range'
        if (min(pair%vcm, pair%vrel) <= 0.0_dp) &
            error stop 'Landau pair speeds are below numerical range'
        if (pair%prefactor < 0.0_dp) error stop 'negative Landau pair weight'
        ! The ordered 1/2 sum counts a cross pair twice and a self pair once.
        if (pair%a == pair%b) pair%prefactor = pair%prefactor/2.0_dp
    end subroutine initialize_pair_scales

    subroutine integrate_pair(pair, q, gradient, matrix)
        type(landau_pair_t), intent(in) :: pair
        type(landau_quadrature_t), intent(in) :: q
        procedure(landau_gradient) :: gradient
        complex(dp), intent(out) :: matrix(:, :)
        type(landau_workspace_t) :: workspace
        real(dp) :: radii(2), weight
        integer :: iy, it

        matrix = cmplx(0.0_dp, 0.0_dp, dp)
        if (pair%prefactor <= 0.0_dp) return
        call initialize_workspace(pair, size(matrix, 1), workspace)
        do it = 1, size(q%t)
            do iy = 1, size(q%y)
                radii = [pair%vcm*sqrt(q%y(iy)), pair%vrel*sqrt(q%t(it))]
                weight = pair%prefactor*q%wy(iy)*q%wt(it)
                call integrate_angles(pair, q, radii, weight, gradient, &
                                      workspace, matrix)
            end do
        end do
    end subroutine integrate_pair

    subroutine initialize_workspace(pair, ndof, workspace)
        type(landau_pair_t), intent(in) :: pair
        integer, intent(in) :: ndof
        type(landau_workspace_t), intent(out) :: workspace

        allocate (workspace%ga(3, pair%last_a - pair%first_a + 1))
        allocate (workspace%gb(3, pair%last_b - pair%first_b + 1))
        allocate (workspace%delta(3, ndof), workspace%transverse(3, ndof))
    end subroutine initialize_workspace

    subroutine integrate_angles(pair, q, radii, radial_weight, gradient, work, matrix)
        type(landau_pair_t), intent(in) :: pair
        type(landau_quadrature_t), intent(in) :: q
        real(dp), intent(in) :: radii(2), radial_weight
        procedure(landau_gradient) :: gradient
        type(landau_workspace_t), intent(inout) :: work
        complex(dp), intent(inout) :: matrix(:, :)
        real(dp) :: angles(3), weight
        integer :: ic, ir, ip

        do ip = 1, size(q%phi)
            do ir = 1, size(q%mu)
                do ic = 1, size(q%mu)
                    angles = [q%mu(ic), q%mu(ir), q%phi(ip)]
                    weight = radial_weight*q%wmu(ic)*q%wmu(ir)*q%wphi(ip)
                    call accumulate_node(pair, radii, angles, weight, gradient, &
                                         work, matrix)
                end do
            end do
        end do
    end subroutine integrate_angles

    subroutine accumulate_node(pair, radii, angles, weight, gradient, work, matrix)
        type(landau_pair_t), intent(in) :: pair
        real(dp), intent(in) :: radii(2), angles(3), weight
        procedure(landau_gradient) :: gradient
        type(landau_workspace_t), intent(inout) :: work
        complex(dp), intent(inout) :: matrix(:, :)
        real(dp) :: va(3), vb(3), direction(3)

        call landau_pair_velocities(radii(1), radii(2), pair%fa, pair%fb, &
            angles(1), angles(2), angles(3), va(1), va(2), va(3), &
            vb(1), vb(2), vb(3), direction(1), direction(2), direction(3))
        if (.not. all(ieee_is_finite(va)) .or. .not. all(ieee_is_finite(vb))) &
            error stop 'Landau quadrature velocities exceed numerical range'
        call gradient(pair%a, va, work%ga)
        call gradient(pair%b, vb, work%gb)
        call validate_gradients(work%ga)
        call validate_gradients(work%gb)
        work%delta = cmplx(0.0_dp, 0.0_dp, dp)
        work%delta(:, pair%first_a:pair%last_a) = work%ga/pair%ma
        work%delta(:, pair%first_b:pair%last_b) = &
            work%delta(:, pair%first_b:pair%last_b) - work%gb/pair%mb
        call transverse_gradients(direction, work%delta, work%transverse)
        call accumulate_gram(weight, work%transverse, matrix)
    end subroutine accumulate_node

    subroutine validate_gradients(gradients)
        complex(dp), intent(in) :: gradients(:, :)

        if (.not. all(ieee_is_finite(real(gradients, dp)))) &
            error stop 'Landau basis gradient real part must be finite'
        if (.not. all(ieee_is_finite(aimag(gradients)))) &
            error stop 'Landau basis gradient imaginary part must be finite'
    end subroutine validate_gradients

    subroutine transverse_gradients(direction, delta, transverse)
        real(dp), intent(in) :: direction(3)
        complex(dp), intent(in) :: delta(:, :)
        complex(dp), intent(out) :: transverse(:, :)
        real(dp) :: re(3), im(3)
        integer :: j

        do j = 1, size(delta, 2)
            call landau_transverse(direction(1), direction(2), direction(3), &
                real(delta(1, j), dp), real(delta(2, j), dp), &
                real(delta(3, j), dp), re(1), re(2), re(3))
            call landau_transverse(direction(1), direction(2), direction(3), &
                aimag(delta(1, j)), aimag(delta(2, j)), aimag(delta(3, j)), &
                im(1), im(2), im(3))
            transverse(:, j) = cmplx(re, im, dp)
        end do
    end subroutine transverse_gradients

    subroutine accumulate_gram(weight, gradients, matrix)
        real(dp), intent(in) :: weight
        complex(dp), intent(in) :: gradients(:, :)
        complex(dp), intent(inout) :: matrix(:, :)
        integer :: i, j

        do j = 1, size(matrix, 2)
            do i = 1, j
                matrix(i, j) = matrix(i, j) + weight* &
                    dot_product(gradients(:, i), gradients(:, j))
            end do
        end do
    end subroutine accumulate_gram

    subroutine combine_blocks(blocks, matrix)
        complex(dp), intent(in) :: blocks(:, :, :)
        complex(dp), intent(out) :: matrix(:, :)
        integer :: k, i, j

        matrix = cmplx(0.0_dp, 0.0_dp, dp)
        do k = 1, size(blocks, 3)
            matrix = matrix + blocks(:, :, k)
        end do
        do j = 1, size(matrix, 2)
            matrix(j, j) = cmplx(real(matrix(j, j), dp), 0.0_dp, dp)
            do i = j + 1, size(matrix, 1)
                matrix(i, j) = conjg(matrix(j, i))
            end do
        end do
        call validate_gradients(matrix)
    end subroutine combine_blocks

    include 'generated/landau_projection_leaf.inc'
end module neort_landau_projection
