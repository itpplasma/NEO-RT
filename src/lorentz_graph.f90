module neort_lorentz_graph
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_quadrature, only: gauss_legendre
    implicit none
    private

    ! Each branch is oriented from the shared junction to its natural endpoint.
    ! A one-branch graph is an interval with two natural endpoints.
    type, public :: lorentz_branch_t
        real(dp), allocatable :: eta(:)
        complex(dp), allocatable :: response(:)
        complex(dp) :: work = (0.0_dp, 0.0_dp)
        real(dp) :: diffusion_work = 0.0_dp
        real(dp) :: reaction_work = 0.0_dp
    end type lorentz_branch_t

    abstract interface
        subroutine lorentz_sample_i(branch, eta, mass, diffusion, omega, nu, drive)
            import :: dp
            integer, intent(in) :: branch
            real(dp), intent(in) :: eta
            real(dp), intent(out) :: mass, diffusion, omega, nu
            complex(dp), intent(out) :: drive
        end subroutine lorentz_sample_i
    end interface

    type :: branch_system_t
        complex(dp), allocatable :: diagonal(:), off_diagonal(:), rhs(:)
        complex(dp), allocatable :: interior(:, :), drive(:, :)
        real(dp), allocatable :: mass(:, :), diffusion(:, :), nu(:, :)
        real(dp), allocatable :: weight(:, :)
        integer :: info = 0
        logical :: has_reaction = .false.
    end type branch_system_t

    public :: lorentz_sample_i, solve_lorentz_graph, neort_lorentz_diffusion

    interface
        subroutine zgtsv(n, nrhs, dl, d, du, b, ldb, info)
            import :: dp
            integer, intent(in) :: n, nrhs, ldb
            complex(dp), intent(inout) :: dl(*), d(*), du(*), b(ldb, *)
            integer, intent(out) :: info
        end subroutine zgtsv
    end interface

contains

    subroutine solve_lorentz_graph(branches, sample, info, parallel_samples)
        ! Galerkin solve of -(a f')' + M (nu + i omega) f = M H.
        ! Sharing the junction DOF enforces continuity and the weak Kirchhoff
        ! condition. The three arms are eliminated by pivoted LAPACK solves;
        ! the remaining Schur complement has one complex unknown.
        type(lorentz_branch_t), intent(inout) :: branches(:)
        procedure(lorentz_sample_i) :: sample
        integer, intent(out) :: info
        logical, intent(in), optional :: parallel_samples
        type(branch_system_t), allocatable :: systems(:)
        real(dp) :: nodes(3), weights(3)
        complex(dp) :: junction_diagonal, junction_rhs, junction
        complex(dp) :: compatibility, moment, mean_response
        real(dp) :: source_scale, volume, total_volume
        integer :: branch
        logical :: parallel, has_reaction

        info = 0
        parallel = .true.
        if (present(parallel_samples)) parallel = parallel_samples
        if (size(branches) < 1) then
            info = -1
            return
        end if
        allocate(systems(size(branches)))
        call gauss_legendre(3, nodes, weights)
        !$omp parallel do default(shared) private(branch) if(parallel)
        do branch = 1, size(branches)
            call assemble_branch(branches(branch), branch, sample, nodes, weights, &
                systems(branch))
        end do
        !$omp end parallel do
        do branch = 1, size(branches)
            if (systems(branch)%info /= 0) then
                info = systems(branch)%info
                return
            end if
        end do
        junction_diagonal = (0.0_dp, 0.0_dp)
        junction_rhs = (0.0_dp, 0.0_dp)
        compatibility = (0.0_dp, 0.0_dp)
        source_scale = 0.0_dp
        has_reaction = .false.
        do branch = 1, size(branches)
            associate(system => systems(branch))
                compatibility = compatibility + sum(system%rhs)
                source_scale = source_scale + sum(abs(system%rhs))
                has_reaction = has_reaction .or. system%has_reaction
                junction_diagonal = junction_diagonal + system%diagonal(1) + &
                    system%off_diagonal(1)*system%interior(1, 2)
                junction_rhs = junction_rhs + system%rhs(1) - &
                    system%off_diagonal(1)*system%interior(1, 1)
            end associate
        end do
        if (.not. has_reaction) then
            ! The constant null mode has no artificial damping. A compatible
            ! source permits a gauge f(junction)=0; otherwise no stationary
            ! response exists in this projected model.
            if (abs(compatibility) > 1.0e-11_dp*max(source_scale, tiny(1.0_dp))) then
                info = 2
                return
            end if
            junction = (0.0_dp, 0.0_dp)
        else
            if (abs(junction_diagonal) <= tiny(1.0_dp)) then
                info = 1
                return
            end if
            junction = junction_rhs/junction_diagonal
        end if
        !$omp parallel do default(shared) private(branch)
        do branch = 1, size(branches)
            call recover_branch(branches(branch), systems(branch), junction, &
                nodes)
        end do
        !$omp end parallel do
        if (.not. has_reaction) then
            mean_response = (0.0_dp, 0.0_dp)
            total_volume = 0.0_dp
            do branch = 1, size(branches)
                call branch_mass_moments(branches(branch), systems(branch), &
                    nodes, moment, volume)
                mean_response = mean_response + moment
                total_volume = total_volume + volume
            end do
            ! The unique gauge is zero orbit-volume-weighted mean. Updating
            ! the shared unknown shifts every arm by its constant null mode.
            junction = -mean_response/total_volume
            !$omp parallel do default(shared) private(branch)
            do branch = 1, size(branches)
                call recover_branch(branches(branch), systems(branch), junction, &
                    nodes)
            end do
            !$omp end parallel do
        end if
    end subroutine solve_lorentz_graph

    subroutine assemble_branch(branch_data, branch, sample, nodes, weights, system)
        type(lorentz_branch_t), intent(in) :: branch_data
        integer, intent(in) :: branch
        procedure(lorentz_sample_i) :: sample
        real(dp), intent(in) :: nodes(3), weights(3)
        type(branch_system_t), intent(out) :: system
        real(dp) :: span, eta, mass, diffusion, omega, nu, shape(2), gradient(2)
        real(dp) :: local_weight, local_diffusion
        complex(dp) :: drive, reaction, local_matrix(2, 2)
        complex(dp), allocatable :: diagonal(:), lower(:), upper(:)
        integer :: n, element, point, i, j

        if (.not. allocated(branch_data%eta)) then
            system%info = -2
            return
        end if
        n = size(branch_data%eta)
        if (n < 2) then
            system%info = -2
            return
        end if
        if (.not. all(ieee_is_finite(branch_data%eta))) then
            system%info = -3
            return
        end if
        span = branch_data%eta(n) - branch_data%eta(1)
        if (span == 0.0_dp) then
            system%info = -3
            return
        end if
        do element = 1, n - 1
            if ((branch_data%eta(element + 1) - branch_data%eta(element))*span &
                <= 0.0_dp) then
                system%info = -3
                return
            end if
        end do
        allocate(system%diagonal(n), system%off_diagonal(n - 1), system%rhs(n))
        allocate(system%interior(n - 1, 2), system%drive(3, n - 1))
        allocate(system%mass(3, n - 1), system%diffusion(3, n - 1))
        allocate(system%nu(3, n - 1), system%weight(3, n - 1))
        system%diagonal = (0.0_dp, 0.0_dp)
        system%off_diagonal = (0.0_dp, 0.0_dp)
        system%rhs = (0.0_dp, 0.0_dp)
        do element = 1, n - 1
            span = branch_data%eta(element + 1) - branch_data%eta(element)
            local_matrix = (0.0_dp, 0.0_dp)
            do point = 1, 3
                call neort_lorentz_shape(nodes(point), span, shape(1), shape(2), &
                    gradient(1), gradient(2))
                eta = sum(branch_data%eta(element:element + 1)*shape)
                call sample(branch, eta, mass, diffusion, omega, nu, drive)
                if (.not. all(ieee_is_finite([mass, diffusion, omega, nu, &
                    real(drive, dp), aimag(drive)]))) then
                    system%info = -4
                    return
                end if
                if (mass <= 0.0_dp .or. diffusion < 0.0_dp .or. nu < 0.0_dp) then
                    system%info = -4
                    return
                end if
                local_weight = weights(point)*abs(span)/2.0_dp
                system%weight(point, element) = local_weight
                system%mass(point, element) = mass
                system%diffusion(point, element) = diffusion
                system%nu(point, element) = nu
                system%drive(point, element) = drive
                reaction = mass*cmplx(nu, omega, dp)
                if (nu /= 0.0_dp .or. omega /= 0.0_dp) system%has_reaction = .true.
                local_diffusion = local_weight*diffusion
                do j = 1, 2
                    do i = 1, 2
                        local_matrix(i, j) = local_matrix(i, j) + &
                            local_diffusion*gradient(i)*gradient(j) + &
                            local_weight*reaction*shape(i)*shape(j)
                    end do
                end do
                system%rhs(element:element + 1) = &
                    system%rhs(element:element + 1) + &
                    local_weight*mass*drive*shape
            end do
            system%diagonal(element:element + 1) = &
                system%diagonal(element:element + 1) + &
                [local_matrix(1, 1), local_matrix(2, 2)]
            system%off_diagonal(element) = local_matrix(1, 2)
        end do
        if (all(abs(system%diagonal) == 0.0_dp) .and. &
                all(abs(system%off_diagonal) == 0.0_dp)) then
            ! No stationary susceptibility is defined by a zero operator.
            ! Distinguish this exact limit from a failed LAPACK factorization.
            system%info = -5
            return
        end if
        ! Allocate a spare off-diagonal for the n=1 LAPACK interior problem.
        allocate(diagonal(n - 1), lower(max(1, n - 2)), upper(max(1, n - 2)))
        diagonal = system%diagonal(2:)
        lower = (0.0_dp, 0.0_dp)
        upper = (0.0_dp, 0.0_dp)
        if (n > 2) then
            lower(1:n - 2) = system%off_diagonal(2:)
            upper(1:n - 2) = system%off_diagonal(2:)
        end if
        system%interior(:, 1) = system%rhs(2:)
        system%interior(:, 2) = (0.0_dp, 0.0_dp)
        system%interior(1, 2) = -system%off_diagonal(1)
        call zgtsv(n - 1, 2, lower, diagonal, upper, system%interior, n - 1, &
            system%info)
    end subroutine assemble_branch

    subroutine recover_branch(branch, system, junction, nodes)
        type(lorentz_branch_t), intent(inout) :: branch
        type(branch_system_t), intent(in) :: system
        complex(dp), intent(in) :: junction
        real(dp), intent(in) :: nodes(3)
        real(dp) :: span, shape(2), gradient(2), weight
        complex(dp) :: value, derivative
        integer :: n, element, point

        n = size(branch%eta)
        if (allocated(branch%response)) deallocate(branch%response)
        allocate(branch%response(n))
        branch%response(1) = junction
        branch%response(2:) = system%interior(:, 1) + &
            system%interior(:, 2)*junction
        branch%work = (0.0_dp, 0.0_dp)
        branch%diffusion_work = 0.0_dp
        branch%reaction_work = 0.0_dp
        do element = 1, n - 1
            span = branch%eta(element + 1) - branch%eta(element)
            do point = 1, 3
                call neort_lorentz_shape(nodes(point), span, shape(1), shape(2), &
                    gradient(1), gradient(2))
                value = sum(branch%response(element:element + 1)*shape)
                derivative = sum(branch%response(element:element + 1)*gradient)
                weight = system%weight(point, element)
                branch%work = branch%work + weight*system%mass(point, element)* &
                    conjg(system%drive(point, element))*value
                branch%diffusion_work = branch%diffusion_work + &
                    weight*system%diffusion(point, element)*abs(derivative)**2
                branch%reaction_work = branch%reaction_work + &
                    weight*system%mass(point, element)*system%nu(point, element)* &
                    abs(value)**2
            end do
        end do
    end subroutine recover_branch

    subroutine branch_mass_moments(branch, system, nodes, moment, volume)
        type(lorentz_branch_t), intent(in) :: branch
        type(branch_system_t), intent(in) :: system
        real(dp), intent(in) :: nodes(3)
        complex(dp), intent(out) :: moment
        real(dp), intent(out) :: volume
        real(dp) :: shape(2), gradient(2), span, measure
        integer :: element, point

        moment = (0.0_dp, 0.0_dp)
        volume = 0.0_dp
        do element = 1, size(branch%eta) - 1
            span = branch%eta(element + 1) - branch%eta(element)
            do point = 1, 3
                call neort_lorentz_shape(nodes(point), span, shape(1), shape(2), &
                    gradient(1), gradient(2))
                measure = system%weight(point, element)*system%mass(point, element)
                moment = moment + measure* &
                    sum(branch%response(element:element + 1)*shape)
                volume = volume + measure
            end do
        end do
    end subroutine branch_mass_moments

    include 'generated/lorentz_leaf.inc'

end module neort_lorentz_graph
