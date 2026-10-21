module neort_fp_tensor
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_quadrature, only: gauss_legendre
    implicit none
    private

    ! Every arm starts at the same pitch junction. Its outer end is natural.
    ! The speed origin has one value, shared by all pitch nodes and all arms.
    type, public :: fp_branch_t
        real(dp), allocatable :: eta(:)
        complex(dp), allocatable :: response(:, :, :)
        complex(dp) :: work(2) = (0.0_dp, 0.0_dp)
        complex(dp) :: self_work(2) = (0.0_dp, 0.0_dp)
        real(dp) :: diffusion_work(2) = 0.0_dp
    end type fp_branch_t

    type, public :: fp_diagnostics_t
        integer :: n_dof = 0
        integer :: bandwidth = 0
        real(dp) :: volume = 0.0_dp
        real(dp) :: constant_residual = 0.0_dp
        real(dp) :: particle_residual = 0.0_dp
        real(dp) :: symmetry_residual = 0.0_dp
        complex(dp) :: source_moment(2) = (0.0_dp, 0.0_dp)
        complex(dp) :: response_mean(2) = (0.0_dp, 0.0_dp)
    end type fp_diagnostics_t

    abstract interface
        subroutine fp_sample_i(branch, u, eta, w, p, e, omega, drive)
            import :: dp
            integer, intent(in) :: branch
            real(dp), intent(in) :: u, eta
            real(dp), intent(out) :: w, p, e, omega
            complex(dp), intent(out) :: drive(2)
        end subroutine fp_sample_i
    end interface

    type :: fp_cache_t
        real(dp), allocatable :: w(:, :, :, :), p(:, :, :, :)
        real(dp), allocatable :: e(:, :, :, :), omega(:, :, :, :)
        complex(dp), allocatable :: drive(:, :, :, :, :)
        integer :: info = 0
    end type fp_cache_t

    public :: fp_sample_i, solve_fp_tensor, neort_fp_q1_shape
    public :: neort_fp_tensor_coeff, neort_fp_phase_scale

    interface
        subroutine zgbsv(n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info)
            import :: dp
            integer, intent(in) :: n, kl, ku, nrhs, ldab, ldb
            complex(dp), intent(inout) :: ab(ldab, *), b(ldb, *)
            integer, intent(out) :: ipiv(*), info
        end subroutine zgbsv
    end interface

contains

    subroutine solve_fp_tensor(u, branches, sample, info, parallel_samples, &
            diagnostics, physical_coefficients, physical_response, source_moments)
        ! Q1 Galerkin solve of K_pitch + K_energy + i M omega. LAPACK owns
        ! the band factorization. Independent speed slabs are assembled in
        ! two colours; callback sampling defaults to the calling thread.
        !
        ! At identically zero omega, the optional physical coefficients combine
        ! the original two RHS before the global compatibility test. Column 1
        ! then contains that physical response, column 2 is zero, and work(1)
        ! continues to observe the ORIGINAL first drive. No RHS is projected.
        real(dp), intent(in) :: u(:)
        type(fp_branch_t), intent(inout) :: branches(:)
        procedure(fp_sample_i) :: sample
        integer, intent(out) :: info
        logical, intent(in), optional :: parallel_samples
        type(fp_diagnostics_t), intent(out), optional :: diagnostics
        complex(dp), intent(in), optional :: physical_coefficients(2)
        logical, intent(out), optional :: physical_response
        complex(dp), intent(out), optional :: source_moments(2)
        type(fp_cache_t), allocatable :: cache(:)
        type(fp_diagnostics_t) :: diagnostic
        complex(dp), allocatable :: band(:, :), rhs(:, :)
        real(dp), allocatable :: mass(:)
        integer, allocatable :: offsets(:), pivots(:)
        real(dp) :: nodes(3), weights(3), source_scale(2), scale, moment_norm
        integer :: branch, npitch, n, k, colour, slab, column, index
        logical :: parallel, has_reaction, physical, zero_source, incompatible

        info = 0
        physical = .false.
        parallel = .false.
        if (present(parallel_samples)) parallel = parallel_samples
        if (present(physical_response)) physical_response = .false.
        if (present(source_moments)) source_moments = (0.0_dp, 0.0_dp)
        if (present(diagnostics)) diagnostics = diagnostic
        do branch = 1, size(branches)
            if (allocated(branches(branch)%response)) &
                deallocate(branches(branch)%response)
            branches(branch)%work = (0.0_dp, 0.0_dp)
            branches(branch)%self_work = (0.0_dp, 0.0_dp)
            branches(branch)%diffusion_work = 0.0_dp
        end do
        call validate_mesh(u, branches, info)
        if (info /= 0) return
        if (present(physical_coefficients)) then
            if (.not. all(ieee_is_finite(real(physical_coefficients, dp))) .or. &
                    .not. all(ieee_is_finite(aimag(physical_coefficients)))) then
                info = -4
                return
            end if
        end if
        allocate(offsets(size(branches)), cache(size(branches)))
        npitch = 1
        do branch = 1, size(branches)
            offsets(branch) = npitch
            npitch = npitch + size(branches(branch)%eta) - 1
        end do
        n = 1 + (size(u) - 1)*npitch
        k = min(n - 1, 2*npitch - 1)
        allocate(band(3*k + 1, n), rhs(n, 2), mass(n), pivots(n))
        band = (0.0_dp, 0.0_dp)
        rhs = (0.0_dp, 0.0_dp)
        mass = 0.0_dp
        call gauss_legendre(3, nodes, weights)
        !$omp parallel do default(shared) private(branch) if(parallel)
        do branch = 1, size(branches)
            call sample_branch(u, branches(branch)%eta, branch, sample, &
                nodes, cache(branch))
        end do
        !$omp end parallel do
        do branch = 1, size(branches)
            if (cache(branch)%info /= 0) then
                info = cache(branch)%info
                return
            end if
        end do
        do colour = 1, 2
            !$omp parallel do default(shared) private(slab)
            do slab = colour, size(u) - 1, 2
                call assemble_slab(u, branches, cache, offsets, npitch, slab, &
                    nodes, weights, k, band, rhs, mass)
            end do
            !$omp end parallel do
        end do
        if (.not. all(ieee_is_finite(real(band, dp))) .or. &
                .not. all(ieee_is_finite(aimag(band))) .or. &
                .not. all(ieee_is_finite(real(rhs, dp))) .or. &
                .not. all(ieee_is_finite(aimag(rhs))) .or. &
                .not. all(ieee_is_finite(mass))) then
            info = -4
            return
        end if
        call inspect_system(band, rhs, mass, k, diagnostic)
        if (.not. ieee_is_finite(diagnostic%volume) .or. &
                diagnostic%volume <= 0.0_dp .or. &
                .not. all(ieee_is_finite(real(diagnostic%source_moment, dp))) .or. &
                .not. all(ieee_is_finite(aimag(diagnostic%source_moment)))) then
            info = -4
            return
        end if
        if (present(source_moments)) source_moments = diagnostic%source_moment
        has_reaction = .false.
        do branch = 1, size(branches)
            has_reaction = has_reaction .or. any(abs(cache(branch)%omega) > 0.0_dp)
        end do
        if (.not. has_reaction) then
            if (present(physical_coefficients)) then
                rhs(:, 1) = physical_coefficients(1)*rhs(:, 1) + &
                    physical_coefficients(2)*rhs(:, 2)
                rhs(:, 2) = (0.0_dp, 0.0_dp)
                physical = .true.
                if (present(physical_response)) physical_response = .true.
            end if
            if (.not. all(ieee_is_finite(real(rhs, dp))) .or. &
                    .not. all(ieee_is_finite(aimag(rhs)))) then
                info = -4
                if (present(diagnostics)) diagnostics = diagnostic
                return
            end if
            do column = 1, 2
                source_scale(column) = sum(hypot(real(rhs(:, column), dp), &
                    aimag(rhs(:, column))))
                moment_norm = hypot(real(sum(rhs(:, column)), dp), &
                    aimag(sum(rhs(:, column))))
                incompatible = moment_norm > source_scale(column)
                if (moment_norm > 0.0_dp .and. source_scale(column) > 0.0_dp) then
                    if (moment_norm <= source_scale(column)) then
                        ! A larger exponent gap implies a ratio below 2**(-40),
                        ! safely below 1e-11. Otherwise division stays normal,
                        ! including for subnormal source norms.
                        if (exponent(moment_norm) >= &
                                exponent(source_scale(column)) - 40) then
                            incompatible = moment_norm/source_scale(column) > 1.0e-11_dp
                        end if
                    end if
                end if
                if (incompatible) then
                    info = 2
                    if (present(diagnostics)) diagnostics = diagnostic
                    return
                end if
            end do
        end if
        scale = maxval(abs(band))
        if (.not. ieee_is_finite(scale)) then
            info = -4
            if (present(diagnostics)) diagnostics = diagnostic
            return
        end if
        zero_source = all(abs(real(rhs, dp)) <= 0.0_dp) .and. &
            all(abs(aimag(rhs)) <= 0.0_dp)
        if (.not. zero_source) then
            if (scale <= tiny(1.0_dp)) then
                info = -5
                if (present(diagnostics)) diagnostics = diagnostic
                return
            end if
            if (.not. has_reaction) then
                ! Fix the single global constant by its value at the origin.
                ! Subtract its W mean after solving; no damping is introduced.
                do index = 1, min(n, k + 1)
                    band(2*k + index, 1) = (0.0_dp, 0.0_dp)
                    band(2*k + 2 - index, index) = (0.0_dp, 0.0_dp)
                end do
                band(2*k + 1, 1) = cmplx(scale, 0.0_dp, dp)
                rhs(1, :) = (0.0_dp, 0.0_dp)
            end if
            call zgbsv(n, k, k, 2, band, size(band, 1), pivots, rhs, n, info)
            if (info /= 0) then
                if (present(diagnostics)) diagnostics = diagnostic
                return
            end if
            if (.not. all(ieee_is_finite(real(rhs, dp))) .or. &
                    .not. all(ieee_is_finite(aimag(rhs)))) then
                info = -6
                if (present(diagnostics)) diagnostics = diagnostic
                return
            end if
            if (.not. has_reaction) then
                do column = 1, 2
                    rhs(:, column) = rhs(:, column) - &
                        sum(mass*rhs(:, column))/diagnostic%volume
                end do
            end if
        end if
        do column = 1, 2
            diagnostic%response_mean(column) = &
                sum(mass*rhs(:, column))/diagnostic%volume
        end do
        !$omp parallel do default(shared) private(branch)
        do branch = 1, size(branches)
            call recover_branch(u, branches(branch), cache(branch), &
                offsets(branch), npitch, rhs, nodes, weights, physical, &
                physical_coefficients)
        end do
        !$omp end parallel do
        if (present(diagnostics)) diagnostics = diagnostic
    end subroutine solve_fp_tensor

    subroutine validate_mesh(u, branches, info)
        real(dp), intent(in) :: u(:)
        type(fp_branch_t), intent(in) :: branches(:)
        integer, intent(out) :: info
        real(dp) :: direction
        integer :: branch, node

        info = 0
        if (size(branches) < 1) then
            info = -1
            return
        end if
        if (size(u) < 2) then
            info = -2
            return
        end if
        if (.not. all(ieee_is_finite(u))) then
            info = -3
            return
        end if
        if (abs(u(1)) > 0.0_dp) then
            info = -3
            return
        end if
        do node = 2, size(u)
            if (u(node) <= u(node - 1)) then
                info = -3
                return
            end if
        end do
        do branch = 1, size(branches)
            if (.not. allocated(branches(branch)%eta)) then
                info = -2
                return
            end if
            if (size(branches(branch)%eta) < 2) then
                info = -2
                return
            end if
            if (.not. all(ieee_is_finite(branches(branch)%eta))) then
                info = -3
                return
            end if
            if (abs(branches(branch)%eta(1) - branches(1)%eta(1)) > 0.0_dp) then
                info = -3
                return
            end if
            direction = branches(branch)%eta(2) - branches(branch)%eta(1)
            if (abs(direction) <= 0.0_dp) then
                info = -3
                return
            end if
            do node = 2, size(branches(branch)%eta)
                if ((branches(branch)%eta(node) - &
                        branches(branch)%eta(node - 1))*direction <= 0.0_dp) then
                    info = -3
                    return
                end if
            end do
        end do
    end subroutine validate_mesh

    subroutine sample_branch(u, eta, branch, sample, nodes, cache)
        real(dp), intent(in) :: u(:), eta(:), nodes(3)
        integer, intent(in) :: branch
        procedure(fp_sample_i) :: sample
        type(fp_cache_t), intent(out) :: cache
        real(dp) :: uq, eq, w, p, e, omega
        complex(dp) :: drive(2)
        integer :: slab, element, qu, qe, ns, ne

        ns = size(u) - 1
        ne = size(eta) - 1
        allocate(cache%w(3, 3, ne, ns), cache%p(3, 3, ne, ns))
        allocate(cache%e(3, 3, ne, ns), cache%omega(3, 3, ne, ns))
        allocate(cache%drive(2, 3, 3, ne, ns))
        do slab = 1, ns
            do element = 1, ne
                do qe = 1, 3
                    eq = 0.5_dp*((1.0_dp - nodes(qe))*eta(element) + &
                        (1.0_dp + nodes(qe))*eta(element + 1))
                    do qu = 1, 3
                        uq = 0.5_dp*((1.0_dp - nodes(qu))*u(slab) + &
                            (1.0_dp + nodes(qu))*u(slab + 1))
                        call sample(branch, uq, eq, w, p, e, omega, drive)
                        if (.not. all(ieee_is_finite([w, p, e, omega, &
                                real(drive, dp), aimag(drive)]))) then
                            cache%info = -4
                            return
                        end if
                        if (w <= 0.0_dp .or. p < 0.0_dp .or. e < 0.0_dp) then
                            cache%info = -4
                            return
                        end if
                        cache%w(qu, qe, element, slab) = w
                        cache%p(qu, qe, element, slab) = p
                        cache%e(qu, qe, element, slab) = e
                        cache%omega(qu, qe, element, slab) = omega
                        cache%drive(:, qu, qe, element, slab) = drive
                    end do
                end do
            end do
        end do
    end subroutine sample_branch

    subroutine assemble_slab(u, branches, cache, offsets, npitch, slab, &
            nodes, weights, k, band, rhs, mass)
        real(dp), intent(in) :: u(:), nodes(3), weights(3)
        type(fp_branch_t), intent(in) :: branches(:)
        type(fp_cache_t), intent(in) :: cache(:)
        integer, intent(in) :: offsets(:), npitch, slab, k
        complex(dp), intent(inout) :: band(:, :), rhs(:, :)
        real(dp), intent(inout) :: mass(:)
        complex(dp) :: local_matrix(4, 4), local_rhs(4, 2), reaction
        real(dp) :: local_mass(4), shape(4), du(4), de(4), span_u, span_eta
        real(dp) :: measure, w, p, e
        integer :: branch, element, qu, qe, i, j, dof(4), row

        span_u = u(slab + 1) - u(slab)
        do branch = 1, size(branches)
            do element = 1, size(branches(branch)%eta) - 1
                span_eta = branches(branch)%eta(element + 1) - &
                    branches(branch)%eta(element)
                local_matrix = (0.0_dp, 0.0_dp)
                local_rhs = (0.0_dp, 0.0_dp)
                local_mass = 0.0_dp
                do qe = 1, 3
                    do qu = 1, 3
                        call neort_fp_q1_shape(nodes(qu), nodes(qe), &
                            span_u, span_eta, shape, du, de)
                        measure = weights(qu)*weights(qe)* &
                            span_u*abs(span_eta)/4.0_dp
                        w = cache(branch)%w(qu, qe, element, slab)
                        p = cache(branch)%p(qu, qe, element, slab)
                        e = cache(branch)%e(qu, qe, element, slab)
                        reaction = cmplx(0.0_dp, &
                            w*cache(branch)%omega(qu, qe, element, slab), dp)
                        do j = 1, 4
                            do i = 1, 4
                                local_matrix(i, j) = local_matrix(i, j) + &
                                    measure*(p*de(i)*de(j) + e*du(i)*du(j) + &
                                    reaction*shape(i)*shape(j))
                            end do
                        end do
                        do j = 1, 2
                            local_rhs(:, j) = local_rhs(:, j) + &
                                measure*w*shape* &
                                cache(branch)%drive(j, qu, qe, element, slab)
                        end do
                        local_mass = local_mass + measure*w*shape
                    end do
                end do
                call element_dofs(slab, element, offsets(branch), npitch, dof)
                do j = 1, 4
                    do i = 1, 4
                        row = 2*k + 1 + dof(i) - dof(j)
                        band(row, dof(j)) = band(row, dof(j)) + local_matrix(i, j)
                    end do
                end do
                do i = 1, 4
                    mass(dof(i)) = mass(dof(i)) + local_mass(i)
                    rhs(dof(i), :) = rhs(dof(i), :) + local_rhs(i, :)
                end do
            end do
        end do
    end subroutine assemble_slab

    pure subroutine element_dofs(slab, element, offset, npitch, dof)
        integer, intent(in) :: slab, element, offset, npitch
        integer, intent(out) :: dof(4)

        dof(1) = node_dof(slab, element, offset, npitch)
        dof(2) = node_dof(slab + 1, element, offset, npitch)
        dof(3) = node_dof(slab + 1, element + 1, offset, npitch)
        dof(4) = node_dof(slab, element + 1, offset, npitch)
    end subroutine element_dofs

    pure integer function node_dof(speed_node, pitch_node, offset, npitch) &
            result(dof)
        integer, intent(in) :: speed_node, pitch_node, offset, npitch
        integer :: pitch_dof

        if (speed_node == 1) then
            dof = 1
            return
        end if
        pitch_dof = 1
        if (pitch_node > 1) pitch_dof = offset + pitch_node - 1
        dof = 1 + (speed_node - 2)*npitch + pitch_dof
    end function node_dof

    subroutine inspect_system(band, rhs, mass, k, diagnostic)
        complex(dp), intent(in) :: band(:, :), rhs(:, :)
        real(dp), intent(in) :: mass(:)
        integer, intent(in) :: k
        type(fp_diagnostics_t), intent(out) :: diagnostic
        real(dp), allocatable :: row_sum(:), column_sum(:)
        real(dp) :: value, scale, asymmetry
        integer :: i, j, n

        n = size(mass)
        allocate(row_sum(n), column_sum(n))
        row_sum = 0.0_dp
        column_sum = 0.0_dp
        scale = 0.0_dp
        asymmetry = 0.0_dp
        do j = 1, n
            do i = max(1, j - k), min(n, j + k)
                value = real(band(2*k + 1 + i - j, j), dp)
                row_sum(i) = row_sum(i) + value
                column_sum(j) = column_sum(j) + value
                scale = max(scale, abs(value))
                asymmetry = max(asymmetry, &
                    abs(band(2*k + 1 + i - j, j) - &
                    band(2*k + 1 + j - i, i)))
            end do
        end do
        diagnostic%n_dof = n
        diagnostic%bandwidth = k
        diagnostic%volume = sum(mass)
        diagnostic%source_moment = sum(rhs, dim=1)
        diagnostic%constant_residual = maxval(abs(row_sum))/ &
            max(scale, tiny(1.0_dp))
        diagnostic%particle_residual = maxval(abs(column_sum))/ &
            max(scale, tiny(1.0_dp))
        diagnostic%symmetry_residual = asymmetry/ &
            max(maxval(abs(band)), tiny(1.0_dp))
    end subroutine inspect_system

    subroutine recover_branch(u, branch, cache, offset, npitch, solution, &
            nodes, weights, physical, physical_coefficients)
        real(dp), intent(in) :: u(:), nodes(3), weights(3)
        type(fp_branch_t), intent(inout) :: branch
        type(fp_cache_t), intent(in) :: cache
        integer, intent(in) :: offset, npitch
        complex(dp), intent(in) :: solution(:, :)
        logical, intent(in) :: physical
        complex(dp), intent(in), optional :: physical_coefficients(2)
        complex(dp) :: value, derivative_u, derivative_eta, drive(2), self_drive
        complex(dp) :: local_values(4)
        real(dp) :: shape(4), du(4), de(4), span_u, span_eta, measure
        integer :: iu, ie, column, qu, qe, dof(4)

        allocate(branch%response(size(u), size(branch%eta), 2))
        do column = 1, 2
            do ie = 1, size(branch%eta)
                do iu = 1, size(u)
                    branch%response(iu, ie, column) = &
                        solution(node_dof(iu, ie, offset, npitch), column)
                end do
            end do
        end do
        do iu = 1, size(u) - 1
            span_u = u(iu + 1) - u(iu)
            do ie = 1, size(branch%eta) - 1
                span_eta = branch%eta(ie + 1) - branch%eta(ie)
                call element_dofs(iu, ie, offset, npitch, dof)
                do qe = 1, 3
                    do qu = 1, 3
                        call neort_fp_q1_shape(nodes(qu), nodes(qe), &
                            span_u, span_eta, shape, du, de)
                        measure = weights(qu)*weights(qe)* &
                            span_u*abs(span_eta)/4.0_dp
                        drive = cache%drive(:, qu, qe, ie, iu)
                        do column = 1, 2
                            local_values = solution(dof, column)
                            value = sum(shape*local_values)
                            derivative_u = sum(du*local_values)
                            derivative_eta = sum(de*local_values)
                            branch%work(column) = branch%work(column) + &
                                measure*cache%w(qu, qe, ie, iu)*conjg(drive(1))*value
                            self_drive = drive(column)
                            if (physical) then
                                self_drive = (0.0_dp, 0.0_dp)
                                if (column == 1) then
                                    if (present(physical_coefficients)) &
                                        self_drive = sum(physical_coefficients*drive)
                                end if
                            end if
                            branch%self_work(column) = branch%self_work(column) + &
                                measure*cache%w(qu, qe, ie, iu)* &
                                conjg(self_drive)*value
                            branch%diffusion_work(column) = &
                                branch%diffusion_work(column) + measure*( &
                                cache%p(qu, qe, ie, iu)*abs(derivative_eta)**2 + &
                                cache%e(qu, qe, ie, iu)*abs(derivative_u)**2)
                        end do
                    end do
                end do
            end do
        end do
    end subroutine recover_branch

    pure subroutine neort_fp_q1_shape(node_u, node_eta, span_u, span_eta, &
            shape, derivative_u, derivative_eta)
        real(dp), intent(in) :: node_u, node_eta, span_u, span_eta
        real(dp), intent(out) :: shape(4), derivative_u(4), derivative_eta(4)

        call neort_fp_q1_shape_scalar(node_u, node_eta, span_u, span_eta, &
            shape(1), shape(2), shape(3), shape(4), &
            derivative_u(1), derivative_u(2), derivative_u(3), derivative_u(4), &
            derivative_eta(1), derivative_eta(2), derivative_eta(3), &
            derivative_eta(4))
    end subroutine neort_fp_q1_shape

    include 'generated/fp_tensor_leaf.inc'

end module neort_fp_tensor
