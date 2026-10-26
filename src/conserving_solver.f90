module neort_conserving_solver
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private

    type, public :: conserving_solver_diagnostics_t
        integer :: rank = 0, nullity = 0
        real(dp) :: mass_condition = 0.0_dp
        real(dp) :: singular_cutoff = 0.0_dp
        real(dp) :: compatibility = 0.0_dp, residual = 0.0_dp
        real(dp) :: entropy = 0.0_dp
        real(dp) :: collision_hermitian_error = 0.0_dp
        real(dp) :: streaming_skew_error = 0.0_dp
        logical, allocatable :: cross_qualified(:)
        real(dp), allocatable :: cross_compatibility(:)
    end type conserving_solver_diagnostics_t

    public :: solve_conserving_system

    interface
        subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
            import :: dp
            character(len=1), intent(in) :: jobz, uplo
            integer, intent(in) :: n, lda, lwork
            real(dp), intent(inout) :: a(lda, *), work(*)
            real(dp), intent(out) :: w(*)
            integer, intent(out) :: info
        end subroutine dsyev
        subroutine zgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, &
                work, lwork, rwork, info)
            import :: dp
            character(len=1), intent(in) :: jobu, jobvt
            integer, intent(in) :: m, n, lda, ldu, ldvt, lwork
            complex(dp), intent(inout) :: a(lda, *), work(*)
            real(dp), intent(out) :: s(*), rwork(*)
            complex(dp), intent(out) :: u(ldu, *), vt(ldvt, *)
            integer, intent(out) :: info
        end subroutine zgesvd
    end interface

contains

    subroutine solve_conserving_system(mass, collision, streaming, precession, &
            source, coefficient, cross_coefficient, diagnostics, info, &
            rank_rtol, compatibility_rtol, singular_values)
        ! Physical Galerkin coefficients: no species row normalization is
        ! discarded. Whiten the positive physical mass Gram, then determine
        ! the nullspace of the ACTUAL complete operator by its SVD. A source
        ! component in its left nullspace is rejected, never subtracted.
        ! The zero-nullspace component is a numerical gauge, not a damping.
        real(dp), intent(in) :: mass(:, :), streaming(:, :), precession(:, :)
        complex(dp), intent(in) :: collision(:, :), source(:, :)
        complex(dp), allocatable, intent(out) :: coefficient(:)
        complex(dp), allocatable, intent(out) :: cross_coefficient(:, :)
        type(conserving_solver_diagnostics_t), intent(out) :: diagnostics
        integer, intent(out) :: info
        real(dp), intent(in) :: rank_rtol, compatibility_rtol
        real(dp), allocatable, intent(out), optional :: singular_values(:)
        real(dp), allocatable :: eigenvectors(:, :), eigenvalues(:)
        real(dp), allocatable :: mass_work(:), inverse_sqrt(:, :), sv(:), rwork(:)
        complex(dp), allocatable :: operator_(:, :), whitened(:, :), saved(:, :)
        complex(dp), allocatable :: left(:, :), right_h(:, :), work(:)
        complex(dp), allocatable :: rhs(:, :), solution(:, :), total(:), y(:)
        real(dp) :: mass_scale, operator_scale, collision_scale, stream_scale
        integer :: n, drivers, i, j, rank, lapack_info

        info = -30
        n = size(mass, 1)
        drivers = size(source, 2)
        if (n < 1 .or. drivers < 1) return
        if (size(mass, 2) /= n .or. size(collision, 1) /= n .or. &
                size(collision, 2) /= n .or. size(streaming, 1) /= n .or. &
                size(streaming, 2) /= n .or. size(precession, 1) /= n .or. &
                size(precession, 2) /= n .or. size(source, 1) /= n) return
        if (.not. all(ieee_is_finite(mass)) .or. &
                .not. all(ieee_is_finite(streaming)) .or. &
                .not. all(ieee_is_finite(precession))) return
        if (.not. finite_complex(collision) .or. .not. finite_complex(source)) return
        if (.not. ieee_is_finite(rank_rtol) .or. &
                .not. ieee_is_finite(compatibility_rtol)) return
        if (rank_rtol <= 0.0_dp .or. rank_rtol >= 1.0_dp .or. &
                compatibility_rtol <= 0.0_dp .or. &
                compatibility_rtol >= 1.0_dp) return
        mass_scale = maxval(abs(mass))
        if (mass_scale <= 0.0_dp) return
        if (maxval(abs(mass - transpose(mass))) > &
                128.0_dp*epsilon(1.0_dp)*mass_scale) return
        allocate(eigenvectors(n, n), eigenvalues(n), mass_work(max(1, 3*n)))
        eigenvectors = mass/mass_scale
        call dsyev('V', 'L', n, eigenvectors, n, eigenvalues, mass_work, &
            size(mass_work), lapack_info)
        info = -31
        if (lapack_info /= 0) return
        ! A nonpositive mass Gram identifies an invalid/redundant basis. It is
        ! not repaired with an eigenvalue floor or by deleting trial functions.
        info = -32
        if (minval(eigenvalues) <= 0.0_dp) return
        diagnostics%mass_condition = maxval(eigenvalues)/minval(eigenvalues)
        allocate(inverse_sqrt(n, n))
        inverse_sqrt = 0.0_dp
        do j = 1, n
            do i = 1, n
                inverse_sqrt(:, i) = inverse_sqrt(:, i) + &
                    eigenvectors(:, j)*eigenvectors(i, j)/ &
                    sqrt(mass_scale*eigenvalues(j))
            end do
        end do
        allocate(operator_(n, n), whitened(n, n), saved(n, n))
        operator_ = collision + cmplx(streaming, precession, dp)
        whitened = matmul(inverse_sqrt, matmul(operator_, inverse_sqrt))
        if (.not. finite_complex(whitened)) return
        saved = whitened
        collision_scale = maxval(abs(collision))
        if (collision_scale > 0.0_dp) diagnostics%collision_hermitian_error = &
            maxval(abs(collision - conjg(transpose(collision))))/collision_scale
        stream_scale = maxval(abs(streaming))
        if (stream_scale > 0.0_dp) diagnostics%streaming_skew_error = &
            maxval(abs(streaming + transpose(streaming)))/stream_scale
        allocate(left(n, n), right_h(n, n), sv(n), rwork(5*n))
        allocate(work(max(1, 8*n*n)))
        call zgesvd('A', 'A', n, n, whitened, n, sv, left, n, right_h, n, &
            work, size(work), rwork, lapack_info)
        info = -33
        if (lapack_info /= 0) return
        if (.not. all(ieee_is_finite(sv))) return
        diagnostics%singular_cutoff = rank_rtol*maxval(sv)
        rank = count(sv > diagnostics%singular_cutoff)
        diagnostics%rank = rank
        diagnostics%nullity = n - rank
        allocate(rhs(n, drivers), solution(n, drivers), total(n), y(n))
        allocate(coefficient(n), cross_coefficient(n, drivers))
        allocate(diagnostics%cross_qualified(drivers))
        allocate(diagnostics%cross_compatibility(drivers))
        rhs = matmul(inverse_sqrt, source)
        total = sum(rhs, dim=2)
        call compatible_solve(total, y, diagnostics%compatibility)
        info = -35
        if (diagnostics%compatibility > compatibility_rtol) return
        coefficient = matmul(inverse_sqrt, y)
        operator_scale = sqrt(sum(abs(saved)**2))
        diagnostics%residual = norm(matmul(saved, y) - total)/ &
            max(norm(total), operator_scale*norm(y), tiny(1.0_dp))
        diagnostics%entropy = real(dot_product(coefficient, &
            matmul(collision, coefficient)), dp)
        do j = 1, drivers
            call compatible_solve(rhs(:, j), solution(:, j), &
                diagnostics%cross_compatibility(j))
            diagnostics%cross_qualified(j) = &
                diagnostics%cross_compatibility(j) <= compatibility_rtol
            cross_coefficient(:, j) = (0.0_dp, 0.0_dp)
            if (diagnostics%cross_qualified(j)) cross_coefficient(:, j) = &
                matmul(inverse_sqrt, solution(:, j))
        end do
        info = -36
        if (.not. finite_complex(reshape(coefficient, [n, 1]))) return
        if (.not. ieee_is_finite(diagnostics%residual)) return
        if (diagnostics%residual > 16.0_dp*compatibility_rtol) return
        if (present(singular_values)) singular_values = sv
        info = 0

    contains

        subroutine compatible_solve(b, x, compatibility)
            complex(dp), intent(in) :: b(:)
            complex(dp), intent(out) :: x(:)
            real(dp), intent(out) :: compatibility
            complex(dp) :: projected(n)
            integer :: k
            projected = matmul(conjg(transpose(left)), b)
            compatibility = 0.0_dp
            if (norm(b) > 0.0_dp) &
                compatibility = norm(projected(rank + 1:n))/norm(b)
            x = (0.0_dp, 0.0_dp)
            do k = 1, rank
                x = x + conjg(right_h(k, :))*projected(k)/sv(k)
            end do
        end subroutine compatible_solve
    end subroutine solve_conserving_system

    pure logical function finite_complex(array)
        complex(dp), intent(in) :: array(:, :)
        finite_complex = all(ieee_is_finite(real(array, dp))) .and. &
            all(ieee_is_finite(aimag(array)))
    end function finite_complex

    pure real(dp) function norm(vector)
        complex(dp), intent(in) :: vector(:)
        norm = sqrt(sum(abs(vector)**2))
    end function norm
end module neort_conserving_solver
