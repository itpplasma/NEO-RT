program test_conserving_solver
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_conserving_solver, only: solve_conserving_system, &
        conserving_solver_diagnostics_t
    implicit none
    real(dp) :: mass(2, 2), streaming(2, 2), precession(2, 2)
    real(dp) :: rotation(2, 2), transformed_mass(2, 2), transformed_d(2, 2)
    real(dp) :: transformed_p(2, 2), dissipation
    complex(dp) :: collision(2, 2), source(2, 2), expected(2), determinant
    complex(dp) :: transformed_c(2, 2), transformed_source(2, 2), actual_work
    complex(dp), allocatable :: coefficient(:), cross(:, :), covariant(:)
    complex(dp), allocatable :: covariant_cross(:, :)
    type(conserving_solver_diagnostics_t) :: diagnostic, covariant_diagnostic
    integer :: info

    ! Two genuinely coupled relaxation coordinates, with different physical
    ! mass measures. The explicit 2x2 inverse is independent of SVD whitening.
    mass = 0.0_dp
    mass(1, 1) = 3.0_dp
    mass(2, 2) = 11.0_dp
    collision = reshape(cmplx([2.0_dp, -2.0_dp, -2.0_dp, 2.0_dp], &
        0.0_dp, dp), [2, 2])
    streaming = 0.0_dp
    precession = 0.7_dp*mass
    source = (0.0_dp, 0.0_dp)
    source(1, 1) = (1.0_dp, -0.3_dp)
    source(2, 2) = (-0.4_dp, 1.7_dp)
    determinant = (2.0_dp + cmplx(0.0_dp, 2.1_dp, dp))* &
        (2.0_dp + cmplx(0.0_dp, 7.7_dp, dp)) - 4.0_dp
    expected(1) = ((2.0_dp + cmplx(0.0_dp, 7.7_dp, dp))* &
        source(1, 1) + 2.0_dp*source(2, 2))/determinant
    expected(2) = (2.0_dp*source(1, 1) + &
        (2.0_dp + cmplx(0.0_dp, 2.1_dp, dp))*source(2, 2))/determinant
    call solve()
    call require(info == 0, 'coupled nonsingular solve')
    call require(diagnostic%rank == 2, 'complete operator rank')
    call require(maxval(abs(coefficient - expected)) < 2.0e-15_dp, &
        'unequal physical mass normalization')
    call require(maxval(abs(sum(cross, dim=2) - coefficient)) < 2.0e-15_dp, &
        'independent sources sum to physical solution')
    call require(abs(cross(2, 1)) > 1.0e-3_dp, 'cross-species response survives')
    call require(all(diagnostic%cross_qualified), 'qualified cross responses')
    dissipation = 2.0_dp*abs(coefficient(1) - coefficient(2))**2
    call require(abs(diagnostic%entropy - dissipation) < 2.0e-15_dp, &
        'actual collision entropy')

    ! Invertible basis mixing changes all Gram/operator/source matrices by
    ! congruence. The physical function, work, and response must be unchanged.
    rotation = reshape([1.0_dp, 0.2_dp, -0.4_dp, 1.3_dp], [2, 2])
    transformed_mass = matmul(transpose(rotation), matmul(mass, rotation))
    transformed_c = matmul(transpose(rotation), matmul(collision, rotation))
    transformed_d = matmul(transpose(rotation), matmul(streaming, rotation))
    transformed_p = matmul(transpose(rotation), matmul(precession, rotation))
    transformed_source = matmul(transpose(rotation), source)
    call solve_conserving_system(transformed_mass, transformed_c, &
        transformed_d, transformed_p, transformed_source, covariant, &
        covariant_cross, covariant_diagnostic, info, 1.0e-12_dp, 1.0e-10_dp)
    call require(info == 0, 'mixed-basis solve')
    call require(maxval(abs(matmul(rotation, covariant) - coefficient)) < &
        2.0e-15_dp, 'physical basis covariance')
    actual_work = dot_product(sum(source, dim=2), coefficient)
    call require(abs(dot_product(sum(transformed_source, dim=2), covariant) - &
        actual_work) < 3.0e-15_dp, 'physical bra covariance')

    ! With no precession the operator's common-flow null is real. Compatible
    ! equal/opposite forcing is solved in the MASS-minimal gauge, with no
    ! alteration of the source or collision matrix.
    precession = 0.0_dp
    source = (0.0_dp, 0.0_dp)
    source(:, 1) = [(1.0_dp, 0.2_dp), (-1.0_dp, -0.2_dp)]
    call solve()
    call require(info == 0, 'compatible physical nullspace')
    call require(diagnostic%rank == 1 .and. diagnostic%nullity == 1, &
        'conserved common-flow null')
    call require(abs(3.0_dp*coefficient(1) + 11.0_dp*coefficient(2)) < &
        2.0e-15_dp, 'physical mass gauge')
    call require(maxval(abs(matmul(collision, coefficient) - &
        sum(source, dim=2))) < 2.0e-15_dp, 'source is not projected')
    source(2, 1) = (0.0_dp, 0.0_dp)
    call solve()
    call require(info == -35, 'incompatible physical source rejected')
    call require(diagnostic%compatibility > 0.1_dp, &
        'left nullspace compatibility measured')

    ! A nonnormal complete operator has a LEFT null distinct from its right
    ! null. Two individually incompatible source columns can cancel to an
    ! admissible physical drive. Only their combined response is qualified.
    collision = (0.0_dp, 0.0_dp)
    collision(1, 1) = (1.0_dp, 0.0_dp)
    collision(1, 2) = (1.0_dp, 0.0_dp)
    source(:, 1) = [(0.0_dp, 0.0_dp), (1.0_dp, 0.0_dp)]
    source(:, 2) = -source(:, 1)
    call solve()
    call require(info == 0, 'combined physical source compatibility')
    call require(.not. any(diagnostic%cross_qualified), &
        'incompatible individual responses remain unqualified')
    call require(maxval(abs(coefficient)) == 0.0_dp, 'zero physical response')

    collision = (0.0_dp, 0.0_dp)
    precession(1, 1) = 2.0_dp
    precession(2, 2) = 3.0e-15_dp
    source = (0.0_dp, 0.0_dp)
    source(1, 1) = (1.0_dp, 0.0_dp)
    call solve()
    call require(info == 0 .and. diagnostic%nullity == 1, &
        'explicit relative rank threshold')
    call require(diagnostic%singular_cutoff > 0.0_dp, 'cutoff reported')
    source(2, 1) = (1.0_dp, 0.0_dp)
    call solve()
    call require(info == -35, 'no inverse floor for unresolved direction')
    mass(2, 2) = -1.0_dp
    call solve()
    call require(info == -32, 'invalid mass basis rejected')
    print '(a)', 'Conserving joint physical-mass/SVD source oracles passed'

contains

    subroutine solve()
        call solve_conserving_system(mass, collision, streaming, precession, &
            source, coefficient, cross, diagnostic, info, 1.0e-12_dp, 1.0e-10_dp)
    end subroutine solve

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message
        if (.not. condition) then
            print '(a)', message
            error stop 'conserving solver oracle failed'
        end if
    end subroutine require
end program test_conserving_solver
