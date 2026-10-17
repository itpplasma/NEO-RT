program test_conserving_basis
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_conserving_basis, only: conserving_basis_t, &
                                  initialize_conserving_basis, evaluate_conserving_basis
    implicit none
    type(conserving_basis_t) :: basis
    real(dp) :: values(17), gradient(3, 17), plus(17), minus(17), ignored(3, 17)
    real(dp) :: vp(3), vm(3), velocity(3), dvalues(17), fd(17), error, step
    real(dp) :: b, xi, u, h, mirror(17), scaled(17), discontinuity
    integer :: info, i, j, point, branch

    call initialize_conserving_basis(2, 2, 3.0_dp, 0.8_dp, 1.2_dp, 4.0_dp, &
                                     basis, info)
    if (info /= 0 .or. basis%ndof /= 17) error stop 'basis initialization'
    b = 1.0_dp
    error = 0.0_dp
    do point = 1, 12
        u = 0.13_dp + 0.31_dp*point
        xi = 0.13_dp + 0.047_dp*point
        if (mod(point, 2) == 0) xi = -xi
        velocity = 4*u*[sqrt(1 - xi**2)*0.6_dp, sqrt(1 - xi**2)*0.8_dp, xi]
        call evaluate_conserving_basis(basis, b, velocity, values, gradient, &
                                       branch, dvalues, 0.023_dp, 0.17_dp)
        if (abs(sum(values(:15)) - 1) > 8.0e-16_dp) &
            error stop 'graph partition of unity'
        if (maxval(abs(sum(gradient(:, :15), dim=2))) > 2.0e-15_dp) &
            error stop 'density gradient is nonzero'
        do j = 1, 3
            step = 1.0e-5_dp
            vp = velocity
            vm = velocity
            vp(j) = vp(j) + step
            vm(j) = vm(j) - step
            call evaluate_conserving_basis(basis, b, vp, plus, ignored)
            call evaluate_conserving_basis(basis, b, vm, minus, ignored)
            fd = (plus - minus)/(2*step)
            error = max(error, maxval(abs(fd - gradient(j, :))))
        end do
        ! Independent physical mirror force in SI: dz/dt=-v_perp² dlogB/(2L).
        h = 1.0e-5_dp
        vp = velocity
        vm = velocity
        vp(3) = vp(3) - h*sum(velocity(:2)**2)*100*0.023_dp*0.17_dp/2
        vm(3) = vm(3) + h*sum(velocity(:2)**2)*100*0.023_dp*0.17_dp/2
        call evaluate_conserving_basis(basis, b, vp, plus, ignored)
        call evaluate_conserving_basis(basis, b, vm, minus, ignored)
        if (abs((plus(16) - minus(16))/(2*h) - dvalues(16)) > 1.0e-8_dp) &
            error stop 'native flow Liouville derivative'
    end do
    if (error > 2.0e-9_dp) error stop 'physical velocity finite differences'
    ! Junction traces join trapped/co/counter and trapped mirror legs.
    discontinuity = 0
    do i = 1, 2
        xi = sqrt(1 - b*0.8_dp)
        if (i == 2) xi = -xi
        velocity = 5*[sqrt(1 - xi**2), 0.0_dp, xi]
        vp = velocity
        vm = velocity
        vp(3) = vp(3) + 1.0e-10_dp
        vm(3) = vm(3) - 1.0e-10_dp
        call evaluate_conserving_basis(basis, b, vp, plus, ignored)
        call evaluate_conserving_basis(basis, b, vm, minus, ignored)
        discontinuity = max(discontinuity, maxval(abs(plus - minus)))
    end do
    velocity = [4.2_dp, 0.0_dp, 1.0e-10_dp]
    call evaluate_conserving_basis(basis, b, velocity, mirror, ignored)
    velocity(3) = -velocity(3)
    call evaluate_conserving_basis(basis, b, velocity, values, ignored)
    if (maxval(abs(values(:15) - mirror(:15))) > 1.0e-14_dp) &
        error stop 'trapped mirror discontinuity'
    ! All angular traces collapse to the same origin value. Outer graph trace
    ! continues at every speed beyond umax, while analytic invariants continue.
    velocity = [1.0e-10_dp, 0.0_dp, 0.0_dp]
    call evaluate_conserving_basis(basis, b, velocity, values, ignored)
    if (abs(values(1) - 1) > 3.0e-11_dp) error stop 'origin continuity'
    velocity = [4.0_dp, 0.0_dp, 2.0_dp]
    velocity = velocity/sqrt(sum(velocity**2))*12.0_dp
    call evaluate_conserving_basis(basis, b, velocity, values, ignored)
    call evaluate_conserving_basis(basis, b, velocity*1.7_dp, scaled, ignored)
    if (maxval(abs(values(:15) - scaled(:15))) > 3.0e-15_dp) &
        error stop 'outer speed trace continuation'
    if (discontinuity > 2.0e-10_dp) error stop 'TP junction discontinuity'
    print *, 'velocity gradient FD, junction trace:', error, discontinuity
    print *, 'test_conserving_basis PASSED'
end program test_conserving_basis
