program test_fp_tensor
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_fp_tensor, only: fp_branch_t, fp_diagnostics_t, solve_fp_tensor
    implicit none
    real(dp), parameter :: lengths(3) = [1.0_dp, 2.0_dp, 3.0_dp]
    real(dp), parameter :: masses(3) = [1.0_dp, 2.0_dp, 4.0_dp]
    real(dp), parameter :: pitch(3) = [1.0_dp, 2.0_dp, 3.0_dp]
    real(dp), parameter :: energy(3) = [1.0_dp, 3.0_dp, 2.0_dp]
    real(dp), parameter :: amplitude(3) = [1.0_dp, 1.0_dp, -1.0_dp]
    real(dp), parameter :: outer = sqrt(5 + 2*sqrt(10.0_dp/7))/3
    real(dp), parameter :: inner = sqrt(5 - 2*sqrt(10.0_dp/7))/3
    real(dp), parameter :: gx(5) = [-outer, -inner, 0.0_dp, inner, outer]
    real(dp), parameter :: gw(5) = [(322 - 13*sqrt(70.0_dp))/900, &
        (322 + 13*sqrt(70.0_dp))/900, 128.0_dp/225, &
        (322 + 13*sqrt(70.0_dp))/900, (322 - 13*sqrt(70.0_dp))/900]
    real(dp), parameter :: diffusivity = 0.7_dp
    real(dp), allocatable :: u(:)
    type(fp_branch_t), allocatable :: graph(:)
    type(fp_diagnostics_t) :: diagnostics
    integer :: mode, info
    real(dp) :: omega, upper, moments(0:10), radial_mean

    call test_manufactured_graph()
    call test_ornstein_uhlenbeck()
    call test_energy_noncommutation()
    call test_density_null()
    call test_physical_cancellation()
    print *, 'test_fp_tensor PASSED'

contains

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message
        if (condition) return
        print *, message
        error stop 'FP tensor independent oracle mismatch'
    end subroutine require

    subroutine mesh(speed_cells, pitch_cells, arms)
        integer, intent(in) :: speed_cells, pitch_cells, arms
        integer :: i, branch
        real(dp) :: left, right
        if (allocated(u)) deallocate(u)
        if (allocated(graph)) deallocate(graph)
        allocate(u(speed_cells + 1), graph(arms))
        do i = 0, speed_cells
            u(i + 1) = upper*real(i, dp)/speed_cells
        end do
        do branch = 1, arms
            left = 0.0_dp
            right = lengths(branch)
            if (arms == 1) then
                left = -1.0_dp
                right = 1.0_dp
            end if
            allocate(graph(branch)%eta(pitch_cells + 1))
            do i = 0, pitch_cells
                graph(branch)%eta(i + 1) = left + (right - left)*real(i, dp)/pitch_cells
            end do
        end do
    end subroutine mesh

    subroutine exact_fields(branch, speed, eta, value, du, deta)
        integer, intent(in) :: branch
        real(dp), intent(in) :: speed, eta
        complex(dp), intent(out) :: value(2), du(2), deta(2)
        real(dp) :: q, dq, r, p
        complex(dp) :: k, denominator
        q = speed**2 - speed**4/(2.0_dp*upper**2)
        dq = 2.0_dp*speed - 2.0_dp*speed**3/upper**2
        select case (mode)
        case (1)
            p = eta - eta**2/(2.0_dp*lengths(branch))
            r = 1.0_dp + amplitude(branch)*p
            value = [cmplx(1.0_dp, 1.0_dp, dp), cmplx(2.0_dp, -1.0_dp, dp)] + &
                q*r*[cmplx(1.0_dp, 0.0_dp, dp), cmplx(0.0_dp, 1.0_dp, dp)]
            if (omega == 0.0_dp) value = (q*r - 28.0_dp/255.0_dp)* &
                [cmplx(1.0_dp, 0.0_dp, dp), cmplx(0.0_dp, 1.0_dp, dp)]
            du = dq*r*[cmplx(1.0_dp, 0.0_dp, dp), cmplx(0.0_dp, 1.0_dp, dp)]
            deta = q*amplitude(branch)*(1.0_dp - eta/lengths(branch))* &
                [cmplx(1.0_dp, 0.0_dp, dp), cmplx(0.0_dp, 1.0_dp, dp)]
        case (2)
            value = [cmplx((speed - speed**3/(3*upper**2))*eta, 0.0_dp, dp), &
                cmplx(q - radial_mean, 0.0_dp, dp)]
            du = [cmplx((1 - speed**2/upper**2)*eta, 0.0_dp, dp), &
                cmplx(dq, 0.0_dp, dp)]
            deta = [cmplx(speed - speed**3/(3*upper**2), 0.0_dp, dp), &
                cmplx(0.0_dp, 0.0_dp, dp)]
        case (3)
            denominator = cmplx(0.0_dp, omega, dp)
            k = sqrt(denominator/diffusivity)
            value(1) = 1.0_dp/denominator
            value(2) = speed**2/denominator + 2*diffusivity/denominator**2 - &
                2*upper*cosh(k*speed)/(denominator*k*sinh(k*upper))
            du = [cmplx(0.0_dp, 0.0_dp, dp), 2*speed/denominator - &
                2*upper*sinh(k*speed)/(denominator*sinh(k*upper))]
            deta = (0.0_dp, 0.0_dp)
        case (7)
            value = [cmplx(q - 7.0_dp/30, 0.0_dp, dp), &
                cmplx(0.0_dp, 0.0_dp, dp)]
            du = [cmplx(dq, 0.0_dp, dp), cmplx(0.0_dp, 0.0_dp, dp)]
            deta = (0.0_dp, 0.0_dp)
        end select
    end subroutine exact_fields

    subroutine sample(branch, speed, eta, weight, p, e, frequency, drive)
        integer, intent(in) :: branch
        real(dp), intent(in) :: speed, eta
        real(dp), intent(out) :: weight, p, e, frequency
        complex(dp), intent(out) :: drive(2)
        complex(dp) :: value(2), du(2), deta(2), differential
        real(dp) :: q, ddq, r, angular, radial
        frequency = omega
        weight = 1.0_dp
        p = 1.0_dp
        e = diffusivity
        drive = (0.0_dp, 0.0_dp)
        select case (mode)
        case (1)
            weight = masses(branch)
            p = pitch(branch)
            e = energy(branch)
            q = speed**2 - speed**4/2.0_dp
            ddq = 2.0_dp - 6.0_dp*speed**2
            r = 1.0_dp + amplitude(branch)*(eta - eta**2/(2*lengths(branch)))
            differential = (-e*ddq*r + p*q*amplitude(branch)/lengths(branch))/weight
            call exact_fields(branch, speed, eta, value, du, deta)
            drive = differential*[cmplx(1.0_dp, 0.0_dp, dp), &
                cmplx(0.0_dp, 1.0_dp, dp)] + cmplx(0.0_dp, omega, dp)*value
        case (2)
            weight = speed**2*exp(-speed**2)
            p = diffusivity*exp(-speed**2)*(1.0_dp - eta**2)
            e = diffusivity*weight
            ! Closed Cartesian OU identities, with explicit reflecting taper.
            angular = diffusivity*((2 + 10/(3*upper**2))*speed - &
                2*speed**3/upper**2)*eta
            radial = diffusivity*(-6 + (4 + 10/upper**2)*speed**2 - &
                4*speed**4/upper**2)
            call exact_fields(branch, speed, eta, value, du, deta)
            drive = [cmplx(angular, 0.0_dp, dp), cmplx(radial, 0.0_dp, dp)] + &
                cmplx(0.0_dp, omega, dp)*value
        case (3)
            drive = [cmplx(1.0_dp, 0.0_dp, dp), cmplx(speed**2, 0.0_dp, dp)]
        case (4)
            drive(1) = (1.0_dp, 0.0_dp)
        case (8)
            drive(1) = cmplx(1.0e-300_dp, -1.0e-300_dp, dp)
        case (6)
            drive = cmplx(0.0_dp, omega, dp)*[1.0_dp, 2.0_dp]
        case (7)
            radial = diffusivity*(-2.0_dp + 6.0_dp*speed**2)
            drive = [cmplx(1.0_dp + 2*radial, 0.0_dp, dp), &
                cmplx(1.0_dp - radial, 0.0_dp, dp)]
        end select
    end subroutine sample

    subroutine solve()
        call solve_fp_tensor(u, graph, sample, info, parallel_samples=.false., &
            diagnostics=diagnostics)
        call require(info == 0, 'Tensor solve failed')
    end subroutine solve

    subroutine check_passivity()
        complex(dp) :: work(2)
        real(dp) :: dissipation(2)
        integer :: branch
        work = (0.0_dp, 0.0_dp)
        dissipation = 0.0_dp
        do branch = 1, size(graph)
            work = work + graph(branch)%self_work
            dissipation = dissipation + graph(branch)%diffusion_work
        end do
        call require(all(dissipation >= 0.0_dp), 'Negative FP dissipation')
        call require(all(abs(real(work, dp) - dissipation) < &
            2.0e-9_dp*max(abs(work), 1.0_dp)), 'FP complex-work dissipation identity')
        call require(diagnostics%constant_residual < 1.0e-11_dp, 'Constant right null')
        call require(diagnostics%particle_residual < 1.0e-11_dp, 'Particle left null')
        call require(diagnostics%symmetry_residual < 1.0e-11_dp, 'Equal-T symmetry')
    end subroutine check_passivity

    subroutine errors(l2, gradient)
        real(dp), intent(out) :: l2, gradient
        complex(dp) :: value(2), du(2), deta(2), fh(2), fu(2), fe(2), node(2)
        real(dp) :: hu, he, speed, eta, su(2), se(2), w, mass, p, e, norm
        integer :: branch, iu, ie, qu, qe, a, b
        l2 = 0.0_dp
        gradient = 0.0_dp
        norm = 0.0_dp
        do branch = 1, size(graph)
            do ie = 1, size(graph(branch)%eta) - 1
                he = graph(branch)%eta(ie + 1) - graph(branch)%eta(ie)
                do iu = 1, size(u) - 1
                    hu = u(iu + 1) - u(iu)
                    do qe = 1, size(gx)
                        se = [(1 - gx(qe))/2, (1 + gx(qe))/2]
                        eta = sum(graph(branch)%eta(ie:ie + 1)*se)
                        do qu = 1, size(gx)
                            su = [(1 - gx(qu))/2, (1 + gx(qu))/2]
                            speed = sum(u(iu:iu + 1)*su)
                            fh = (0.0_dp, 0.0_dp)
                            fu = (0.0_dp, 0.0_dp)
                            fe = (0.0_dp, 0.0_dp)
                            do b = 0, 1
                                do a = 0, 1
                                    node = graph(branch)%response(iu + a, ie + b, :)
                                    fh = fh + su(a + 1)*se(b + 1)*node
                                    fu = fu + real(2*a - 1, dp)/hu*se(b + 1)*node
                                    fe = fe + su(a + 1)*real(2*b - 1, dp)/he*node
                                end do
                            end do
                            call exact_fields(branch, speed, eta, value, du, deta)
                            mass = 1.0_dp
                            p = 1.0_dp
                            e = diffusivity
                            if (mode == 1) then
                                mass = masses(branch)
                                p = pitch(branch)
                                e = energy(branch)
                            else if (mode == 2) then
                                mass = speed**2*exp(-speed**2)
                                p = diffusivity*exp(-speed**2)*(1 - eta**2)
                                e = diffusivity*mass
                            end if
                            w = hu*abs(he)*gw(qu)*gw(qe)/4
                            l2 = l2 + w*mass*sum(abs(fh - value)**2)
                            norm = norm + w*mass*sum(abs(value)**2)
                            gradient = gradient + w*(e*sum(abs(fu - du)**2) + &
                                p*sum(abs(fe - deta)**2))
                        end do
                    end do
                end do
            end do
        end do
        l2 = sqrt(l2/norm)
        gradient = sqrt(gradient)
    end subroutine errors

    subroutine test_manufactured_graph()
        complex(dp), parameter :: overlap = cmplx(113.0_dp/5, -187751.0_dp/3780, dp)
        real(dp), parameter :: norm1 = 147661.0_dp/3780, norm2 = 312217.0_dp/3780
        real(dp), parameter :: dissipation = 61921.0_dp/9450
        real(dp) :: l2(2), gradient(2), work_error(2), ratio
        complex(dp) :: reference(2), work(2), self_work(2)
        integer :: refinement, branch
        mode = 1
        upper = 1.0_dp
        omega = 1.3_dp
        reference = [cmplx(dissipation, -omega*norm1, dp), &
            cmplx(0.0_dp, dissipation, dp) - cmplx(0.0_dp, omega, dp)*overlap]
        do refinement = 1, 2
            call mesh(4*2**refinement, 4*2**refinement, 3)
            call solve()
            call check_passivity()
            work = (0.0_dp, 0.0_dp)
            self_work = (0.0_dp, 0.0_dp)
            do branch = 1, 3
                work = work + graph(branch)%work
                self_work = self_work + graph(branch)%self_work
                call require(maxval(abs(graph(branch)%response(1, :, 1) - &
                    graph(1)%response(1, 1, 1))) < 1.0e-12_dp, 'Anisotropic u=0')
                call require(maxval(abs(graph(branch)%response(:, 1, :) - &
                    graph(1)%response(:, 1, :))) < 1.0e-12_dp, 'Broken graph junction')
            end do
            call errors(l2(refinement), gradient(refinement))
            work_error(refinement) = maxval(abs((work - reference)/reference))
            call require(abs(self_work(2) - cmplx(dissipation, -omega*norm2, dp)) < &
                0.03_dp*abs(self_work(2)), 'Second source self-work moment')
        end do
        ratio = l2(1)/l2(2)
        call require(ratio > 3.5_dp .and. ratio < 4.5_dp, '2D response refinement')
        ratio = gradient(1)/gradient(2)
        call require(ratio > 1.8_dp .and. ratio < 2.2_dp, '2D gradient refinement')
        ratio = work_error(1)/work_error(2)
        call require(ratio > 3.4_dp .and. ratio < 4.6_dp, '2D cross-work refinement')
        call require(l2(2) < 0.004_dp, '2D manufactured response accuracy')
        omega = 0.0_dp
        call solve()
        call check_passivity()
        call errors(l2(1), gradient(1))
        call require(l2(1) < 0.02_dp, 'One global density gauge, not one per speed')
        call require(maxval(abs(diagnostics%response_mean)) < 1.0e-10_dp, &
            'Global weighted-mean gauge')
    end subroutine test_manufactured_graph

    subroutine test_ornstein_uhlenbeck()
        real(dp) :: norm1, norm2, d1, d2, l2, gradient
        complex(dp) :: self_work(2), reference(2)
        integer :: n
        mode = 2
        upper = 3.0_dp
        omega = 1.1_dp
        moments(0) = sqrt(acos(-1.0_dp))*erf(upper)/2
        moments(1) = (1 - exp(-upper**2))/2
        do n = 2, 10
            moments(n) = (n - 1)*moments(n - 2)/2 - &
                upper**(n - 1)*exp(-upper**2)/2
        end do
        radial_mean = (moments(4) - moments(6)/(2*upper**2))/moments(2)
        norm1 = 2.0_dp/3*(moments(4) - 2*moments(6)/(3*upper**2) + &
            moments(8)/(9*upper**4))
        norm2 = 2*(moments(6) - moments(8)/upper**2 + moments(10)/(4*upper**4) - &
            2*radial_mean*(moments(4) - moments(6)/(2*upper**2)) + &
            radial_mean**2*moments(2))
        d1 = 2*diffusivity/3*(moments(2) - 2*moments(4)/upper**2 + &
            moments(6)/upper**4) + 4*diffusivity/3*(moments(2) - &
            2*moments(4)/(3*upper**2) + moments(6)/(9*upper**4))
        d2 = 8*diffusivity*(moments(4) - 2*moments(6)/upper**2 + &
            moments(8)/upper**4)
        reference = [cmplx(d1, -omega*norm1, dp), cmplx(d2, -omega*norm2, dp)]
        call mesh(32, 16, 1)
        call solve()
        call check_passivity()
        self_work = graph(1)%self_work
        call require(maxval(abs((self_work - reference)/reference)) < 0.006_dp, &
            'OU drag, pitch factor two, or reflecting-boundary normalization')
        call require(abs(graph(1)%work(2)) < 1.0e-9_dp*abs(self_work(1)), &
            'Odd/even OU source columns mixed')
        call errors(l2, gradient)
        call require(l2 < 0.015_dp, 'Cartesian OU manufactured response')
    end subroutine test_ornstein_uhlenbeck

    subroutine test_energy_noncommutation()
        real(dp) :: l2, gradient
        complex(dp) :: expected, denominator, k, integral
        mode = 3
        upper = 1.0_dp
        omega = 1.2_dp
        denominator = cmplx(0.0_dp, omega, dp)
        k = sqrt(denominator/diffusivity)
        integral = upper**2*sinh(k*upper)/k - 2*upper*cosh(k*upper)/k**2 + &
            2*sinh(k*upper)/k**3
        expected = 2*(upper**5/(5*denominator) + &
            2*diffusivity*upper**3/(3*denominator**2) - &
            2*upper*integral/(denominator*k*sinh(k*upper)))
        call mesh(24, 4, 1)
        call solve()
        call check_passivity()
        call errors(l2, gradient)
        call require(l2 < 5.0e-4_dp, 'Neumann energy-resolvent response')
        call require(abs(graph(1)%self_work(2)/expected - 1) < 0.004_dp, &
            'Independent u2 source self-work')
        call require(abs(graph(1)%work(2) - 2*upper**3/(3*denominator)) < &
            1.0e-9_dp, 'Cross-column particle moment')
        call require(abs(graph(1)%response(1, 1, 2)) > 0.1_dp, &
            'Energy response illegally computed as u2 times first response')
    end subroutine test_energy_noncommutation

    subroutine test_density_null()
        mode = 6
        omega = 1.0_dp
        upper = 1.0_dp
        call mesh(4, 3, 1)
        call solve()
        call require(maxval(abs(graph(1)%response(:, :, 1) - 1.0_dp)) < &
            1.0e-10_dp, 'Constant Maxwellian right null')
        call require(maxval(abs(graph(1)%response(:, :, 2) - 2.0_dp)) < &
            1.0e-10_dp, 'Second constant Maxwellian right null')
        mode = 4
        omega = 0.0_dp
        call solve_fp_tensor(u, graph, sample, info, parallel_samples=.false.)
        call require(info == 2, 'Global density source must be incompatible')
        mode = 8
        call solve_fp_tensor(u, graph, sample, info, parallel_samples=.false.)
        call require(info == 2, 'Tiny density source must remain incompatible')
        mode = 5
        call solve()
        call require(maxval(abs(graph(1)%response)) == 0.0_dp, 'Zero source, no floor')
    end subroutine test_density_null

    subroutine test_physical_cancellation()
        complex(dp), parameter :: coefficients(2) = &
            [cmplx(1.0_dp/3, 0.0_dp, dp), cmplx(-1.0_dp/3, 0.0_dp, dp)]
        complex(dp) :: original_moments(2)
        real(dp) :: l2, gradient, expected
        logical :: physical
        mode = 7
        upper = 1.0_dp
        omega = 0.0_dp
        call mesh(24, 4, 1)
        call solve_fp_tensor(u, graph, sample, info, parallel_samples=.false.)
        call require(info == 2, 'Each original source has nonzero density moment')
        call solve_fp_tensor(u, graph, sample, info, parallel_samples=.false., &
            diagnostics=diagnostics, physical_coefficients=coefficients, &
            physical_response=physical, source_moments=original_moments)
        call require(info == 0 .and. physical, 'Compatible physical-source combination')
        call require(maxval(abs(original_moments - 2.0_dp)) < 1.0e-11_dp, &
            'Original source moments were overwritten by physical combination')
        call require(maxval(abs(graph(1)%response(:, :, 2))) == 0.0_dp, &
            'Physical response must have one active column')
        call errors(l2, gradient)
        call require(l2 < 0.006_dp, 'Physical combination global-gauge response')
        expected = 64*diffusivity/105
        call require(abs(graph(1)%self_work(1)/expected - 1) < 0.006_dp, &
            'Combined-source dissipative work')
        call require(abs(graph(1)%work(1)/(2*expected) - 1) < 0.006_dp, &
            'Original receiving dual must survive physical combination')
        call check_passivity()
    end subroutine test_physical_cancellation

end program test_fp_tensor
