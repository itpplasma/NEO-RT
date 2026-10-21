program test_lorentz_graph
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_lorentz_graph, only: lorentz_branch_t, solve_lorentz_graph
    implicit none

    real(dp), parameter :: lengths(3) = [1.0_dp, 2.0_dp, 3.0_dp]
    real(dp), parameter :: masses(3) = [1.0_dp, 2.0_dp, 4.0_dp]
    real(dp), parameter :: coefficients(3) = [1.0_dp, 2.0_dp, 3.0_dp]
    real(dp), parameter :: amplitudes(3) = [1.0_dp, 1.0_dp, -1.0_dp]
    real(dp), parameter :: gx(3) = [-sqrt(0.6_dp), 0.0_dp, sqrt(0.6_dp)]
    real(dp), parameter :: gw(3) = [5.0_dp/9, 8.0_dp/9, 5.0_dp/9]
    integer :: mode, degree
    real(dp) :: frequency

    call test_legendre()
    call test_junction()
    call test_null_mode()
    call test_invalid_input()
    print *, 'test_lorentz_graph PASSED'

contains

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message
        if (condition) return
        print *, message
        error stop 'Lorentz graph oracle mismatch'
    end subroutine require

    subroutine interval(branch, cells, left, right)
        type(lorentz_branch_t), intent(out) :: branch
        integer, intent(in) :: cells
        real(dp), intent(in) :: left, right
        integer :: i
        allocate(branch%eta(cells + 1))
        do i = 0, cells
            branch%eta(i + 1) = left + (right - left)*real(i, dp)/cells
        end do
    end subroutine interval

    real(dp) function polynomial(ell, x) result(value)
        integer, intent(in) :: ell
        real(dp), intent(in) :: x
        select case (ell)
        case (0)
            value = 1.0_dp
        case (1)
            value = x
        case default
            value = (3.0_dp*x*x - 1.0_dp)/2.0_dp
        end select
    end function polynomial

    subroutine sample(branch, eta, mass, diffusion, omega, nu, drive)
        integer, intent(in) :: branch
        real(dp), intent(in) :: eta
        real(dp), intent(out) :: mass, diffusion, omega, nu
        complex(dp), intent(out) :: drive
        real(dp) :: p
        mass = 1.0_dp
        diffusion = 1.0_dp
        omega = 0.0_dp
        nu = 0.0_dp
        drive = (0.0_dp, 0.0_dp)
        select case (mode)
        case (1)
            ! Physical nu_D=2: eigenvalues nu_D*ell*(ell+1)/2.
            diffusion = 1.0_dp - eta*eta
            nu = 0.5_dp
            omega = frequency
            drive = cmplx(polynomial(degree, eta), 0.0_dp, dp)
        case (2)
            mass = masses(branch)
            diffusion = coefficients(branch)
            omega = 3.0_dp
            nu = 2.0_dp
            p = eta - eta*eta/(2.0_dp*lengths(branch))
            drive = diffusion*amplitudes(branch)/(mass*lengths(branch)) + &
                cmplx(2.0_dp, 3.0_dp, dp)* &
                (cmplx(1.0_dp, 1.0_dp, dp) + amplitudes(branch)*p)
        case (3)
            diffusion = -1.0_dp
        case (4)
            drive = (1.0_dp, 0.0_dp)
        case (8)
            drive = cmplx(1.0e-300_dp, -1.0e-300_dp, dp)
        case (6)
            drive = cmplx(eta, 0.0_dp, dp)
        case (7)
            diffusion = 0.0_dp
        end select
    end subroutine sample

    subroutine test_legendre()
        type(lorentz_branch_t) :: graph(1)
        complex(dp), allocatable :: positive(:)
        complex(dp) :: denominator, exact_work, positive_work
        real(dp) :: error2, norm2, x, h, shape(2), weight, tolerance
        complex(dp) :: value, exact
        integer :: info, element, point

        mode = 1
        frequency = 1.25_dp
        call interval(graph(1), 128, -1.0_dp, 1.0_dp)
        do degree = 0, 2
            denominator = cmplx(0.5_dp + degree*(degree + 1), frequency, dp)
            exact_work = 2.0_dp/((2*degree + 1)*denominator)
            call solve_lorentz_graph(graph, sample, info, parallel_samples=.false.)
            call require(info == 0, 'Legendre solve failed')
            tolerance = 1.0e-10_dp
            if (degree == 2) tolerance = 2.0e-4_dp
            call require(abs(graph(1)%work/exact_work - 1.0_dp) < tolerance, &
                'Legendre work or physical collision normalization')
            error2 = 0.0_dp
            do element = 1, size(graph(1)%eta) - 1
                h = graph(1)%eta(element + 1) - graph(1)%eta(element)
                do point = 1, 3
                    shape = [(1.0_dp - gx(point))/2, (1.0_dp + gx(point))/2]
                    x = sum(graph(1)%eta(element:element + 1)*shape)
                    value = sum(graph(1)%response(element:element + 1)*shape)
                    exact = polynomial(degree, x)/denominator
                    weight = h*gw(point)/2.0_dp
                    error2 = error2 + weight*abs(value - exact)**2
                end do
            end do
            norm2 = 2.0_dp/((2*degree + 1)*abs(denominator)**2)
            call require(sqrt(error2/norm2) < tolerance, 'Legendre response')
            call passivity(graph)
        end do
        degree = 1
        call solve_lorentz_graph(graph, sample, info, parallel_samples=.false.)
        call require(info == 0, 'Positive-frequency solve failed')
        positive = graph(1)%response
        positive_work = graph(1)%work
        frequency = -frequency
        call solve_lorentz_graph(graph, sample, info, parallel_samples=.false.)
        call require(info == 0, 'Negative-frequency solve failed')
        call require(maxval(abs(graph(1)%response - conjg(positive))) < 1.0e-10_dp, &
            'Frequency reversal must conjugate response')
        call require(abs(graph(1)%work - conjg(positive_work)) < 1.0e-10_dp, &
            'Frequency reversal must conjugate work')
    end subroutine test_legendre

    subroutine passivity(graph)
        type(lorentz_branch_t), intent(in) :: graph(:)
        complex(dp) :: work
        real(dp) :: dissipation
        integer :: i
        work = (0.0_dp, 0.0_dp)
        dissipation = 0.0_dp
        do i = 1, size(graph)
            work = work + graph(i)%work
            call require(graph(i)%diffusion_work >= 0.0_dp, 'Negative diffusion work')
            call require(graph(i)%reaction_work >= 0.0_dp, 'Negative reaction work')
            dissipation = dissipation + graph(i)%diffusion_work + graph(i)%reaction_work
        end do
        call require(abs(real(work, dp) - dissipation) < &
            1.0e-10_dp*max(abs(work), 1.0_dp), 'Discrete passivity identity')
    end subroutine passivity

    subroutine test_junction()
        ! Exact junction fluxes (1,2,-3) cancel; exterior fluxes vanish.
        ! Continuous mass norm=98/3, diffusion work=14/3, work=70-98i.
        type(lorentz_branch_t) :: graph(3)
        real(dp) :: errors(3, 2), error2, gradient2, norm2, h, x, weight, shape(2)
        complex(dp) :: value, derivative, exact, work
        integer :: refinement, cells, branch, info, element, point

        mode = 2
        do refinement = 1, 2
            cells = 8*2**refinement
            do branch = 1, 3
                call interval(graph(branch), int(lengths(branch))*cells, &
                    0.0_dp, lengths(branch))
            end do
            call solve_lorentz_graph(graph, sample, info, parallel_samples=.false.)
            call require(info == 0, 'Weighted junction solve failed')
            work = sum(graph%work)
            call passivity(graph)
            error2 = 0.0_dp
            gradient2 = 0.0_dp
            norm2 = 0.0_dp
            do branch = 1, 3
                call require(abs(graph(branch)%response(1) - &
                    graph(1)%response(1)) < 1.0e-12_dp, 'Junction discontinuity')
                do element = 1, size(graph(branch)%eta) - 1
                    h = graph(branch)%eta(element + 1) - graph(branch)%eta(element)
                    derivative = (graph(branch)%response(element + 1) - &
                        graph(branch)%response(element))/h
                    do point = 1, 3
                        shape = [(1.0_dp - gx(point))/2, (1.0_dp + gx(point))/2]
                        x = sum(graph(branch)%eta(element:element + 1)*shape)
                        value = sum(graph(branch)%response(element:element + 1)*shape)
                        exact = cmplx(1.0_dp, 1.0_dp, dp) + amplitudes(branch)* &
                            (x - x*x/(2.0_dp*lengths(branch)))
                        weight = h*gw(point)/2.0_dp
                        error2 = error2 + weight*masses(branch)*abs(value - exact)**2
                        norm2 = norm2 + weight*masses(branch)*abs(value)**2
                        gradient2 = gradient2 + weight*coefficients(branch)* &
                            abs(derivative - amplitudes(branch)* &
                            (1.0_dp - x/lengths(branch)))**2
                    end do
                end do
            end do
            errors(:, refinement) = [sqrt(error2/(98.0_dp/3)), &
                sqrt(gradient2/(14.0_dp/3)), &
                abs(work - cmplx(70.0_dp, -98.0_dp, dp))/abs(work)]
            call require(abs(aimag(work) + 3.0_dp*norm2) < 1.0e-8_dp, &
                'Complex work sign or mass measure')
        end do
        call require(errors(1, 2) < 2.0e-5_dp, 'Weighted graph L2 accuracy')
        call require(errors(3, 2) < 2.0e-6_dp, 'Exact continuous graph work')
        call require(all(errors([1, 3], 1)/errors([1, 3], 2) > 3.6_dp), &
            'Graph L2/work convergence too slow')
        call require(all(errors([1, 3], 1)/errors([1, 3], 2) < 4.4_dp), &
            'Graph L2/work refinement outside asymptotic range')
        call require(errors(2, 1)/errors(2, 2) > 1.9_dp .and. &
            errors(2, 1)/errors(2, 2) < 2.1_dp, 'Graph energy convergence')
    end subroutine test_junction

    subroutine test_null_mode()
        type(lorentz_branch_t) :: graph(1)
        complex(dp) :: moment
        real(dp) :: exact, error, span
        integer :: info, element
        call interval(graph(1), 8, 0.0_dp, 1.0_dp)
        mode = 4
        call solve_lorentz_graph(graph, sample, info, parallel_samples=.false.)
        call require(info == 2, 'Undamped isotropic source must be unsolvable')
        mode = 8
        call solve_lorentz_graph(graph, sample, info, parallel_samples=.false.)
        call require(info == 2, 'Tiny isotropic source must remain incompatible')
        mode = 5
        call solve_lorentz_graph(graph, sample, info, parallel_samples=.false.)
        call require(info == 0, 'Compatible zero source must solve without a floor')
        call require(maxval(abs(graph(1)%response)) == 0.0_dp, 'Zero-source response')
        call require(abs(graph(1)%work) == 0.0_dp, 'Zero-source work')
        mode = 7
        call solve_lorentz_graph(graph, sample, info, parallel_samples=.false.)
        call require(info == -5, 'A zero operator cannot define susceptibility')
        call interval(graph(1), 64, -1.0_dp, 1.0_dp)
        mode = 6
        call solve_lorentz_graph(graph, sample, info, parallel_samples=.false.)
        call require(info == 0, 'Nonzero compatible source must solve')
        moment = (0.0_dp, 0.0_dp)
        error = 0.0_dp
        do element = 1, size(graph(1)%eta) - 1
            span = graph(1)%eta(element + 1) - graph(1)%eta(element)
            moment = moment + span*sum(graph(1)%response(element:element + 1))/2
            exact = graph(1)%eta(element)/2 - graph(1)%eta(element)**3/6
            error = max(error, abs(graph(1)%response(element) - exact))
        end do
        call require(abs(moment) < 1.0e-10_dp, &
            'Constant null mode needs mean-zero gauge')
        call require(error < 1.0e-4_dp, 'Compatible polynomial response')
        call passivity(graph)
    end subroutine test_null_mode

    subroutine test_invalid_input()
        type(lorentz_branch_t) :: graph(1)
        integer :: info
        call interval(graph(1), 2, 0.0_dp, 1.0_dp)
        mode = 3
        call solve_lorentz_graph(graph, sample, info, parallel_samples=.false.)
        call require(info == -4, 'Negative diffusion must be rejected')
        graph(1)%eta = [0.0_dp, 1.0_dp, 0.5_dp]
        mode = 5
        call solve_lorentz_graph(graph, sample, info, parallel_samples=.false.)
        call require(info == -3, 'Nonmonotonic nodes must be rejected')
    end subroutine test_invalid_input

end program test_lorentz_graph
