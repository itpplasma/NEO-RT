program test_krook_passing_weight
    ! Independent angular quadrature of the passing Hamiltonian average.
    ! Time evolution of v_parallel loses the energy invariant close to the
    ! separatrix, making this weight too noisy for the Krook pitch quadrature.
    use iso_fortran_env, only: dp => real64
    use fortnum_quadrature, only: gauss_legendre
    use neort_lib, only: neort_init, neort_prepare_splines, neort_setup_at_s
    use neort_transport, only: transport_bounce
    use neort_orbit, only: nvar
    use do_magfie_mod, only: do_magfie, s
    use driftorbit, only: etatp, vth, mth, mph, m0, epsmn, sign_vpar, &
        nopassing, comptorque
    implicit none
    real(dp) :: eta, coarse(2), fine(2), average(nvar), relative_error, worst
    integer :: k, direction, kstate

    call neort_init('driftorbit.in', 'in_file')
    nopassing = .false.
    comptorque = .false.
    call neort_prepare_splines('plasma.in', 'profile.in')
    call neort_setup_at_s(0.5_dp)
    ! No phase: H = 2 - eta*B. Its time average and transit time follow
    ! directly from dt/dtheta = 1/(v*sqrt(1-eta*B)*|h^theta|).
    mth = 0
    mph = 0
    m0 = 0
    epsmn = 1.0_dp
    worst = 0.0_dp
    do k = 1, 3
        eta = etatp*(1.0_dp - real(k, dp)*1.0e-8_dp)
        call angular_reference(512, coarse)
        call angular_reference(1024, fine)
        if (any(abs(fine-coarse) > 2.0e-8_dp*abs(fine))) &
            error stop 'passing angular reference has not converged'
        do direction = -1, 1, 2
            sign_vpar = real(direction, dp)
            call transport_bounce(vth, eta, fine(1), average, kstate)
            if (kstate /= 2) error stop 'passing transport bounce failed'
            relative_error = abs(average(3)-fine(2))/abs(fine(2))
            worst = max(worst, relative_error)
            print *, 'gap, direction, weight relative error:', &
                k, direction, relative_error
            if (relative_error > 2.0e-7_dp) &
                error stop 'passing Hamiltonian average disagrees with angular oracle'
        end do
    end do
    print *, 'passing weight maximum relative error:', worst
    print *, 'test_krook_passing_weight PASSED'
contains
    subroutine angular_reference(panels, result_)
        integer, intent(in) :: panels
        real(dp), intent(out) :: result_(2)
        real(dp), parameter :: pi = acos(-1.0_dp)
        real(dp) :: nodes(32), weights(32), x(3), bmod, sqrtg
        real(dp) :: hder(3), hcovar(3), hctrvr(3), hcurl(3), dt, width
        integer :: panel, node
        call gauss_legendre(32, nodes, weights)
        result_ = 0.0_dp
        width = 2.0_dp*pi/real(panels, dp)
        x(1:2) = [s, 0.0_dp]
        do panel = 1, panels
            do node = 1, size(nodes)
                x(3) = width*(real(panel, dp) - 0.5_dp + nodes(node)/2.0_dp)
                call do_magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
                dt = width/2.0_dp*weights(node)/ &
                    (vth*sqrt(1.0_dp - eta*bmod)*abs(hctrvr(3)))
                result_ = result_ + dt*[1.0_dp, 2.0_dp - eta*bmod]
            end do
        end do
        result_(2) = result_(2)/result_(1)
    end subroutine angular_reference
end program test_krook_passing_weight
