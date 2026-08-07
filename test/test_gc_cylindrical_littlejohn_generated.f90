program test_gc_cylindrical_littlejohn_generated
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_cylindrical_littlejohn_symbolic, only: &
        evaluate_neort_cylindrical_littlejohn
    use util_for_test, only: pass_test
    implicit none

    real(dp), parameter :: mass = 2.0_dp, charge = -3.0_dp
    real(dp), parameter :: c_light = 5.0_dp, mu = 0.7_dp, bmod = 4.0_dp
    real(dp), parameter :: electrostatic_potential = 0.2_dp
    real(dp), parameter :: radius = 2.3_dp, p_parallel = -0.8_dp
    real(dp), parameter :: psi = 1.1_dp
    real(dp), parameter :: b(3) = [0.8_dp, 2.4_dp, 3.2_dp]
    real(dp), parameter :: bhat(3) = [0.2_dp, 0.6_dp, 0.8_dp]
    real(dp), parameter :: curl_bhat(3) = [0.03_dp, -0.02_dp, 0.04_dp]
    real(dp), parameter :: grad_b(3) = [0.4_dp, -0.1_dp, 0.2_dp]
    real(dp), parameter :: grad_phi(3) = [0.05_dp, 0.3_dp, -0.2_dp]
    real(dp) :: b_star(3), b_parallel_star, force(3), cross_force(3)
    real(dp) :: velocity(3), derivative(5), measure, d_hamiltonian_dt
    real(dp) :: hamiltonian, canonical_p_phi
    real(dp) :: expected_b_star(3), expected_b_parallel_star
    real(dp) :: expected_force(3), expected_cross(3), expected_velocity(3)
    real(dp) :: expected_hamiltonian, expected_p_phi

    call evaluate_neort_cylindrical_littlejohn(mass, charge, c_light, mu, bmod, &
        electrostatic_potential, radius, p_parallel, psi, b(1), b(2), b(3), &
        bhat(1), bhat(2), bhat(3), curl_bhat(1), curl_bhat(2), curl_bhat(3), &
        grad_b(1), grad_b(2), grad_b(3), grad_phi(1), grad_phi(2), grad_phi(3), &
        b_star(1), b_star(2), b_star(3), b_parallel_star, force(1), force(2), &
        force(3), cross_force(1), cross_force(2), cross_force(3), velocity(1), &
        velocity(2), velocity(3), derivative(1), derivative(2), derivative(3), &
        derivative(4), derivative(5), measure, d_hamiltonian_dt, hamiltonian, &
        canonical_p_phi)

    expected_b_star = b + c_light/charge*p_parallel*curl_bhat
    expected_b_parallel_star = dot_product(bhat, expected_b_star)
    expected_force = mu*grad_b + charge*grad_phi
    expected_cross = [bhat(2)*expected_force(3) - bhat(3)*expected_force(2), &
        bhat(3)*expected_force(1) - bhat(1)*expected_force(3), &
        bhat(1)*expected_force(2) - bhat(2)*expected_force(1)]
    expected_velocity = (p_parallel/mass*expected_b_star + &
        c_light/charge*expected_cross)/expected_b_parallel_star
    expected_hamiltonian = p_parallel**2/(2.0_dp*mass) + mu*bmod &
        + charge*electrostatic_potential
    expected_p_phi = charge/c_light*psi + p_parallel*radius*bhat(2)

    call require_vector('B_star', b_star, expected_b_star)
    call require_close('B_parallel_star', b_parallel_star, &
        expected_b_parallel_star)
    call require_vector('force', force, expected_force)
    call require_vector('b cross force', cross_force, expected_cross)
    call require_vector('velocity', velocity, expected_velocity)
    call require_close('dot R', derivative(1), expected_velocity(1))
    call require_close('dot Z', derivative(2), expected_velocity(3))
    call require_close('dot phi', derivative(3), expected_velocity(2)/radius)
    call require_close('dot p_parallel', derivative(4), &
        -dot_product(expected_b_star, expected_force)/expected_b_parallel_star)
    call require_close('dot mu', derivative(5), 0.0_dp)
    call require_close('R B_parallel_star measure', measure, &
        radius*expected_b_parallel_star)
    call require_close('dH/dt', d_hamiltonian_dt, 0.0_dp)
    call require_close('Hamiltonian', hamiltonian, expected_hamiltonian)
    call require_close('canonical p_phi', canonical_p_phi, expected_p_phi)

    write (*, '(A)') 'test_gc_cylindrical_littlejohn_generated OK'
    call pass_test

contains

    subroutine require_vector(label, actual, expected)
        character(*), intent(in) :: label
        real(dp), intent(in) :: actual(3), expected(3)
        integer :: k

        do k = 1, 3
            call require_close(label, actual(k), expected(k))
        end do
    end subroutine require_vector

    subroutine require_close(label, actual, expected)
        character(*), intent(in) :: label
        real(dp), intent(in) :: actual, expected
        real(dp) :: scale

        scale = max(1.0_dp, max(abs(actual), abs(expected)))
        if (abs(actual - expected) > 2.0e-12_dp*scale) then
            write (*, '(A,2(1X,ES24.16))') trim(label), actual, expected
            error stop 'cylindrical Littlejohn oracle failed'
        end if
    end subroutine require_close

end program test_gc_cylindrical_littlejohn_generated
