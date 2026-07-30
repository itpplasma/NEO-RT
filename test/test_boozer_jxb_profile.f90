program test_boozer_jxb_profile
    ! Independent full-profile oracle for the manufactured vector field in the
    ! analytic circular-torus chartmap.  NEO-RT obtains K=Jcal*J by curling the
    ! NetCDF covariant field and contracts it through the generated coordinate
    ! kernel.  The oracle instead:
    !
    !   1. evaluates the known analytic curl(curl(A)) amplitudes;
    !   2. converts B and J separately to Cartesian vectors using the analytic
    !      circular-torus basis;
    !   3. samples the real fields over their common toroidal phase; and
    !   4. directly evaluates Jcal*e_phi dot (J cross B).
    !
    ! Thus neither the implementation's Ampere derivative, metric cofactors,
    ! coordinate cross-product identity, nor complex 1/2 contraction is reused.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use boozer_chartmap_geometry, only: boozer_chartmap_geometry_t, &
        read_boozer_chartmap_geometry
    use boozer_jxb_profile, only: compute_boozer_jxb_profile
    use boozer_vector_io, only: boozer_vector_harmonics_t, &
        read_boozer_vector_harmonics
    implicit none

    integer, parameter :: phase_count = 256
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: mu0 = 2.5_dp
    real(dp), parameter :: tolerance = 3.0e-11_dp
    character(1024) :: geometry_path, vector_path
    type(boozer_chartmap_geometry_t) :: geometry
    type(boozer_vector_harmonics_t) :: field
    real(dp), allocatable :: actual(:)
    real(dp) :: expected
    integer :: radial_index

    if (command_argument_count() /= 2) then
        error stop "chartmap and vector fixture paths required"
    end if
    call get_command_argument(1, geometry_path)
    call get_command_argument(2, vector_path)
    call read_boozer_chartmap_geometry(trim(geometry_path), geometry)
    call read_boozer_vector_harmonics(trim(vector_path), field)
    call compute_boozer_jxb_profile(geometry, field, mu0, actual)

    do radial_index = 1, field%radial_count
        expected = cartesian_surface_oracle( &
            field, radial_index, geometry%theta)
        if (abs(actual(radial_index) - expected) > &
            tolerance*max(1.0_dp, abs(expected))) then
            print *, "surface", radial_index, "expected", expected, &
                "got", actual(radial_index)
            error stop "Boozer jxb full-profile oracle failed"
        end if
    end do
    print *, "test_boozer_jxb_profile: all checks passed"

contains

    function cartesian_surface_oracle(field_data, radial, theta_grid) result(profile)
        type(boozer_vector_harmonics_t), intent(in) :: field_data
        integer, intent(in) :: radial
        real(dp), intent(in) :: theta_grid(:)
        real(dp) :: profile
        complex(dp) :: b_covariant(3), k_contravariant(3), phase
        complex(dp) :: b_physical(3), j_physical(3)
        real(dp) :: basis(3, 3), metric_diagonal(3), jacobian
        real(dp) :: b_real(3), j_real(3), common_phase
        integer :: mode_index, phase_index, theta_index

        profile = 0.0_dp
        do theta_index = 1, size(theta_grid)
            call circular_torus_basis( &
                field_data%s(radial), theta_grid(theta_index), &
                basis, metric_diagonal, jacobian)
            b_covariant = cmplx(0.0_dp, 0.0_dp, dp)
            k_contravariant = cmplx(0.0_dp, 0.0_dp, dp)
            do mode_index = 1, field_data%mode_count
                phase = exp(cmplx( &
                    0.0_dp, real(field_data%modes(mode_index), dp)* &
                    theta_grid(theta_index), dp))
                b_covariant = b_covariant + &
                    field_data%b_covariant(:, radial, mode_index)*phase
                k_contravariant = k_contravariant + &
                    analytic_weighted_current( &
                    mode_index, field_data%s(radial), &
                    field_data%modes(mode_index), field_data%toroidal_mode)*phase
            end do
            b_physical = matmul(basis, b_covariant/metric_diagonal)
            j_physical = matmul(basis, k_contravariant/jacobian)
            do phase_index = 0, phase_count - 1
                common_phase = 2.0_dp*pi* &
                    real(phase_index, dp)/real(phase_count, dp)
                phase = cmplx(cos(common_phase), sin(common_phase), dp)
                b_real = real(b_physical*phase, dp)
                j_real = real(j_physical*phase, dp)
                profile = profile + jacobian*dot_product( &
                    basis(:, 2), cross(j_real, b_real))
            end do
        end do
        profile = 4.0_dp*pi**2*profile/ &
            real(size(theta_grid)*phase_count, dp)
    end function cartesian_surface_oracle

    subroutine circular_torus_basis(s, theta, basis, metric_diagonal, jacobian)
        real(dp), intent(in) :: s, theta
        real(dp), intent(out) :: basis(3, 3), metric_diagonal(3), jacobian
        real(dp), parameter :: major_radius = 3.0_dp
        real(dp), parameter :: minor_scale = 0.5_dp
        real(dp) :: rho, minor_radius, cylindrical_radius

        rho = sqrt(s)
        minor_radius = minor_scale*rho
        cylindrical_radius = major_radius + minor_radius*cos(theta)
        basis(:, 1) = minor_scale/(2.0_dp*rho)* &
            [cos(theta), 0.0_dp, sin(theta)]
        basis(:, 2) = [0.0_dp, cylindrical_radius, 0.0_dp]
        basis(:, 3) = minor_radius* &
            [-sin(theta), 0.0_dp, cos(theta)]
        metric_diagonal = [ &
            minor_scale**2/(4.0_dp*rho**2), &
            cylindrical_radius**2, minor_radius**2]
        jacobian = minor_scale**2*cylindrical_radius/2.0_dp
    end subroutine circular_torus_basis

    function analytic_weighted_current(mode_index, s, m, n) result(current)
        integer, intent(in) :: mode_index, m, n
        real(dp), intent(in) :: s
        complex(dp) :: current(3)
        complex(dp) :: p, q, p_prime, q_prime, p_second, q_second
        complex(dp), parameter :: imaginary = cmplx(0.0_dp, 1.0_dp, dp)

        call potential_values( &
            mode_index, s, p, q, p_prime, q_prime, p_second, q_second)
        current(1) = imaginary*(real(n, dp)*p_prime + real(m, dp)*q_prime)/mu0
        current(2) = ( &
            real(m*m, dp)*p - real(m*n, dp)*q - p_second)/mu0
        current(3) = ( &
            real(n*n, dp)*q - real(m*n, dp)*p - q_second)/mu0
    end function analytic_weighted_current

    subroutine potential_values( &
            mode_index, s, p, q, p_prime, q_prime, p_second, q_second)
        integer, intent(in) :: mode_index
        real(dp), intent(in) :: s
        complex(dp), intent(out) :: p, q, p_prime, q_prime, p_second, q_second
        complex(dp) :: p0, p1, p2, q0, q1, q2

        if (mode_index == 1) then
            p0 = cmplx(1.0_dp, 0.5_dp, dp)
            p1 = cmplx(2.0_dp, -0.25_dp, dp)
            p2 = cmplx(0.3_dp, 0.2_dp, dp)
            q0 = cmplx(-0.2_dp, 0.3_dp, dp)
            q1 = cmplx(0.4_dp, 0.6_dp, dp)
            q2 = cmplx(-0.1_dp, 0.05_dp, dp)
        else
            p0 = cmplx(-0.4_dp, 0.7_dp, dp)
            p1 = cmplx(0.8_dp, 0.1_dp, dp)
            p2 = cmplx(-0.2_dp, 0.15_dp, dp)
            q0 = cmplx(0.9_dp, -0.2_dp, dp)
            q1 = cmplx(-0.3_dp, 0.4_dp, dp)
            q2 = cmplx(0.25_dp, -0.1_dp, dp)
        end if
        p = p0 + p1*s + p2*s**2
        q = q0 + q1*s + q2*s**2
        p_prime = p1 + 2.0_dp*p2*s
        q_prime = q1 + 2.0_dp*q2*s
        p_second = 2.0_dp*p2
        q_second = 2.0_dp*q2
    end subroutine potential_values

    pure function cross(left, right) result(product)
        real(dp), intent(in) :: left(3), right(3)
        real(dp) :: product(3)

        product = [ &
            left(2)*right(3) - left(3)*right(2), &
            left(3)*right(1) - left(1)*right(3), &
            left(1)*right(2) - left(2)*right(1)]
    end function cross

end program test_boozer_jxb_profile
