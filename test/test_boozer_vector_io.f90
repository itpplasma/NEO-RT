program test_boozer_vector_io
    ! Independent oracle: the fixture defines B=curl(A) in an identity chart
    ! from two quadratic radial vector-potential profiles.  Hence
    !
    ! J_s     = i*(n*p' + m*q')
    ! J_phi   = m^2*p - m*n*q - p''
    ! J_theta = n^2*q - m*n*p - q''
    !
    ! exactly.  This checks NetCDF layout, Fourier signs, complex quadratures,
    ! and the nonuniform three-point radial derivative against curl(curl(A)).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use boozer_vector_io, only: boozer_vector_harmonics_t, &
        compute_boozer_weighted_current, read_boozer_vector_harmonics
    implicit none

    character(1024) :: path
    type(boozer_vector_harmonics_t) :: field
    complex(dp), allocatable :: current(:, :, :)
    real(dp), parameter :: mu0 = 2.5_dp
    real(dp), parameter :: tolerance = 2.0e-13_dp
    integer :: k, radial
    complex(dp) :: p, q, p_prime, q_prime, p_second, q_second
    complex(dp) :: expected(3)

    if (command_argument_count() /= 1) error stop "fixture path required"
    call get_command_argument(1, path)
    call read_boozer_vector_harmonics(trim(path), field)
    call compute_boozer_weighted_current(field, mu0, current)
    if (field%radial_count /= 4 .or. field%mode_count /= 2) then
        error stop "fixture dimensions were not preserved"
    end if
    if (field%toroidal_mode /= -1) error stop "toroidal mode was not preserved"

    do k = 1, field%mode_count
        do radial = 1, field%radial_count
            call potential_values(k, field%s(radial), &
                p, q, p_prime, q_prime, p_second, q_second)
            expected(1) = cmplx(0.0_dp, 1.0_dp, dp)*( &
                real(field%toroidal_mode, dp)*p_prime + &
                real(field%modes(k), dp)*q_prime)/mu0
            expected(2) = (real(field%modes(k)**2, dp)*p - &
                real(field%modes(k)*field%toroidal_mode, dp)*q - p_second)/mu0
            expected(3) = (real(field%toroidal_mode**2, dp)*q - &
                real(field%modes(k)*field%toroidal_mode, dp)*p - q_second)/mu0
            if (maxval(abs(current(:, radial, k) - expected)) > tolerance) then
                print *, "mode/radial: ", k, radial
                print *, "expected: ", expected
                print *, "actual:   ", current(:, radial, k)
                error stop "Boozer vector current oracle failed"
            end if
        end do
    end do
    print *, "test_boozer_vector_io: all checks passed"

contains

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

end program test_boozer_vector_io
