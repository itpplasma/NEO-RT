program test_boozer_ampere
    ! Independent oracle: analytic curl in cylindrical coordinates (r,phi,z).
    ! For B_z = Re[f(r)*exp(i*(n*phi + m*z))],
    !
    !   (curl B)^r   = i*n*f/r,
    !   (curl B)^phi = -f'/r,
    !   (curl B)^z   = 0
    !
    ! in contravariant components.  The cylindrical signed Jacobian is r.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use boozer_ampere, only: boozer_harmonic_current
    implicit none

    real(dp), parameter :: tolerance = 2.0e-14_dp
    real(dp), parameter :: radius = 2.5_dp
    real(dp), parameter :: mu0 = 3.0_dp
    real(dp), parameter :: amplitude = 7.0_dp
    real(dp), parameter :: radial_derivative = -4.0_dp
    integer, parameter :: m = 2, n = -3
    complex(dp) :: b_covariant(3), db_covariant_ds(3), expected(3), current(3)

    b_covariant = cmplx(0.0_dp, 0.0_dp, dp)
    b_covariant(3) = cmplx(amplitude, 0.0_dp, dp)
    db_covariant_ds = cmplx(0.0_dp, 0.0_dp, dp)
    db_covariant_ds(3) = cmplx(radial_derivative, 0.0_dp, dp)
    expected = [ &
        cmplx(0.0_dp, real(n, dp)*amplitude/(mu0*radius), dp), &
        cmplx(-radial_derivative/(mu0*radius), 0.0_dp, dp), &
        cmplx(0.0_dp, 0.0_dp, dp)]

    call boozer_harmonic_current( &
        radius, mu0, m, n, b_covariant, db_covariant_ds, current)
    if (maxval(abs(current - expected)) > tolerance) then
        print *, "expected: ", expected
        print *, "actual:   ", current
        error stop "cylindrical Ampere oracle failed"
    end if

    call test_general_complex_cylindrical_field()
    print *, "test_boozer_ampere: all checks passed"

contains

    subroutine test_general_complex_cylindrical_field()
        ! Independent physical-component formula for
        ! (B_r,B_phi,B_z)=(A,C,D)*exp(i*(n*phi+m*z)):
        !
        ! curl_r   = i*n*D/r - i*m*C
        ! curl_phi = i*m*A - D'
        ! curl_z   = (C+r*C' - i*n*A)/r.
        !
        ! Coordinate inputs are covariant (A,r*C,D), while the phi output is
        ! contravariant curl_phi/r.
        complex(dp), parameter :: imaginary = cmplx(0.0_dp, 1.0_dp, dp)
        complex(dp), parameter :: a = cmplx(2.0_dp, -5.0_dp, dp)
        complex(dp), parameter :: c = cmplx(-1.5_dp, 0.75_dp, dp)
        complex(dp), parameter :: d = cmplx(3.25_dp, 1.25_dp, dp)
        complex(dp), parameter :: c_prime = cmplx(-0.4_dp, 0.2_dp, dp)
        complex(dp), parameter :: d_prime = cmplx(0.6_dp, -0.8_dp, dp)
        complex(dp) :: field(3), derivatives(3), oracle(3), got(3)

        field = [a, radius*c, d]
        derivatives = [cmplx(9.0_dp, -7.0_dp, dp), c + radius*c_prime, d_prime]
        oracle(1) = (imaginary*real(n, dp)*d/radius - &
            imaginary*real(m, dp)*c)/mu0
        oracle(2) = (imaginary*real(m, dp)*a - d_prime)/(mu0*radius)
        oracle(3) = (c + radius*c_prime - imaginary*real(n, dp)*a)/(mu0*radius)
        call boozer_harmonic_current( &
            radius, mu0, m, n, field, derivatives, got)
        if (maxval(abs(got - oracle)) > tolerance) then
            print *, "expected: ", oracle
            print *, "actual:   ", got
            error stop "general cylindrical Ampere oracle failed"
        end if
        call boozer_harmonic_current( &
            -radius, mu0, m, n, field, derivatives, got)
        if (maxval(abs(got + oracle)) > tolerance) then
            error stop "signed-Jacobian orientation oracle failed"
        end if
    end subroutine test_general_complex_cylindrical_field

end program test_boozer_ampere
