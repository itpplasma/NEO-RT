program test_gc_generated_geometry_interp_endpoints
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_cylindrical_geometry_symbolic, only: &
        evaluate_neort_cylindrical_geometry
    use neort_cylindrical_bilinear_complex_symbolic, only: &
        evaluate_neort_cylindrical_bilinear_complex
    use neort_profile_endpoint_symbolic, only: evaluate_neort_profile_endpoints
    use util_for_test, only: pass_test
    implicit none

    real(dp) :: geometry(19), interp(9), endpoints(8)
    real(dp), parameter :: radius = 2.0_dp
    real(dp), parameter :: b(3) = [1.0_dp, 2.0_dp, 2.0_dp]
    real(dp), parameter :: db(3,3) = reshape([ &
        0.1_dp, 0.2_dp, 0.3_dp, 0.4_dp, 0.5_dp, 0.6_dp, &
        0.7_dp, 0.8_dp, 0.9_dp], [3,3])
    real(dp) :: bmod, grad_b(3), expected_curl(3), bhat(3), dbhat(3,3)
    real(dp), parameter :: u = 0.25_dp, v = 0.75_dp
    real(dp), parameter :: cell_r0 = 1.0_dp, cell_r1 = 3.0_dp
    real(dp), parameter :: cell_z0 = -2.0_dp, cell_z1 = 2.0_dp
    real(dp), parameter :: coordinate_r = 1.5_dp, coordinate_z = 1.0_dp
    real(dp), parameter :: amplitude_scale = 2.5_dp
    real(dp), parameter :: r00 = 1.0_dp, r10 = 2.0_dp, r01 = 3.0_dp, r11 = 4.0_dp
    real(dp), parameter :: i00 = -1.0_dp, i10 = -2.0_dp, i01 = -3.0_dp, i11 = -4.0_dp
    real(dp), parameter :: s0 = 0.2_dp, s1 = 0.8_dp, f0 = 3.0_dp, f1 = 6.0_dp
    real(dp), parameter :: rho = 0.0_dp

    call evaluate_neort_cylindrical_geometry(radius, b(1), b(2), b(3), &
        db(1,1), db(2,1), db(3,1), db(1,2), db(2,2), db(3,2), db(1,3), &
        db(2,3), db(3,3), geometry(1), geometry(2), geometry(3), geometry(4), &
        geometry(5), geometry(6), geometry(7), geometry(8), geometry(9), &
        geometry(10), geometry(11), geometry(12), geometry(13), geometry(14), &
        geometry(15), geometry(16), geometry(17), geometry(18), geometry(19))
    bmod = sqrt(sum(b**2))
    bhat = b/bmod
    grad_b = [dot_product(bhat, db(:,1)), dot_product(bhat, db(:,2)), &
        dot_product(bhat, db(:,3))]
    dbhat = 0.0_dp
    dbhat(:,1) = (db(:,1)-bhat*grad_b(1))/bmod
    dbhat(:,2) = (db(:,2)-bhat*grad_b(2))/bmod
    dbhat(:,3) = (db(:,3)-bhat*grad_b(3))/bmod
    expected_curl = [dbhat(3,2)-dbhat(2,3), dbhat(1,3)-dbhat(3,1), &
        dbhat(2,1)+bhat(2)/radius-dbhat(1,2)]
    call require_close('geometry bmod', geometry(1), bmod)
    call require_vector('geometry grad B', geometry(5:7), grad_b)
    call require_vector('cylindrical curl connection', geometry(17:19), expected_curl)

    call evaluate_neort_cylindrical_bilinear_complex(coordinate_r, &
        coordinate_z, cell_r0, cell_r1, cell_z0, cell_z1, amplitude_scale, &
        r00, i00, r10, i10, r01, i01, r11, i11, interp(1), interp(2), &
        interp(3), interp(4), interp(5), interp(6), interp(7), interp(8), &
        interp(9))
    call require_close('bilinear normalized R', interp(1), u)
    call require_close('bilinear normalized Z', interp(2), v)
    call require_close('bilinear partition', &
        interp(3)+interp(4)+interp(5)+interp(6), 1.0_dp)
    call require_close('bilinear real value', interp(7), &
        amplitude_scale*((1-u)*(1-v)*r00 + u*(1-v)*r10 + &
        (1-u)*v*r01 + u*v*r11))
    call require_close('bilinear imaginary value', interp(8), &
        amplitude_scale*((1-u)*(1-v)*i00 + u*(1-v)*i10 + &
        (1-u)*v*i01 + u*v*i11))
    call require_close('bilinear partition residual', interp(9), 0.0_dp)

    call evaluate_neort_profile_endpoints(s0, s1, f0, f1, rho, endpoints(1), &
        endpoints(2), endpoints(3), endpoints(4), endpoints(5), endpoints(6), &
        endpoints(7), endpoints(8))
    call require_close('profile endpoint s=0', endpoints(1), 2.0_dp)
    call require_close('profile endpoint s=1', endpoints(2), 7.0_dp)
    call require_close('profile endpoint derivative', endpoints(3), 5.0_dp)
    call require_close('profile endpoint slope', endpoints(5), 5.0_dp)
    call require_close('profile axis regularity', endpoints(7), 0.0_dp)
    call require_close('profile axis residual', endpoints(8), 0.0_dp)

    write (*, '(A)') 'test_gc_generated_geometry_interp_endpoints OK'
    call pass_test

contains

    subroutine require_vector(label, actual, expected)
        character(*), intent(in) :: label
        real(dp), intent(in) :: actual(:), expected(:)
        integer :: k
        do k = 1, size(actual)
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
            error stop 'generated geometry/interpolation oracle failed'
        end if
    end subroutine require_close

end program test_gc_generated_geometry_interp_endpoints
