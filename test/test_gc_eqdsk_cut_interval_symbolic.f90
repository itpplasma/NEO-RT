program test_gc_eqdsk_cut_interval_symbolic
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_eqdsk_cut_numerator_interval_symbolic, only: &
        evaluate_neort_eqdsk_cut_numerator_interval
    use neort_eqdsk_quintic_cell_jet_interval_symbolic, only: &
        evaluate_neort_eqdsk_quintic_cell_jet_interval
    use neort_eqdsk_quintic_profile_jet_interval_symbolic, only: &
        evaluate_neort_eqdsk_quintic_profile_jet_interval
    use neort_gc_outward_interval, only: gc_outward_interval, &
        gc_outward_interval_is_valid, gc_outward_interval_t, operator(/)
    implicit none

    type(gc_outward_interval_t) :: coefficient(0:5,0:5)
    type(gc_outward_interval_t) :: profile_coefficient(0:5)
    type(gc_outward_interval_t) :: jet(10), profile(4), numerator(3)
    type(gc_outward_interval_t) :: delta_r, delta_z, profile_delta
    type(gc_outward_interval_t) :: radius, psi_sep
    real(dp) :: coefficient_value(0:5,0:5), profile_value(0:5)
    real(dp) :: point_jet(10), point_profile(3), point_numerator(3)
    real(dp) :: r_value, z_value, delta_value, radius_value
    integer :: i, j, ir, iz

    do i = 0, 5
        do j = 0, 5
            coefficient_value(i,j) = real((-1)**(i+j)*(i+1)*(2*j+1), dp) &
                /real(64*(i+j+1), dp)
            coefficient(i,j) = gc_outward_interval(coefficient_value(i,j), &
                coefficient_value(i,j))
        end do
        profile_value(i) = real((-1)**i*(i+2), dp)/128.0_dp
        profile_coefficient(i) = gc_outward_interval(profile_value(i), &
            profile_value(i))
    end do
    delta_r = gc_outward_interval(0.125_dp, 0.375_dp)
    delta_z = gc_outward_interval(-0.25_dp, 0.25_dp)
    radius = gc_outward_interval(1.625_dp, 1.875_dp)
    psi_sep = gc_outward_interval(8.0_dp, 8.0_dp)

    call interval_cell_jet(delta_r, delta_z, coefficient, jet)
    call require(all_valid(jet), 'cell-jet interval is invalid')
    profile_delta = jet(1)/psi_sep
    call interval_profile_jet(profile_delta, profile_coefficient, profile)
    call require(all_valid(profile), 'profile-jet interval is invalid')
    call evaluate_neort_eqdsk_cut_numerator_interval(radius, jet(2), jet(3), &
        jet(4), jet(5), jet(6), jet(7), jet(8), jet(9), jet(10), profile(1), &
        profile(2), psi_sep, numerator(1), numerator(2), numerator(3))
    call require(all_valid(numerator), 'Eq.13 numerator interval is invalid')

    do ir = 0, 8
        r_value = delta_r%lo + real(ir,dp)*(delta_r%hi-delta_r%lo)/8.0_dp
        radius_value = 1.5_dp+r_value
        do iz = 0, 8
            z_value = delta_z%lo + real(iz,dp)*(delta_z%hi-delta_z%lo)/8.0_dp
            call independent_cell_jet(r_value, z_value, coefficient_value, &
                point_jet)
            delta_value = point_jet(1)/8.0_dp
            call independent_profile_jet(delta_value, profile_value, point_profile)
            call independent_cut_numerator(radius_value, point_jet, &
                point_profile(1), point_profile(2), 8.0_dp, point_numerator)
            do i = 1, 10
                call require(encloses(jet(i), point_jet(i)), &
                    'cell-jet enclosure missed an independent point')
            end do
            do i = 1, 3
                call require(encloses(profile(i), point_profile(i)), &
                    'profile enclosure missed an independent point')
                call require(encloses(numerator(i), point_numerator(i)), &
                    'Eq.13 enclosure missed an independent point')
            end do
        end do
    end do

    write (*, '(a)') 'test_gc_eqdsk_cut_interval_symbolic OK'

contains

    subroutine interval_cell_jet(r, z, c, value)
        type(gc_outward_interval_t), intent(in) :: r, z, c(0:5,0:5)
        type(gc_outward_interval_t), intent(out) :: value(10)

        call evaluate_neort_eqdsk_quintic_cell_jet_interval(r, z, &
            c(0,0), c(0,1), c(0,2), c(0,3), c(0,4), c(0,5), &
            c(1,0), c(1,1), c(1,2), c(1,3), c(1,4), c(1,5), &
            c(2,0), c(2,1), c(2,2), c(2,3), c(2,4), c(2,5), &
            c(3,0), c(3,1), c(3,2), c(3,3), c(3,4), c(3,5), &
            c(4,0), c(4,1), c(4,2), c(4,3), c(4,4), c(4,5), &
            c(5,0), c(5,1), c(5,2), c(5,3), c(5,4), c(5,5), &
            value(1), value(2), value(3), value(4), value(5), value(6), &
            value(7), value(8), value(9), value(10))
    end subroutine interval_cell_jet

    subroutine interval_profile_jet(delta, c, value)
        type(gc_outward_interval_t), intent(in) :: delta, c(0:5)
        type(gc_outward_interval_t), intent(out) :: value(4)
        type(gc_outward_interval_t) :: vacuum_b, vacuum_r

        vacuum_b = gc_outward_interval(2.0_dp, 2.0_dp)
        vacuum_r = gc_outward_interval(3.0_dp, 3.0_dp)
        call evaluate_neort_eqdsk_quintic_profile_jet_interval(delta, &
            c(0), c(1), c(2), c(3), c(4), c(5), vacuum_b, vacuum_r, &
            value(1), value(2), value(3), value(4))
    end subroutine interval_profile_jet

    subroutine independent_cell_jet(r, z, c, value)
        real(dp), intent(in) :: r, z, c(0:5,0:5)
        real(dp), intent(out) :: value(10)
        integer :: p, q

        value = 0.0_dp
        do p = 0, 5
            do q = 0, 5
                value(1) = value(1)+c(p,q)*r**p*z**q
                if (p >= 1) value(2) = value(2)+p*c(p,q)*r**(p-1)*z**q
                if (q >= 1) value(3) = value(3)+q*c(p,q)*r**p*z**(q-1)
                if (p >= 2) value(4) = value(4)+p*(p-1)*c(p,q) &
                    *r**(p-2)*z**q
                if (p >= 1 .and. q >= 1) value(5) = value(5)+p*q*c(p,q) &
                    *r**(p-1)*z**(q-1)
                if (q >= 2) value(6) = value(6)+q*(q-1)*c(p,q) &
                    *r**p*z**(q-2)
                if (p >= 3) value(7) = value(7)+p*(p-1)*(p-2)*c(p,q) &
                    *r**(p-3)*z**q
                if (p >= 2 .and. q >= 1) value(8) = value(8) &
                    +p*(p-1)*q*c(p,q)*r**(p-2)*z**(q-1)
                if (p >= 1 .and. q >= 2) value(9) = value(9) &
                    +p*q*(q-1)*c(p,q)*r**(p-1)*z**(q-2)
                if (q >= 3) value(10) = value(10)+q*(q-1)*(q-2)*c(p,q) &
                    *r**p*z**(q-3)
            end do
        end do
    end subroutine independent_cell_jet

    subroutine independent_profile_jet(delta, c, value)
        real(dp), intent(in) :: delta, c(0:5)
        real(dp), intent(out) :: value(3)
        integer :: p

        value = 0.0_dp
        do p = 0, 5
            value(1) = value(1)+c(p)*delta**p
            if (p >= 1) value(2) = value(2)+p*c(p)*delta**(p-1)
            if (p >= 2) value(3) = value(3)+p*(p-1)*c(p)*delta**(p-2)
        end do
    end subroutine independent_profile_jet

    subroutine independent_cut_numerator(r, p, f, f_prime, separatrix, value)
        real(dp), intent(in) :: r, p(10), f, f_prime, separatrix
        real(dp), intent(out) :: value(3)
        real(dp) :: pr2, pz2, g, difference, mixed, square_difference

        pr2 = p(2)**2
        pz2 = p(3)**2
        g = f**2+pr2+pz2
        difference = p(6)-p(4)
        mixed = p(2)*p(3)*difference
        square_difference = pr2-pz2
        value(1) = p(3)*g+r*(mixed+p(5)*square_difference)
        value(2) = mixed+p(5)*g+p(5)*square_difference &
            +p(3)*(2.0_dp*f*f_prime*p(2)/separatrix &
            +2.0_dp*p(2)*p(4)+2.0_dp*p(5)*p(3)) &
            +r*(p(2)*p(5)*difference+p(2)*p(3)*(p(9)-p(7)) &
            +p(4)*p(3)*difference+p(8)*square_difference &
            +p(5)*(2.0_dp*p(2)*p(4)-2.0_dp*p(5)*p(3)))
        value(3) = p(3)*(2.0_dp*f*f_prime*p(3)/separatrix &
            +2.0_dp*p(2)*p(5)+2.0_dp*p(3)*p(6))+p(6)*g &
            +r*(p(2)*p(3)*(p(10)-p(8))+p(2)*p(6)*difference &
            +p(5)*p(3)*difference+p(5)*(2.0_dp*p(2)*p(5) &
            -2.0_dp*p(3)*p(6))+p(9)*square_difference)
    end subroutine independent_cut_numerator

    logical function all_valid(values)
        type(gc_outward_interval_t), intent(in) :: values(:)
        integer :: k

        all_valid = .true.
        do k = 1, size(values)
            all_valid = all_valid .and. gc_outward_interval_is_valid(values(k))
        end do
    end function all_valid

    logical function encloses(interval, value)
        type(gc_outward_interval_t), intent(in) :: interval
        real(dp), intent(in) :: value

        encloses = value >= interval%lo .and. value <= interval%hi
    end function encloses

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message

        if (.not. condition) error stop message
    end subroutine require

end program test_gc_eqdsk_cut_interval_symbolic
