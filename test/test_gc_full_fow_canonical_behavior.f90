program test_gc_full_fow_canonical_behavior
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_eqdsk_cut_flux_coordinate_interval_symbolic, only: &
        evaluate_neort_eqdsk_cut_flux_coordinate_interval
    use neort_eqdsk_cut_flux_coordinate_symbolic, only: &
        evaluate_neort_eqdsk_cut_flux_coordinate
    use neort_full_fow_canonical_numerator_interval_symbolic, only: &
        evaluate_neort_full_fow_canonical_numerator_interval
    use neort_full_fow_canonical_numerator_symbolic, only: &
        evaluate_neort_full_fow_canonical_numerator
    use neort_full_fow_canonical_symmetry_symbolic, only: &
        evaluate_neort_full_fow_canonical_symmetry
    use neort_full_fow_canonical_turning_interval_symbolic, only: &
        evaluate_neort_full_fow_canonical_turning_interval
    use neort_full_fow_canonical_turning_symbolic, only: &
        evaluate_neort_full_fow_canonical_turning
    use neort_gc_outward_interval, only: gc_outward_interval, &
        gc_outward_interval_t
    use util_for_test, only: pass_test
    implicit none

    real(dp), parameter :: mass = 1.8_dp, c_light = 2.6_dp
    real(dp), parameter :: charge = -1.2_dp
    real(dp), parameter :: regular_radius = 0.73_dp
    real(dp), parameter :: turning_radius = 1.1_dp
    real(dp), parameter :: turning_side = 1.0_dp
    real(dp), parameter :: turning_slope = 1.25_dp
    real(dp), parameter :: tolerance = 2.0e-6_dp
    real(dp) :: sigma, n_value, n_r_value, n_rr_value
    real(dp) :: n_fd(3), n_r_fd(3), n_rr_fd(3)
    real(dp) :: h_values(3), h, scale
    real(dp) :: n_turning(2), dn_dy_turning(2)
    real(dp) :: y_values(3), regular_turning_n, turning_slope_fd
    real(dp) :: symmetry(21), n_reference, n_r_reference, n_rr_reference
    integer :: i, branch

    h_values = [0.16_dp, 0.08_dp, 0.04_dp]
    do branch = 1, 2
        sigma = real(2*branch-3, dp)
        call evaluate_regular_numerator(sigma, n_value, n_r_value, &
            n_rr_value)
        do i = 1, 3
            h = h_values(i)
            n_fd(i) = direct_n(regular_radius, h/8.0_dp, sigma)
            n_r_fd(i) = direct_n_first(regular_radius, h, sigma)
            n_rr_fd(i) = direct_n_second(regular_radius, h, sigma)
            call require_close('N versus direct psi-star oracle', n_value, &
                n_fd(i), 5.0e-5_dp)
            call require_close('N_R versus direct oracle', n_r_value, &
                n_r_fd(i), 5.0e-4_dp)
            call require_close('N_RR versus direct oracle', n_rr_value, &
                n_rr_fd(i), 5.0e-3_dp)
        end do
        call require_convergence('N finite difference', n_fd)
        call require_convergence('N_R finite difference', n_r_fd)
        call require_convergence('N_RR finite difference', n_rr_fd)
    end do

    do branch = 1, 2
        sigma = real(2*branch-3, dp)
        call evaluate_neort_full_fow_canonical_turning(turning_side, &
            turning_slope, turning_g(turning_radius), &
            turning_g_first(turning_radius), regular_psi_first(turning_radius), &
            mass, c_light, sigma, charge, n_turning(branch), &
            dn_dy_turning(branch))
        y_values = [0.02_dp, 0.01_dp, 0.005_dp]
        do i = 1, 3
            regular_turning_n = direct_turning_n(y_values(i), sigma)
            turning_slope_fd = (regular_turning_n-n_turning(branch))/ &
                y_values(i)
            call require_close('regular N reaches turning limit', &
                regular_turning_n, n_turning(branch), 2.0e-2_dp)
            call require_close('regular N reaches y derivative limit', &
                turning_slope_fd, dn_dy_turning(branch), 2.0e-2_dp)
        end do
    end do
    call require_close('turning y derivative is sigma independent', &
        dn_dy_turning(1), dn_dy_turning(2), tolerance)
    call require_close('turning y derivative direct oracle', &
        dn_dy_turning(1), 2.0_dp*regular_psi_first(turning_radius)*sqrt( &
        turning_slope), tolerance)
    call require_close('turning function vanishes at root', &
        turning_F(turning_radius), 0.0_dp, tolerance)
    call require_nonzero('turning samples are on allowed side', &
        turning_F(turning_radius+turning_side*y_values(3)**2))
    call require_nonzero('turning N retains sigma sign', &
        n_turning(1)-n_turning(2))

    call check_canonical_interval_enclosure
    call check_turning_interval_enclosure
    call check_cut_flux_interval_enclosure

    call evaluate_regular_numerator(1.0_dp, n_reference, n_r_reference, &
        n_rr_reference)
    call evaluate_neort_full_fow_canonical_symmetry( &
        regular_F(regular_radius), regular_F_first(regular_radius), &
        regular_F_second(regular_radius), regular_F_third(regular_radius), &
        regular_G(regular_radius), regular_G_first(regular_radius), &
        regular_G_second(regular_radius), regular_G_third(regular_radius), &
        regular_psi_first(regular_radius), &
        regular_psi_second(regular_radius), regular_psi_third(regular_radius), &
        mass, c_light, 1.0_dp, charge, 2.0_dp, -1.5_dp, 1.7_dp, -0.6_dp, &
        0.83_dp, -0.37_dp, 1.4_dp, 0.61_dp, -0.29_dp, 0.8_dp, -1.1_dp, &
        0.47_dp, -0.23_dp, 0.8_dp, -1.1_dp, &
        symmetry(1), symmetry(2), symmetry(3), symmetry(4), symmetry(5), &
        symmetry(6), symmetry(7), symmetry(8), symmetry(9), symmetry(10), &
        symmetry(11), symmetry(12), symmetry(13), symmetry(14), symmetry(15), &
        symmetry(16), symmetry(17), symmetry(18), symmetry(19), symmetry(20), &
        symmetry(21))
    call require_close('coordinate reversal N', symmetry(1), -n_reference, &
        tolerance)
    call require_close('coordinate reversal N_R', symmetry(2), n_r_reference, &
        tolerance)
    call require_close('coordinate reversal N_RR', symmetry(3), -n_rr_reference, &
        tolerance)
    call require_close('time phase reversal', symmetry(4), &
        -0.83_dp*symmetry(5), tolerance)
    call require_close('time resonance reversal', symmetry(5), &
        -(2.0_dp*1.7_dp-1.5_dp*(-0.6_dp)), tolerance)
    call require_close('Fourier pair resonance reversal', symmetry(6), &
        symmetry(5), tolerance)
    call require_close('combined resonance', symmetry(7), &
        -symmetry(5), tolerance)
    call require_close('toroidal phase relabelling', symmetry(8), &
        -1.5_dp*(-0.37_dp)-0.83_dp*symmetry(7), tolerance)
    call require_close('relabelled P_phi', symmetry(9), &
        -regular_G(regular_radius)*0.61_dp-charge*1.4_dp/c_light, tolerance)
    call require_nonzero('sigma reversal changes canonical offset', symmetry(10))
    do i = 11, 20
        call require_close('symmetry residual', symmetry(i), 0.0_dp, tolerance)
    end do
    call require_close('toroidal torque component changes sign', symmetry(21), &
        -2.0_dp*0.8_dp, tolerance)

    write (*, '(A)') 'test_gc_full_fow_canonical_behavior OK'
    call pass_test

contains

    subroutine evaluate_regular_numerator(sigma_value, n_out, nr_out, nrr_out)
        real(dp), intent(in) :: sigma_value
        real(dp), intent(out) :: n_out, nr_out, nrr_out

        call evaluate_neort_full_fow_canonical_numerator( &
            regular_F(regular_radius), regular_F_first(regular_radius), &
            regular_F_second(regular_radius), regular_F_third(regular_radius), &
            regular_G(regular_radius), regular_G_first(regular_radius), &
            regular_G_second(regular_radius), regular_G_third(regular_radius), &
            regular_psi_first(regular_radius), &
            regular_psi_second(regular_radius), regular_psi_third(regular_radius), &
            mass, c_light, sigma_value, charge, n_out, nr_out, nrr_out)
    end subroutine evaluate_regular_numerator

    pure real(dp) function regular_F(radius)
        real(dp), intent(in) :: radius

        regular_F = 1.5_dp + 0.2_dp*sin(0.7_dp*radius) + &
            0.05_dp*cos(1.3_dp*radius)
    end function regular_F

    pure real(dp) function regular_F_first(radius)
        real(dp), intent(in) :: radius

        regular_F_first = 0.14_dp*cos(0.7_dp*radius) - &
            0.065_dp*sin(1.3_dp*radius)
    end function regular_F_first

    pure real(dp) function regular_F_second(radius)
        real(dp), intent(in) :: radius

        regular_F_second = -0.098_dp*sin(0.7_dp*radius) - &
            0.0845_dp*cos(1.3_dp*radius)
    end function regular_F_second

    pure real(dp) function regular_F_third(radius)
        real(dp), intent(in) :: radius

        regular_F_third = -0.0686_dp*cos(0.7_dp*radius) + &
            0.10985_dp*sin(1.3_dp*radius)
    end function regular_F_third

    pure real(dp) function regular_G(radius)
        real(dp), intent(in) :: radius

        regular_G = 0.9_dp + 0.1_dp*cos(0.5_dp*radius) + &
            0.03_dp*sin(1.1_dp*radius)
    end function regular_G

    pure real(dp) function regular_G_first(radius)
        real(dp), intent(in) :: radius

        regular_G_first = -0.05_dp*sin(0.5_dp*radius) + &
            0.033_dp*cos(1.1_dp*radius)
    end function regular_G_first

    pure real(dp) function regular_G_second(radius)
        real(dp), intent(in) :: radius

        regular_G_second = -0.025_dp*cos(0.5_dp*radius) - &
            0.0363_dp*sin(1.1_dp*radius)
    end function regular_G_second

    pure real(dp) function regular_G_third(radius)
        real(dp), intent(in) :: radius

        regular_G_third = 0.0125_dp*sin(0.5_dp*radius) - &
            0.03993_dp*cos(1.1_dp*radius)
    end function regular_G_third

    pure real(dp) function regular_psi(radius)
        real(dp), intent(in) :: radius

        regular_psi = 0.4_dp*radius + 0.07_dp*sin(0.9_dp*radius) + &
            0.02_dp*cos(1.7_dp*radius)
    end function regular_psi

    pure real(dp) function regular_psi_first(radius)
        real(dp), intent(in) :: radius

        regular_psi_first = 0.4_dp + 0.063_dp*cos(0.9_dp*radius) - &
            0.034_dp*sin(1.7_dp*radius)
    end function regular_psi_first

    pure real(dp) function regular_psi_second(radius)
        real(dp), intent(in) :: radius

        regular_psi_second = -0.0567_dp*sin(0.9_dp*radius) - &
            0.0578_dp*cos(1.7_dp*radius)
    end function regular_psi_second

    pure real(dp) function regular_psi_third(radius)
        real(dp), intent(in) :: radius

        regular_psi_third = -0.05103_dp*cos(0.9_dp*radius) + &
            0.09826_dp*sin(1.7_dp*radius)
    end function regular_psi_third

    pure real(dp) function psi_star_direct(radius, sigma_value)
        real(dp), intent(in) :: radius, sigma_value

        psi_star_direct = regular_psi(radius) + mass*c_light*sigma_value/ &
            charge*sqrt(regular_F(radius))*regular_G(radius)
    end function psi_star_direct

    pure real(dp) function direct_first(radius, h_value, sigma_value)
        real(dp), intent(in) :: radius, h_value, sigma_value

        direct_first = (psi_star_direct(radius-2.0_dp*h_value, sigma_value) &
            - 8.0_dp*psi_star_direct(radius-h_value, sigma_value) + &
            8.0_dp*psi_star_direct(radius+h_value, sigma_value) - &
            psi_star_direct(radius+2.0_dp*h_value, sigma_value))/ &
            (12.0_dp*h_value)
    end function direct_first

    pure real(dp) function direct_second(radius, h_value, sigma_value)
        real(dp), intent(in) :: radius, h_value, sigma_value

        direct_second = (-psi_star_direct(radius+2.0_dp*h_value, sigma_value) &
            + 16.0_dp*psi_star_direct(radius+h_value, sigma_value) - &
            30.0_dp*psi_star_direct(radius, sigma_value) + &
            16.0_dp*psi_star_direct(radius-h_value, sigma_value) - &
            psi_star_direct(radius-2.0_dp*h_value, sigma_value))/ &
            (12.0_dp*h_value**2)
    end function direct_second

    pure real(dp) function direct_n(radius, h_value, sigma_value)
        real(dp), intent(in) :: radius, h_value, sigma_value

        direct_n = 2.0_dp*sqrt(regular_F(radius))* &
            direct_first(radius, h_value, sigma_value)
    end function direct_n

    pure real(dp) function direct_n_first(radius, h_value, sigma_value)
        real(dp), intent(in) :: radius, h_value, sigma_value
        real(dp) :: inner_h

        inner_h = h_value/8.0_dp
        direct_n_first = (direct_n(radius-2.0_dp*h_value, inner_h, sigma_value) &
            - 8.0_dp*direct_n(radius-h_value, inner_h, sigma_value) + &
            8.0_dp*direct_n(radius+h_value, inner_h, sigma_value) - &
            direct_n(radius+2.0_dp*h_value, inner_h, sigma_value))/ &
            (12.0_dp*h_value)
    end function direct_n_first

    pure real(dp) function direct_n_second(radius, h_value, sigma_value)
        real(dp), intent(in) :: radius, h_value, sigma_value
        real(dp) :: inner_h

        inner_h = h_value/8.0_dp
        direct_n_second = (-direct_n(radius+2.0_dp*h_value, inner_h, sigma_value) &
            + 16.0_dp*direct_n(radius+h_value, inner_h, sigma_value) - &
            30.0_dp*direct_n(radius, inner_h, sigma_value) + &
            16.0_dp*direct_n(radius-h_value, inner_h, sigma_value) - &
            direct_n(radius-2.0_dp*h_value, inner_h, sigma_value))/ &
            (12.0_dp*h_value**2)
    end function direct_n_second

    pure real(dp) function turning_F(radius)
        real(dp), intent(in) :: radius

        turning_F = turning_slope*(radius-turning_radius)
    end function turning_F

    pure real(dp) function turning_g(radius)
        real(dp), intent(in) :: radius

        turning_g = regular_G(radius)
    end function turning_g

    pure real(dp) function turning_g_first(radius)
        real(dp), intent(in) :: radius

        turning_g_first = regular_G_first(radius)
    end function turning_g_first

    pure real(dp) function turning_psi_star(y_value, sigma_value)
        real(dp), intent(in) :: y_value, sigma_value
        real(dp) :: radius

        radius = turning_radius+turning_side*y_value**2
        turning_psi_star = regular_psi(radius) + &
            mass*c_light*sigma_value/charge*sqrt(turning_F(radius))* &
            turning_g(radius)
    end function turning_psi_star

    pure real(dp) function turning_psi_y_first(y_value, h_value, sigma_value)
        real(dp), intent(in) :: y_value, h_value, sigma_value

        turning_psi_y_first = (turning_psi_star(y_value-2.0_dp*h_value, &
            sigma_value) - 8.0_dp*turning_psi_star(y_value-h_value, &
            sigma_value) + 8.0_dp*turning_psi_star(y_value+h_value, &
            sigma_value) - turning_psi_star(y_value+2.0_dp*h_value, &
            sigma_value))/(12.0_dp*h_value)
    end function turning_psi_y_first

    pure real(dp) function direct_turning_n(y_value, sigma_value)
        real(dp), intent(in) :: y_value, sigma_value
        real(dp) :: h_value, radius

        h_value = y_value/32.0_dp
        radius = turning_radius+turning_side*y_value**2
        direct_turning_n = sqrt(turning_F(radius))/(turning_side*y_value)* &
            turning_psi_y_first(y_value, h_value, sigma_value)
    end function direct_turning_n

    subroutine check_canonical_interval_enclosure
        type(gc_outward_interval_t) :: input(15), output(3)
        real(dp) :: radius, scalar_output(3)
        integer :: j

        input(1) = gc_outward_interval(1.35_dp, 1.7_dp)
        input(2) = gc_outward_interval(-0.2_dp, 0.2_dp)
        input(3) = gc_outward_interval(-0.2_dp, 0.2_dp)
        input(4) = gc_outward_interval(-0.2_dp, 0.2_dp)
        input(5) = gc_outward_interval(0.85_dp, 1.05_dp)
        input(6) = gc_outward_interval(-0.1_dp, 0.1_dp)
        input(7) = gc_outward_interval(-0.1_dp, 0.1_dp)
        input(8) = gc_outward_interval(-0.1_dp, 0.1_dp)
        input(9) = gc_outward_interval(0.35_dp, 0.5_dp)
        input(10) = gc_outward_interval(-0.12_dp, 0.12_dp)
        input(11) = gc_outward_interval(-0.15_dp, 0.15_dp)
        input(12) = gc_outward_interval(mass, mass)
        input(13) = gc_outward_interval(c_light, c_light)
        input(14) = gc_outward_interval(1.0_dp, 1.0_dp)
        input(15) = gc_outward_interval(charge, charge)
        call evaluate_neort_full_fow_canonical_numerator_interval( &
            input(1), input(2), input(3), input(4), input(5), input(6), &
            input(7), input(8), input(9), input(10), input(11), input(12), &
            input(13), input(14), input(15), output(1), output(2), output(3))
        do j = 0, 8
            radius = 0.2_dp+0.1_dp*real(j, dp)
            call evaluate_neort_full_fow_canonical_numerator( &
                regular_F(radius), regular_F_first(radius), &
                regular_F_second(radius), regular_F_third(radius), &
                regular_G(radius), regular_G_first(radius), &
                regular_G_second(radius), regular_G_third(radius), &
                regular_psi_first(radius), regular_psi_second(radius), &
                regular_psi_third(radius), mass, c_light, 1.0_dp, charge, &
                scalar_output(1), scalar_output(2), scalar_output(3))
            call require_contains(output(1), scalar_output(1), &
                'canonical interval N')
            call require_contains(output(2), scalar_output(2), &
                'canonical interval N_R')
            call require_contains(output(3), scalar_output(3), &
                'canonical interval N_RR')
        end do
    end subroutine check_canonical_interval_enclosure

    subroutine check_turning_interval_enclosure
        type(gc_outward_interval_t) :: output(2)
        real(dp) :: scalar_output(2)
        integer :: j

        call evaluate_neort_full_fow_canonical_turning_interval( &
            gc_outward_interval(1.0_dp, 1.0_dp), &
            gc_outward_interval(1.1_dp, 1.4_dp), &
            gc_outward_interval(0.85_dp, 1.0_dp), &
            gc_outward_interval(-0.1_dp, 0.1_dp), &
            gc_outward_interval(0.35_dp, 0.5_dp), &
            gc_outward_interval(mass, mass), gc_outward_interval(c_light, &
            c_light), gc_outward_interval(1.0_dp, 1.0_dp), &
            gc_outward_interval(charge, charge), output(1), output(2))
        do j = 0, 6
            call evaluate_neort_full_fow_canonical_turning(1.0_dp, &
                1.1_dp+0.05_dp*real(j, dp), 0.85_dp+0.02_dp*real(j, dp), &
                0.0_dp, 0.35_dp+0.02_dp*real(j, dp), mass, c_light, &
                1.0_dp, charge, scalar_output(1), scalar_output(2))
            call require_contains(output(1), scalar_output(1), &
                'turning interval N')
            call require_contains(output(2), scalar_output(2), &
                'turning interval dN/dy')
        end do
    end subroutine check_turning_interval_enclosure

    subroutine check_cut_flux_interval_enclosure
        type(gc_outward_interval_t) :: output(3)
        real(dp) :: scalar_output(3), a, b, c, d
        integer :: ia, ib, ic, id

        call evaluate_neort_eqdsk_cut_flux_coordinate_interval( &
            gc_outward_interval(1.8_dp, 2.2_dp), &
            gc_outward_interval(0.9_dp, 1.1_dp), &
            gc_outward_interval(2.0_dp, 2.4_dp), &
            gc_outward_interval(-0.3_dp, 0.2_dp), output(1), output(2), &
            output(3))
        do ia = 0, 2
            a = 1.8_dp+0.2_dp*real(ia, dp)
            do ib = 0, 2
                b = 0.9_dp+0.1_dp*real(ib, dp)
                do ic = 0, 2
                    c = 2.0_dp+0.2_dp*real(ic, dp)
                    do id = 0, 2
                        d = -0.3_dp+0.25_dp*real(id, dp)
                        call evaluate_neort_eqdsk_cut_flux_coordinate(a, b, &
                            c, d, scalar_output(1), scalar_output(2), &
                            scalar_output(3))
                        call require_contains(output(1), scalar_output(1), &
                            'cut flux interval d psi/d rho')
                        call require_contains(output(2), scalar_output(2), &
                            'cut flux interval dR/d rho')
                        call require_contains(output(3), scalar_output(3), &
                            'cut flux interval dZ/d rho')
                    end do
                end do
            end do
        end do
    end subroutine check_cut_flux_interval_enclosure

    subroutine require_convergence(label, values)
        character(*), intent(in) :: label
        real(dp), intent(in) :: values(3)
        real(dp) :: coarse_change, fine_change

        coarse_change = abs(values(2)-values(1))
        fine_change = abs(values(3)-values(2))
        if (fine_change > 0.2_dp*coarse_change+1.0e-10_dp) then
            write (*, '(A,2(1X,ES24.16))') trim(label), coarse_change, &
                fine_change
            error stop 'finite-difference oracle did not converge'
        end if
    end subroutine require_convergence

    subroutine require_close(label, actual, expected, local_tolerance)
        character(*), intent(in) :: label
        real(dp), intent(in) :: actual, expected, local_tolerance

        scale = max(1.0_dp, max(abs(actual), abs(expected)))
        if (abs(actual-expected) > local_tolerance*scale) then
            write (*, '(A,2(1X,ES24.16))') trim(label), actual, expected
            error stop 'canonical behavioral oracle failed'
        end if
    end subroutine require_close

    subroutine require_nonzero(label, value)
        character(*), intent(in) :: label
        real(dp), intent(in) :: value

        if (abs(value) <= tolerance) then
            write (*, '(A,1X,ES24.16)') trim(label), value
            error stop 'canonical sign distinction was lost'
        end if
    end subroutine require_nonzero

    subroutine require_contains(interval, value, label)
        type(gc_outward_interval_t), intent(in) :: interval
        real(dp), intent(in) :: value
        character(*), intent(in) :: label

        if (value < interval%lo .or. value > interval%hi) then
            write (*, '(A,3(1X,ES24.16))') trim(label), interval%lo, &
                value, interval%hi
            error stop 'interval behavioral oracle failed'
        end if
    end subroutine require_contains

end program test_gc_full_fow_canonical_behavior
