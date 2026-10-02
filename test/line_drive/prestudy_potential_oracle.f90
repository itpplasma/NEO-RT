module prestudy_potential_oracle
    !! Independent m_b=0 trapped-orbit electrostatic toggle for the circular
    !! pre-study. Geometry and radial displacement are specified physically;
    !! no production field, perturbation, drift or gauge routines are called.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_quadrature, only: gauss_legendre_ab

    implicit none
    private
    public :: prestudy_potential_harmonic

    real(dp), parameter :: pi = acos(-1.0_dp)
    complex(dp), parameter :: imaginary_unit = (0.0_dp, 1.0_dp)

contains

    subroutine prestudy_potential_harmonic( &
            r, edge_flux, amp, eta, speed, theta_orientation, nq, &
            potential_prefactor, expected, convergence)
        !! expected = (qi Phi0'/E) <xi^s exp(i n q theta_B)>_time.
        !! r is minor radius / R0 (R0=1 m); edge_flux is the dimensionless
        !! toroidal edge flux 1-sqrt(1-a^2), amp is the physical displacement
        !! amplitude / R0, eta has units 1/G, and speed has units cm/s.
        !! theta_orientation is +/-1; nq is the signed product n*q.
        !! potential_prefactor = qi*Phi0prime/(mi*speed^2/2), supplied by host.
        !! The full trapped bounce contains two equal crossing contributions;
        !! their factor of two cancels in the normalized time average.
        real(dp), intent(in) :: r, edge_flux, amp, eta, speed
        real(dp), intent(in) :: theta_orientation, nq, potential_prefactor
        complex(dp), intent(out) :: expected
        real(dp), intent(out) :: convergence
        complex(dp) :: coarse, fine
        real(dp) :: time_coarse, time_fine, norm

        if (r <= 0.0_dp) error stop 'potential oracle: radius must be positive'
        if (r >= 1.0_dp) error stop 'potential oracle: radius must be below R0'
        if (edge_flux <= 0.0_dp) error stop 'potential oracle: invalid edge flux'
        if (eta <= 0.0_dp) error stop 'potential oracle: eta must be positive'
        if (speed <= 0.0_dp) error stop 'potential oracle: speed must be positive'
        if (abs(abs(theta_orientation) - 1.0_dp) > 1.0e-12_dp) then
            error stop 'potential oracle: theta orientation must be +/-1'
        end if
        call integrate_crossing(128, r, edge_flux, amp, eta, speed, &
            theta_orientation, nq, coarse, time_coarse)
        call integrate_crossing(256, r, edge_flux, amp, eta, speed, &
            theta_orientation, nq, fine, time_fine)
        expected = potential_prefactor*fine
        norm = max(abs(fine), abs(coarse), tiny(1.0_dp))
        convergence = max(abs(fine - coarse)/norm, &
            abs(time_fine - time_coarse)/time_fine)
    end subroutine prestudy_potential_harmonic

    subroutine integrate_crossing(order, r, edge_flux, amp, eta, speed, &
            theta_orientation, nq, average, crossing_time)
        integer, intent(in) :: order
        real(dp), intent(in) :: r, edge_flux, amp, eta, speed
        real(dp), intent(in) :: theta_orientation, nq
        complex(dp), intent(out) :: average
        real(dp), intent(out) :: crossing_time
        real(dp) :: nodes(order), weights(order)
        real(dp) :: qbar, bscale, htheta, turn_cosine, turning_angle, k, dsdr
        real(dp) :: g1, g2, u, sin_u, cos_u, theta_geom, theta_b, mirror
        real(dp) :: small_angle, large_angle, cosine_difference, time_weight
        complex(dp) :: radial_displacement, xi_s, numerator
        integer :: j

        qbar = 1.5_dp + 8.0_dp*r**2
        bscale = 1.0e4_dp*sqrt(1.0_dp + (r/qbar)**2)
        ! B(theta_geom) = bscale/(1+r*cos(theta_geom)).
        htheta = 1.0_dp/(100.0_dp*qbar*sqrt(1.0_dp + (r/qbar)**2))
        turn_cosine = (eta*bscale - 1.0_dp)/r
        if (turn_cosine <= -1.0_dp) then
            error stop 'potential oracle: passing orbit is outside its contract'
        end if
        if (turn_cosine >= 1.0_dp) then
            error stop 'potential oracle: trapped endpoint is outside its contract'
        end if
        turning_angle = acos(turn_cosine)
        k = sqrt((1.0_dp - r)/(1.0_dp + r))
        dsdr = r/(sqrt(1.0_dp - r**2)*edge_flux)
        g1 = exp(-((r - 0.2_dp)/0.15_dp)**2)
        g2 = 0.5_dp*exp(-((r - 0.25_dp)/0.2_dp)**2)
        call gauss_legendre_ab(order, -0.5_dp*pi, 0.5_dp*pi, nodes, weights)
        numerator = 0
        crossing_time = 0
        do j = 1, order
            u = nodes(j)
            sin_u = sin(u)
            cos_u = cos(u)
            theta_geom = turning_angle*sin_u
            theta_b = theta_orientation*2.0_dp*atan2( &
                k*sin(0.5_dp*theta_geom), cos(0.5_dp*theta_geom))
            ! cos(theta_geom)-cos(turning_angle), evaluated without cancellation.
            ! The small factor uses 1-|sin(u)|=cos(u)^2/(1+|sin(u)|).
            small_angle = turning_angle*cos_u**2/(1.0_dp + abs(sin_u))
            large_angle = turning_angle*(1.0_dp + abs(sin_u))
            cosine_difference = 2.0_dp*sin(0.5_dp*small_angle) &
                *sin(0.5_dp*large_angle)
            mirror = r*cosine_difference/(1.0_dp + r*cos(theta_geom))
            if (mirror <= 0.0_dp) error stop 'potential oracle: invalid mirror factor'
            time_weight = weights(j)*turning_angle*cos_u &
                /(speed*abs(htheta)*sqrt(mirror))
            radial_displacement = amp*(g1*exp(-2.0_dp*imaginary_unit*theta_geom) &
                + g2*exp(-3.0_dp*imaginary_unit*theta_geom))
            xi_s = dsdr*radial_displacement
            numerator = numerator + time_weight*xi_s &
                *exp(imaginary_unit*nq*theta_b)
            crossing_time = crossing_time + time_weight
        end do
        if (crossing_time <= 0.0_dp) error stop 'potential oracle: zero crossing time'
        average = numerator/crossing_time
    end subroutine integrate_crossing

end module prestudy_potential_oracle
