module line_period_oracle
    !! Independent thin-orbit period from the conserved parallel energy.
    !! No orbit integration, shooting, frequency spline, or line-drive call.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use do_magfie_mod, only: do_magfie
    use fortnum_quadrature, only: gauss_legendre_ab
    implicit none
    private
    real(dp), parameter :: pi = acos(-1.0_dp)
    public :: thin_period_oracle
contains
    subroutine thin_period_oracle(s0, theta0, v, eta, trapped, period, &
            relative_error, order)
        real(dp), intent(in) :: s0, theta0, v, eta
        logical, intent(in) :: trapped
        real(dp), intent(out) :: period, relative_error
        integer, intent(out) :: order
        real(dp) :: left, right, previous, current, b, htheta, dbdtheta
        integer :: n

        if (v <= 0.0_dp) error stop "period oracle requires positive speed"
        call field(theta0, b, htheta, dbdtheta)
        if (1.0_dp - eta*b <= 0.0_dp) then
            error stop "period oracle requires an allowed initial point"
        end if
        left = theta0
        right = theta0 + 2.0_dp*pi
        if (trapped) then
            call turning_point(-1.0_dp, left)
            call turning_point(1.0_dp, right)
        end if
        previous = 0.0_dp
        relative_error = huge(1.0_dp)
        n = 32
        do
            call integrate_period(n, current)
            if (n > 32) then
                relative_error = abs(current - previous)/current
                if (relative_error < 5.0e-12_dp) exit
            end if
            if (n >= 512) error stop "period oracle Gaussian rule did not converge"
            previous = current
            n = 2*n
        end do
        period = current
        order = n
    contains
        subroutine field(theta, bmod, hth, dbth)
            real(dp), intent(in) :: theta
            real(dp), intent(out) :: bmod, hth, dbth
            real(dp) :: sqrtg, bder(3), hcov(3), hctr(3), hcurl(3)

            call do_magfie([s0, 0.0_dp, theta], bmod, sqrtg, &
                bder, hcov, hctr, hcurl)
            hth = hctr(3)
            dbth = bmod*bder(3)
        end subroutine field

        subroutine turning_point(direction, root)
            !! Nearest physical mirror point in this initial-point well.
            real(dp), intent(in) :: direction
            real(dp), intent(out) :: root
            real(dp) :: allowed, forbidden, midpoint, bmod, hth, dbth, f
            integer :: j

            allowed = theta0
            do j = 1, 4096
                forbidden = theta0 + direction*j*2.0_dp*pi/4096.0_dp
                call field(forbidden, bmod, hth, dbth)
                if (1.0_dp - eta*bmod <= 0.0_dp) exit
                allowed = forbidden
            end do
            if (j > 4096) error stop "period oracle failed to bracket mirror point"
            do j = 1, 80
                midpoint = 0.5_dp*(allowed + forbidden)
                if (midpoint == allowed) exit
                if (midpoint == forbidden) exit
                call field(midpoint, bmod, hth, dbth)
                f = 1.0_dp - eta*bmod
                if (f > 0.0_dp) then
                    allowed = midpoint
                else
                    forbidden = midpoint
                end if
            end do
            root = 0.5_dp*(allowed + forbidden)
        end subroutine turning_point

        function mirror_gap(theta, distance, lower_endpoint) result(gap)
            !! Integrate dB/dtheta near a mirror point to avoid subtracting
            !! nearly equal B values at high-order Gaussian endpoint nodes.
            real(dp), intent(in) :: theta, distance
            logical, intent(in) :: lower_endpoint
            real(dp) :: gap, nodes(4), weights(4), bmod, hth, dbth, integral
            real(dp) :: a, z
            integer :: k

            a = theta
            z = right
            if (lower_endpoint) then
                a = left
                z = theta
            end if
            call gauss_legendre_ab(4, 0.0_dp, 1.0_dp, nodes, weights)
            integral = 0.0_dp
            do k = 1, 4
                call field(a + nodes(k)*(z - a), bmod, hth, dbth)
                integral = integral + weights(k)*dbth
            end do
            gap = eta*distance*integral
            if (lower_endpoint) gap = -gap
        end function mirror_gap

        subroutine integrate_period(nrule, value)
            integer, intent(in) :: nrule
            real(dp), intent(out) :: value
            real(dp) :: nodes(nrule), weights(nrule), theta, bmod, hth, dbth
            real(dp) :: center, radius, u, jacobian, gap, dl, dr
            integer :: k

            center = 0.5_dp*(left + right)
            radius = 0.5_dp*(right - left)
            if (trapped) then
                call gauss_legendre_ab(nrule, -0.5_dp*pi, 0.5_dp*pi, nodes, weights)
            else
                call gauss_legendre_ab(nrule, left, right, nodes, weights)
            end if
            value = 0.0_dp
            do k = 1, nrule
                theta = nodes(k)
                jacobian = 1.0_dp
                if (trapped) then
                    u = nodes(k)
                    dl = 2.0_dp*radius*sin(0.5_dp*(u + 0.5_dp*pi))**2
                    dr = 2.0_dp*radius*sin(0.5_dp*(0.5_dp*pi - u))**2
                    theta = center + radius*sin(u)
                    if (dl < dr) theta = left + dl
                    if (dr <= dl) theta = right - dr
                    jacobian = radius*cos(u)
                end if
                call field(theta, bmod, hth, dbth)
                gap = 1.0_dp - eta*bmod
                if (trapped) then
                    if (dl < 1.0e-3_dp) gap = mirror_gap(theta, dl, .true.)
                    if (dr < 1.0e-3_dp) gap = mirror_gap(theta, dr, .false.)
                end if
                if (gap <= 0.0_dp) error stop "period oracle sampled forbidden point"
                if (abs(hth) <= tiny(1.0_dp)) error stop "period oracle htheta is zero"
                value = value + weights(k)*jacobian/(abs(hth)*sqrt(gap)*v)
            end do
            if (trapped) value = 2.0_dp*value
        end subroutine integrate_period
    end subroutine thin_period_oracle
end module line_period_oracle
