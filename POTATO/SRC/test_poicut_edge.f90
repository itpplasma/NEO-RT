module field_sub
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp) :: psif, dpsidr, dpsidz, flux_sign
end module field_sub

module poicut_mod
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp) :: rmagaxis = 3.0_dp, zmagaxis = 0.4_dp
end module poicut_mod

! Independent analytic field: psi=(R-3)^2+(Z-.4)^2, B=1/R.
! The rectangular field provider returns zero beyond R=[2,4], just as a
! convex-wall field reader does. The requested cut endpoint is still INSIDE.
subroutine gpsigb_and_ders(r, z, g, gr, gz)
    use field_sub, only: dp, psif, dpsidr, dpsidz, flux_sign
    use poicut_mod, only: rmagaxis, zmagaxis
    implicit none
    real(dp), intent(in) :: r, z
    real(dp), intent(out) :: g, gr, gz

    if (r < 2.0_dp .or. r > 4.0_dp) then
        psif = 0.0_dp
        dpsidr = 0.0_dp
        dpsidz = 0.0_dp
        g = 0.0_dp
        gr = 0.0_dp
        gz = 0.0_dp
        return
    end if
    psif = flux_sign*((r - rmagaxis)**2 + (z - zmagaxis)**2)
    dpsidr = flux_sign*2.0_dp*(r - rmagaxis)
    dpsidz = flux_sign*2.0_dp*(z - zmagaxis)
    g = flux_sign*2.0_dp*(z - zmagaxis)/r**2
    gr = -flux_sign*4.0_dp*(z - zmagaxis)/r**3
    gz = flux_sign*2.0_dp/r**2
end subroutine gpsigb_and_ders

program test_poicut_edge
    use field_sub, only: dp, psif, flux_sign
    use poicut_mod, only: rmagaxis, zmagaxis
    implicit none
    real(dp), parameter :: radii(3) = [0.95_dp, 0.995_dp, 0.9995_dp]
    real(dp) :: rend, zend, radius
    integer :: side, k, orientation

    do orientation = -1, 1, 2
        flux_sign = real(orientation, dp)
        do k = 1, 3
            radius = radii(k)
            do side = -1, 1, 2
                call trace_poicut_end(side*rmagaxis/100.0_dp, flux_sign*radius**2, &
                                      flux_sign, 1.0e-12_dp, rend, zend)
                if (abs(rend - (rmagaxis + side*radius)) > 1.0e-10_dp) &
                    error stop 'Poincare endpoint moved from its analytic flux surface'
                if (abs(zend - zmagaxis) > 1.0e-12_dp) &
                    error stop 'Poincare endpoint left the magnetic-axis branch'
                ! The final field state is consumed immediately by the caller.
                if (abs(psif - flux_sign*radius**2) > 1.0e-10_dp) &
                    error stop 'Poincare endpoint left stale field state'
            end do
        end do
    end do
end program test_poicut_edge
