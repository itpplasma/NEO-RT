module poicut_edge_fixture
    ! Independent analytic field: psi=(R-3)^2+(Z-.4)^2, B=1/R. The rectangular
    ! field provider returns zero beyond R=[2,4], just as a convex-wall field
    ! reader does. The requested cut endpoint is still INSIDE that rectangle, so
    ! a fixed continuation step that overshoots must be retried, not followed.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), parameter :: r_axis = 3.0_dp, z_axis = 0.4_dp
    real(dp) :: flux_sign = 1.0_dp
    real(dp) :: psi_last = 0.0_dp
contains

    subroutine field_at_rz(r, z, psi, dpsi_dr, dpsi_dz, g, gr, gz)
        real(dp), intent(in) :: r, z
        real(dp), intent(out) :: psi, dpsi_dr, dpsi_dz, g, gr, gz

        if (r < 2.0_dp .or. r > 4.0_dp) then
            psi = 0.0_dp
            dpsi_dr = 0.0_dp
            dpsi_dz = 0.0_dp
            g = 0.0_dp
            gr = 0.0_dp
            gz = 0.0_dp
        else
            psi = flux_sign*((r - r_axis)**2 + (z - z_axis)**2)
            dpsi_dr = flux_sign*2.0_dp*(r - r_axis)
            dpsi_dz = flux_sign*2.0_dp*(z - z_axis)
            g = flux_sign*2.0_dp*(z - z_axis)/r**2
            gr = -flux_sign*4.0_dp*(z - z_axis)/r**3
            gz = flux_sign*2.0_dp/r**2
        end if
        psi_last = psi
    end subroutine field_at_rz

end module poicut_edge_fixture

program test_poicut_edge
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use poicut_tracing_mod, only: trace_poicut_end
    use poicut_edge_fixture, only: r_axis, z_axis, flux_sign, psi_last, field_at_rz
    implicit none
    real(dp), parameter :: radii(3) = [0.95_dp, 0.995_dp, 0.9995_dp]
    real(dp) :: rend, zend, radius
    integer :: side, k, orientation

    do orientation = -1, 1, 2
        flux_sign = real(orientation, dp)
        do k = 1, 3
            radius = radii(k)
            do side = -1, 1, 2
                call trace_poicut_end(r_axis, z_axis, side*r_axis/100.0_dp, &
                                      flux_sign*radius**2, flux_sign, 1.0e-12_dp, &
                                      field_at_rz, rend, zend)
                if (abs(rend - (r_axis + side*radius)) > 1.0e-10_dp) &
                    error stop 'Poincare endpoint moved from its analytic flux surface'
                if (abs(zend - z_axis) > 1.0e-12_dp) &
                    error stop 'Poincare endpoint left the magnetic-axis branch'
                ! The final field state is consumed immediately by the caller.
                if (abs(psi_last - flux_sign*radius**2) > 1.0e-10_dp) &
                    error stop 'Poincare endpoint left stale field state'
            end do
        end do
    end do
end program test_poicut_edge
