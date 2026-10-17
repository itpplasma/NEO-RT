program test_krook_work_oracle
    use iso_fortran_env, only: dp => real64
    use util, only: pi, c, qi, mi, imun
    use neort_transport, only: Tphi_int
    use neort_krook, only: integrate_krook
    use neort_profiles, only: ns, vth, A1, A2
    use do_magfie_mod, only: q, psi_pr, sign_theta
    use driftorbit, only: mph
    implicit none
    real(dp) :: u, taub, nu, detuning, jacobian, f0, fh, fp, force, expected, result(1)
    real(dp) :: max_error, tolerance
    complex(dp) :: h, response
    integer :: charge_case, geometry_case, source_case, status
    ! Manufactured action-space fixture; use the actual native c, with qi=+/-c.
    ! m=vth=ns=1 implies T=1/2 and f0=exp(-u**2)/pi**(3/2).
    mi = 1.0_dp
    vth = 1.0_dp
    ns = 1.0_dp
    taub = 0.7_dp
    h = cmplx(0.4_dp, 0.7_dp, dp)
    nu = 2.0_dp
    max_error = 0.0_dp
    do charge_case = -1, 1, 2
        qi = charge_case*c
        do geometry_case = 0, 3
            q = 1.0_dp
            psi_pr = 1.0_dp
            mph = 2
            detuning = 1.0_dp
            if (geometry_case == 1) q = -1.0_dp
            if (geometry_case == 2) psi_pr = -1.0_dp
            if (geometry_case == 3) then
                q = -1.0_dp
                psi_pr = -1.0_dp
                mph = -2
                detuning = -3.0_dp
            end if
            do source_case = 0, 1
                u = 0.8_dp + 0.5_dp*source_case
                A1 = real(source_case, dp)
                A2 = 0.25_dp*source_case
                force = A1 + A2*u**2
                f0 = exp(-u**2)/pi**1.5_dp
                fh = -2.0_dp*f0
                fp = -c*q/(qi*sign_theta*psi_pr)*f0*force
                ! Canonical phase-space Jacobian per du deta ds, independent of Tphi_int:
                ! (2pi)^3 * [m^2*c*v^2/(2|qi|)] * [taub/(2pi)]
                ! * [m*vth^2*u] * |qi*psi_pr/(c*q)|.
                jacobian = (2*pi)**3*(mi**2*c*(u*vth)**2/(2*abs(qi))) &
                    *(taub/(2*pi))*(mi*vth**2*u)*abs(qi*psi_pr/(c*q))
                response = fh*h + imun*real(mph,dp)*fp*h/(nu+imun*detuning)
                expected = jacobian*phase_work(response)
                call integrate_krook(residual, weight, 0.0_dp, 1.0_dp, nu, 4, &
                    1.0e-11_dp, 1.0e-14_dp, result, status)
                if (status /= 0) error stop 'manufactured work quadrature failed'
                tolerance = 2.0e-12_dp*max(1.0_dp, abs(expected))
                if (abs(result(1)-expected) > tolerance) then
                    print *, 'native torque and independent perturbation work:', result(1), expected
                    error stop 'Krook work normalization or nonadiabatic source mismatch'
                end if
                max_error = max(max_error, abs(result(1)-expected))
                if (source_case == 0 .and. result(1) /= 0.0_dp) &
                    error stop 'static zero-drive Maxwellian has nonzero torque'
            end do
        end do
    end do
    print *, 'Krook direct-work absolute maximum error =', max_error
    print *, 'test_krook_work_oracle PASSED'
contains
    real(dp) function phase_work(f) result(average)
        complex(dp), intent(in) :: f
        complex(dp) :: phase
        integer :: k
        average = 0.0_dp
        do k = 0, 3
            phase = exp(imun*real(k,dp)*pi/2)
            average = average-real(imun*real(mph,dp)*h*phase,dp)*real(f*phase,dp)/4
        end do
    end function phase_work
    subroutine residual(eta, value, derivative)
        real(dp), intent(in) :: eta
        real(dp), intent(out) :: value, derivative
        value = detuning
        derivative = 0.0_dp
    end subroutine residual
    subroutine weight(eta, values)
        real(dp), intent(in) :: eta
        real(dp), intent(out) :: values(:)
        values(1) = Tphi_int(u, taub, abs(h)**2)
    end subroutine weight
end program test_krook_work_oracle
