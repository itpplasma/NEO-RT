module neort_action_trace_contract
    use iso_fortran_env, only: dp => real64
    implicit none

contains

    pure function complex_orbit_action(bounce_real, bounce_imag, kinetic_energy) &
            result(action)
        real(dp), intent(in) :: bounce_real, bounce_imag, kinetic_energy
        complex(dp) :: action

        action = kinetic_energy*cmplx(bounce_real, bounce_imag, dp)
    end function complex_orbit_action

    pure function action_modulus_squared(action) result(norm2)
        complex(dp), intent(in) :: action
        real(dp) :: norm2

        norm2 = real(action*conjg(action), dp)
    end function action_modulus_squared

    pure function action_trace_weight(du, tphi, jacobian, attenuation) result(weight)
        real(dp), intent(in) :: du, tphi, jacobian, attenuation
        real(dp) :: weight

        weight = du * tphi / abs(jacobian) * attenuation
    end function action_trace_weight

end module neort_action_trace_contract
