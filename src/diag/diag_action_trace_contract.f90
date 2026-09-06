module neort_action_trace_contract
    use iso_fortran_env, only: dp => real64
    implicit none

contains

    pure function action_trace_weight(du, tphi, jacobian, attenuation) result(weight)
        real(dp), intent(in) :: du, tphi, jacobian, attenuation
        real(dp) :: weight

        weight = du * tphi / abs(jacobian) * attenuation
    end function action_trace_weight

end module neort_action_trace_contract
