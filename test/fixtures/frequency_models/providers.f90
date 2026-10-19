! Analytic orbit providers for testing the actual frequency/root implementation.
module logger
    implicit none
    integer, parameter :: LOG_TRACE=5
contains
    subroutine debug(msg)
        character(*), intent(in) :: msg
    end subroutine
    subroutine trace(msg)
        character(*), intent(in) :: msg
    end subroutine
    subroutine warning(msg)
        character(*), intent(in) :: msg
        print '(a)', trim(msg)
    end subroutine
    integer function get_log_level()
        get_log_level=0
    end function
end module
module util
    use iso_fortran_env, only: dp=>real64
    implicit none
    real(dp), parameter :: pi=acos(-1.0_dp)
end module
module driftorbit
    use iso_fortran_env, only: dp=>real64
    implicit none
    real(dp) :: etamin=0, etamax=0, etatp=1, etadt=1.5_dp
    real(dp) :: epsst_spl=1e-8_dp, epst_spl=1e-6_dp, epst=1e-8_dp
    real(dp) :: epssp_spl=1e-8_dp, epsp_spl=1e-6_dp, sign_vpar=1, sign_vpar_htheta=1, vth=1
    integer :: mph=1, mth=0, nlev=100, magdrift_passing=1
    logical :: magdrift=.true., nonlin=.false., supban=.false.
end module
module neort_profiles
    use driftorbit, only: vth
    use iso_fortran_env, only: dp=>real64
    implicit none
    real(dp) :: Om_tE=.3_dp, dOm_tEds=0
end module
module do_magfie_mod
    use iso_fortran_env, only: dp=>real64
    implicit none
    real(dp) :: iota=.5_dp, s=.09_dp, Bthcov=1, Bphcov=1, q=2
end module
module shaing
    use iso_fortran_env, only: dp=>real64
    implicit none
contains
    real(dp) function omph_shaing(v,eta)
        real(dp), intent(in) :: v,eta
        omph_shaing=0
    end function
end module
module neort_orbit
    use iso_fortran_env, only: dp=>real64
    use driftorbit, only: etatp
    use util, only: pi
    implicit none
    integer, parameter :: nvar=7
    abstract interface
        subroutine timestep_i(v,eta,neq,t,y,ydot)
            import dp
            integer, intent(in) :: neq
            real(dp), intent(in) :: v,eta,t,y(neq)
            real(dp), intent(out) :: ydot(neq)
        end subroutine
    end interface
contains
    real(dp) function bounce_time(v,eta,taub_estimate) result(taub)
        real(dp), intent(in) :: v,eta
        real(dp), intent(in), optional :: taub_estimate
        taub=(5.0_dp-.125_dp*log(abs(eta-etatp)))/v
    end function
    subroutine bounce_fast(v,eta,taub,bounceavg,ts,istate_out)
        real(dp), intent(in) :: v,eta,taub
        real(dp), intent(out) :: bounceavg(nvar)
        procedure(timestep_i) :: ts
        integer, intent(out), optional :: istate_out
        bounceavg=0
        bounceavg(3)=(1.25_dp+.05_dp*log(abs(eta-etatp)))*2*pi/(v*taub)
        if(present(istate_out)) istate_out=2
    end subroutine
    subroutine timestep(v,eta,neq,t,y,ydot)
        integer, intent(in) :: neq
        real(dp), intent(in) :: v,eta,t,y(neq)
        real(dp), intent(out) :: ydot(neq)
        ydot=0
    end subroutine
end module
