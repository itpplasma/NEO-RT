module thetadata_mod
    use iso_fortran_env, only: dp => real64

    implicit none

    logical :: attenuation_data_initialized = .false.
    integer :: npoiarg
    real(dp) :: argmin, argstepinv
    real(dp), dimension(:), allocatable :: argvals, thetavals

    ! All variables are shared (read-only after init_attenuation_data)

contains

    subroutine init_attenuation_data()
        ! Read attenuation factor data from file once
        ! MUST be called BEFORE any parallel region
        integer, parameter :: iunitarg = 741
        integer :: i
        real(dp) :: arglog, theta_val

        if (attenuation_data_initialized) return

        open (iunitarg, file='thetafun_inp.dat', action='read', status='old')
        npoiarg = 0
        do
            read (iunitarg, *, end=1) arglog
            npoiarg = npoiarg + 1
        end do
        1       close (iunitarg)

        if (allocated(argvals)) deallocate(argvals)
        if (allocated(thetavals)) deallocate(thetavals)
        allocate (argvals(npoiarg), thetavals(npoiarg))

        open (iunitarg, file='thetafun_inp.dat', action='read', status='old')
        do i = 1, npoiarg
            read (iunitarg, *) arglog, theta_val
            argvals(i) = log(arglog)
            thetavals(i) = log(theta_val)
        end do
        close (iunitarg)

        argmin = argvals(1)
        argstepinv = 1.0_dp / (argvals(2) - argvals(1))

        attenuation_data_initialized = .true.
    end subroutine init_attenuation_data

end module thetadata_mod

subroutine attenuation_factor(D, Theta)
    use iso_fortran_env, only: dp => real64
    use thetadata_mod, only: npoiarg, argmin, argstepinv, argvals, thetavals
    use polylag_3, only: mp, indef, plag1d

    implicit none

    real(dp), intent(in) :: D
    real(dp), intent(out) :: Theta
    real(dp) :: arglog, fun, der
    integer, dimension(mp) :: indu
    real(dp), dimension(mp) :: xp, fp

    if (D <= 0.0_dp) then
        Theta = 0.0_dp
        return
    end if

    ! The table resolves the transition, not the asymptotic tails. Avoid
    ! unconstrained cubic extrapolation in log space outside its domain.
    if (D <= exp(argvals(1))) then
        Theta = 1.7555267_dp * D
        return
    else if (D >= exp(argvals(npoiarg))) then
        Theta = 1.0_dp
        return
    end if

    arglog = log(D)

    call indef(arglog, argmin, argstepinv, npoiarg, indu)

    xp = argvals(indu)
    fp = thetavals(indu)

    call plag1d(arglog, fp, argstepinv, xp, fun, der)

    Theta = exp(fun)

end subroutine attenuation_factor
