program potato_drift_probe
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use parmot_mod, only: rmu, ro0
    use phielec_of_psi_mod, only: polyphi

    implicit none

    integer, parameter :: nstate = 6
    character(len=16), parameter :: labels(nstate) = [character(len=16) :: &
        'parallel', 'qplus_Ezero', 'qplus_Eplus', 'qplus_Eminus', &
        'qminus_Ezero', 'qminus_Eplus']
    real(dp), parameter :: rho0_probe = 1.0e-3_dp
    real(dp), parameter :: potential_slope = 1.0e10_dp
    real(dp), parameter :: p_probe = 1.0_dp
    real(dp), parameter :: lambda_probe = 0.2_dp

    character(len=1024) :: points_file
    integer :: input_unit, ios, sample, state, status
    real(dp) :: z(5), vz(5), vcart(3), state_ro0(nstate), state_slope(nstate)
    real(dp) :: cphi, sphi

    interface
        subroutine velo(tau, z, vz)
            import dp
            real(dp), intent(in) :: tau, z(5)
            real(dp), intent(out) :: vz(5)
        end subroutine velo
    end interface

    call get_command_argument(1, points_file, status=status)
    if (status /= 0 .or. len_trim(points_file) == 0) then
        error stop 'usage: potato_drift_probe.x POINTS_FILE'
    end if

    state_ro0 = [0.0_dp, rho0_probe, rho0_probe, rho0_probe, &
        -rho0_probe, -rho0_probe]
    state_slope = [0.0_dp, 0.0_dp, potential_slope, -potential_slope, &
        0.0_dp, -potential_slope]
    rmu = 1.0e30_dp

    open(newunit=input_unit, file=trim(points_file), status='old', action='read')
    sample = 0
    do
        read(input_unit, *, iostat=ios) z(1:3)
        if (ios < 0) exit
        if (ios > 0) error stop 'invalid drift-probe point'
        sample = sample + 1
        z(4) = p_probe
        z(5) = lambda_probe
        cphi = cos(z(2))
        sphi = sin(z(2))

        do state = 1, nstate
            ro0 = state_ro0(state)
            polyphi = 0.0_dp
            ! profile_input evaluates polyphi(8)*s_pol + polyphi(9).
            polyphi(8) = state_slope(state)
            call velo(0.0_dp, z, vz)
            vcart = [ &
                vz(1)*cphi - z(1)*vz(2)*sphi, &
                vz(1)*sphi + z(1)*vz(2)*cphi, &
                vz(3)]
            write(*, '(a,1x,i0,1x,a,1x,*(es24.16e3,1x))') &
                'POTATO_DRIFT_STATE', sample, trim(labels(state)), &
                z, ro0, state_slope(state), vz, vcart
        end do
    end do
    close(input_unit)

end program potato_drift_probe
