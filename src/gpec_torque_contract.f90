module gpec_torque_contract
    ! Strict access to the integrated torque reported by GPEC_CONTROL.
    !
    ! GPEC_CONTROL does not contain a volume perturbed-current vector or a
    ! radial JxB profile. Its scalar torque can therefore check only an
    ! independently integrated result, not the electromagnetic-torque kernel.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

    implicit none
    private

    type, public :: gpec_control_torque
        integer :: toroidal_mode = 0
        real(dp) :: torque_nm = 0.0_dp
    end type gpec_control_torque

    public :: read_gpec_control_torque
    public :: require_gpec_volume_current

contains

    subroutine read_gpec_control_torque(filename, toroidal_mode, result)
        character(*), intent(in) :: filename
        integer, intent(in) :: toroidal_mode
        type(gpec_control_torque), intent(out) :: result
        character(512) :: line
        integer :: unit, io_status, separator, torque_count
        logical :: header_found

        if (toroidal_mode == 0) error stop "GPEC toroidal mode must be nonzero"

        open(newunit=unit, file=filename, status="old", action="read", iostat=io_status)
        if (io_status /= 0) error stop "cannot open GPEC_CONTROL file"

        header_found = .false.
        torque_count = 0
        do
            read(unit, "(a)", iostat=io_status) line
            if (io_status < 0) exit
            if (io_status > 0) error stop "cannot read GPEC_CONTROL file"

            if (index(adjustl(line), "GPEC_CONTROL:") == 1) header_found = .true.
            if (index(adjustl(line), "toroidal torque") == 1) then
                separator = index(line, "=")
                if (separator == 0) error stop "malformed GPEC toroidal torque"
                torque_count = torque_count + 1
                read(line(separator + 1:), *, iostat=io_status) result%torque_nm
                if (io_status /= 0) error stop "malformed GPEC toroidal torque"
            end if
        end do
        close(unit)

        if (.not. header_found) error stop "not a GPEC_CONTROL file"
        if (torque_count /= 1) error stop "GPEC_CONTROL must contain one toroidal torque"
        if (.not. ieee_is_finite(result%torque_nm)) then
            error stop "GPEC toroidal torque must be finite"
        end if
        result%toroidal_mode = toroidal_mode
    end subroutine read_gpec_control_torque

    subroutine require_gpec_volume_current()
        error stop "GPEC output has no volume perturbed-current vector for JxB"
    end subroutine require_gpec_volume_current

end module gpec_torque_contract
