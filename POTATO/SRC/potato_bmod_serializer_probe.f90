program potato_bmod_serializer_probe
    implicit none

    character(len=512) :: points_path
    complex(8) :: amplitude
    double precision :: radius, height
    integer :: index, input_unit, status

    call get_command_argument(1, points_path)
    open(newunit=input_unit, file=trim(points_path), status="old", action="read")

    index = 0
    do
        read(input_unit, *, iostat=status) radius, height
        if (status < 0) exit
        if (status > 0) error stop "invalid POTATO bmod probe point"
        index = index + 1
        call bmod_pert(radius, height, amplitude)
        write(*, '(A,I0,4ES25.16)') "POTATO_BMOD_SAMPLE ", index, radius, height, &
            real(amplitude), aimag(amplitude)
    end do
    close(input_unit)
end program potato_bmod_serializer_probe
