program potato_profile_readback
    use potato_profile_data_mod, only: load_potato_profiles
    use field_eq_mod, only: psi_axis, psi_sep
    implicit none
    double precision, parameter :: e_charge = 4.8032d-10, ev = 1.6022d-12
    double precision :: energy, charge, energy_scale, field_scale, psi
    double precision, allocatable :: x(:), values(:, :)
    character(len=1024) :: filename, argument
    integer :: unit, count, k

    if (command_argument_count() /= 8) error stop &
        'usage: profile_readback FILE Ealpha Z escale fscale psi_axis psi_sep QUERY_FILE'
    call get_command_argument(1, filename)
    call get_command_argument(2, argument)
    read(argument, *) energy
    call get_command_argument(3, argument)
    read(argument, *) charge
    call get_command_argument(4, argument)
    read(argument, *) energy_scale
    call get_command_argument(5, argument)
    read(argument, *) field_scale
    call get_command_argument(6, argument)
    read(argument, *) psi_axis
    call get_command_argument(7, argument)
    read(argument, *) psi_sep
    if (psi_axis == psi_sep) error stop 'zero poloidal flux span'
    call load_potato_profiles(trim(filename), energy, charge, energy_scale, field_scale, e_charge, ev)
    call get_command_argument(8, filename)
    open(newunit=unit, file=trim(filename), status='old')
    read(unit, *) count
    allocate(x(count), values(7, count))
    do k = 1, count
        read(unit, *) x(k)
    end do
    close(unit)
    !$omp parallel do private(psi)
    do k = 1, count
        psi = psi_axis + x(k)*(psi_sep-psi_axis)
        values(1, k) = x(k)
        call denstemp_of_psi(psi, values(2, k), values(3, k), values(5, k), values(6, k))
        call phielec_of_psi(psi, values(4, k), values(7, k))
    end do
    !$omp end parallel do
    do k = 1, count
        write(*, '(7ES25.16)') values(:, k)
    end do
end program potato_profile_readback
