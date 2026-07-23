program potato_profile_serializer_probe
    use field_eq_mod, only: psi_axis, psi_sep
    use phielec_of_psi_mod, only: npolyphi, polydens, polyphi, polytemp
    implicit none

    double precision, parameter :: e_charge = 4.8032d-10
    double precision, parameter :: ev = 1.6022d-12
    character(len=512) :: argument, profile_path
    double precision :: charge, ddens, dphi_dpsi, dtemp, energy_ev
    double precision :: density, phi_elec, psi, s_pol, temperature
    integer :: index, input_unit

    call get_command_argument(1, profile_path)
    call get_command_argument(2, argument)
    read(argument, *) psi_axis
    call get_command_argument(3, argument)
    read(argument, *) psi_sep
    call get_command_argument(4, argument)
    read(argument, *) energy_ev
    call get_command_argument(5, argument)
    read(argument, *) charge

    open(newunit=input_unit, file=trim(profile_path), status="old", action="read")
    read(input_unit, *)
    read(input_unit, *)
    read(input_unit, *) polydens
    read(input_unit, *) polytemp
    read(input_unit, *) polytemp
    read(input_unit, *) polyphi
    close(input_unit)

    ! Match the production ingestion in SRC/tt.f90 for unit scaling.
    polytemp = polytemp/energy_ev
    polyphi = polyphi*charge*e_charge/(energy_ev*ev)

    write(*, '(A,2ES25.16)') "POTATO_PROFILE_FLUX ", psi_axis, psi_sep
    write(*, '(A,10ES25.16)') "POTATO_PROFILE_DENSITY_COEFFICIENTS ", polydens
    write(*, '(A,10ES25.16)') "POTATO_PROFILE_TEMPERATURE_COEFFICIENTS ", polytemp
    write(*, '(A,10ES25.16)') "POTATO_PROFILE_POTENTIAL_COEFFICIENTS ", polyphi

    do index = 0, 8
        s_pol = dble(index)/8.0d0
        psi = psi_axis + s_pol*(psi_sep - psi_axis)
        call phielec_of_psi(psi, phi_elec, dphi_dpsi)
        call denstemp_of_psi(psi, density, temperature, ddens, dtemp)
        write(*, '(A,I0,7ES25.16)') "POTATO_PROFILE_SAMPLE ", index, s_pol, psi, &
            density, temperature, phi_elec, dphi_dpsi, ddens
    end do
end program potato_profile_serializer_probe
