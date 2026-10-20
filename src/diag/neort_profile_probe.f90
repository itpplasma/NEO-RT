program neort_profile_probe
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use collis_alp, only: collision_ion_mass_g, collision_electron_mass_g, &
        collision_charge_statcoulomb, collision_ev_erg, efcolf, velrat, enrat
    use do_magfie_mod, only: bfac, inp_swi, R0, set_s, s_outside_data, &
        magfie_thread_init, read_boozer_file, init_magfie_at_s
    use neort_profiles, only: read_plasma_input, prepare_plasma_splines, &
        prepare_profile_splines, init_plasma_at_s, init_profile_at_s, nonlin, &
        ni1, dni1ds, ni2, dni2ds, Ti1, dTi1ds, Ti2, dTi2ds, Te, dTeds, &
        M_t, dM_tds, vth, dvthds, Om_tE, dOm_tEds
    use util, only: readdata, mi, qi, mu, ev
    use logger, only: set_log_level
    implicit none

    character(len=4096) :: bc_path, plasma_path, profile_path, query_path
    character(len=4096) :: output_path, argument, line
    real(dp), allocatable :: plasma(:, :), rotation(:, :)
    real(dp) :: am1, am2, z1, z2, efac, lower, upper, coordinate
    integer :: n, input_unit, output_unit, status, rows

    call arguments()
    call load_inputs()
    call scan()

contains

    subroutine scan()
        open (newunit=input_unit, file=trim(query_path), status='old', action='read')
        open (newunit=output_unit, file=trim(output_path), status='replace')
        call header()
        rows = 0
        do
            read (input_unit, '(a)', iostat=status) line
            if (status < 0) exit
            if (status /= 0) error stop 'cannot read query'
            line = adjustl(line)
            if (len_trim(line) == 0) cycle
            if (line(1:1) == '#') cycle
            read (line, *, iostat=status) coordinate
            if (status /= 0) error stop 'invalid radial query'
            call point(coordinate)
            rows = rows + 1
        end do
        close (input_unit)
        close (output_unit)
        if (rows == 0) error stop 'empty query'
    end subroutine scan

    subroutine arguments()
        if (command_argument_count() /= 8) then
            error stop 'usage: probe BC SWITCH PLASMA MACH EFAC BFAC POINTS OUTPUT'
        end if
        call get_command_argument(1, bc_path)
        call get_command_argument(2, argument)
        read (argument, *, iostat=status) inp_swi
        if (status /= 0) error stop 'invalid input switch'
        if (inp_swi /= 8 .and. inp_swi /= 9) error stop 'requires Boozer switch 8/9'
        call get_command_argument(3, plasma_path)
        call get_command_argument(4, profile_path)
        call scale(5, efac)
        call scale(6, bfac)
        if (abs(bfac) <= tiny(bfac)) error stop 'zero BFAC'
        call get_command_argument(7, query_path)
        call get_command_argument(8, output_path)
        if (any([trim(output_path) == trim(bc_path), &
                 trim(output_path) == trim(plasma_path), &
                 trim(output_path) == trim(profile_path), &
                 trim(output_path) == trim(query_path)])) then
            error stop 'output must differ from inputs'
        end if
    end subroutine arguments

    subroutine scale(index, value)
        integer, intent(in) :: index
        real(dp), intent(out) :: value
        call get_command_argument(index, argument)
        read (argument, *, iostat=status) value
        if (status /= 0) error stop 'invalid scale'
        if (.not. ieee_is_finite(value)) error stop 'nonfinite scale'
    end subroutine scale

    subroutine load_inputs()
        call set_log_level(-1)
        call magfie_thread_init()
        call read_boozer_file(trim(bc_path))
        if (.not. ieee_is_finite(R0)) error stop 'nonfinite carrier R0'
        if (R0 <= 0.0_dp) error stop 'nonpositive carrier R0'
        call read_plasma_input(trim(plasma_path), n, am1, am2, z1, z2, plasma)
        if (.not. all(ieee_is_finite([am1, am2, z1, z2]))) then
            error stop 'nonfinite species header'
        end if
        if (min(am1, am2, z1, z2) <= 0.0_dp) then
            error stop 'requires positive background ion masses and charges'
        end if
        call grid(plasma)
        if (any(plasma(:, 2) <= 0.0_dp)) error stop 'nonpositive first density'
        if (any(plasma(:, 3) < 0.0_dp)) error stop 'negative second density'
        if (any(plasma(:, 4:6) <= 0.0_dp)) error stop 'nonpositive temperature'
        call readdata(trim(profile_path), 2, rotation)
        call grid(rotation)
        lower = max(plasma(1, 1), rotation(1, 1))
        upper = min(plasma(n, 1), rotation(size(rotation, 1), 1))
        call prepare_plasma_splines(n, am1, am2, z1, z2, plasma)
        call prepare_profile_splines(rotation)
        nonlin = .false.
    end subroutine load_inputs

    subroutine grid(table)
        real(dp), intent(in) :: table(:, :)
        integer :: count
        count = size(table, 1)
        if (count < 4) error stop 'profile requires four points'
        if (.not. all(ieee_is_finite(table))) error stop 'nonfinite profile'
        if (any(table(2:, 1) <= table(:count - 1, 1))) then
            error stop 'profile grid must strictly increase'
        end if
    end subroutine grid

    subroutine point(s)
        real(dp), intent(in) :: s
        if (.not. ieee_is_finite(s)) error stop 'nonfinite query'
        if (s < lower .or. s > upper) error stop 'outside common profile support'
        if (s_outside_data(s)) error stop 'outside carrier support'
        call set_s(s)
        call init_magfie_at_s()
        call init_plasma_at_s()
        call init_profile_at_s(R0, efac, bfac)
        call row(s)
    end subroutine point

    subroutine scalar(name, value)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: value
        write (output_unit, '(a,es25.16e3)') '# '//name//'=', value
    end subroutine scalar

    subroutine header()
        write (output_unit, '(a)') &
            '# Native profile/collision readback; nonlin=false; no orbits/torque'
        call scalar('am1', am1)
        call scalar('am2', am2)
        call scalar('Z1', z1)
        call scalar('Z2', z2)
        call scalar('efac', efac)
        call scalar('bfac', bfac)
        call scalar('orbit_mass_unit_g', mu)
        call scalar('orbit_ev_erg', ev)
        call scalar('kernel_mass_unit_g', collision_ion_mass_g)
        call scalar('kernel_electron_mass_g', collision_electron_mass_g)
        call scalar('kernel_ev_erg', collision_ev_erg)
        call scalar('kernel_charge_statC', collision_charge_statcoulomb)
        write (output_unit, '(a)') &
            's_tor,n1_cm3,dn1ds,n2_cm3,dn2ds,T1_eV,dT1ds,T2_eV,dT2ds,'// &
            'Te_eV,dTeds,Mach,dMachds,vth_cm_s,dvthds,Omega_rad_s,dOmegads,'// &
            'orbit_m1_g,orbit_q1_statC,ne_cm3,Zeff,pressure_erg_cm3,'// &
            'kernel_m1_g,kernel_m2_g,kernel_me_g,efcolf1,efcolf2,efcolf3,'// &
            'velrat1,velrat2,velrat3,enrat1,enrat2,enrat3,R0_cm'
    end subroutine header

    subroutine row(s)
        real(dp), intent(in) :: s
        real(dp) :: values(35), ne, zeff, pressure
        ne = z1 * ni1 + z2 * ni2
        zeff = (z1**2 * ni1 + z2**2 * ni2) / ne
        pressure = (ne * Te + ni1 * Ti1 + ni2 * Ti2) * ev
        values(1:17) = [s, ni1, dni1ds, ni2, dni2ds, Ti1, dTi1ds, Ti2, &
                        dTi2ds, Te, dTeds, M_t, dM_tds, vth, dvthds, &
                        Om_tE, dOm_tEds]
        values(18:25) = [mi, qi, ne, zeff, pressure, &
                        am1 * collision_ion_mass_g, &
                        am2 * collision_ion_mass_g, collision_electron_mass_g]
        values(26:35) = [efcolf, velrat, enrat, R0]
        if (.not. all(ieee_is_finite(values))) error stop 'nonfinite native result'
        write (output_unit, '(34(es25.16e3,","),es25.16e3)') values
    end subroutine row
end program neort_profile_probe
