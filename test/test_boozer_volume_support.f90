program test_boozer_volume_support
    ! Synthetic circular geometry with constant |B|, iota and currents.
    use iso_fortran_env, only: dp => real64
    use do_magfie_mod, only: inp_swi, read_boozer_file, set_s, init_magfie_at_s
    implicit none
    character(len=32) :: mode
    integer :: unit, k
    real(dp) :: radius, dvds, surface
    real(dp), parameter :: pi = acos(-1.0_dp)

    call get_command_argument(1, mode)
    inp_swi = 8
    open (newunit=unit, file='synthetic.bc', status='replace')
    do k = 1, 5
        write (unit, '(a)') 'synthetic circular Boozer data'
    end do
    write (unit, *) 1, 0, 5, 1, -1.0_dp, 1.0_dp, 3.0_dp
    do k = 1, 5
        surface = real(k, dp)/5.0_dp
        radius = sqrt(surface)
        ! 4*pi**2 * sign_theta * flux/(2*pi) * ItoB * I_phi / B**2
        ! SI volume: -0.4*pi*1e-6 for flux=-1, I_phi=B=1.
        dvds = -0.4_dp*pi*1.0e-6_dp
        if (k == 5) then
            dvds = dvds*1.2_dp
            if (trim(mode) == 'sign') dvds = -dvds
        end if
        write (unit, '(a)') 's iota I_phi I_theta pprime dVds'
        write (unit, '(a)') 'synthetic units'
        write (unit, *) surface, 0.5_dp, 1.0_dp, 0.0_dp, 0.0_dp, dvds
        write (unit, '(a)') 'm n R Z shift B'
        write (unit, *) 0, 0, 3.0_dp, 0.0_dp, 0.0_dp, 1.0_dp
        write (unit, *) 1, 0, radius, radius, 0.0_dp, 0.0_dp
    end do
    close (unit)
    call read_boozer_file('synthetic.bc')
    select case (trim(mode))
    case ('outside', 'sign')
        call set_s(0.3_dp)
    case ('inside')
        call set_s(1.0_dp)
    case ('support')
        call set_s(0.9_dp)
    case ('later')
        call set_s(0.3_dp)
        call init_magfie_at_s()
        call set_s(0.9_dp)
    case default
        error stop 'unknown test mode'
    end select
    call init_magfie_at_s()
    print *, 'accepted outside mismatch'
end program test_boozer_volume_support
