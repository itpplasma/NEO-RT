module profile_exit_oracle_mod
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use potato_profile_data_mod, only: evaluate_profile_table
    implicit none
contains
    subroutine circular_velocity(tau, z, vz)
        real(dp), intent(in) :: tau, z(5)
        real(dp), intent(out) :: vz(5)
        real(dp) :: value, derivative

        ! Exact circular orbit; its profile label spans [-0.2,1.2] and returns
        ! inside support. Sampling a finite edge must not erase that excursion.
        call evaluate_profile_table(0.5_dp+0.7_dp*z(3), 3, value, derivative)
        vz = 0.0_dp
        vz(1) = -z(3)
        vz(2) = 0.25_dp
        vz(3) = z(1)-3.0_dp
    end subroutine circular_velocity
end module profile_exit_oracle_mod

program test_profile_exit
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use potato_profile_data_mod, only: load_potato_profiles, evaluate_profile_table, &
        reset_profile_support, profile_support_left
    use potato_input_mod, only: profile_exit_absorb, potato_nml, &
        canonical_axis_constant_extension, zero_canonical_axis_drive
    use profile_exit_oracle_mod, only: circular_velocity
    use poicut_mod, only: npc, rpc_beg, h_rpc, rpc_arr, zpc_arr
    use field_eq_mod, only: ierrfield
    implicit none
    real(dp) :: value, derivative, coordinate, expected
    real(dp) :: z(5), extra(0), taub, delphi
    integer :: unit, i, column, side, failures, ierr
    character(len=16) :: mode

    call get_command_argument(1, mode)
    if (profile_exit_absorb) error stop 'absorption must default off'
    open(newunit=unit, file='profile-exit.table', status='replace')
    write(unit, '(A)') 'POTATO_PROFILE_TABLE_V1'
    write(unit, *) 5
    do i = 0, 4
        write(unit, *) real(i,dp)/4, 2.0_dp+real(i,dp)/4, &
            1.0_dp+real(i,dp)/4, 3.0_dp+real(i,dp)/4
    end do
    close(unit)
    call load_potato_profiles('profile-exit.table', 1.0_dp, 1.0_dp, &
        1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp)
    if (trim(mode) == 'off') then
        call evaluate_profile_table(1.1_dp, 3, value, derivative)
        error stop 'default accepted outside support'
    endif
    ! Exercise the actual namelist registration as well as the default.
    open(newunit=unit, status='scratch')
    write(unit, '(A)') '&potato_nml profile_exit_absorb=.true. /'
    rewind(unit)
    read(unit, nml=potato_nml)
    close(unit)
    if (.not. profile_exit_absorb) error stop 'namelist did not enable absorption'
    do column = 1, 3
        do side = 0, 1
            call reset_profile_support()
            coordinate = -0.1_dp+1.2_dp*real(side,dp)
            call evaluate_profile_table(coordinate, column, value, derivative)
            expected = real(column,dp)+real(side,dp)
            if (column == 1) expected = 2.0_dp+real(side,dp)
            if (column == 2) expected = 1.0_dp+real(side,dp)
            if (abs(value-expected) > 1.0e-12_dp) error stop 'wrong edge value'
            if (.not. ieee_is_finite(value)) error stop 'nonfinite edge value'
            if (derivative /= 0.0_dp) error stop 'edge must be constant'
            if (.not. profile_support_left) error stop 'profile exit not recorded'
            call evaluate_profile_table(0.5_dp, column, value, derivative)
            if (.not. profile_support_left) error stop 'exit flag was not sticky'
        enddo
    enddo
    call reset_profile_support()
    call evaluate_profile_table(1.0_dp+0.5e-12_dp, 3, value, derivative)
    if (profile_support_left) error stop 'roundoff tolerance changed'
    if (abs(derivative-1.0_dp) > 1.0e-12_dp) error stop 'interior slope changed'

    failures = 0
    !$omp parallel num_threads(4) private(value,derivative) reduction(+:failures)
    ! A nested inactive region matches the energy -> class/mode call context.
    !$omp parallel if(.false.) private(value,derivative)
    call reset_profile_support()
    !$omp end parallel
    !$omp do schedule(static,1)
    do i = 1, 4
        if (mod(i,2) == 0) then
            call evaluate_profile_table(1.1_dp, 3, value, derivative)
        else
            call evaluate_profile_table(0.5_dp, 3, value, derivative)
        endif
    enddo
    !$omp end do
    ! The barrier ensures a shared loss flag cannot accidentally pass.
    !$omp do schedule(static,1)
    do i = 1, 4
        if (profile_support_left .neqv. (mod(i,2) == 0)) failures=failures+1
        call reset_profile_support()
        if (profile_support_left) failures=failures+1
    enddo
    !$omp end do
    !$omp end parallel
    if (failures /= 0) error stop 'profile loss leaked between OpenMP workers'

    npc = 4
    rpc_beg = 1.0_dp
    h_rpc = 1.0_dp
    allocate(rpc_arr(0:npc), zpc_arr(0:npc))
    rpc_arr = [(real(i,dp)+1.0_dp, i=0,npc)]
    zpc_arr = 0.0_dp
    ierrfield = 0
    z = [4.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.7_dp]
    call reset_profile_support()
    call find_bounce(0, circular_velocity, 0.02_dp, z, taub, delphi, extra, ierr)
    if (ierr /= 0) error stop 'manufactured orbit did not close'
    if (abs(taub-2*acos(-1.0_dp)) > 2.0e-6_dp) error stop 'wrong orbit period'
    if (.not. profile_support_left) error stop 'orbit excursion was missed'
    ! A lost physical orbit must not qualify for the canonical-axis shortcut.
    canonical_axis_constant_extension = .true.
    if (zero_canonical_axis_drive(-0.1_dp,0.0_dp,1.0_dp, &
                                  .not. profile_support_left)) &
        error stop 'canonical closure admitted a lost physical orbit'
    call reset_profile_support()
    call evaluate_profile_table(0.5_dp, 3, value, derivative)
    if (profile_support_left) error stop 'next confined orbit inherited loss'
    print *, 'profile absorption: edges, orbit, namelist and OpenMP isolation PASS'
end program test_profile_exit
