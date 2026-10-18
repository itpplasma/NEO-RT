module bounce_oracle_mod
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp) :: direction
    logical :: trapped = .false.
contains
    subroutine circular_velocity(tau, z, vz)
        real(dp), intent(in) :: tau, z(5)
        real(dp), intent(out) :: vz(5)

        vz = 0.0_dp
        if (trapped) then
            ! A trapped oscillation visits the same RZ after half a cycle,
            ! with opposite pitch. Its full phase-space period is 2*pi.
            vz(1) = -direction*z(5)
            vz(2) = 0.25_dp
            vz(5) = direction*(z(1) - 3.0_dp)
            return
        end if
        vz(1) = -direction*z(3)
        vz(2) = 0.25_dp
        vz(3) = direction*(z(1) - 3.0_dp)
        vz(5) = -0.7_dp*direction*z(3)
    end subroutine circular_velocity
end module bounce_oracle_mod

program test_bounce_offcut
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use bounce_oracle_mod, only: dp, direction, circular_velocity, trapped
    use orbit_dim_mod, only: write_orb, iunit1
    use poicut_mod, only: npc, rpc_beg, h_rpc, rpc_arr, zpc_arr
    use field_eq_mod, only: ierrfield
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: steps(3) = [0.1_dp, 0.01_dp, 0.02_dp]
    real(dp), parameter :: phases(3) = [0.0_dp, 0.03_dp, -0.03_dp]
    real(dp) :: z(5), last(6), sample(6), extra(0), taub, delphi
    real(dp) :: theta, previous, winding, jump, phase, closure, expected_winding
    integer :: i, j, k, n, ierr, ios, failures, orbit_kind
    logical :: valid
    external :: find_bounce

    npc = 4
    rpc_beg = 1.0_dp
    h_rpc = 1.0_dp
    allocate (rpc_arr(0:npc), zpc_arr(0:npc))
    rpc_arr = [(real(i, dp) + 1.0_dp, i=0, npc)]
    zpc_arr = 0.0_dp
    ierrfield = 0
    write_orb = .true.
    failures = 0
    do orbit_kind = 1, 2
        trapped = orbit_kind == 2
        do k = -1, 1, 2
            direction = real(k, dp)
            expected_winding = merge(0.0_dp, direction, trapped)
            do j = 1, merge(1, size(phases), trapped)
                phase = phases(j)
                z = [3.0_dp + cos(phase), 0.0_dp, sin(phase), &
                     1.0_dp, 0.7_dp*cos(phase)]
                if (trapped) z = [3.0_dp, 0.0_dp, 0.1_dp, 1.0_dp, 1.0_dp]
                do i = 1, size(steps)
                    open (newunit=iunit1, status='scratch', form='formatted')
                    call find_bounce(0, circular_velocity, steps(i), z, &
                                     taub, delphi, extra, ierr)
                    rewind (iunit1)
                    n = 0
                    winding = 0.0_dp
                    do
                        read (iunit1, *, iostat=ios) sample
                        if (ios /= 0) exit
                        if (ieee_is_nan(sample(1))) exit
                        theta = atan2(sample(3), sample(1) - 3.0_dp)
                        if (n > 0) then
                            jump = modulo(theta - previous + pi, 2*pi) - pi
                            winding = winding + jump/(2*pi)
                        end if
                        previous = theta
                        last = sample
                        n = n + 1
                    end do
                    close (iunit1)
                    ! The exact solution returns after 2*pi. Check the whole state,
                    ! including the sign-changing pitch, rather than ierr alone.
                    closure = maxval(abs(last([1, 3, 4, 5]) - z([1, 3, 4, 5])))
                    valid = ierr == 0 .and. n > 2
                    valid = valid .and. abs(taub - 2*pi) < 2.0e-6_dp
                    valid = valid .and. abs(delphi - 0.5_dp*pi) < 5.0e-7_dp
                    valid = valid .and. closure < 2.0e-6_dp
                    valid = valid .and. &
                            abs(winding - expected_winding) < 5.0e-7_dp
                    write (*, '(a,l1,a,i2,3(a,es13.5),a,l1)') &
                        'trapped=', trapped, ' direction=', k, &
                        ' phase=', phase, ' step=', steps(i), &
                        ' period=', taub, ' full_return=', valid
                    if (.not. valid) failures = failures + 1
                end do
            end do
        end do
    end do
    write_orb = .false.
    if (failures /= 0) error stop 'find_bounce did not return to its starting state'
    print *, 'test_bounce_offcut: PASS (24 independent full-return oracles)'
end program test_bounce_offcut
