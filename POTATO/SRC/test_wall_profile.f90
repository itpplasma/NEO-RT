program test_wall_profile
    use libneo_kinds, only: dp
    use wall_loss_mod, only: load_wall, wall_loaded, outside_wall, &
        outside_wall_exact, fastpath_resolves_circle, fastpath_resolves_ellipse
    implicit none

    integer, parameter :: ng = 400, nrep = 200
    real(dp), parameter :: rlo = 80.0_dp, rhi = 240.0_dp
    real(dp), parameter :: zlo = -130.0_dp, zhi = 110.0_dp
    character(len=512) :: wallfile
    integer :: i, k, rep, npts, ncirc, nell, mism
    integer(8) :: acc_fast, acc_exact
    real(dp) :: r, z, t0, t1, t_fast, t_exact

    if (command_argument_count() < 1) then
        print *, 'usage: test_wall_profile.x <convexwall.dat>'
        error stop 1
    end if
    call get_command_argument(1, wallfile)
    call load_wall(trim(wallfile))
    if (.not. wall_loaded) then
        print *, 'could not load wall: ', trim(wallfile)
        error stop 1
    end if

    npts = 0; ncirc = 0; nell = 0; mism = 0
    do i = 1, ng
        r = rlo + (rhi - rlo) * (i - 0.5_dp) / ng
        do k = 1, ng
            z = zlo + (zhi - zlo) * (k - 0.5_dp) / ng
            npts = npts + 1
            if (fastpath_resolves_circle(r, z)) ncirc = ncirc + 1
            if (fastpath_resolves_ellipse(r, z)) nell = nell + 1
            if (outside_wall(r, z) .neqv. outside_wall_exact(r, z)) mism = mism + 1
        end do
    end do

    acc_fast = 0
    call cpu_time(t0)
    do rep = 1, nrep
        do i = 1, ng
            r = rlo + (rhi - rlo) * (i - 0.5_dp) / ng
            do k = 1, ng
                z = zlo + (zhi - zlo) * (k - 0.5_dp) / ng
                if (outside_wall(r, z)) acc_fast = acc_fast + 1
            end do
        end do
    end do
    call cpu_time(t1)
    t_fast = t1 - t0

    acc_exact = 0
    call cpu_time(t0)
    do rep = 1, nrep
        do i = 1, ng
            r = rlo + (rhi - rlo) * (i - 0.5_dp) / ng
            do k = 1, ng
                z = zlo + (zhi - zlo) * (k - 0.5_dp) / ng
                if (outside_wall_exact(r, z)) acc_exact = acc_exact + 1
            end do
        end do
    end do
    call cpu_time(t1)
    t_exact = t1 - t0
    if (acc_fast /= acc_exact) error stop 'profile accumulators differ'

    print '(A,I0)', 'grid points:            ', npts
    print '(A,F6.2,A)', 'circle fast-path:       ', 100.0_dp*ncirc/npts, '% O(1)'
    print '(A,F6.2,A)', 'ellipse fast-path:      ', 100.0_dp*nell/npts, '% O(1)'
    print '(A,F6.2,A)', 'polygon tests, circle:  ', 100.0_dp*(npts-ncirc)/npts, '% of calls'
    print '(A,F6.2,A)', 'polygon tests, ellipse: ', 100.0_dp*(npts-nell)/npts, '% of calls'
    print '(A,ES10.3,A,ES10.3,A,F5.2,A)', 'time fast=', t_fast, 's exact=', t_exact, &
        's  speedup=', t_exact/t_fast, 'x'
    if (mism /= 0) then
        print '(A,I0)', 'FAIL fast-path disagrees with polygon: ', mism
        error stop 1
    end if
    print *, 'fast-path matches exact on all points'
end program test_wall_profile
