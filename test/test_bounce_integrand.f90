program test_bounce_integrand
    ! Dump the Om_tB bounce-average integrand along one thin orbit.
    !
    ! Usage: test_bounce_integrand.x <runname> <ux> <dist>
    !   runname : like neo_rt.x / test_frequencies.x
    !   ux      : v = ux*vth
    !   dist    : (eta - etatp)/etatp of the traced pitch
    !
    ! Writes <runname>_integrand.dat with columns
    !   t/taub  theta  vpar/v  OmtB_dens*v^2 [rad/s]  running_avg*v^2 [rad/s]
    ! and prints taub, Om_tE, Ti1 and the integrated vs splined Om_tB.
    use iso_fortran_env, only: dp => real64
    use do_magfie_mod, only: s, do_magfie
    use neort_lib, only: neort_init, neort_prepare_splines, neort_setup_at_s
    use neort_orbit, only: bounce_time, timestep
    use neort_freq, only: Om_th, Om_tB
    use neort_profiles, only: vth, Om_tE, Ti1
    use driftorbit, only: etatp, etadt, sign_vpar
    use util, only: files_exist, pi

    implicit none

    integer, parameter :: neq = 3, nsteps = 400000, nout = 4000
    character(1024) :: runname, arg
    real(dp) :: ux, dist, v, eta, taub, dt, t
    real(dp) :: y(neq), k1(neq), k2(neq), k3(neq), k4(neq), yt(neq)
    real(dp) :: bmod, sqrtg, x(3), hder(3), hcovar(3), hctrvr(3), hcurl(3)
    real(dp) :: OmtB_spl, Omth_, d1, d2, d3, d4
    integer :: i, fid

    call get_command_argument(1, runname)
    call get_command_argument(2, arg)
    read (arg, *) ux
    call get_command_argument(3, arg)
    read (arg, *) dist

    call neort_init(trim(runname)//".in", "in_file", "in_file_pert")
    if (files_exist("plasma.in", "profile.in")) then
        call neort_prepare_splines("plasma.in", "profile.in")
    end if
    call neort_setup_at_s(s)
    sign_vpar = 1

    v = ux*vth
    eta = etatp*(1.0_dp + dist)
    taub = bounce_time(v, eta)

    x = [s, 0.0_dp, 0.0_dp]
    call do_magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
    y = [0.0_dp, v*sqrt(max(1.0_dp - eta*bmod, 0.0_dp)), 0.0_dp]

    call Om_tB(v, eta, OmtB_spl, d1, d2)
    call Om_th(v, eta, Omth_, d3, d4)
    print *, 'Ti1[eV]=', Ti1, ' vth=', vth, ' E_kin[eV]=', ux**2*Ti1
    print *, 'etatp=', etatp, ' eta=', eta, ' dist=', dist
    print *, 'taub[s]=', taub, ' omega_b=', 2.0_dp*pi/taub, ' (Om_th spl=', Omth_, ')'
    print *, 'Om_tE=', Om_tE, ' Om_tB(spline)=', OmtB_spl, &
        ' Om_ph=', Om_tE + OmtB_spl

    dt = taub/nsteps
    open (newunit=fid, file=trim(runname)//'_integrand.dat')
    write (fid, '(A)') '# t/taub theta vpar/v OmtB_dens*v2[rad/s] runavg*v2[rad/s]'
    t = 0.0_dp
    do i = 1, nsteps
        call timestep(v, eta, neq, t, y, k1)
        if (mod(i, nsteps/nout) == 0) then
            write (fid, *) t/taub, y(1), y(2)/v, k1(3)*v**2, &
                y(3)/max(t, 1.0e-300_dp)*v**2
        end if
        yt = y + 0.5_dp*dt*k1
        call timestep(v, eta, neq, t, yt, k2)
        yt = y + 0.5_dp*dt*k2
        call timestep(v, eta, neq, t, yt, k3)
        yt = y + dt*k3
        call timestep(v, eta, neq, t, yt, k4)
        y = y + dt/6.0_dp*(k1 + 2.0_dp*k2 + 2.0_dp*k3 + k4)
        t = t + dt
    end do
    close (fid)
    print *, 'final theta=', y(1), ' vpar/v=', y(2)/v
    print *, 'Om_tB(integrated)=', y(3)/taub*v**2

end program test_bounce_integrand
