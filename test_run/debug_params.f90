program debug_params
    use driftorbit, only: vth, B0, s, eps, R0, iota
    use do_magfie_mod, only: init_magfie => init

    implicit none
    real(8) :: v, eta, bmod, taub_estimate, vperp_val
    
    s = 0.3d0
    vth = 4.0d7
    
    call init_magfie()
    
    \! Check critical parameters
    write(*,*) 'R0 =', R0
    write(*,*) 'eps =', eps  
    write(*,*) 'iota =', iota
    write(*,*) 'B0 =', B0
    write(*,*) 'vth =', vth
    
    \! Test bounce time estimate
    v = 0.5d0
    eta = 0.5d0
    bmod = B0
    vperp_val = vth * v * sqrt(1d0 - eta**2)
    
    if (eps > 0d0 .and. R0 > 0d0 .and. iota /= 0d0) then
        taub_estimate = 2.0*3.14159/abs(vperp_val*iota/R0*sqrt(eps/2d0))
        write(*,*) 'Bounce time estimate:', taub_estimate
    else
        write(*,*) 'Invalid parameters\!'
    end if
    
end program debug_params
EOF < /dev/null
