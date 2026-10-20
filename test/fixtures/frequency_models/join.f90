program test_frequency_join_manufactured
    use iso_fortran_env, only: dp=>real64
    use ieee_arithmetic, only: ieee_next_after
    use neort_freq, only: init_canon_freq_trapped_spline, init_canon_freq_passing_spline, Om_th, Om_tB, Om_ph
    use driftorbit, only: etatp, epst_spl, epsp_spl, sign_vpar, magdrift, magdrift_passing
    use neort_profiles, only: Om_tE
    use do_magfie_mod, only: iota
    use util, only: pi
    implicit none
    integer :: branch,sgn,side
    real(dp) :: join,eta,w,b,p,dv,de,tau,oracle_w,oracle_b,v
    call init_canon_freq_trapped_spline()
    call init_canon_freq_passing_spline()
    v=1.7_dp
    do branch=-1,1,2
        join=etatp*(1+merge(epst_spl,-epsp_spl,branch==1))
        do sgn=-1,1,2
            sign_vpar=real(sgn,dp)
            do side=-1,1,2
                eta=ieee_next_after(join,join+real(side,dp))
                tau=5-.125_dp*log(abs(eta-etatp))
                oracle_w=sgn*2*pi*v/tau
                oracle_b=(1.25_dp+.05_dp*log(abs(eta-etatp)))*2*pi*v**2/tau
                call Om_th(v,eta,w,dv,de)
                call Om_tB(v,eta,b,dv,de)
                call Om_ph(v,eta,p,dv,de)
                print *, branch,sgn,side,w,oracle_w,b,oracle_b
                if(abs(w-oracle_w)>1e-10_dp) error stop 'poloidal logarithmic join'
                if(abs(b-oracle_b)>1e-10_dp) error stop 'drift logarithmic join'
                if(branch==-1) oracle_b=oracle_b+oracle_w/iota
                if(abs(p-(Om_tE+oracle_b))>1e-10_dp) error stop 'toroidal logarithmic join'
            end do
        end do
    end do
    ! The passing switch is independent of the trapped magnetic-drift switch.
    magdrift=.false.
    magdrift_passing=1
    sign_vpar=1
    eta=etatp*(1-.5_dp*epsp_spl)
    call Om_th(v,eta,w,dv,de)
    call Om_tB(v,eta,b,dv,de)
    call Om_ph(v,eta,p,dv,de)
    if(abs(b)<1e-3_dp) error stop 'passing drift control is vacuous'
    if(abs(p-(Om_tE+w/iota+b))>1e-10_dp) error stop 'independent passing drift switch'
end program
