program test_resonance_discontinuity
    use iso_fortran_env, only: dp=>real64
    use neort_freq, only: init_canon_freq_trapped_spline, OmtB_spl_coeff, &
        k_taub_t, d_taub_t, k_OmtB_t, d_OmtB_t
    use neort_resonance, only: driftorbit_root, resonance_value
    use neort_profiles, only: Om_tE
    use driftorbit, only: etatp, epst_spl
    use spline, only: spline_coeff
    implicit none
    real(dp) :: join,x(3),y(3),root(2),res,deriv
    call init_canon_freq_trapped_spline()
    Om_tE=0
    join=etatp*(1+epst_spl)
    k_taub_t=-.125_dp
    d_taub_t=5
    k_OmtB_t=0
    d_OmtB_t=-1
    x=[join,join+.125_dp,join+.25_dp]
    ! Lower branch is strictly negative; upper is >=1: no zero exists.
    y=1+x-join
    OmtB_spl_coeff=spline_coeff(x,y)
    root=driftorbit_root(1.0_dp,1e-12_dp,etatp+.5_dp*epst_spl,join+.125_dp)
    print *, 'no-root result:',root
    if(root(1)>=0 .or. root(2)/=0) error stop 'discontinuous bracket accepted as root'
    ! Analytic continuous root, with the same actual native refinement path.
    y=x-join-.05_dp
    OmtB_spl_coeff=spline_coeff(x,y)
    root=driftorbit_root(1.0_dp,1e-12_dp,join+.01_dp,join+.10_dp)
    call resonance_value(1.0_dp,root(1),res,deriv)
    print *, 'continuous root:',root,res
    if(abs(root(1)-(join+.05_dp))>1e-12_dp) error stop 'continuous root lost'
    if(abs(res)>1e-9_dp .or. abs(deriv-1)>1e-12_dp) error stop 'continuous residual/Jacobian'
end program
