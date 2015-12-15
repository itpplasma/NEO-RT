! Resonant regimes in tokamaks
! in the action-angle formalism
! Christopher Albert, 2015

MODULE driftorbit
  USE do_magfie_mod
  USE neo_magfie_perturbation, ONLY: neo_read_pert_control, neo_read_pert,&
       neo_init_spline_pert, neo_magfie_pert_amp, m_phi
  USE spline
  
  IMPLICIT NONE
  SAVE

  INTEGER, PARAMETER :: nvar = 6

  ! Orbit parameters
  REAL(8) :: v, eta
  REAL(8) :: taub, bounceavg(nvar)

  ! Plasma parameters
  REAL(8) :: Om_tE, vth, M_t, n0

  ! For splining in the trapped eta region
  INTEGER, PARAMETER :: netaspl = 100
  REAL(8) :: OmtB_spl_coeff(netaspl-1, 5)
  REAL(8) :: Omth_spl_coeff(netaspl-1, 5)
  REAL(8) :: vres_spl_coeff(netaspl-1, 5)

  ! For splining in the passing eta region
  INTEGER, PARAMETER :: netaspl_pass = 100
  REAL(8) :: OmtB_pass_spl_coeff(netaspl_pass-1, 5)
  REAL(8) :: Omth_pass_spl_coeff(netaspl_pass-1, 5)
  REAL(8) :: vres_pass_spl_coeff(netaspl-1, 5)

  ! Harmonics TODO: make dynamic, multiple harmonics
  REAL(8) :: epsmn        ! perturbation amplitude B1/B0
  INTEGER :: m0, mph      ! Boozer toroidal and poloidal perturbation mode
  INTEGER :: mth          ! canonical poloidal mode
  LOGICAL :: supban       ! calculate superbanana plateau
  LOGICAL :: magdrift     ! consider magnetic drift 
  LOGICAL :: nopassing    ! neglect passing particles
  LOGICAL :: noshear      ! neglect magnetic shear
  LOGICAL :: calcflux     ! calculation based on flux instead of torque
  LOGICAL :: pertfile     ! read perturbation from file with neo_magfie_pert

  ! Flux surface TODO: make a dynamic, multiple flux surfaces support
  REAL(8) :: dVds, etadt, etatp
  REAL(8) :: etamin, etamax
  REAL(8) :: vmin2, vmax2 ! permanent vmin vmax for integration bracketing

  ! Check if init is done
  LOGICAL :: init_done

  ! TODO: better B0 calculation (average magnetic field on flux surface)
  REAL(8) :: B0
  REAL(8) :: Bmin, Bmax, th0

  REAL(8), PARAMETER :: epst_spl=1d-6, epsp_spl=1d-6   ! dist to tpb for spline
  REAL(8), PARAMETER :: epsst_spl=1d-3, epssp_spl=1d-3 ! dist to deep for spline
  REAL(8), PARAMETER :: epst=1d-8, epsp=1d-8 ! smallest eta distance to tp bound
  REAL(8) :: k_taub_p, d_taub_p, k_taub_t, d_taub_t ! extrapolation at tp bound
  REAL(8) :: k_OmtB_p, d_Omtb_p, k_Omtb_t, d_Omtb_t ! extrapolation at tp bound

  ! Number of levels for coarse root finding
  INTEGER, PARAMETER :: nlev = 100
CONTAINS

  SUBROUTINE init
    init_done = .FALSE.
    CALL do_magfie_init
    IF (pertfile) THEN
       CALL neo_read_pert_control
       CALL neo_read_pert
       CALL neo_init_spline_pert
       mph = m_phi
    END IF
    CALL init_fsa
    CALL init_misc
    CALL init_Om_spl
    CALL init_Om_pass_spl
    init_done = .TRUE.
  END SUBROUTINE init

  SUBROUTINE init_Om_spl
    ! Initialise splines for canonical frequencies of trapped orbits
    
    REAL(8) :: etarange(netaspl), Om_tB_v(netaspl), Omth_v(netaspl)
    INTEGER :: k
    REAL(8) :: aa, b
    REAL(8) :: taub0, taub1, leta0, leta1, OmtB0, OmtB1
    
    taub0 = 0d0
    taub1 = 0d0
    leta0 = 0d0
    leta1 = 0d0
    OmtB0 = 0d0
    OmtB1 = 0d0
    
    v = vth    
    etamin = etatp
    etamax = etatp + (etadt-etatp)*(1d0-epsst_spl)

    b = LOG(epst_spl)
    aa = 1d0/(netaspl-1d0)*(LOG(etamax/etamin - 1d0) - b)
    
    DO k = netaspl-1, 0, -1
       !eta = etamin * (1d0 + (k/(netaspl-1d0))**4*(etamax/etamin-1d0))
       eta = etamin*(1d0 + EXP(aa*k+b))
       etarange(k+1) = eta
       IF (k == netaspl-1) THEN
          CALL bounce
       ELSE
          CALL bounce(taub)
       END IF
       IF (magdrift) Om_tB_v(k+1) = bounceavg(3)
       Omth_v(k+1) = 2*pi/(v*taub)
       IF (k==0) THEN
          leta0 = LOG(eta-etatp)
          taub0 = v*taub
          IF (magdrift) OmtB0 = Om_tB_v(k+1)/Omth_v(k+1)
       END IF
       IF (k==1) THEN
          leta1 = LOG(eta-etatp)
          taub1 = v*taub
          IF (magdrift) OmtB1 = Om_tB_v(k+1)/Omth_v(k+1)
       END IF
    END DO
    
    k_taub_t = (taub1-taub0)/(leta1-leta0)
    d_taub_t = taub0 - k_taub_t*leta0
    Omth_spl_coeff = spline_coeff(etarange, Omth_v)
        
    IF (magdrift) THEN
       k_OmtB_t = (OmtB1-OmtB0)/(leta1-leta0)
       d_OmtB_t = OmtB0 - k_OmtB_t*leta0
       OmtB_spl_coeff = spline_coeff(etarange, Om_tB_v)
    END IF

    
  END SUBROUTINE init_Om_spl
  
  SUBROUTINE init_Om_pass_spl
    ! Initialise splines for canonical frequencies of passing orbits
    
    REAL(8) :: etarange(netaspl_pass), Om_tB_v(netaspl_pass), Omth_v(netaspl_pass)
    REAL(8) :: aa, b
    INTEGER :: k
    REAL(8) :: leta0, leta1, taub0, taub1, OmtB0, OmtB1
    taub0 = 0d0
    taub1 = 0d0
    leta0 = 0d0
    leta1 = 0d0
    OmtB0 = 0d0
    OmtB1 = 0d0
    
    v = vth
    !etamin = etatp*1d-9
    !etamax = etatp*(1d0-1d-9)
    etamin = etatp*epssp_spl
    etamax = etatp
    
    b = LOG((etamax-etamin)/etamax)
    aa = 1d0/(netaspl_pass-1d0)*(LOG(epsp_spl) - b)
    
    DO k = netaspl_pass-1, 0, -1
       !eta = etamin * (1d0 + k/(netaspl_pass-1d0)*(etamax/etamin-1d0))
       eta = etamax * (1d0 - EXP(aa*k+b))
       etarange(k+1) = eta
       IF (k == netaspl_pass-1) THEN
          CALL bounce
       ELSE
          CALL bounce(taub)
       END IF
       IF (magdrift) Om_tB_v(k+1) = bounceavg(3)
       Omth_v(k+1) = 2*pi/(v*taub)
       IF (k==netaspl_pass-2) THEN
          leta0 = LOG(etatp-eta)
          taub0 = v*taub
          IF (magdrift) OmtB0 = Om_tB_v(k+1)/Omth_v(k+1)
       END IF
       IF (k==netaspl_pass-1) THEN
          leta1 = LOG(etatp-eta)
          taub1 = v*taub
          IF (magdrift) OmtB1 = Om_tB_v(k+1)/Omth_v(k+1)
       END IF
    END DO

    k_taub_p = (taub1-taub0)/(leta1-leta0)
    d_taub_p = taub0 - k_taub_p*leta0
    Omth_pass_spl_coeff = spline_coeff(etarange, Omth_v)
    
    IF (magdrift) THEN
       k_OmtB_p = (OmtB1-OmtB0)/(leta1-leta0)
       d_OmtB_p = OmtB0 - k_OmtB_p*leta0
       OmtB_pass_spl_coeff = spline_coeff(etarange, Om_tB_v)
    END IF
  END SUBROUTINE init_Om_pass_spl

  SUBROUTINE init_misc
    ! TODO: fine search for minima and maxima
    etatp = 1d0/Bmax
    etadt = 1d0/Bmin
  END SUBROUTINE init_misc
  
  FUNCTION Jperp()
    REAL(8) :: Jperp
    Jperp = mi**2*c/(2d0*qi)*eta*v**2
  END FUNCTION Jperp
  
  FUNCTION vpar(bmod)
    !   parallel velocity
    REAL(8) :: bmod, vpar
    vpar = v*SQRT(1 - eta*bmod)
    IF (isnan(vpar)) THEN
       vpar = 0d0
    END IF
  END FUNCTION vpar
  
  FUNCTION vperp(bmod)
    !   perpendicular velocity
    REAL(8) :: bmod, vperp
    vperp = v*SQRT(eta*bmod)
    IF (isnan(vperp)) THEN
       vperp = 0d0
    END IF
  END FUNCTION vperp
  
  SUBROUTINE jac 
  END SUBROUTINE jac

  SUBROUTINE timestep2(neq, t, y, ydot)
    INTEGER, INTENT (in) :: neq
    REAL(8), INTENT (in) :: t
    REAL(8), INTENT (in) :: y(neq)
    REAL(8), INTENT (out) :: ydot(neq)

    CALL timestep(t, y, ydot)
  END SUBROUTINE timestep2
  
  SUBROUTINE timestep(t, y, ydot)
    !
    !  Timestep function for orbit integration.
    !  Includes poloidal angle theta and parallel velocity.
    !  More integrands may be added starting from y(3)
    !

    !integer, intent (in) :: neq
    REAL(8), INTENT (in) :: t
    REAL(8), INTENT (in) :: y(*)
    REAL(8), INTENT (out) :: ydot(*)
    
    REAL(8) :: bmod, sqrtg, x(3), hder(3), hcovar(3), hctrvr(3), hcurl(3)
    REAL(8) :: Om_tB_v
    REAL(8) :: Omth, dOmthdv, dOmthdeta
    REAL(8) :: t0
    REAL(8) :: shearterm
    COMPLEX(8) :: epsn, Hn ! relative amplitude of perturbation field epsn=Bn/B0
                           ! and Hamiltonian Hn = (H - H0)_n
    
    x(1) = s
    x(2) = 0d0
    x(3) = y(1)
    CALL do_magfie( x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl )
    CALL Om_th(Omth, dOmthdv, dOmthdeta)
    
    IF (pertfile) THEN
       CALL neo_magfie_pert_amp( x, epsn )
       epsn = epsn/bmod
    ELSE
       epsn = epsmn*EXP(imun*m0*y(1))
    END IF

    shearterm = Bphcov*dqds
    IF (noshear) THEN
       shearterm = 0
    END IF

!!$    Om_tB_v = mi*c/(2d0*qi*sqrtg*hctrvr(3)*bmod**2)*(&      ! Om_tB/v**2
!!$         -(2d0-eta*bmod)*bmod*hder(1)&
!!$         +2d0*(1d0-eta*bmod)*hctrvr(3)*&
!!$         (dBthcovds+q*dBphcovds+shearterm))

    Om_tB_v = mi*c*q/(2d0*qi*psi_pr*bmod)*(&      ! Om_tB/v**2
         -(2d0-eta*bmod)*bmod*hder(1)&
         +2d0*(1d0-eta*bmod)*hctrvr(3)*&
         (dBthcovds+q*dBphcovds+shearterm))

    !print *, c*mi*vth**2/(2*qi*psi_pr), Om_tB_v*vth**2
   
    ydot(1) = y(2)*hctrvr(3)                                    ! theta
    ydot(2) = -v**2*eta/2d0*hctrvr(3)*hder(3)*bmod              ! v_par
    ydot(3) = Om_tB_v                                           ! v_ph

    IF (init_done) THEN
       IF (eta > etatp) THEN
          !t0 = 0.25*2*pi/Omth ! Quarter of bounce time to set phase 0 at tips
          t0 = 0d0
          Hn = (2d0-eta*bmod)*epsn*EXP(imun*(q*mph*(y(1))-mth*(t-t0)*Omth))
          !ydot(4) = (2d0 - eta*bmod)*epsmn*&
          !     cos((m0+q*mph)*y(1)-mth*(t-t0)*Omth) ! Re Hn
          !ydot(5) = (2d0 - eta*bmod)*epsmn*&
          !     sin((m0+q*mph)*y(1)-mth*(t-t0)*Omth) ! Im Hn
       ELSE
          Hn = (2d0-eta*bmod)*epsn*EXP(imun*(q*mph*(y(1))-(mth+q*mph)*t*Omth))
          !ydot(4) = (2d0 - eta*bmod)*epsmn*&
          !     cos((m0+q*mph)*y(1)-(mth+q*mph)*t*Omth) ! Re Hn
          !ydot(5) = (2d0 - eta*bmod)*epsmn*&
          !     sin((m0+q*mph)*y(1)-(mth+q*mph)*t*Omth) ! Im Hn
       END IF
       ydot(4) = REAL(Hn)
       ydot(5) = AIMAG(Hn)
    ELSE
       ydot(4:6) = 0d0
    END IF
  END SUBROUTINE timestep
  
  FUNCTION findroot2(y0, dt, tol)
    !
    !  Finds the root of an orbit after the first turn
    !
    USE dvode_f90_m
    
    REAL(8) :: y0(nvar), dt, tol, findroot2(nvar+1)
    INTEGER :: n

    INTEGER :: k, state, rootstate
    REAL(8) :: ti, told
    REAL(8) :: y(nvar), yold(nvar)

    LOGICAL :: passing
    
    REAL(8) :: atol(nvar), rtol, tout
    INTEGER :: neq, itask, istate
    TYPE (vode_opts) :: options

    neq = nvar
    rtol = 1d-9
    atol = 1d-10
    itask = 1
    istate = 1
    options = set_normal_opts(abserr_vector=atol, relerr=rtol, nevents=2)
    
    ! check for passing orbit
    passing = .FALSE.
    IF (eta < etatp) THEN
       passing = .TRUE.
    END IF

    n = 500
    rootstate = -1
    
    y = y0
    yold = y0
    ti = 0d0
    state = 1
    DO k = 2,n
       yold = y
       told = ti

       tout = ti+dt
       CALL dvode_f90(timestep2, neq, y, ti, tout, itask, istate, options,&
            g_fcn = bounceroots)
       !print *, ti
       IF (istate == 3) THEN
          IF (passing .OR. (yold(1)-th0)<0) THEN
             EXIT
          END IF
          
       END IF

       istate = 2
    END DO
    IF (istate /= 3) THEN
       PRINT *, "ERROR: findroot2 did not converge after 500 iterations"
       PRINT *, eta, etamin, etamax
    END IF
    
    findroot2(1)  = ti
    findroot2(2:) = y
  END FUNCTION findroot2
  
  SUBROUTINE bounceroots (NEQ, T, Y, NG, GOUT)
     INTEGER, INTENT(in) :: NEQ, NG
     REAL(8), INTENT(in) :: T, Y(neq)
     REAL(8), INTENT(out) :: GOUT(ng)
     GOUT(1) = Y(1) - th0
     GOUT(2) = 2d0*pi - (Y(1) - th0)
     RETURN
   END SUBROUTINE bounceroots

  SUBROUTINE bounce(taub_estimate)
    ! calculate all bounce averages
    REAL(8) :: findroot_res(nvar+1)
    REAL(8), OPTIONAL :: taub_estimate
    REAL(8) :: taub_est
    REAL(8) :: bmod, sqrtg, x(3), hder(3), hcovar(3), hctrvr(3), hcurl(3)
    REAL(8) :: y0(nvar)

    x(1) = s
    x(2) = 0d0
    x(3) = th0
    CALL do_magfie( x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl )
    y0 = 1d-15
    y0(1) = th0
    y0(2) = vpar(bmod)
    y0(3) = 0d0

    IF (PRESENT(taub_estimate)) THEN
       taub_est = taub_estimate
    ELSE
       taub_est = 2.0*pi/(vperp(bmod)*iota/R0*SQRT(eps/2d0))
    END IF
       
    findroot_res = findroot2(y0, taub_est/5d0, 1d-10)
    !print *, taub
    taub = findroot_res(1)
    bounceavg = findroot_res(2:)/taub
  END SUBROUTINE bounce
  
  SUBROUTINE Om_tB(OmtB, dOmtBdv, dOmtBdeta)
    ! returns bounce averaged toroidal magnetic drift frequency
    ! and derivatives w.r.t. v and eta
    REAL(8), INTENT(out) :: OmtB, dOmtBdv, dOmtBdeta
    REAL(8) :: splineval(3)
    REAL(8) :: Omth, dOmthdv, dOmthdeta
    IF (eta > etatp) THEN
       IF (eta > etatp*(1+epst_spl)) THEN
          splineval = spline_val_0(OmtB_spl_coeff, eta)
       ELSE ! extrapolation
          CALL Om_th(Omth, dOmthdv, dOmthdeta)
          splineval(1) = (k_OmtB_t*LOG(eta-etatp) + d_OmtB_t)*Omth/v
          splineval(2) = Omth/v*k_OmtB_t/(eta-etatp) +&
               dOmthdeta/v*(k_OmtB_t*LOG(eta-etatp) + d_OmtB_t)
       END IF
    ELSE
       IF (eta < etatp*(1-epsp_spl)) THEN
          splineval = spline_val_0(OmtB_pass_spl_coeff, eta)
       ELSE ! extrapolation
          CALL Om_th(Omth, dOmthdv, dOmthdeta)
          splineval(1) = (k_OmtB_p*LOG(etatp-eta) + d_OmtB_p)*Omth/v
          splineval(2) = Omth/v*k_OmtB_p/(eta-etatp) +&
               dOmthdeta/v*(k_OmtB_p*LOG(etatp-eta) + d_OmtB_p)
       END IF
    END IF
    OmtB = splineval(1)*v**2
    dOmtBdv = 2d0*splineval(1)*v
    dOmtBdeta = splineval(2)*v**2
  END SUBROUTINE Om_tB
    
  SUBROUTINE Om_ph(Omph, dOmphdv, dOmphdeta)
    ! returns canonical toroidal frequency
    ! and derivatives w.r.t. v and eta
    REAL(8), INTENT(out) :: Omph, dOmphdv, dOmphdeta
    REAL(8) :: Omth, dOmthdv, dOmthdeta
    REAL(8) :: OmtB, dOmtBdv, dOmtBdeta

    IF (eta > etatp) THEN
       Omph = Om_tE
       dOmphdv = 0d0
       dOmphdeta = 0d0
       IF (magdrift) THEN
          CALL Om_tB(OmtB, dOmtBdv, dOmtBdeta)
          Omph = Omph + OmtB
          dOmphdv = dOmphdv + dOmtBdv
          dOmphdeta = dOmphdeta + dOmtBdeta
       END IF
    ELSE
       CALL Om_th(Omth, dOmthdv, dOmthdeta)
       Omph = Om_tE + Omth/iota
       dOmphdv = dOmthdv/iota
       dOmphdeta = dOmthdeta/iota
       IF (magdrift) THEN
          CALL Om_tB(OmtB, dOmtBdv, dOmtBdeta)
          Omph = Omph + OmtB
          dOmphdv = dOmphdv + dOmtBdv
          dOmphdeta = dOmphdeta + dOmtBdeta
       END IF
    END IF
  END SUBROUTINE Om_ph

  SUBROUTINE Om_th(Omth, dOmthdv, dOmthdeta)
    ! returns canonical poloidal frequency
    ! and derivatives w.r.t. v and eta
    REAL(8), INTENT(out) :: Omth, dOmthdv, dOmthdeta
    REAL(8) :: splineval(3)

    IF (eta > etatp) THEN
       IF (eta > etatp*(1+epst_spl)) THEN
          splineval = spline_val_0(Omth_spl_coeff, eta)
       ELSE ! extrapolation
          splineval(1) = 2d0*pi/(k_taub_t*LOG(eta-etatp) + d_taub_t)
          splineval(2) = -splineval(1)**2/(2d0*pi) * k_taub_t/(eta-etatp)
       END IF
    ELSE
       IF (eta < etatp*(1-epsp_spl)) THEN
          splineval = spline_val_0(Omth_pass_spl_coeff, eta)
       ELSE ! extrapolation
          splineval(1) = 2d0*pi/(k_taub_p*LOG(etatp-eta) + d_taub_p)
          splineval(2) = -splineval(1)**2/(2d0*pi) * k_taub_p/(eta-etatp)
       END IF
    END IF
    Omth = splineval(1)*v
    dOmthdv = splineval(1)
    dOmthdeta = splineval(2)*v
  END SUBROUTINE Om_th

  SUBROUTINE driftorbit_coarse(eta_min, eta_max, roots, nroots)
    REAL(8), INTENT(in) :: eta_min, eta_max
    REAL(8), INTENT(out) :: roots(:,:)
    REAL(8) :: deta
    REAL(8) :: Omph, dOmphdv, dOmphdeta
    REAL(8) :: Omth, dOmthdv, dOmthdeta
    REAL(8) :: res, dresdv, dresdeta
    REAL(8) :: resold, dresdvold, dresdetaold
    INTEGER :: k, ninterv, nroots

    ninterv = SIZE(roots,1)
    
    deta = (eta_max-eta_min)*1d0/ninterv
    nroots = 0
    
    DO k = 0,ninterv
       eta = eta_min + k*deta
       CALL Om_th(Omth, dOmthdv, dOmthdeta)
       CALL Om_ph(Omph, dOmphdv, dOmphdeta)
       res = mth*Omth + mph*Omph
       dresdv = mth*dOmthdv + mph*dOmphdv
       dresdeta = mth*dOmthdeta + mph*dOmphdeta
       IF (k>0) THEN
          IF (SIGN(1d0,res) /= SIGN(1d0,resold)) THEN
             nroots = nroots+1
             roots(nroots, 1) = eta - deta
             roots(nroots, 2) = eta
          END IF
       END IF
       resold = res
       dresdvold = dresdv
       dresdetaold = dresdeta
    END DO
  END SUBROUTINE driftorbit_coarse
  
  FUNCTION driftorbit_nroot(nlev, eta_min, eta_max)
    INTEGER :: driftorbit_nroot, nlev
    REAL(8), OPTIONAL :: eta_min, eta_max
    REAL(8) :: etamin2, etamax2
    REAL(8) :: Omph_etamin, Omph_etamax, dummy, dummy2, &
         Omth_etamin, Omth_etamax, res_etamin, res_etamax
    
    IF (PRESENT(eta_min) .AND. PRESENT(eta_max)) THEN
       etamin2 = eta_min
       etamax2 = eta_max
    ELSE
       ! default behavior for trapped particles
       etamin2 = etatp*(1d0+epst)
       etamax2 = etadt*(1d0-epst)
    END IF

    driftorbit_nroot = 0
    ! TODO: return number of possible roots instead of 0 and 1
    eta = etamax2
    CALL Om_ph(Omph_etamin, dummy, dummy2)
    CALL Om_th(Omth_etamin, dummy, dummy2)
    eta = etamin2
    CALL Om_ph(Omph_etamax, dummy, dummy2)
    CALL Om_th(Omth_etamax, dummy, dummy2)

    res_etamin = mph*Omph_etamin+mth*Omth_etamin
    res_etamax = mph*Omph_etamax+mth*Omth_etamax
    IF(SIGN(1d0,res_etamin) /= SIGN(1d0,res_etamax)) THEN
       driftorbit_nroot = 1
    !else
       !print *, "driftorbit_nroot: ", v/vth,(etamin-etatp)/(etadt-etatp),&
       !     (etamax-etatp)/(etadt-etatp), res_etamin, res_etamax
    END IF
    IF(isnan(res_etamin) .OR. isnan(res_etamax)) THEN
       PRINT *, "ERROR: driftorbit_nroot found NaN value in Om_ph_ba"
       RETURN
    END IF
  END FUNCTION driftorbit_nroot

  FUNCTION driftorbit_root(tol, eta_min, eta_max)
    REAL(8) :: driftorbit_root(2)
    REAL(8) :: tol, res, res_old, eta0, eta_old
    REAL(8) :: Omph, dOmphdv, dOmphdeta
    REAL(8) :: Omth, dOmthdv, dOmthdeta
    INTEGER :: maxit, k, state
    REAL(8), OPTIONAL:: eta_min, eta_max
    REAL(8) :: etamin2, etamax2
    maxit = 100
    state = -2
    eta0 = eta
    eta_old = 0d0
    res = 0d0
    IF (PRESENT(eta_min) .AND. PRESENT(eta_max)) THEN
       etamin2 = eta_min
       etamax2 = eta_max
    ELSE
       ! default behavior for trapped particles
       etamin2 = etatp*(1d0+epst)
       etamax2 = etadt*(1d0-epst)
    END IF

    eta = etamin2
   
    IF(driftorbit_nroot(1, etamin2, etamax2) == 0) THEN
       PRINT *, "ERROR: driftorbit_root couldn't bracket 0 for v/vth = ", v/vth
       PRINT *, "ERROR: etamin = ", etamin2, " etamax = ", etamax2
       eta = etamin2
       CALL Om_ph(Omph, dOmphdv, dOmphdeta)
       CALL Om_th(Omth, dOmthdv, dOmthdeta)
       res = mph*Omph + mth*Omth
       PRINT *, res
       eta = etamax2
       CALL Om_ph(Omph, dOmphdv, dOmphdeta)
       CALL Om_th(Omth, dOmthdv, dOmthdeta)
       res = mph*Omph + mth*Omth
       PRINT *, res
       RETURN
    END IF
    DO k = 1,maxit
       res_old = res
       CALL Om_ph(Omph, dOmphdv, dOmphdeta)
       CALL Om_th(Omth, dOmthdv, dOmthdeta)
       res = mph*Omph + mth*Omth
       
       driftorbit_root(1) = eta

       IF (ABS(res) < tol) THEN
          state = 1
          driftorbit_root(2) = mph*dOmphdeta + mth*dOmthdeta
          !eta = eta*(1+1d-10)
          !call Om_ph(Omph, dOmphdv, dOmphdeta)
          !call Om_th(Omth, dOmthdv, dOmthdeta)
          !driftorbit_root(2) = (mph*Omph + mth*Omth - res)/(eta*1d-10)
          EXIT
       ! TODO: better condition based on slope of function
       ELSEIF ((mth >= 0 .AND. res > 0 .AND. eta > etatp)      .OR.&
               (mth >= -mph*q .AND. res < 0 .AND. eta < etatp) .OR.&
               (mth < 0 .AND. res < 0 .AND. eta > etatp)       .OR.&
               (mth < -mph*q .AND. res > 0 .AND. eta < etatp)) THEN
          etamax2 = eta
          eta_old = eta
          eta = (eta+etamin2)/2d0
       ELSE
          etamin2 = eta
          eta_old = eta
          eta = (eta+etamax2)/2d0
       END IF
    END DO
    IF (state < 0) THEN
       driftorbit_root(1) = 0
       driftorbit_root(2) = 0
       PRINT *, "ERROR: driftorbit_root did not converge in 100 iterations"
       PRINT *, "v/vth = ", v/vth, "mth = ", mth, "mph = ", mph
    END IF
    eta = eta0
  END FUNCTION driftorbit_root

  FUNCTION driftorbit_root2(tol, eta_min, eta_max)
    REAL(8) :: driftorbit_root2(2)
    REAL(8) :: tol, res, res_old, eta0, eta_old
    REAL(8) :: Omph, dOmphdv, dOmphdeta
    REAL(8) :: Omth, dOmthdv, dOmthdeta
    INTEGER :: maxit, k, state
    REAL(8), INTENT(in) :: eta_min, eta_max
    REAL(8) :: etamin2, etamax2
    LOGICAL :: slope_pos
    
    maxit = 100
    state = -2
    eta0 = eta
    eta_old = 0d0
    res = 0d0

    etamin2 = eta_min
    etamax2 = eta_max

    eta = etamin2
    CALL Om_ph(Omph, dOmphdv, dOmphdeta)
    CALL Om_th(Omth, dOmthdv, dOmthdeta)
    res = mph*Omph + mth*Omth
    
    eta = etamax2
    CALL Om_ph(Omph, dOmphdv, dOmphdeta)
    CALL Om_th(Omth, dOmthdv, dOmthdeta)
    IF(mph*Omph + mth*Omth - res > 0) THEN
       slope_pos = .TRUE.
    ELSE
       slope_pos = .FALSE.
    END IF
   
    IF(driftorbit_nroot(1, etamin2, etamax2) == 0) THEN
       PRINT *, "ERROR: driftorbit_root2 couldn't bracket 0 for v/vth = ", v/vth
       PRINT *, "ERROR: etamin = ", etamin2, " etamax = ", etamax2
       
       RETURN
    END IF
    
    DO k = 1,maxit
       res_old = res
       CALL Om_ph(Omph, dOmphdv, dOmphdeta)
       CALL Om_th(Omth, dOmthdv, dOmthdeta)
       res = mph*Omph + mth*Omth
       
       driftorbit_root2(1) = eta

       IF (ABS(res) < tol) THEN
          state = 1
          driftorbit_root2(2) = mph*dOmphdeta + mth*dOmthdeta
          EXIT
       ELSEIF ((slope_pos .AND. res > 0) .OR. &
               ((.NOT. slope_pos) .AND. res < 0)) THEN
          etamax2 = eta
          eta_old = eta
          eta = (eta+etamin2)/2d0
       ELSE
          etamin2 = eta
          eta_old = eta
          eta = (eta+etamax2)/2d0
       END IF
    END DO
    IF (state < 0) THEN
       driftorbit_root2(1) = 0
       driftorbit_root2(2) = 0
       PRINT *, "ERROR: driftorbit_root did not converge in 100 iterations"
       PRINT *, "v/vth = ", v/vth, "mth = ", mth, "mph = ", mph
    END IF
    eta = eta0
  END FUNCTION driftorbit_root2
  
  FUNCTION find_vmin(vmin0, vmax0)
    REAL(8) :: find_vmin, tol, vmin0, vmax0, vmin, vmax, v0
    INTEGER, PARAMETER :: nit = 100
    INTEGER :: k
    ! Bisection search for smallest possible v
    v0 = v
    tol = 1d-12*vth
    vmin = vmin0
    vmax = vmax0
    
    DO k=1,nit
       IF(ABS(vmax-vmin) < tol) THEN
          EXIT
       END IF
       v = (vmax + vmin)/2d0
       IF (driftorbit_nroot(1) /= 0) THEN
          vmax = v
       ELSE
          vmin = v
       END IF
    END DO
    v = v0
    find_vmin = vmax
  END FUNCTION find_vmin
  
  SUBROUTINE find_vlim(vmin, vmax)
    REAL(8) :: vmin, vmax, eta0
    REAL(8) :: Omth, dOmthdv, dOmthdeta

    IF (supban) THEN
       vmin = find_vmin(vmin, vmax)
       RETURN
    END IF

    ! if doing drift-orbit resonances,
    ! there is both a minimum and a maximum frequency
    eta0 = eta
    
    etamin = etatp*(1d0+epst)
    etamax = etadt*(1d0-epst)

    ! trapped orbits
    eta = etamax
    CALL Om_th(Omth, dOmthdv, dOmthdeta)
    vmin = MAX(vmin,-(mph*Om_tE)/(mth*Omth/v))

    eta = etamin
    CALL Om_th(Omth, dOmthdv, dOmthdeta)
    vmax = MIN(vmax,-(mph*Om_tE)/(mth*Omth/v))

    ! passing orbits
    etamin = etatp*epsp
    etamax = etatp*(1d0-epsp)
    
    eta = etamin
    CALL Om_th(Omth, dOmthdv, dOmthdeta)
    IF (-(mph*Om_tE)/((q*mph+mth)*Omth/v) > 0) THEN
       vmin = MIN(vmin,-(mph*Om_tE)/((q*mph+mth)*Omth/v))
    END IF
    
    eta = etamax
    CALL Om_th(Omth, dOmthdv, dOmthdeta)
    vmax = MAX(vmax,-(mph*Om_tE)/((q*mph+mth)*Omth/v))
        
    eta = eta0
  END SUBROUTINE find_vlim

  SUBROUTINE find_mthlim_t
  END SUBROUTINE find_mthlim_t

  SUBROUTINE find_mthlim_p
  END SUBROUTINE find_mthlim_p
  
!!$  subroutine find_vlim_p_magdrift(vmin, vmax)
!!$    real(8) :: vmin, vmax, eta0
!!$    real(8) :: Omth, dOmthdv, dOmthdeta
!!$    real(8) :: OmtB, dOmtBdv, dOmtBdeta
!!$    real(8) :: splineval(3)
!!$    real(8) :: k
!!$
!!$    do k = 1:
!!$       
!!$    
!!$    splineval = spline_val_0(OmtB_spl_coeff, eta)*v**2
!!$    OmtB = splineval(1)
!!$    dOmtBdv = 2*splineval(1)/v
!!$    dOmtBdeta = splineval(2)
!!$           
!!$  end subroutine find_vlim_p_magdrift

  
  SUBROUTINE find_vlim_p(vmin, vmax)
    REAL(8) :: vmin, vmax, eta0
    REAL(8) :: Omth, dOmthdv, dOmthdeta

    IF (supban) THEN
       vmin = find_vmin(vmin, vmax)
       RETURN
    END IF

    !if (magdrift) then
    !   find_vlim_p_magdrift(vmin, vmax)
    !   return
    !end if

    ! for drift-orbit resonances,
    ! there is both a minimum and a maximum frequency
    eta0 = eta

    ! passing orbits
    etamin = etatp*epsp
    etamax = etatp*(1d0-epsp)
    
    eta = etamin
    CALL Om_th(Omth, dOmthdv, dOmthdeta)
    vmin = MAX(vmin,-Om_tE/((mth*1d0/mph+q)*Omth/v))
    
    eta = etamax
    CALL Om_th(Omth, dOmthdv, dOmthdeta)
    vmax = MIN(vmax,-(mph*Om_tE)/((q*mph+mth)*Omth/v))
        
    eta = eta0
  END SUBROUTINE find_vlim_p
  
  
  SUBROUTINE find_vlim_t(vmin, vmax)
    REAL(8) :: vmin, vmax, eta0
    REAL(8) :: Omth, dOmthdv, dOmthdeta

    IF (mth == 0) THEN
       RETURN
    END IF
    
    IF (supban) THEN
       vmin = find_vmin(vmin, vmax)
       RETURN
    END IF

    ! if not supban (drift-orbit resonances),
    ! there is both a minimum and a maximum frequency
    eta0 = eta
    
    etamin = etatp*(1d0+epst)
    etamax = etadt*(1d0-epst)

    ! trapped orbits
    eta = etamax
    CALL Om_th(Omth, dOmthdv, dOmthdeta)
    vmin = MAX(vmin,-(mph*Om_tE)/(mth*Omth/v))

    eta = etamin
    CALL Om_th(Omth, dOmthdv, dOmthdeta)
    vmax = MIN(vmax,-(mph*Om_tE)/(mth*Omth/v))
        
    eta = eta0
  END SUBROUTINE find_vlim_t
  
  SUBROUTINE init_fsa
    ! Calculate the flux surface areas
    INTEGER, PARAMETER :: nth = 1000
    INTEGER :: k
    REAL(8) :: thrange(nth), dth
    REAL(8) :: bmod, sqrtg, x(3), hder(3), hcovar(3), hctrvr(3), hcurl(3)

    thrange = -pi + (/(k*2*pi/nth,k=1,nth)/)
    
    dth = thrange(2) - thrange(1) 
    x(1) = s
    x(2) = 0d0
    x(3) = 0d0
    
    dVds = 0d0
    B0  = 0d0
    PRINT *, "eps orig: ", eps
    eps = 0d0

    Bmin = -1d0
    Bmax = 0d0
    
    DO k = 1, nth
       x(3) = thrange(k)
       CALL do_magfie( x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl )
       dVds = dVds + sqrtg*dth
       B0   = B0  + bmod*dth
       eps  = eps - COS(x(3))*bmod*dth

       ! TODO: do fine search
       IF ((Bmin < 0) .OR. (bmod < Bmin)) THEN
          Bmin = bmod
          th0 = x(3)
       END IF
       IF (bmod > Bmax) Bmax = bmod
    END DO

    PRINT *, th0
    !th0 = 0d0 ! TODO remove this

    dVds = 2d0*pi*dVds
    B0   = B0/(2d0*pi)
    eps  = eps/(B0*pi)
    PRINT *, "eps calc:  ", eps
    PRINT *, "Bmin,Bmax: ", Bmin, Bmax

   ! call disp('init_fsa: iota       = ', iota)
   ! call disp('init_fsa: fsa/psi_pr = ', fsa/psi_pr)

  END SUBROUTINE init_fsa
  
  FUNCTION D11int(ux, etax)
    REAL(8) :: D11int
    REAL(8) :: ux, etax
    REAL(8) :: Omth, dummy, dummy2
    REAL(8) :: Hmn2
    
    v = ux*vth
    eta = etax
    CALL Om_th(Omth, dummy, dummy2)
    CALL bounce(2d0*pi/Omth)
    Hmn2 = (bounceavg(4)**2 + bounceavg(5)**2)*(mi*(ux*vth)**2/2d0)**2

    IF (calcflux) THEN
       D11int = pi**(3d0/2d0)*mph*c**2*q*vth/&
            (qi**2*dVds*psi_pr)*ux**3*EXP(-ux**2)*&
            taub*Hmn2*(mph-(mth+q*mph)*Bphcov*Omth*bounceavg(6))
    ELSE          
       D11int = pi**(3d0/2d0)*mph**2*c**2*q*vth/&
         (qi**2*dVds*psi_pr)*ux**3*EXP(-ux**2)*&
         taub*Hmn2
    END IF
  END FUNCTION D11int
  
  FUNCTION D12int(ux, etax)
    REAL(8) :: D12int
    REAL(8) :: ux, etax
    REAL(8) :: Omth, dummy, dummy2
    REAL(8) :: Hmn2
    
    v = ux*vth
    eta = etax
    CALL Om_th(Omth, dummy, dummy2)
    CALL bounce(2d0*pi/Omth)
    Hmn2 = (bounceavg(4)**2 + bounceavg(5)**2)*(mi*(ux*vth)**2/2d0)**2

    IF (calcflux) THEN
       D12int = pi**(3d0/2d0)*mph*c**2*q*vth/&
            (qi**2*dVds*psi_pr)*ux**3*EXP(-ux**2)*&
            taub*Hmn2*(mph-(mth+q*mph)*Bphcov*Omth*bounceavg(6))*ux**2
    ELSE          
       D12int = pi**(3d0/2d0)*mph**2*c**2*q*vth/&
         (qi**2*dVds*psi_pr)*ux**3*EXP(-ux**2)*&
         taub*Hmn2*ux**2
    END IF
  END FUNCTION D12int
  
  FUNCTION D11int_u(ux)
    REAL(8) :: ux
    REAL(8) :: D11int_u
    REAL(8) :: eta_res(2)
    REAL(8) :: roots(nlev, 3)
    INTEGER :: nroots, kr
    
    v = ux*vth
    D11int_u = 0d0

    CALL driftorbit_coarse(etamin, etamax, roots, nroots)
    IF(nroots == 0) RETURN
    
    !if (abs(M_t) > 1d-12) then
    !   eta_res = driftorbit_root(1d-8*abs(Om_tE), etamin, etamax)
    !else
    !   eta_res=driftorbit_root(1d-8*abs(c*mi*vth**2/(2*qi*psi_pr)),&
    !        etamin, etamax)
    !end if
    
    DO kr = 1,nroots
       eta_res = driftorbit_root2(1d-8*ABS(Om_tE), roots(kr,1), roots(kr,2))
       eta = eta_res(1)
       D11int_u = D11int_u + D11int(ux, eta_res(1))*1d0/ABS(eta_res(2))
    END DO
    
  END FUNCTION D11int_u
  
  FUNCTION D12int_u(ux)
    REAL(8) :: ux
    REAL(8) :: D12int_u
    REAL(8) :: eta_res(2)
    REAL(8) :: roots(nlev, 3)
    INTEGER :: nroots, kr
    v = ux*vth
    D12int_u = 0d0
    
    CALL driftorbit_coarse(etamin, etamax, roots, nroots)
    IF(nroots == 0) RETURN
    
    DO kr = 1,nroots       
       eta_res = driftorbit_root2(1d-8*ABS(Om_tE), roots(kr,1), roots(kr,2))
       eta = eta_res(1)
       D12int_u = D12int_u + D12int(ux, eta_res(1))*1d0/ABS(eta_res(2))
    END DO
  END FUNCTION D12int_u
  
  FUNCTION flux_integral(vmin, vmax)
    REAL(8) :: vmin, vmax
    REAL(8) :: flux_integral(2), err
    REAL(8) :: Dp, dsdreff
    INTEGER :: neval, ier
      
    flux_integral = 0d0
      
    CALL qag(D11int_u ,&
         (vmin+(vmax-vmin)*1d-10)/vth, (vmax-(vmax-vmin)*1d-10)/vth,&
         1d-9, 1d-3, 6, flux_integral(1), err, neval, ier)       

    CALL qag(D12int_u ,&
         (vmin+(vmax-vmin)*1d-10)/vth, (vmax-(vmax-vmin)*1d-10)/vth,&
         1d-9, 1d-3, 6, flux_integral(2), err, neval, ier)   
      
    Dp = pi*vth**3/(16d0*R0*iota*(qi*B0/(mi*c))**2)
    dsdreff = 2d0/a*SQRT(s)
    flux_integral = dsdreff**(-2)*flux_integral/Dp

  END FUNCTION flux_integral

!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
! WIP: integration with ODEs                                                   !
!------------------------------------------------------------------------------!
END MODULE driftorbit
MODULE lineint
  USE driftorbit

  CONTAINS
  
  SUBROUTINE intstep(neq, t, y, ydot)
    !
    !  Timestep function for orbit integration.
    !  Includes poloidal angle theta and parallel velocity.
    !  More integrands may be added starting from y(3)
    !
    INTEGER, INTENT (in) :: neq
    REAL(8), INTENT (in) :: t
    REAL(8), INTENT (in) :: y(neq)
    REAL(8), INTENT (out) :: ydot(neq)

    REAL(8) :: Omph, dOmphdv, dOmphdeta
    REAL(8) :: Omth, dOmthdv, dOmthdeta
    REAL(8) :: G, dGdv, dGdeta, absgradG
    REAL(8) :: deltatp
    REAL(8) :: Sveta

    ! check if trapped or passing region
    IF (etamin > etatp) THEN
       deltatp = 0d0
    ELSE
       deltatp = 1d0
    END IF

    Sveta = (vmax2-vmin2)*(etamax-etamin)

    v = vmin2 + y(1)*(vmax2-vmin2)
    eta = etamin + y(2)*(etamax-etamin)
    
    CALL Om_ph(Omph, dOmphdv, dOmphdeta)
    CALL Om_th(Omth, dOmthdv, dOmthdeta)
    G = (mph*q*deltatp + mth)*Omth + mph*Omph
    dGdv = ((mph*q*deltatp + mth)*dOmthdv + mph*dOmphdv)*(vmax2-vmin2)
    dGdeta = ((mph*q*deltatp + mth)*dOmthdeta + mph*dOmphdeta)*(etamax-etamin)
    absgradG = SQRT(dGdv**2 + dGdeta**2)
    
    ydot(1) = dGdeta           ! vbar
    ydot(2) = -dGdv            ! etabar
    CALL bounce(2d0*pi/Omth)
    ydot(3) = (vmax2-vmin2)*(etamax-etamin)*D11_ode()/vth            ! D11
    ydot(4) = (vmax2-vmin2)*(etamax-etamin)*D11_ode()*(v/vth)**2/vth ! D12
    
    ydot = ydot/absgradG

    ! always go in positive velocity direction
    IF (etamin > etatp) THEN
       ydot(1) = -ydot(1)
       ydot(2) = -ydot(2)
    END IF
  END SUBROUTINE intstep

  FUNCTION flux_integral_ode(vmin, vmax)
    USE dvode_f90_m2
    REAL(8) :: flux_integral_ode(2)
    REAL(8) :: vmin, vmax
    
    REAL(8) :: sk, ds
    REAL(8) :: Dp, dsdreff, eta_res(2)
    INTEGER :: ks, ns

    REAL(8) :: y(4), t, tout, rtol(4), atol(4), rwork(1024)
    INTEGER :: neq, itask, istate, iopt, itol, iwork(1024), ng, jt
    TYPE(VODE_OPTS) :: OPTIONS

    vmin2 = vmin + 1d-8*(vmax-vmin)
    vmax2 = vmax - 1d-8*(vmax-vmin)

    neq = 4
    ng = 4
    y(1) = 1d0-1d-15
    y(2) = 1d-15
    y(3) = 0d0
    y(4) = 0d0
    itol = 4
    t = 0.0d0
    rtol(1) = 1d-9
    rtol(2) = 1d-9
    rtol(3) = 1d-9
    rtol(4) = 1d-9
    atol(1) = 1d-9
    atol(2) = 1d-9
    atol(3) = 1d-9
    atol(4) = 1d-9
    itask = 1
    istate = 1
    iopt = 1
    jt = 2
    rwork = 0d0
    iwork = 0
    iwork(6) = 10000

    ns = 1000
    ds = 4d0/ns
    sk = 0d0

    dsdreff = 2d0/a*SQRT(s)
    Dp = pi*vth**3/(16d0*R0*iota*(qi*B0/(mi*c))**2)
    
    OPTIONS = SET_OPTS(DENSE_J=.TRUE.,RELERR_VECTOR=RTOL, &
         ABSERR_VECTOR=ATOL,NEVENTS=NG)

    v = vmax2
    IF (mth /= 0) THEN
       IF (etamin < etatp) THEN
          y(2) = 1d0 - 1d-15
       ELSE
          y(2) = 1d-15
       END IF
    ELSE
       eta_res = driftorbit_root(MAX(1d-9*ABS(Om_tE),1d-12), etamin, etamax)
       y(2) = (eta_res(1)-etamin)/(etamax-etamin)
    END IF    
    
    DO ks = 1, ns
       PRINT *, ks, y(2)
       tout = t+ds
       CALL dvode_f90(intstep,neq,y,t,tout,itask,istate,options,g_fcn=bbox)
       IF (istate == 3) THEN
          EXIT
       END IF
       PRINT *, ks
    END DO
    
    flux_integral_ode(1) = dsdreff**(-2)*y(3)/Dp 
    flux_integral_ode(2) = dsdreff**(-2)*y(4)/Dp
    PRINT *, flux_integral_ode
  END FUNCTION flux_integral_ode

  
  SUBROUTINE dummyjac
  END SUBROUTINE dummyjac

  SUBROUTINE bbox (NEQ, T, Y, NG, GOUT)
    INTEGER, INTENT(in) :: NEQ, NG
    REAL(8), INTENT(in) :: T, Y(neq)
    REAL(8), INTENT(out) :: GOUT(ng)
    !if (M_t > 0) then
    !   GOUT(1) = Y(1)
    !else
       GOUT(1) = 1d0
    !end if
    GOUT(2) = 1d0 - Y(1)
    GOUT(3) = Y(2)
    GOUT(4) = 1d0 - Y(2)
    RETURN
  END SUBROUTINE bbox
  
  FUNCTION D11_ode()
    REAL(8) :: D11_ode
    REAL(8) :: ux
    REAL(8) :: Hmn2
    
    ux = v/vth
    
    Hmn2 = (bounceavg(4)**2 + bounceavg(5)**2)*(mi*(ux*vth)**2/2d0)**2
    
    D11_ode = pi**(3d0/2d0)*mph**2*c**2*q*vth/&
         (qi**2*dVds*psi_pr)*ux**3*EXP(-ux**2)*&
         taub*Hmn2
  END FUNCTION D11_ode

END MODULE lineint
