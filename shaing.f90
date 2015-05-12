

! TODO: adapt this to neo2 magfie

  function kappa2()
    real(8) :: kappa2
    kappa2 = (1d0 - eta*B0mn(fsind,1)*(1-eps(fsind)))/(2d0*eta*(B0mn(fsind,1)*eps(fsind)))
    !kappa2 = sqrt((mi*v**2/2d0-J_perp()*om_c(0d0))/(2*J_perp()*&
    !     om_c(pi/2)*eps(fsind)))
    !print *, mi*v**2/2d0, E() - qi*Phi_E(fsind)
  end function kappa2
  
  ! Shaing2009-035009 - (8) for bounce averaged toroidal drift 
  function Omph_shaing()
    real(8) :: Omph_shaing, kappa, Eell, Kell         
    
    ! Elliptic integrals for Shaing formula
    !kappa = sqrt((E()-qi*Phi_E(fsind)-J_perp()*om_c(0.0d0))/(2*J_perp()*&
    !     om_c(pi/2)*eps(fsind)))
    kappa = sqrt((mi*v**2/2d0-J_perp()*om_c(0d0))/(2*J_perp()*&
         om_c(pi/2)*eps))
    !kappa = sqrt(0.83)               ! test case so that 2*E/K = 1
    Kell = ellf(pi/2d0, kappa)
    Eell = elle(pi/2d0, kappa)

    ! Shaing2009-035009 - (8) for bounce averaged toroidal drift 
    Omph_shaing = Om_tE +&
         J_perp()*B0mn(fsind,1)/mi*depsdr(fsind)*(2*Eell/Kell-1)
  end function Omph_shaing
