module orbit
  
  
contains

  function vecprod(u,v)
    real(8) :: vecprod(3)
    real(8), intent(in) :: u(3), v(3)

    vecprod(1) = u(2)*v(3) - u(3)*v(2)
    vecprod(2) = u(3)*v(1) - u(1)*v(3)
    vecprod(3) = u(1)*v(2) - u(2)*v(1)
  end function vecprod
  
  subroutine step_zeroorder(neq, t, y, ydot)
    integer, intent (in) :: neq
    real(8), intent (in) :: t
    real(8), intent (in) :: y(neq)
    real(8), intent (out) :: ydot(neq)
    
    real(8) :: bmod, sqrtg, hder(3), hcovar(3), hctrvr(3), hcurl(3)

    call do_magfie(y(1:3), bmod, sqrtg, bder, hcovar, hctrvr, hcurl)

    shearterm = Bphcov*dqds
    if (noshear) then
       shearterm = 0
    end if

    Om_tB_v = mi*c*q/(2d0*qi*psi_pr*bmod)*(&      ! Om_tB/v**2
         -(2d0-eta*bmod)*bmod*hder(1)&
         +2d0*(1d0-eta*bmod)*hctrvr(3)*&
         (dBthcovds+q*dBphcovds+shearterm))

    ydot(1) = 0.0                                               ! r
    ydot(2) = Om_tB_v*v**2 + Om_tE                              ! phi
    ydot(3) = y(4)*hctrvr(3)                                    ! theta
    ydot(4) = -v**2*eta/2d0*hctrvr(3)*hder(3)*bmod              ! v_par     
    
  end subroutine step_zeroorder
  
  
  subroutine step_full(neq, t, y, ydot)
    ! TODO: add gradPhi for fast rotation
    
    integer, intent (in) :: neq
    real(8), intent (in) :: t
    real(8), intent (in) :: y(neq)
    real(8), intent (out) :: ydot(neq)

    real(8) :: bmod, sqrtg, hder(3), hcovar(3), hctrvr(3), hcurl(3)
    
    real(8) :: hstar(3)       ! hstar = h + vpar/(omc)*curl(h)
    real(8) :: hcrossgradH(3) ! h x gradH = h x (mu*gradB + e*gradPhi) (contravariant comp)
    real(8) :: hstarpar       ! hstar_par = h*hstar
    real(8) :: Er             ! Er = dphi/dr = psi_pr/(q*c)*Om_tE, no fast rotation yet

    call do_magfie(y(1:3), bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
    
    hstar = hctrvr + mi*c/(qi*bmod)*y(4)*hcurl
    hstarpar = 1.0 + mi*c/(qi*bmod)*y(4)*sum(hcovar*hcurl)

    Er = psi_pr/(q*c)*Om_tE
    hcrossgradH = sqrtg*vecprod(hcovar, mu*bmod*hder + (/qi*Er,0,0/)) 
    
    ydot(1:3) = 1.0d0/hstarpar*(y(4)*hstar + c/(qi*bmod**2)*hcrossgradH)
    ydot(4) = -1.0d0/(mi*hstarpar)*mu*bmod*sum(hstar*hder)
    
  end subroutine step_full
end module orbit
