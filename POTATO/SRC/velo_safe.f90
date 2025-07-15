!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
! Safe version of velo routine that prevents floating point exceptions
!
subroutine velo_safe(tau,z,vz)
    use parmot_mod, only : rmu,ro0,gradpsiast,dpsiast_dR,dpsiast_dZ
    use field_eq_mod, only : ierrfield
    
    implicit none
    
    integer :: i
    
    double precision, intent(in)  :: tau, z(5)
    double precision, intent(out) :: vz(5)
    double precision :: x(3),bmod,sqrtg,bder(3),hcovar(3),hctrvr(3),hcurl(3)
    double precision :: derphi(3),phi_elec
    double precision :: p,alambd,p2,ovmu,gamma2,gamma,ppar,vpa,coala
    double precision :: rmumag,rovsqg,rosqgb,rovbm
    double precision :: a_phi(3),a_b(3),a_c(3),hstar(3)
    double precision :: s_hc,hpstar,phidot,blodot,bra
    
    ! Safety parameters to prevent division by zero
    double precision, parameter :: SMALL_P = 1.0d-10
    double precision, parameter :: SMALL_BMOD = 1.0d-10
    double precision, parameter :: SMALL_SQRTG = 1.0d-10
    double precision, parameter :: SMALL_GAMMA = 1.0d-10
    double precision, parameter :: SMALL_HPSTAR = 1.0d-10
    double precision, parameter :: SMALL_VPA = 1.0d-10
    
    x(1:3) = z(1:3)
    
    call magfie(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
    
    if(ierrfield.ne.0) then
        ! Return zero velocities on field error
        vz = 0.0d0
        return
    endif
    
    ! Ensure minimum values to prevent division by zero
    bmod = max(bmod, SMALL_BMOD)
    sqrtg = max(sqrtg, SMALL_SQRTG)
    
    call elefie(x,phi_elec,derphi)
    
    p = max(z(4), SMALL_P)  ! Ensure minimum momentum
    alambd = z(5)
    
    ! Ensure alambd is in valid range [-1, 1]
    alambd = max(-1.0d0, min(1.0d0, alambd))
    
    p2 = p*p
    ovmu = 2.0d0/rmu
    gamma2 = p2*ovmu + 1.0d0
    gamma = sqrt(gamma2)
    gamma = max(gamma, SMALL_GAMMA)  ! Ensure minimum gamma
    
    ppar = p*alambd
    vpa = ppar/gamma
    
    ! Prevent vpa from being exactly zero for gradpsiast calculation
    if (abs(vpa) < SMALL_VPA .and. gradpsiast) then
        vpa = sign(SMALL_VPA, vpa)
        if (vpa == 0.0d0) vpa = SMALL_VPA
    end if
    
    coala = (1.0d0 - alambd**2)
    rmumag = 0.5d0*p2*coala/bmod
    
    rovsqg = ro0/sqrtg
    rosqgb = 0.5d0*rovsqg/bmod
    rovbm = ro0/bmod
    
    a_phi(1) = (hcovar(2)*derphi(3) - hcovar(3)*derphi(2))*rosqgb
    a_b(1) = (hcovar(2)*bder(3) - hcovar(3)*bder(2))*rovsqg
    a_phi(2) = (hcovar(3)*derphi(1) - hcovar(1)*derphi(3))*rosqgb
    a_b(2) = (hcovar(3)*bder(1) - hcovar(1)*bder(3))*rovsqg
    a_phi(3) = (hcovar(1)*derphi(2) - hcovar(2)*derphi(1))*rosqgb
    a_b(3) = (hcovar(1)*bder(2) - hcovar(2)*bder(1))*rovsqg
    
    s_hc = 0.0d0
    do i = 1, 3
        a_c(i) = hcurl(i)*rovbm
        s_hc = s_hc + a_c(i)*hcovar(i)
        hstar(i) = hctrvr(i) + ppar*a_c(i)
    enddo
    hpstar = 1.0d0 + ppar*s_hc
    
    ! Ensure minimum hpstar to prevent division by zero
    if (abs(hpstar) < SMALL_HPSTAR) then
        hpstar = sign(SMALL_HPSTAR, hpstar)
    end if
    
    ! Velocities in coordinate space
    phidot = 0.0d0
    blodot = 0.0d0
    do i = 1, 3
        bra = vpa*hstar(i) + a_phi(i) + a_b(i)*rmumag/gamma
        vz(i) = bra/hpstar
        phidot = phidot + vz(i)*derphi(i)
        blodot = blodot + vz(i)*bder(i)
    enddo
    
    ! Velocities in phase space
    vz(4) = -0.5d0*gamma*phidot/p
    
    ! Safe computation of vz(5) avoiding division issues
    if (abs(hpstar) > SMALL_HPSTAR .and. p > SMALL_P .and. gamma > SMALL_GAMMA) then
        vz(5) = -(0.5d0*coala/hpstar)*(sum(hstar*derphi)/p &
                + p*sum(hstar*bder)/gamma + alambd*sum(a_phi*bder))
    else
        vz(5) = 0.0d0
    end if
    
    if(gradpsiast) then
        if (abs(vpa) > SMALL_VPA) then
            dpsiast_dR = bmod*hpstar*x(1)/vpa
            dpsiast_dZ = dpsiast_dR
            dpsiast_dR = dpsiast_dR*vz(3)
            dpsiast_dZ = -dpsiast_dZ*vz(1)
        else
            dpsiast_dR = 0.0d0
            dpsiast_dZ = 0.0d0
        end if
    endif
    
    return
end subroutine velo_safe
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc