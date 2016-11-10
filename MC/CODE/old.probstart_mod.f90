MODULE probstart_mod
  USE neo_precision
  USE polylag_3,  ONLY : mp, indef, plag1d
  IMPLICIT NONE
  !
  PRIVATE c, e_charge, e_mass, p_mass, ev
  REAL(kind=dp),PARAMETER  :: c=2.9979d10
  REAL(kind=dp),PARAMETER  :: e_charge=4.8032d-10
  REAL(kind=dp),PARAMETER  :: e_mass=9.1094d-28
  REAL(kind=dp),PARAMETER  :: p_mass=1.6726d-24
  REAL(kind=dp),PARAMETER  :: ev=1.6022d-12
  !
CONTAINS
  INCLUDE 'cx_ion.f90'
  !
  SUBROUTINE interp_prof_vec(es_prof_vec, prof_vec, es_lvec, interp_prof)
    REAL(kind=dp), DIMENSION(:), INTENT(in)  :: es_prof_vec, prof_vec, es_lvec
    REAL(kind=dp), DIMENSION(:), INTENT(out) :: interp_prof
    ! local definitions
    INTEGER :: ind, num_prof, num_l
    INTEGER, DIMENSION(mp) :: indu
    REAL(kind=dp), DIMENSION(mp) :: xp,fp
    REAL(kind=dp) :: s, der, dxm1
    !
    num_l = SIZE(es_lvec,1)
    num_prof = SIZE(es_prof_vec,1)
    dxm1 = 1.0_dp/(es_prof_vec(2)-es_prof_vec(1))    
    DO ind=1,num_l
       ! prepare interpolation
       s = es_lvec(ind)
       CALL indef(s,es_prof_vec(1),dxm1,num_prof,indu)
       ! interpolate
       xp=es_prof_vec(indu)
       fp=prof_vec(indu)
       CALL plag1d(s,fp,dxm1,xp,interp_prof(ind),der)
    END DO
  END SUBROUTINE interp_prof_vec
  !
  SUBROUTINE calc_probstart(lvec, es_lvec, plasma_mat, E_beam, amb, probstart_vec)
    REAL(kind=dp), DIMENSION(:), INTENT(in)   :: lvec, es_lvec
    REAL(kind=dp), DIMENSION(:,:), INTENT(in) :: plasma_mat
    REAL(kind=dp), INTENT(in)                 :: E_beam, amb
    REAL(kind=dp), DIMENSION(:), INTENT(out)  :: probstart_vec
    ! local definitions
    INTEGER :: ind, num_prof, num_l
    REAL(kind=dp) :: v0, delta_l, nu_ov_v0, sigma_cx
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: ni_prof, Te_prof, Ti_prof
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: nui_ov_v0, nucx_ov_v0, tau_vec
    !
    num_l = SIZE(es_lvec,1)
    ALLOCATE(ni_prof(num_l),Te_prof(num_l),Ti_prof(num_l))
    CALL interp_prof_vec(plasma_mat(:,1), plasma_mat(:,2), es_lvec, ni_prof)
    CALL interp_prof_vec(plasma_mat(:,1), plasma_mat(:,4), es_lvec, Ti_prof)
    CALL interp_prof_vec(plasma_mat(:,1), plasma_mat(:,6), es_lvec, Te_prof)
!!$    DO ind = 1,num_l
!!$       PRINT *, es_lvec(ind), ni_prof(ind), Ti_prof(ind), Te_prof(ind)
!!$    END DO
    !
    ALLOCATE(nui_ov_v0(num_l),nucx_ov_v0(num_l))
    v0=sqrt(2.d0*E_beam*ev/(amb*p_mass))
    DO ind = 1,num_l
       nui_ov_v0(ind) = sigmav_ion(Te_prof(ind))*ni_prof(ind)/v0
       sigma_cx =  0.6937e-14_dp*((1.0_dp-0.155*LOG10(E_beam))**2)/&
            (1.0_dp+0.1112e-14_dp*(E_beam**3.3_dp))
       nucx_ov_v0(ind) = ni_prof(ind)*sigma_cx
    END DO
    !
    ALLOCATE(tau_vec(num_l))
    tau_vec(1) = 0.0_dp
    delta_l = lvec(2)-lvec(1)
    DO ind = 2,num_l
       nu_ov_v0 = nui_ov_v0(ind) + nucx_ov_v0(ind)
       tau_vec(ind) = tau_vec(ind-1)+delta_l*nu_ov_v0
    END DO
    !
    probstart_vec = 1.0_dp - EXP(-tau_vec)
    !
    DEALLOCATE(ni_prof,Te_prof,Ti_prof)
    DEALLOCATE(nui_ov_v0,nucx_ov_v0)
    DEALLOCATE(tau_vec)
  END SUBROUTINE calc_probstart
END MODULE probstart_mod
