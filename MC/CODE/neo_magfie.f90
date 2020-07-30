MODULE neo_magfie_mod

  USE neo_precision
  USE neo_input,                                                       &
       ONLY: es, ixm, ixn, mnmax, psi_pr, pixm, pixn, nfp
  USE neo_control,                                                     &
       ONLY: fluxs_interp, write_progress, phi_n, theta_n, lab_swi,    &
       inp_swi
  USE neo_sub_mod,                                                     &
       ONLY: neo_read_control, neo_init, neo_init_spline
  USE neo_spline_data,                                                 &
       ONLY: r_mhalf,                                                  &
       a_bmnc, b_bmnc, c_bmnc, d_bmnc,                                 &
       a_bmns, b_bmns, c_bmns, d_bmns,                                 &
       a_rmnc, b_rmnc, c_rmnc, d_rmnc,                                 &
       a_rmns, b_rmns, c_rmns, d_rmns,                                 &
       a_zmnc, b_zmnc, c_zmnc, d_zmnc,                                 &
       a_zmns, b_zmns, c_zmns, d_zmns,                                 &
       a_lmnc, b_lmnc, c_lmnc, d_lmnc,                                 &
       a_lmns, b_lmns, c_lmns, d_lmns,                                 &
       a_iota, b_iota, c_iota, d_iota,                                 &
       a_curr_tor, b_curr_tor, c_curr_tor, d_curr_tor,                 &
       a_curr_pol, b_curr_pol, c_curr_pol, d_curr_pol,                 &
       a_pprime,   b_pprime,   c_pprime,   d_pprime,                   &
       a_sqrtg00,  b_sqrtg00,  c_sqrtg00,  d_sqrtg00
  USE inter_interfaces,                                                &
       ONLY: splint_horner3,                                           &
       tf, tfp, tfpp, tfppp,                                           &
       tfone, tfzero
  USE neo_work,                                                        &
       ONLY: cosmth, cosnph, sinmth, sinnph, theta_int, phi_int,       &
       theta_start, theta_end, phi_start, phi_end
  USE neo_actual_fluxs, ONLY : s_sqrtg00
  USE spline_mod, ONLY: spl2d, poi2d, eva2d, eva2d_fd
  !! Modifications by Andreas F. Martitsch (12.03.2014)
  ! Use this quantity for normalization. Note:
  ! Variable is computed in mag_interface.f90 ("boozer_bmod0").
  ! It is available for the first time after 1st call
  ! of "make_magnetics". Therefore, within the first two calls
  ! of "neo_magfie" this variable is zero, but these calls are
  ! not used for the computation of physical quantities.
  USE partpa_mod,  ONLY : bmod0
  !! End Modifications by Andreas F. Martitsch (12.03.2014)
  !! Modifications by Andreas F. Martitsch (26.01.2016)
  !--> Normalization Bref for stand-alone version
  USE neo_exchange, ONLY : bmref
  !! Modifications by Andreas F. Martitsch (26.01.2016)
  USE polylag_3, ONLY : indef, plag1d , mp_lag => mp
  !! Modifications by Andreas F. Martitsch (27.01.2016)
  ! -> switch on/off use of splined Fourier coefficients within neo_magfie
  USE spline_settings, ONLY : isw_spl_fourier_cof, isw_eval_spl2d_der, &
       isw_eval_bcovars, flux_surf_dist
  !! End Modifications by Andreas F. Martitsch (27.01.2016)

  !---------------------------------------------------------------------------
  !USE var_sub_misc, ONLY: fac_c,iota_m ! fac_m
  !---------------------------------------------------------------------------

  IMPLICIT NONE
  REAL(dp), DIMENSION(:), ALLOCATABLE              :: magfie_sarray
  INTEGER                                          :: magfie_spline    = 0
  INTEGER                                          :: magfie_newspline = 1
  INTEGER                                          :: magfie_newspline_cyl = 1
  INTEGER                                          :: magfie_newspline_bN  = 1
  ! switch for different output:
  !  0:  output for SMT and BMC
  !  1:  output for NEO
  INTEGER                                          :: magfie_result    = 0
  INTEGER                                          :: magfie_sarray_len

  REAL(dp), DIMENSION(:), ALLOCATABLE              :: curr_tor_array
  REAL(dp), DIMENSION(:), ALLOCATABLE              :: curr_tor_s_array
  REAL(dp), DIMENSION(:), ALLOCATABLE              :: curr_pol_array
  REAL(dp), DIMENSION(:), ALLOCATABLE              :: curr_pol_s_array
  REAL(dp), DIMENSION(:), ALLOCATABLE              :: iota_array
  REAL(dp), DIMENSION(:), ALLOCATABLE              :: pprime_array
  REAL(dp), DIMENSION(:), ALLOCATABLE              :: sqrtg00_array

  REAL(dp)                                         :: s_pprime

  REAL(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE, TARGET      :: bmod_spl
  REAL(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE, TARGET      :: bb_s_spl
  REAL(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE, TARGET      :: bb_tb_spl
  REAL(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE, TARGET      :: bb_pb_spl
  REAL(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE, TARGET      :: gval_spl
  !! Modifications by Andreas F. Martitsch (11.03.2014)
  ! Storage arrays for the 2d splines (over the flux-surface) of the additionally
  ! needed metric tensor elements (used to compute the B-field components,
  ! which are necessary for modeling the magnetic rotation).
  ! Once computed these arrays can be used to reconstruct the desired
  ! quantities at intermediate (theta,phi)-values.
  REAL(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE, TARGET      :: gstb_spl
  REAL(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE, TARGET      :: gspb_spl
  REAL(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE, TARGET      :: gstb_tb_spl
  REAL(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE, TARGET      :: gspb_tb_spl
  REAL(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE, TARGET      :: gstb_pb_spl
  REAL(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE, TARGET      :: gspb_pb_spl
  !! End Modifications by Andreas F. Martitsch (11.03.2014)
  !! Modifications by Andreas F. Martitsch (13.11.2014)
  ! Storage arrays for the 2d splines (over the flux-surface) of the additionally
  ! needed quantities for NTV output
  REAL(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE, TARGET      :: R_spl
  REAL(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE, TARGET      :: Z_spl
  REAL(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE, TARGET      :: L_spl
  !! End Modifications by Andreas F. Martitsch (13.11.2014)
  REAL(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE, TARGET      :: R_axi_spl
  REAL(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE, TARGET      :: Z_axi_spl
  REAL(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE, TARGET      :: L_axi_spl

  REAL(dp) :: boozer_iota
  REAL(dp) :: boozer_sqrtg00
  !! Modifications by Andreas F. Martitsch (12.03.2014)
  ! boozer_curr_tor, boozer_curr_pol, boozer_psi_pr,
  ! boozer_sqrtg11 and boozer_isqrg are now converted
  ! to cgs-units.
  ! This step requires changes within rhs_kin.f90 and
  ! ripple_solver.f90!
  REAL(dp) :: boozer_curr_pol_hat
  REAL(dp) :: boozer_curr_tor_hat
  ! Radial derivatives of toroidal / poloidal currents
  ! (In fact these currents are already the respective
  ! covariant B-field components; conversion done within
  ! neo_read)
  REAL(dp) :: boozer_curr_pol_hat_s
  REAL(dp) :: boozer_curr_tor_hat_s
  REAL(dp) :: boozer_psi_pr_hat
  REAL(dp) :: boozer_sqrtg11 ! Test
  REAL(dp) :: boozer_isqrg
  !! End Modifications by Andreas F. Martitsch (12.03.2014)

  REAL(dp), PRIVATE :: av_b2_m ! Klaus

  !! Modifications by Andreas F. Martitsch (11.03.2014)
  ! Transfer the computed values of the additionally needed
  ! B-field components between neo_magfie_a and neo_magfie_b
  ! (best solution at the moment, since keyword optional
  ! does not seem to work for a module procedure)
  REAL(dp), PRIVATE :: dbcovar_s_dtheta
  REAL(dp), PRIVATE :: dbcovar_s_dphi
  !! End Modifications by Andreas F. Martitsch (11.03.2014)

  !! Modifications by Andreas F. Martitsch (13.11.2014)
  ! Local variables for the additionally needed quantities for NTV output
  REAL(dp), PRIVATE :: r_val, z_val, p_val
  REAL(dp), PRIVATE :: r_s_val, z_s_val, p_s_val
  REAL(dp), PRIVATE :: r_tb_val, z_tb_val, p_tb_val
  REAL(dp), PRIVATE :: r_pb_val, z_pb_val, p_pb_val
  !! End Modifications by Andreas F. Martitsch (13.11.2014)

  INTERFACE neo_magfie
     MODULE PROCEDURE neo_magfie_a
  END INTERFACE neo_magfie

!!$  INTERFACE neo_magfie
!!$     MODULE PROCEDURE neo_magfie_a, neo_magfie_b, neo_magfie_c
!!$  END INTERFACE

  INTERFACE cyl_coord
     MODULE PROCEDURE cyl_coord_a
  END INTERFACE cyl_coord

  INTERFACE plot_Bfield
     MODULE PROCEDURE plot_Bfield_a
  END INTERFACE plot_Bfield

CONTAINS

  SUBROUTINE neo_magfie_a( x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl )
    ! input / output
    REAL(dp), DIMENSION(:),       INTENT(in)            :: x
    REAL(dp),                     INTENT(out)           :: bmod
    REAL(dp),                     INTENT(out)           :: sqrtg
    REAL(dp), DIMENSION(SIZE(x)), INTENT(out)           :: bder
    REAL(dp), DIMENSION(SIZE(x)), INTENT(out)           :: hcovar
    REAL(dp), DIMENSION(SIZE(x)), INTENT(out)           :: hctrvr
    REAL(dp), DIMENSION(SIZE(x)), INTENT(out)           :: hcurl
    ! local definitions
    !! Modifications by Andreas F. Martitsch (11.03.2014)
    ! Locally computed values of the additionally needed
    ! metric tensor elements
    REAL(dp)                                         :: gstb, gspb
    REAL(dp)                                         :: gstb_tb, gspb_tb
    REAL(dp)                                         :: gstb_pb, gspb_pb
    ! Locally computed value of the additionally needed
    ! B-field component
    REAL(dp)                                         :: bcovar_s
    !! End Modifications by Andreas F. Martitsch (11.03.2014)
    INTEGER(i4b)                                     :: swd = 1
    INTEGER                                          :: i, m, n
    INTEGER                                          :: npsi
    REAL(dp)                                         :: m0  = 0.0_dp
    REAL(dp)                                         :: yp, ypp, yppp

    REAL(dp)                                         :: bmnc, bmnc_s
    REAL(dp)                                         :: sinv, cosv
    REAL(dp)                                         :: iota
    REAL(dp)                                         :: curr_tor, curr_tor_s
    REAL(dp)                                         :: curr_pol, curr_pol_s
    REAL(dp)                                         :: bb_s, bb_tb, bb_pb
    REAL(dp)                                         :: fac, fac1

    INTEGER                                          :: k_es
    INTEGER                                          :: s_detected
    INTEGER                                          :: imn
    INTEGER                                          :: it, ip, im, in
    INTEGER                                          :: mt = 1
    INTEGER                                          :: mp = 1
    INTEGER                                          :: theta_ind, phi_ind
    INTEGER                                          :: ierr
    REAL(dp)                                         :: s
    REAL(dp)                                         :: magfie_epsi = 1.e-9
    REAL(dp)                                         :: bi, bi_s, ri, zi, li
    !! Modifications by Andreas F. Martitsch (07.03.2014)
    ! Auxiliary variables for the Fourier summation
    REAL(dp)                                         :: ri_s, zi_s, li_s
    REAL(dp)                                         :: bis, bis_s, ris, zis, lis
    REAL(dp)                                         :: ris_s, zis_s, lis_s
    !! End Modifications by Andreas F. Martitsch (07.03.2014)
    REAL(dp)                                         :: theta_d, phi_d

    REAL(dp), DIMENSION(:), ALLOCATABLE              :: s_bmnc, s_bmnc_s
    REAL(dp), DIMENSION(:), ALLOCATABLE              :: s_rmnc, s_zmnc, s_lmnc
    !! Modifications by Andreas F. Martitsch (06.03.2014)
    ! Radial derivatives of (R,Z,phi)-components obtained from the 1d spline
    REAL(dp), DIMENSION(:), ALLOCATABLE              :: s_rmnc_s, s_zmnc_s, s_lmnc_s
    !! End Modifications by Andreas F. Martitsch (06.03.2014)
    !! Modifications by Andreas F. Martitsch (06.08.2014)
    ! Additional data from Boozer files without Stellarator symmetry
    REAL(kind=dp),    DIMENSION(:),     ALLOCATABLE :: s_bmns, s_bmns_s
    REAL(kind=dp),    DIMENSION(:),     ALLOCATABLE :: s_rmns, s_zmns, s_lmns
    REAL(kind=dp),    DIMENSION(:),     ALLOCATABLE :: s_rmns_s, s_zmns_s, s_lmns_s
    !! End Modifications by Andreas F. Martitsch (06.08.2014)
    REAL(dp), DIMENSION(:,:), ALLOCATABLE            :: bmod_a
    REAL(dp), DIMENSION(:,:), ALLOCATABLE            :: bb_s_a
    REAL(dp), DIMENSION(:,:), ALLOCATABLE            :: bb_tb_a
    REAL(dp), DIMENSION(:,:), ALLOCATABLE            :: bb_pb_a
    REAL(dp), DIMENSION(:,:), ALLOCATABLE            :: r,z,l
    REAL(dp), DIMENSION(:,:), ALLOCATABLE            :: r_tb,z_tb,p_tb
    REAL(dp), DIMENSION(:,:), ALLOCATABLE            :: r_pb,z_pb,p_pb
    REAL(dp), DIMENSION(:,:), ALLOCATABLE            :: gtbtb,gpbpb,gtbpb
    !! Modifications by Andreas F. Martitsch (11.03.2014)
    ! Temporary storage arrays for the Fourier summations related to
    ! the radial derivatives of (R,Z,phi)-components and the
    ! additionally needed metric tensor elements
    REAL(dp), DIMENSION(:,:), ALLOCATABLE            :: r_s,z_s,p_s
    REAL(dp), DIMENSION(:,:), ALLOCATABLE            :: r_stb,z_stb,p_stb
    REAL(dp), DIMENSION(:,:), ALLOCATABLE            :: r_tbtb,z_tbtb,p_tbtb
    REAL(dp), DIMENSION(:,:), ALLOCATABLE            :: r_pbtb,z_pbtb,p_pbtb
    REAL(dp), DIMENSION(:,:), ALLOCATABLE            :: r_spb,z_spb,p_spb
    REAL(dp), DIMENSION(:,:), ALLOCATABLE            :: r_pbpb,z_pbpb,p_pbpb
    REAL(dp), DIMENSION(:,:), ALLOCATABLE            :: gstb_a,gspb_a
    REAL(dp), DIMENSION(:,:), ALLOCATABLE            :: gstb_tb_a,gspb_tb_a
    REAL(dp), DIMENSION(:,:), ALLOCATABLE            :: gstb_pb_a,gspb_pb_a
    !! End Modifications by Andreas F. Martitsch (11.03.2014)
    REAL(dp), DIMENSION(:,:), ALLOCATABLE            :: sqrg11_met

    REAL(dp), DIMENSION(:,:,:,:), POINTER            :: p_spl

    REAL(dp) :: isqrg, sqrg11

    !REAL(dp) :: s_sqrtg00_m Klaus

    ! Temporary variables for Lagrange interpolation
    INTEGER :: ind1
    INTEGER, DIMENSION(mp_lag) :: indu
    REAL(dp), DIMENSION(mp_lag) :: xp, fp, fp1
    REAL(dp) :: fun, der, dxm1
    REAL(dp), DIMENSION(2) :: bder_thph

    !*******************************************************************
    ! Initialisation if necessary
    !*******************************************************************
    IF ( .NOT. ALLOCATED(es) ) THEN
       CALL neo_read_control()
       fluxs_interp = 1
       CALL neo_init(npsi)
       PRINT *, 'theta_start,theta_end,phi_start,phi_end'
       PRINT *, theta_start,theta_end,phi_start,phi_end
       !! Modifications by Andreas F. Martitsch (26.01.2016)
       !--> Normalization Bref for stand-alone version
       bmod0 = bmref
       !PRINT *,bmod0
       !STOP
       !PRINT *,es(991)
       !STOP
       !! End Modifications by Andreas F. Martitsch (26.01.2016)
       !! Modifications by Andreas F. Martitsch (27.01.2016)
       ! -> switch on/off use of splined Fourier coefficients within neo_magfie
       IF (isw_spl_fourier_cof .NE. 1) THEN
          ALLOCATE(magfie_sarray(SIZE(es)/flux_surf_dist))
          magfie_sarray = es(1:SIZE(es):flux_surf_dist)
       END IF
       !! End Modifications by Andreas F. Martitsch (27.01.2016)
    END IF
    !*******************************************************************
    ! Spline of surfaces in magfie_sarray
    !*******************************************************************
    IF (magfie_spline .EQ. 1 .AND. magfie_newspline .EQ. 1) THEN
       magfie_sarray_len =  SIZE(magfie_sarray)
       !****************************************************************
       ! Allocation
       !****************************************************************
       ALLOCATE( bmod_spl(4,4,theta_n,phi_n,magfie_sarray_len) )
       IF (isw_eval_spl2d_der .NE. 1) THEN
          ALLOCATE( bb_s_spl(4,4,theta_n,phi_n,magfie_sarray_len) )
          ALLOCATE( bb_tb_spl(4,4,theta_n,phi_n,magfie_sarray_len) )
          ALLOCATE( bb_pb_spl(4,4,theta_n,phi_n,magfie_sarray_len) )
       END IF
       ALLOCATE( gval_spl(4,4,theta_n,phi_n,magfie_sarray_len) )
       !! Modifications by Andreas F. Martitsch (11.03.2014)
       ! Allocate the storage arrays for the 2d spline interpolation
       ! (over the flux-surface) of the additionally needed metric
       ! tensor elements
       IF (isw_eval_bcovars .EQ. 1) THEN
          ALLOCATE( gstb_spl(4,4,theta_n,phi_n,magfie_sarray_len) )
          ALLOCATE( gspb_spl(4,4,theta_n,phi_n,magfie_sarray_len) )
          ALLOCATE( gstb_tb_spl(4,4,theta_n,phi_n,magfie_sarray_len) )
          ALLOCATE( gspb_tb_spl(4,4,theta_n,phi_n,magfie_sarray_len) )
          ALLOCATE( gstb_pb_spl(4,4,theta_n,phi_n,magfie_sarray_len) )
          ALLOCATE( gspb_pb_spl(4,4,theta_n,phi_n,magfie_sarray_len) )
       END IF
       !! End Modifications by Andreas F. Martitsch (11.03.2014)
       !! Modifications by Andreas F. Martitsch (13.11.2014)
!!$       ! Allocate storage arrays for the 2d periodic splines
!!$       ! of the additionally needed quantities for NTV output
!!$       ALLOCATE( R_spl(4,4,theta_n,phi_n,magfie_sarray_len) )
!!$       ALLOCATE( Z_spl(4,4,theta_n,phi_n,magfie_sarray_len) )
!!$       !! End Modifications by Andreas F. Martitsch (13.11.2014)

       ALLOCATE( curr_tor_array(magfie_sarray_len) )
       ALLOCATE( curr_tor_s_array(magfie_sarray_len) )
       ALLOCATE( curr_pol_array(magfie_sarray_len) )
       ALLOCATE( curr_pol_s_array(magfie_sarray_len) )
       ALLOCATE( iota_array(magfie_sarray_len) )
       ALLOCATE( pprime_array(magfie_sarray_len) )
       ALLOCATE( sqrtg00_array(magfie_sarray_len) )
       !****************************************************************
       ! Loop over predefined s-values
       !****************************************************************
       DO k_es = 1, magfie_sarray_len
          s = magfie_sarray(k_es)
          !*************************************************************
          ! Surface
          !*************************************************************
          IF (write_progress .EQ. 1) THEN
             PRINT *, 'Predefined: Initialize Surface, k_es = ',k_es
          END IF
          ALLOCATE ( s_bmnc(mnmax) )
          ALLOCATE ( s_bmnc_s(mnmax) )
          ALLOCATE ( s_rmnc(mnmax) )
          ALLOCATE ( s_zmnc(mnmax) )
          ALLOCATE ( s_lmnc(mnmax) )
          !! Modifications by Andreas F. Martitsch (06.03.2014)
          ! Compute the necessary radial derivatives for the
          ! (R,Z,phi)-components obtained from the 1d spline
          ALLOCATE ( s_rmnc_s(mnmax) ) ! Allocate arrays for additional
          ALLOCATE ( s_zmnc_s(mnmax) ) ! radial derivatives
          ALLOCATE ( s_lmnc_s(mnmax) )
          !! End Modifications by Andreas F. Martitsch (06.03.2014)
          !
          !! Modifications by Andreas F. Martitsch (06.08.2014)
          ! Additional data from Boozer files without Stellarator symmetry
          IF (inp_swi .EQ. 9) THEN        ! ASDEX-U (E. Strumberger)
             ALLOCATE ( s_bmns(mnmax) )
             ALLOCATE ( s_bmns_s(mnmax) )
             ALLOCATE ( s_rmns(mnmax) )
             ALLOCATE ( s_zmns(mnmax) )
             ALLOCATE ( s_lmns(mnmax) )
             ALLOCATE ( s_rmns_s(mnmax) )
             ALLOCATE ( s_zmns_s(mnmax) )
             ALLOCATE ( s_lmns_s(mnmax) )
          END IF
          !! End Modifications by Andreas F. Martitsch (06.08.2014)
          !
          IF (isw_spl_fourier_cof .EQ. 1) THEN
             DO imn = 1, mnmax
                ! Switch swd turns on (1) / off (0) the computation of the
                ! radial derivatives within splint_horner3
                swd = 1
                CALL splint_horner3(es,                                   &
                     a_bmnc(:,imn), b_bmnc(:,imn),                        &
                     c_bmnc(:,imn), d_bmnc(:,imn),                        &
                     swd, r_mhalf(imn),                                   &
                     s, tf, tfp, tfpp, tfppp,                             &
                     s_bmnc(imn), s_bmnc_s(imn), ypp, yppp)
                !swd = 0 ! Now we want to compute the radial derivatives
                swd = 1
                CALL splint_horner3(es,                                   &
                     a_rmnc(:,imn), b_rmnc(:,imn),                        &
                     c_rmnc(:,imn), d_rmnc(:,imn),                        &
                     swd, r_mhalf(imn),                                   &
                     s, tf, tfp, tfpp, tfppp,                             &
                     s_rmnc(imn), s_rmnc_s(imn), ypp, yppp)
                swd = 1
                CALL splint_horner3(es,                                   &
                     a_zmnc(:,imn), b_zmnc(:,imn),                        &
                     c_zmnc(:,imn), d_zmnc(:,imn),                        &
                     swd, r_mhalf(imn),                                   &
                     s, tf, tfp, tfpp, tfppp,                             &
                     s_zmnc(imn), s_zmnc_s(imn), ypp, yppp)
                swd = 1
                CALL splint_horner3(es,                                   &
                     a_lmnc(:,imn), b_lmnc(:,imn),                        &
                     c_lmnc(:,imn), d_lmnc(:,imn),                        &
                     swd, r_mhalf(imn),                                   &
                     s, tf, tfp, tfpp, tfppp,                             &
                     s_lmnc(imn), s_lmnc_s(imn), ypp, yppp)
                !
                !! Modifications by Andreas F. Martitsch (06.08.2014)
                ! Additional data from Boozer files without Stellarator symmetry
                IF (inp_swi .EQ. 9) THEN        ! ASDEX-U (E. Strumberger)
                   swd = 1
                   CALL splint_horner3(es,                                   &
                        a_bmns(:,imn), b_bmns(:,imn),                        &
                        c_bmns(:,imn), d_bmns(:,imn),                        &
                        swd, r_mhalf(imn),                                   &
                        s, tf, tfp, tfpp, tfppp,                             &
                        s_bmns(imn), s_bmns_s(imn), ypp, yppp)
                   swd = 1
                   CALL splint_horner3(es,                                   &
                        a_rmns(:,imn), b_rmns(:,imn),                        &
                        c_rmns(:,imn), d_rmns(:,imn),                        &
                        swd, r_mhalf(imn),                                   &
                        s, tf, tfp, tfpp, tfppp,                             &
                        s_rmns(imn), s_rmns_s(imn), ypp, yppp)
                   swd = 1
                   CALL splint_horner3(es,                                   &
                        a_zmns(:,imn), b_zmns(:,imn),                        &
                        c_zmns(:,imn), d_zmns(:,imn),                        &
                        swd, r_mhalf(imn),                                   &
                        s, tf, tfp, tfpp, tfppp,                             &
                        s_zmns(imn), s_zmns_s(imn), ypp, yppp)
                   swd = 1
                   CALL splint_horner3(es,                                   &
                        a_lmns(:,imn), b_lmns(:,imn),                        &
                        c_lmns(:,imn), d_lmns(:,imn),                        &
                        swd, r_mhalf(imn),                                   &
                        s, tf, tfp, tfpp, tfppp,                             &
                        s_lmns(imn), s_lmns_s(imn), ypp, yppp)
                END IF
                !! End Modifications by Andreas F. Martitsch (06.08.2014)
                !
             END DO
          ELSE
             s_bmnc = a_bmnc(1+(k_es-1)*flux_surf_dist,:)
             s_bmnc_s = 0.0d0
             s_rmnc = a_rmnc(1+(k_es-1)*flux_surf_dist,:)
             s_rmnc_s = 0.0d0
             s_zmnc = a_zmnc(1+(k_es-1)*flux_surf_dist,:)
             s_zmnc_s = 0.0d0
             s_lmnc = a_lmnc(1+(k_es-1)*flux_surf_dist,:)
             s_lmnc_s = 0.0d0
             IF (inp_swi .EQ. 9) THEN        ! ASDEX-U (E. Strumberger)
                s_bmns = a_bmns(1+(k_es-1)*flux_surf_dist,:)
                s_bmns_s = 0.0d0
                s_rmns = a_bmns(1+(k_es-1)*flux_surf_dist,:)
                s_rmns_s = 0.0d0
                s_zmns = a_bmns(1+(k_es-1)*flux_surf_dist,:)
                s_zmns_s = 0.0d0
                s_lmns = a_bmns(1+(k_es-1)*flux_surf_dist,:)
                s_lmns_s = 0.0d0
             END IF
          END IF
          !*************************************************************
          ! Fourier summation for the full theta-phi array
          !*************************************************************
          IF (write_progress .EQ. 1) THEN
             PRINT *, 'Do Fourier'
          END IF
          ALLOCATE( bmod_a(theta_n,phi_n) )
          ALLOCATE( bb_s_a(theta_n,phi_n) )
          ALLOCATE( bb_tb_a(theta_n,phi_n) )
          ALLOCATE( bb_pb_a(theta_n,phi_n) )
          bmod_a  = 0.0_dp
          bb_s_a  = 0.0_dp
          bb_tb_a = 0.0_dp
          bb_pb_a = 0.0_dp

          ALLOCATE( r(theta_n,phi_n) )  ! NEW
          ALLOCATE( z(theta_n,phi_n) )
          ALLOCATE( l(theta_n,phi_n) )
          ALLOCATE( r_tb(theta_n,phi_n) )
          ALLOCATE( z_tb(theta_n,phi_n) )
          ALLOCATE( p_tb(theta_n,phi_n) )
          ALLOCATE( r_pb(theta_n,phi_n) )
          ALLOCATE( z_pb(theta_n,phi_n) )
          ALLOCATE( p_pb(theta_n,phi_n) )
          r = 0.0d0
          z = 0.0d0
          l = 0.0d0
          r_tb = 0.0d0
          z_tb = 0.0d0
          p_tb = 0.0d0
          r_pb = 0.0d0
          z_pb = 0.0d0
          p_pb = 0.0d0
          !! Modifications by Andreas F. Martitsch (11.03.2014)
          ! Allocate temporary storage arrays for the Fourier summations
          ! related to the radial derivatives of (R,Z,phi)-components
          ALLOCATE( r_s(theta_n,phi_n) )
          ALLOCATE( z_s(theta_n,phi_n) )
          ALLOCATE( p_s(theta_n,phi_n) )
          ALLOCATE( r_stb(theta_n,phi_n) )
          ALLOCATE( z_stb(theta_n,phi_n) )
          ALLOCATE( p_stb(theta_n,phi_n) )
          ALLOCATE( r_tbtb(theta_n,phi_n) )
          ALLOCATE( z_tbtb(theta_n,phi_n) )
          ALLOCATE( p_tbtb(theta_n,phi_n) )
          ALLOCATE( r_pbtb(theta_n,phi_n) )
          ALLOCATE( z_pbtb(theta_n,phi_n) )
          ALLOCATE( p_pbtb(theta_n,phi_n) )
          ALLOCATE( r_spb(theta_n,phi_n) )
          ALLOCATE( z_spb(theta_n,phi_n) )
          ALLOCATE( p_spb(theta_n,phi_n) )
          ALLOCATE( r_pbpb(theta_n,phi_n) )
          ALLOCATE( z_pbpb(theta_n,phi_n) )
          ALLOCATE( p_pbpb(theta_n,phi_n) )
          r_s = 0.0d0
          z_s = 0.0d0
          p_s = 0.0d0
          r_stb = 0.0d0
          z_stb = 0.0d0
          p_stb = 0.0d0
          r_tbtb = 0.0d0
          z_tbtb = 0.0d0
          p_tbtb = 0.0d0
          r_pbtb = 0.0d0
          z_pbtb = 0.0d0
          p_pbtb = 0.0d0
          r_spb = 0.0d0
          z_spb = 0.0d0
          p_spb = 0.0d0
          r_pbpb = 0.0d0
          z_pbpb = 0.0d0
          p_pbpb = 0.0d0
          !! End Modifications by Andreas F. Martitsch (11.03.2014)

          DO imn=1,mnmax
             ri = s_rmnc(imn) ! NEW
             zi = s_zmnc(imn)
             li = s_lmnc(imn)
             !! Modifications by Andreas F. Martitsch (07.03.2014)
             ! Auxiliary variables for the Fourier summation
             ri_s = s_rmnc_s(imn)
             zi_s = s_zmnc_s(imn)
             li_s = s_lmnc_s(imn)
             !! End Modifications by Andreas F. Martitsch (07.03.2014)
             bi   = s_bmnc(imn)
             bi_s = s_bmnc_s(imn)
             !! Modifications by Andreas F. Martitsch (06.08.2014)
             ! Additional data from Boozer files without Stellarator symmetry
             IF (inp_swi .EQ. 9) THEN        ! ASDEX-U (E. Strumberger)
                !
                ris = s_rmns(imn)
                zis = s_zmns(imn)
                lis = s_lmns(imn)
                bis = s_bmns(imn)
                !
                ris_s = s_rmns_s(imn)
                zis_s = s_zmns_s(imn)
                lis_s = s_lmns_s(imn)
                bis_s = s_bmns_s(imn)
                !
             END IF
             !! End Modifications by Andreas F. Martitsch (06.08.2014)
             m = ixm(imn)
             n = ixn(imn)
             !IF (n .EQ. 0) CYCLE
             !IF (n .NE. 0) bi=bi*1000.0d0
             !IF (n .NE. 0 .and. inp_swi .EQ. 9) bis=bis*1000.0d0
             im = pixm(imn)
             in = pixn(imn)
             DO ip=1,phi_n
                DO it=1,theta_n
                   !! Modifications by Andreas F. Martitsch (06.08.2014)
                   ! Additional data from Boozer files without Stellarator symmetry
                   IF (inp_swi .EQ. 9) THEN        ! ASDEX-U (E. Strumberger)
                      cosv = cosmth(it,im) * cosnph(ip,in) - sinmth(it,im) * sinnph(ip,in)
                      sinv = sinmth(it,im) * cosnph(ip,in) + cosmth(it,im) * sinnph(ip,in)

                      bmod_a(it,ip) = bmod_a(it,ip)     + bi*cosv   + bis*sinv
                      bb_s_a(it,ip) = bb_s_a(it,ip)     + bi_s*cosv + bis_s*sinv
                      bb_tb_a(it,ip)  = bb_tb_a(it,ip)  - m*bi*sinv + m*bis*cosv
                      bb_pb_a(it,ip)  = bb_pb_a(it,ip)  - n*bi*sinv + n*bis*cosv

                      r(it,ip) = r(it,ip) + ri*cosv + ris*sinv
                      z(it,ip) = z(it,ip) + zi*cosv + zis*sinv
                      l(it,ip) = l(it,ip) + li*cosv + lis*sinv

                      r_tb(it,ip) = r_tb(it,ip) - m*ri*sinv + m*ris*cosv
                      r_pb(it,ip) = r_pb(it,ip) - n*ri*sinv + n*ris*cosv
                      z_tb(it,ip) = z_tb(it,ip) - m*zi*sinv + m*zis*cosv
                      z_pb(it,ip) = z_pb(it,ip) - n*zi*sinv + n*zis*cosv
                      !! Modifications by Andreas F. Martitsch (12.11.2015)
                      !According to Erika Strumberger (Email 11.10.2015)
                      !the conversion from phi_b to phi is given by
                      !"\phi-phi_b = 2\pi/N_p \sum ( c \cos(2\pi (m u + n v) ) + s \sin(2\pi (m u+n v) ) )"
                      !where  \phi=2\pi/N_p v.
                      !This expression differs by a minus sign from the
                      !expression used by J. Geiger ( phi_b-\phi = ... )!
                      !-> previous versions used this definition:
                      !p_tb(it,ip) = p_tb(it,ip) + m*li*sinv - m*lis*cosv ! -l_tb
                      !p_pb(it,ip) = p_pb(it,ip) + n*li*sinv - n*lis*cosv ! -l_pb
                      !-> corrected formulas:
                      p_tb(it,ip) = p_tb(it,ip) - m*li*sinv + m*lis*cosv ! +l_tb
                      p_pb(it,ip) = p_pb(it,ip) - n*li*sinv + n*lis*cosv ! +l_pb
                      !! End Modifications by Andreas F. Martitsch (12.11.2015)

                      r_s(it,ip) = r_s(it,ip) + ri_s*cosv + ris_s*sinv
                      z_s(it,ip) = z_s(it,ip) + zi_s*cosv + zis_s*sinv
                      !! Modifications by Andreas F. Martitsch (12.11.2015)
                      !According to Erika Strumberger (Email 11.10.2015)
                      !the conversion from phi_b to phi is given by
                      !"\phi-phi_b = 2\pi/N_p \sum ( c \cos(2\pi (m u + n v) ) + s \sin(2\pi (m u+n v) ) )"
                      !where  \phi=2\pi/N_p v.
                      !This expression differs by a minus sign from the
                      !expression used by J. Geiger ( phi_b-\phi = ... )!
                      !-> previous versions used this definition:
                      !p_s(it,ip) = p_s(it,ip) - li_s*cosv - lis_s*sinv ! -l_s
                      !-> corrected formulas:
                      p_s(it,ip) = p_s(it,ip) + li_s*cosv + lis_s*sinv ! +l_s
                      !! End Modifications by Andreas F. Martitsch (12.11.2015)

                      r_stb(it,ip) = r_stb(it,ip) - m*ri_s*sinv + m*ris_s*cosv
                      z_stb(it,ip) = z_stb(it,ip) - m*zi_s*sinv + m*zis_s*cosv
                      !! Modifications by Andreas F. Martitsch (12.11.2015)
                      !According to Erika Strumberger (Email 11.10.2015)
                      !the conversion from phi_b to phi is given by
                      !"\phi-phi_b = 2\pi/N_p \sum ( c \cos(2\pi (m u + n v) ) + s \sin(2\pi (m u+n v) ) )"
                      !where  \phi=2\pi/N_p v.
                      !This expression differs by a minus sign from the
                      !expression used by J. Geiger ( phi_b-\phi = ... )!
                      !-> previous versions used this definition:
                      !p_stb(it,ip) = p_stb(it,ip) + m*li_s*sinv - m*lis_s*cosv ! -l_stb
                      !-> corrected formulas:
                      p_stb(it,ip) = p_stb(it,ip) - m*li_s*sinv + m*lis_s*cosv ! +l_stb
                      !! End Modifications by Andreas F. Martitsch (12.11.2015)

                      r_spb(it,ip) = r_spb(it,ip) - n*ri_s*sinv + n*ris_s*cosv
                      z_spb(it,ip) = z_spb(it,ip) - n*zi_s*sinv + n*zis_s*cosv
                      !! Modifications by Andreas F. Martitsch (12.11.2015)
                      !According to Erika Strumberger (Email 11.10.2015)
                      !the conversion from phi_b to phi is given by
                      !"\phi-phi_b = 2\pi/N_p \sum ( c \cos(2\pi (m u + n v) ) + s \sin(2\pi (m u+n v) ) )"
                      !where  \phi=2\pi/N_p v.
                      !This expression differs by a minus sign from the
                      !expression used by J. Geiger ( phi_b-\phi = ... )!
                      !-> previous versions used this definition:
                      !p_spb(it,ip) = p_spb(it,ip) + n*li_s*sinv - n*lis_s*cosv ! -l_spb
                      !-> corrected formulas:
                      p_spb(it,ip) = p_spb(it,ip) - n*li_s*sinv + n*lis_s*cosv ! +l_spb
                      !! End Modifications by Andreas F. Martitsch (12.11.2015)

                      r_tbtb(it,ip) = r_tbtb(it,ip) - m*m*ri*cosv - m*m*ris*sinv
                      z_tbtb(it,ip) = z_tbtb(it,ip) - m*m*zi*cosv - m*m*zis*sinv
                      !! Modifications by Andreas F. Martitsch (12.11.2015)
                      !According to Erika Strumberger (Email 11.10.2015)
                      !the conversion from phi_b to phi is given by
                      !"\phi-phi_b = 2\pi/N_p \sum ( c \cos(2\pi (m u + n v) ) + s \sin(2\pi (m u+n v) ) )"
                      !where  \phi=2\pi/N_p v.
                      !This expression differs by a minus sign from the
                      !expression used by J. Geiger ( phi_b-\phi = ... )!
                      !-> previous versions used this definition:
                      !p_tbtb(it,ip) = p_tbtb(it,ip) + m*m*li*cosv + m*m*lis*sinv ! -l_tbtb
                      !-> corrected formulas:
                      p_tbtb(it,ip) = p_tbtb(it,ip) - m*m*li*cosv - m*m*lis*sinv ! +l_tbtb
                      !! End Modifications by Andreas F. Martitsch (12.11.2015)

                      r_pbtb(it,ip) = r_pbtb(it,ip) - m*n*ri*cosv - m*n*ris*sinv
                      z_pbtb(it,ip) = z_pbtb(it,ip) - m*n*zi*cosv - m*n*zis*sinv
                      !! Modifications by Andreas F. Martitsch (12.11.2015)
                      !According to Erika Strumberger (Email 11.10.2015)
                      !the conversion from phi_b to phi is given by
                      !"\phi-phi_b = 2\pi/N_p \sum ( c \cos(2\pi (m u + n v) ) + s \sin(2\pi (m u+n v) ) )"
                      !where  \phi=2\pi/N_p v.
                      !This expression differs by a minus sign from the
                      !expression used by J. Geiger ( phi_b-\phi = ... )!
                      !-> previous versions used this definition:
                      !p_pbtb(it,ip) = p_pbtb(it,ip) + m*n*li*cosv + m*n*lis*sinv ! -l_pbtb
                      !-> corrected formulas:
                      p_pbtb(it,ip) = p_pbtb(it,ip) - m*n*li*cosv - m*n*lis*sinv ! +l_pbtb
                      !! End Modifications by Andreas F. Martitsch (12.11.2015)

                      r_pbpb(it,ip) = r_pbpb(it,ip) - n*n*ri*cosv - n*n*ris*sinv
                      z_pbpb(it,ip) = z_pbpb(it,ip) - n*n*zi*cosv - n*n*zis*sinv
                      !! Modifications by Andreas F. Martitsch (12.11.2015)
                      !According to Erika Strumberger (Email 11.10.2015)
                      !the conversion from phi_b to phi is given by
                      !"\phi-phi_b = 2\pi/N_p \sum ( c \cos(2\pi (m u + n v) ) + s \sin(2\pi (m u+n v) ) )"
                      !where  \phi=2\pi/N_p v.
                      !This expression differs by a minus sign from the
                      !expression used by J. Geiger ( phi_b-\phi = ... )!
                      !-> previous versions used this definition:
                      !p_pbpb(it,ip) = p_pbpb(it,ip) + n*n*li*cosv + n*n*lis*sinv ! -l_pbpb
                      !-> corrected formulas:
                      p_pbpb(it,ip) = p_pbpb(it,ip) - n*n*li*cosv - n*n*lis*sinv ! +l_pbpb
                      !! End Modifications by Andreas F. Martitsch (12.11.2015)
                   ELSE
                      cosv = cosmth(it,im) * cosnph(ip,in) + sinmth(it,im) * sinnph(ip,in)
                      sinv = sinmth(it,im) * cosnph(ip,in) - cosmth(it,im) * sinnph(ip,in)

                      bmod_a(it,ip)   = bmod_a(it,ip)   +     bi   * cosv
                      bb_s_a(it,ip)   = bb_s_a(it,ip)   +     bi_s * cosv
                      bb_tb_a(it,ip)  = bb_tb_a(it,ip)  - m * bi   * sinv
                      bb_pb_a(it,ip)  = bb_pb_a(it,ip)  + n * bi   * sinv

                      r(it,ip) = r(it,ip) + ri*cosv
                      z(it,ip) = z(it,ip) + zi*sinv
                      l(it,ip) = l(it,ip) + li*sinv

                      r_tb(it,ip) = r_tb(it,ip) - m*ri*sinv
                      r_pb(it,ip) = r_pb(it,ip) + n*ri*sinv
                      z_tb(it,ip) = z_tb(it,ip) + m*zi*cosv
                      z_pb(it,ip) = z_pb(it,ip) - n*zi*cosv
                      p_tb(it,ip) = p_tb(it,ip) - m*li*cosv ! -l_tb
                      p_pb(it,ip) = p_pb(it,ip) + n*li*cosv ! -l_pb
                      !! Modifications by Andreas F. Martitsch (07.03.2014)
                      ! Temporary storage arrays for the Fourier summations
                      ! related to the radial derivatives of (R,Z,phi)-components
                      r_s(it,ip) = r_s(it,ip) + ri_s*cosv
                      z_s(it,ip) = z_s(it,ip) + zi_s*sinv
                      p_s(it,ip) = p_s(it,ip) - li_s*sinv ! -l_s
                      r_stb(it,ip) = r_stb(it,ip) - m*ri_s*sinv
                      z_stb(it,ip) = z_stb(it,ip) + m*zi_s*cosv
                      p_stb(it,ip) = p_stb(it,ip) - m*li_s*cosv ! -l_stb
                      r_tbtb(it,ip) = r_tbtb(it,ip) - m*m*ri*cosv
                      z_tbtb(it,ip) = z_tbtb(it,ip) - m*m*zi*sinv
                      p_tbtb(it,ip) = p_tbtb(it,ip) + m*m*li*sinv ! -l_tbtb
                      r_pbtb(it,ip) = r_pbtb(it,ip) + m*n*ri*cosv
                      z_pbtb(it,ip) = z_pbtb(it,ip) + m*n*zi*sinv
                      p_pbtb(it,ip) = p_pbtb(it,ip) - m*n*li*sinv ! -l_pbtb
                      r_spb(it,ip) = r_spb(it,ip) + n*ri_s*sinv
                      z_spb(it,ip) = z_spb(it,ip) - n*zi_s*cosv
                      p_spb(it,ip) = p_spb(it,ip) + n*li_s*cosv ! -l_spb
                      r_pbpb(it,ip) = r_pbpb(it,ip) - n*n*ri*cosv
                      z_pbpb(it,ip) = z_pbpb(it,ip) - n*n*zi*sinv
                      p_pbpb(it,ip) = p_pbpb(it,ip) + n*n*li*sinv ! -l_pbpb
                      !! End Modifications by Andreas F. Martitsch (07.03.2014)
                   END IF
                   !! End Modifications by Andreas F. Martitsch (06.08.2014)
                END DO
             END DO
          END DO
          DEALLOCATE( s_bmnc )
          DEALLOCATE( s_bmnc_s )
          DEALLOCATE( s_rmnc )
          DEALLOCATE( s_zmnc )
          DEALLOCATE( s_lmnc )
          !! Modifications by Andreas F. Martitsch (07.03.2014)
          ! Deallocate arrays for the additional radial derivatives
          DEALLOCATE( s_rmnc_s )
          DEALLOCATE( s_zmnc_s )
          DEALLOCATE( s_lmnc_s )
          !! End Modifications by Andreas F. Martitsch (07.03.2014)
          !
          !! Modifications by Andreas F. Martitsch (06.08.2014)
          ! Additional data from Boozer files without Stellarator symmetry
          IF (inp_swi .EQ. 9) THEN        ! ASDEX-U (E. Strumberger)
             DEALLOCATE ( s_bmns )
             DEALLOCATE ( s_bmns_s )
             DEALLOCATE ( s_rmns )
             DEALLOCATE ( s_zmns )
             DEALLOCATE ( s_lmns )
             DEALLOCATE ( s_rmns_s )
             DEALLOCATE ( s_zmns_s )
             DEALLOCATE ( s_lmns_s )
          END IF
          !! End Modifications by Andreas F. Martitsch (06.08.2014)
          !
          IF (lab_swi .EQ. 5 .OR. lab_swi .EQ. 3) THEN ! CHS, LHD
             p_tb = - p_tb
             p_pb = 1 - p_pb
             !! Modifications by Andreas F. Martitsch (07.03.2014)
             ! ToDo: Implement conversion for p_s
             PRINT *,'WARNING FROM NEO_MAGFIE: CONVERSION FOR RADIAL DERIVATIVE OF BOOZER-PHI NOT IMPLEMENTED!'
             !! End Modifications by Andreas F. Martitsch (07.03.2014)
          ELSE
             p_tb = p_tb * twopi / nfp
             p_pb = 1.0_dp + p_pb * twopi / nfp
             !! Modifications by Andreas F. Martitsch (11.03.2014)
             ! Conversion factor between Boozer-phi and L (see Boozer-file)
             p_s = p_s * twopi / nfp
             p_stb = p_stb * twopi / nfp
             p_tbtb = p_tbtb * twopi / nfp
             p_pbtb = p_pbtb * twopi / nfp
             p_spb = p_spb * twopi / nfp
             p_pbpb = p_pbpb * twopi / nfp
             !! End Modifications by Andreas F. Martitsch (11.03.2014)
          END IF

          ! **********************************************************************
          ! Ensure periodicity boundaries to be the same
          ! **********************************************************************
          r(theta_n,:) = r(1,:)
          r(:,phi_n)   = r(:,1)
          z(theta_n,:) = z(1,:)
          z(:,phi_n)   = z(:,1)
          l(theta_n,:) = l(1,:)
          l(:,phi_n)   = l(:,1)
          bmod_a(theta_n,:) = bmod_a(1,:)
          bmod_a(:,phi_n)   = bmod_a(:,1)
          r_tb(theta_n,:) = r_tb(1,:)
          r_tb(:,phi_n)   = r_tb(:,1)
          r_pb(theta_n,:) = r_pb(1,:)
          r_pb(:,phi_n)   = r_pb(:,1)
          z_tb(theta_n,:) = z_tb(1,:)
          z_tb(:,phi_n)   = z_tb(:,1)
          z_pb(theta_n,:) = z_pb(1,:)
          z_pb(:,phi_n)   = z_pb(:,1)
          p_tb(theta_n,:) = p_tb(1,:)
          p_tb(:,phi_n)   = p_tb(:,1)
          p_pb(theta_n,:) = p_pb(1,:)
          p_pb(:,phi_n)   = p_pb(:,1)
          !! Modifications by Andreas F. Martitsch (07.03.2014)
          ! Temporary storage arrays for the Fourier summations
          ! related to the radial derivatives of (R,Z,phi)-components
          r_s(theta_n,:) = r_s(1,:)
          r_s(:,phi_n)   = r_s(:,1)
          z_s(theta_n,:) = z_s(1,:)
          z_s(:,phi_n)   = z_s(:,1)
          p_s(theta_n,:) = p_s(1,:)
          p_s(:,phi_n)   = p_s(:,1)
          r_stb(theta_n,:) = r_stb(1,:)
          r_stb(:,phi_n)   = r_stb(:,1)
          z_stb(theta_n,:) = z_stb(1,:)
          z_stb(:,phi_n)   = z_stb(:,1)
          p_stb(theta_n,:) = p_stb(1,:)
          p_stb(:,phi_n)   = p_stb(:,1)
          r_tbtb(theta_n,:) = r_tbtb(1,:)
          r_tbtb(:,phi_n)   = r_tbtb(:,1)
          z_tbtb(theta_n,:) = z_tbtb(1,:)
          z_tbtb(:,phi_n)   = z_tbtb(:,1)
          p_tbtb(theta_n,:) = p_tbtb(1,:)
          p_tbtb(:,phi_n)   = p_tbtb(:,1)
          r_pbtb(theta_n,:) = r_pbtb(1,:)
          r_pbtb(:,phi_n)   = r_pbtb(:,1)
          z_pbtb(theta_n,:) = z_pbtb(1,:)
          z_pbtb(:,phi_n)   = z_pbtb(:,1)
          p_pbtb(theta_n,:) = p_pbtb(1,:)
          p_pbtb(:,phi_n)   = p_pbtb(:,1)
          r_spb(theta_n,:) = r_spb(1,:)
          r_spb(:,phi_n)   = r_spb(:,1)
          z_spb(theta_n,:) = z_spb(1,:)
          z_spb(:,phi_n)   = z_spb(:,1)
          p_spb(theta_n,:) = p_spb(1,:)
          p_spb(:,phi_n)   = p_spb(:,1)
          r_pbpb(theta_n,:) = r_pbpb(1,:)
          r_pbpb(:,phi_n)   = r_pbpb(:,1)
          z_pbpb(theta_n,:) = z_pbpb(1,:)
          z_pbpb(:,phi_n)   = z_pbpb(:,1)
          p_pbpb(theta_n,:) = p_pbpb(1,:)
          p_pbpb(:,phi_n)   = p_pbpb(:,1)
          !! End Modifications by Andreas F. Martitsch (07.03.2014)
          bb_tb_a(theta_n,:) = bb_tb_a(1,:)
          bb_tb_a(:,phi_n)   = bb_tb_a(:,1)
          !! Modifications by Andreas F. Martitsch (25.08.2014)
          ! This seems to be a copy-paste, which moved one version
          ! to another:
          !bb_s_a(theta_n,:)  = bb_tb_a(1,:)
          !bb_s_a(:,phi_n)    = bb_tb_a(:,1)
          ! this should be correct now:
          bb_s_a(theta_n,:)  = bb_s_a(1,:)
          bb_s_a(:,phi_n)    = bb_s_a(:,1)
          !! End Modifications by Andreas F. Martitsch (25.08.2014)
          bb_pb_a(theta_n,:) = bb_pb_a(1,:)
          bb_pb_a(:,phi_n)   = bb_pb_a(:,1)

          ! **********************************************************************
          ! Derived quantities
          ! **********************************************************************
          ALLOCATE( gtbtb(theta_n,phi_n) )
          ALLOCATE( gpbpb(theta_n,phi_n) )
          ALLOCATE( gtbpb(theta_n,phi_n) )
          ALLOCATE( sqrg11_met(theta_n,phi_n) )
          !! Modifications by Andreas F. Martitsch (11.03.2014)
          ! Allocate temporary storage arrays for the Fourier summations
          ! related to the additionally needed metric tensor elements
          ALLOCATE( gstb_a(theta_n,phi_n) )
          ALLOCATE( gspb_a(theta_n,phi_n) )
          ALLOCATE( gstb_tb_a(theta_n,phi_n) )
          ALLOCATE( gspb_tb_a(theta_n,phi_n) )
          ALLOCATE( gstb_pb_a(theta_n,phi_n) )
          ALLOCATE( gspb_pb_a(theta_n,phi_n) )
          !! End Modifications by Andreas F. Martitsch (11.03.2014)
          ! metric tensor
          gtbtb = r_tb*r_tb + z_tb*z_tb + r*r*p_tb*p_tb
          gpbpb = r_pb*r_pb + z_pb*z_pb + r*r*p_pb*p_pb
          gtbpb = r_tb*r_pb + z_tb*z_pb + r*r*p_tb*p_pb
          !! Modifications by Andreas F. Martitsch (11.03.2014)
          ! Compute the additionally needed metric tensor elements
          gstb_a = r_s*r_tb + z_s*z_tb + r*r*p_s*p_tb
          gspb_a = r_s*r_pb + z_s*z_pb + r*r*p_s*p_pb
          gstb_tb_a = r_stb*r_tb + r_s*r_tbtb  + z_stb*z_tb + z_s*z_tbtb + &
               2.0d0*r*r_tb*p_s*p_tb + r*r*(p_stb*p_tb + p_s*p_tbtb)
          gspb_tb_a = r_stb*r_pb + r_s*r_pbtb  + z_stb*z_pb + z_s*z_pbtb + &
               2.0d0*r*r_tb*p_s*p_pb + r*r*(p_stb*p_pb + p_s*p_pbtb)
          gstb_pb_a = r_spb*r_tb + r_s*r_pbtb  + z_spb*z_tb + z_s*z_pbtb + &
               2.0d0*r*r_pb*p_s*p_tb + r*r*(p_spb*p_tb + p_s*p_pbtb)
          gspb_pb_a = r_spb*r_pb + r_s*r_pbpb  + z_spb*z_pb + z_s*z_pbpb + &
               2.0d0*r*r_pb*p_s*p_pb + r*r*(p_spb*p_pb + p_s*p_pbpb)
          !! End Modifications by Andreas F. Martitsch (11.03.2014)

          ! Winny for Klaus
          av_b2_m = theta_n * phi_n / SUM(1 / (bmod_a*bmod_a))
          !PRINT *, 'theta_n,phi_n ',theta_n,phi_n
          !PRINT *, 'av_b2_m ',av_b2_m
          ! Winny for Klaus - Ende


          ! $1/sqrt(g)$
          !! fac = s_curr_pol + s_iota * s_curr_tor  ! (J + iota I)
          !! isqrg  = b*b / fac
          ! $sqrt(g^{11})$
          !IF (g11_swi .EQ. 0) THEN
          !   sqrg11 = 1.0_dp
          !ELSE
             !sqrg11 = SQRT( (gtbtb*gpbpb - gtbpb*gtbpb ) * isqrg**2)
             sqrg11_met = SQRT( (gtbtb*gpbpb - gtbpb*gtbpb ) )
          !END IF

             !PRINT *, 'max_gtbtb = ',maxval(gtbtb)
             !PRINT *, 'min_gtbtb = ',MINVAL(gtbtb)
             !PRINT *, 'max_gpbpb = ',maxval(gpbpb)
             !PRINT *, 'min_gpbpb = ',MINVAL(gpbpb)
             !PRINT *, 'max_gtbpb = ',MAXVAL(gtbpb)
             !PRINT *, 'min_gtbpb = ',MINVAL(gtbpb)

          !*************************************************************
          ! Do the 2-D periodic spline
          !*************************************************************
          IF (write_progress .EQ. 1) THEN
             PRINT *, 'Do 2-D spline of bmod'
          END IF
          p_spl => bmod_spl(:,:,:,:,k_es)
          CALL spl2d(theta_n,phi_n,theta_int,phi_int,mt,mp,            &
               bmod_a,p_spl)
          IF (isw_eval_spl2d_der .NE. 1) THEN
             IF (write_progress .EQ. 1) THEN
                PRINT *, 'Do 2-D spline of bb_s'
             END IF
             p_spl => bb_s_spl(:,:,:,:,k_es)
             CALL spl2d(theta_n,phi_n,theta_int,phi_int,mt,mp,            &
                  bb_s_a,p_spl)
             IF (write_progress .EQ. 1) THEN
                PRINT *, 'Do 2-D spline of bb_tb'
             END IF
             p_spl => bb_tb_spl(:,:,:,:,k_es)
             CALL spl2d(theta_n,phi_n,theta_int,phi_int,mt,mp,            &
                  bb_tb_a,p_spl)
             IF (write_progress .EQ. 1) THEN
                PRINT *, 'Do 2-D spline of bb_pb'
             END IF
             p_spl => bb_pb_spl(:,:,:,:,k_es)
             CALL spl2d(theta_n,phi_n,theta_int,phi_int,mt,mp,            &
                  bb_pb_a,p_spl)
          END IF
          IF (write_progress .EQ. 1) THEN
             PRINT *, 'Do 2-D spline of sqrg11'
          END IF
          p_spl => gval_spl(:,:,:,:,k_es)
          CALL spl2d(theta_n,phi_n,theta_int,phi_int,mt,mp,            &
               sqrg11_met,p_spl)
          !! Modifications by Andreas F. Martitsch (11.03.2014)
          ! Compute the 2d periodic splines (over the flux-surface)
          ! of the additionally needed metric tensor elements
          IF (isw_eval_bcovars .EQ. 1) THEN
             p_spl => gstb_spl(:,:,:,:,k_es)
             CALL spl2d(theta_n,phi_n,theta_int,phi_int,mt,mp,            &
                  gstb_a,p_spl)
             p_spl => gspb_spl(:,:,:,:,k_es)
             CALL spl2d(theta_n,phi_n,theta_int,phi_int,mt,mp,            &
                  gspb_a,p_spl)
             p_spl => gstb_tb_spl(:,:,:,:,k_es)
             CALL spl2d(theta_n,phi_n,theta_int,phi_int,mt,mp,            &
                  gstb_tb_a,p_spl)
             p_spl => gspb_tb_spl(:,:,:,:,k_es)
             CALL spl2d(theta_n,phi_n,theta_int,phi_int,mt,mp,            &
                  gspb_tb_a,p_spl)
             p_spl => gstb_pb_spl(:,:,:,:,k_es)
             CALL spl2d(theta_n,phi_n,theta_int,phi_int,mt,mp,            &
                  gstb_pb_a,p_spl)
             p_spl => gspb_pb_spl(:,:,:,:,k_es)
             CALL spl2d(theta_n,phi_n,theta_int,phi_int,mt,mp,            &
                  gspb_pb_a,p_spl)
          END IF
          !! End Modifications by Andreas F. Martitsch (11.03.2014)
!!$          !! Modifications by Andreas F. Martitsch (13.11.2014)
!!$          ! Compute the 2d periodic splines (over the flux-surface)
!!$          ! of the additionally needed quantities for NTV output
!!$          p_spl => R_spl(:,:,:,:,k_es)
!!$          CALL spl2d(theta_n,phi_n,theta_int,phi_int,mt,mp,            &
!!$               r,p_spl)
!!$          p_spl => Z_spl(:,:,:,:,k_es)
!!$          CALL spl2d(theta_n,phi_n,theta_int,phi_int,mt,mp,            &
!!$               z,p_spl)
!!$          !! End Modifications by Andreas F. Martitsch (13.11.2014)

          !PRINT *, 'max_g11 = ', maxval(sqrg11_met)
          !PRINT *, 'min_g11 = ', minval(sqrg11_met)

          !stop "Compute magnetic field components via a spline (and not via direct Fourier summation)"

          DEALLOCATE( bmod_a )
          DEALLOCATE( bb_s_a )
          DEALLOCATE( bb_tb_a )
          DEALLOCATE( bb_pb_a )

          DEALLOCATE( r )  ! NEW
          DEALLOCATE( z )
          DEALLOCATE( l )
          DEALLOCATE( r_tb )
          DEALLOCATE( z_tb )
          DEALLOCATE( p_tb )
          DEALLOCATE( r_pb )
          DEALLOCATE( z_pb )
          DEALLOCATE( p_pb )

          DEALLOCATE( gtbtb )
          DEALLOCATE( gpbpb )
          DEALLOCATE( gtbpb )
          DEALLOCATE( sqrg11_met )

          !! Modifications by Andreas F. Martitsch (11.03.2014)
          ! Deallocate temporary storage arrays for the Fourier summations
          ! related to the radial derivatives of (R,Z,phi)-components
          ! and the additionally needed metric tensor elements
          DEALLOCATE( r_s )
          DEALLOCATE( z_s )
          DEALLOCATE( p_s )
          DEALLOCATE( r_stb )
          DEALLOCATE( z_stb )
          DEALLOCATE( p_stb )
          DEALLOCATE( r_tbtb )
          DEALLOCATE( z_tbtb )
          DEALLOCATE( p_tbtb )
          DEALLOCATE( r_pbtb )
          DEALLOCATE( z_pbtb )
          DEALLOCATE( p_pbtb )
          DEALLOCATE( gstb_a )
          DEALLOCATE( gspb_a )
          DEALLOCATE( gstb_tb_a )
          DEALLOCATE( gspb_tb_a )
          DEALLOCATE( r_spb )
          DEALLOCATE( z_spb )
          DEALLOCATE( p_spb )
          DEALLOCATE( r_pbpb )
          DEALLOCATE( z_pbpb )
          DEALLOCATE( p_pbpb )
          DEALLOCATE( gstb_pb_a )
          DEALLOCATE( gspb_pb_a )
          !! End Modifications by Andreas F. Martitsch (11.03.2014)

          !*************************************************************
          ! Provide curr_tor, curr_tor_s, curr_pol, curr_pol_s, iota
          !*************************************************************
          IF (write_progress .EQ. 1) THEN
             PRINT *, 'Prep of currents: ',s
          END IF
          swd = 1 ! derivative
          CALL splint_horner3(es,                                      &
               a_curr_tor, b_curr_tor, c_curr_tor, d_curr_tor,         &
               swd, m0,                                                &
               s, tfone, tfzero, tfzero, tfzero,                       &
               curr_tor_array(k_es), curr_tor_s_array(k_es), ypp, yppp)
          swd = 1 ! derivative
          CALL splint_horner3(es,                                      &
               a_curr_pol, b_curr_pol, c_curr_pol, d_curr_pol,         &
               swd, m0,                                                &
               s, tfone, tfzero, tfzero, tfzero,                       &
               curr_pol_array(k_es), curr_pol_s_array(k_es) ,ypp, yppp)
          swd = 0 ! no derivative
          CALL splint_horner3(es,                                      &
               a_iota, b_iota, c_iota, d_iota, swd, m0,                &
               s, tfone, tfzero, tfzero, tfzero,                       &
               iota_array(k_es), yp, ypp, yppp)
          CALL splint_horner3(es,                                      &
               a_pprime, b_pprime, c_pprime, d_pprime, swd, m0,                &
               s, tfone, tfzero, tfzero, tfzero,                       &
               pprime_array(k_es), yp, ypp, yppp)
          CALL splint_horner3(es,                                      &
               a_sqrtg00, b_sqrtg00, c_sqrtg00, d_sqrtg00, swd, m0,                &
               s, tfone, tfzero, tfzero, tfzero,                       &
               sqrtg00_array(k_es), yp, ypp, yppp)
       END DO
       magfie_newspline = 0
    END IF

    s_detected = 0
    IF (magfie_spline .EQ. 1) THEN
       s = x(1)
       !****************************************************************
       ! Detection of index
       !****************************************************************
       dxm1 = 1.0d0/(magfie_sarray(2)-magfie_sarray(1))
       CALL indef(s,magfie_sarray(1),dxm1,magfie_sarray_len,indu)
       xp = magfie_sarray(indu)
       s_detected = 1
       k_es = indu(1)
!!$       DO k_es = 1, magfie_sarray_len
!!$          !IF ( ABS(s-magfie_sarray(k_es)) .LT. magfie_epsi) THEN
!!$          IF ( ABS(s-magfie_sarray(k_es)) .LT. (magfie_sarray(2)-magfie_sarray(1))) THEN
!!$             s_detected = 1
!!$             EXIT
!!$          END IF
!!$       END DO
       IF (s_detected .EQ. 1) THEN
          !PRINT *,magfie_sarray(k_es)
          !STOP "flux surface s detected"
          fp = curr_tor_array(indu)
          CALL plag1d(s,fp,dxm1,xp,fun,der)
          curr_tor   = fun
          curr_tor_s = der
          fp = curr_pol_array(indu)
          CALL plag1d(s,fp,dxm1,xp,fun,der)
          curr_pol   = fun
          curr_pol_s = der
          fp = iota_array(indu)
          CALL plag1d(s,fp,dxm1,xp,fun,der)
          iota       = fun
          fp = pprime_array(indu)
          CALL plag1d(s,fp,dxm1,xp,fun,der)
          s_pprime   = fun ! only local
          fp = sqrtg00_array(indu)
          CALL plag1d(s,fp,dxm1,xp,fun,der)
          s_sqrtg00  = fun
          ! ************************************************************
          ! Evaluation of 2d-splines
          ! ************************************************************
          CALL poi2d(theta_int,phi_int,mt,mp,                          &
               theta_start,theta_end,phi_start,phi_end,                &
               x(3),x(2),theta_ind,phi_ind,theta_d,phi_d,ierr)
          !
          DO ind1=1,4
             p_spl => bmod_spl(:,:,:,:,indu(ind1))
             CALL eva2d(theta_n,phi_n,theta_ind,phi_ind,theta_d,phi_d,    &
                  p_spl,bmod)
             fp(ind1) = bmod
          END DO
          CALL plag1d(s,fp,dxm1,xp,fun,der)
          bmod = fun
          bb_s = der
          !PRINT *, bmod, bb_s
          !PAUSE
          !
          IF (isw_eval_spl2d_der .NE. 1) THEN
             DO ind1=1,4
                p_spl => bb_s_spl(:,:,:,:,indu(ind1))
                CALL eva2d(theta_n,phi_n,theta_ind,phi_ind,theta_d,phi_d,    &
                     p_spl,bb_s)
                fp(ind1) = bb_s
             END DO
             CALL plag1d(s,fp,dxm1,xp,fun,der)
             bb_s = fun
             !
             DO ind1=1,4
                p_spl => bb_tb_spl(:,:,:,:,indu(ind1))
                CALL eva2d(theta_n,phi_n,theta_ind,phi_ind,theta_d,phi_d,    &
                     p_spl,bb_tb)
                fp(ind1) = bb_tb
             END DO
             CALL plag1d(s,fp,dxm1,xp,fun,der)
             bb_tb = fun
             !
             DO ind1=1,4
                p_spl => bb_pb_spl(:,:,:,:,indu(ind1))
                CALL eva2d(theta_n,phi_n,theta_ind,phi_ind,theta_d,phi_d,    &
                     p_spl,bb_pb)
                fp(ind1) = bb_pb
             END DO
             CALL plag1d(s,fp,dxm1,xp,fun,der)
             bb_pb = fun
          ELSE
             !
             DO ind1=1,4
                p_spl => bmod_spl(:,:,:,:,indu(ind1))
                CALL eva2d_fd(theta_n,phi_n,theta_ind,phi_ind,theta_d,phi_d,    &
                     p_spl,bder_thph)
                fp(ind1) = bder_thph(1)
                fp1(ind1) = bder_thph(2)
             END DO
             CALL plag1d(s,fp,dxm1,xp,fun,der)
             bb_tb = fun
             CALL plag1d(s,fp1,dxm1,xp,fun,der)
             bb_pb = fun
             !
          END IF
          !
          DO ind1=1,4
             p_spl => gval_spl(:,:,:,:,indu(ind1))
             CALL eva2d(theta_n,phi_n,theta_ind,phi_ind,theta_d,phi_d,    &
                  p_spl,sqrg11)
             fp(ind1) = sqrg11
          END DO
          CALL plag1d(s,fp,dxm1,xp,fun,der)
          sqrg11 = fun
          !
          !! Modifications by Andreas F. Martitsch (11.03.2014)
          ! Evaluate the 2d periodic splines (over the flux-surface)
          ! of the additionally needed metric tensor elements
          !PRINT *,'x: ',x(1),x(2),x(3)
          IF (isw_eval_bcovars .EQ. 1) THEN
             !
             DO ind1=1,4
                p_spl => gstb_spl(:,:,:,:,indu(ind1))
                CALL eva2d(theta_n,phi_n,theta_ind,phi_ind,theta_d,phi_d,    &
                     p_spl,gstb)
                fp(ind1) = gstb
             END DO
             CALL plag1d(s,fp,dxm1,xp,fun,der)
             gstb = fun
             !PRINT *,'gstb: ', gstb
             !
             DO ind1=1,4
                p_spl => gspb_spl(:,:,:,:,indu(ind1))
                CALL eva2d(theta_n,phi_n,theta_ind,phi_ind,theta_d,phi_d,    &
                     p_spl,gspb)
                fp(ind1) = gspb
             END DO
             CALL plag1d(s,fp,dxm1,xp,fun,der)
             gspb = fun
             !PRINT *,'gspb: ', gspb
             !
             DO ind1=1,4
                p_spl => gstb_tb_spl(:,:,:,:,indu(ind1))
                CALL eva2d(theta_n,phi_n,theta_ind,phi_ind,theta_d,phi_d,    &
                     p_spl,gstb_tb)
                fp(ind1) = gstb_tb
             END DO
             CALL plag1d(s,fp,dxm1,xp,fun,der)
             gstb_tb = fun
             !PRINT *,'gstb_tb: ', gstb_tb
             !
             DO ind1=1,4
                p_spl => gspb_tb_spl(:,:,:,:,indu(ind1))
                CALL eva2d(theta_n,phi_n,theta_ind,phi_ind,theta_d,phi_d,    &
                     p_spl,gspb_tb)
                fp(ind1) = gspb_tb
             END DO
             CALL plag1d(s,fp,dxm1,xp,fun,der)
             gspb_tb = fun
             !PRINT *,'gspb_tb: ', gspb_tb
             !
             DO ind1=1,4
                p_spl => gstb_pb_spl(:,:,:,:,indu(ind1))
                CALL eva2d(theta_n,phi_n,theta_ind,phi_ind,theta_d,phi_d,    &
                     p_spl,gstb_pb)
                fp(ind1) = gstb_pb
             END DO
             CALL plag1d(s,fp,dxm1,xp,fun,der)
             gstb_pb = fun
             !PRINT *,'gstb_pb: ', gstb_pb
             !
             DO ind1=1,4
                p_spl => gspb_pb_spl(:,:,:,:,indu(ind1))
                CALL eva2d(theta_n,phi_n,theta_ind,phi_ind,theta_d,phi_d,    &
                     p_spl,gspb_pb)
                fp(ind1) = gspb_pb
             END DO
             CALL plag1d(s,fp,dxm1,xp,fun,der)
             gspb_pb = fun
             !PRINT *,'gspb_pb: ', gspb_pb
             !
          END IF
          !! End Modifications by Andreas F. Martitsch (11.03.2014)
!!$          !! Modifications by Andreas F. Martitsch (13.11.2014)
!!$          ! Evaluate the 2d periodic splines (over the flux-surface)
!!$          ! of the additionally needed quantities for NTV output
!!$          p_spl => R_spl(:,:,:,:,k_es)
!!$          CALL eva2d(theta_n,phi_n,theta_ind,phi_ind,theta_d,phi_d,    &
!!$               p_spl,r_val)
!!$          p_spl => Z_spl(:,:,:,:,k_es)
!!$          CALL eva2d(theta_n,phi_n,theta_ind,phi_ind,theta_d,phi_d,    &
!!$               p_spl,z_val)
!!$          !! End Modifications by Andreas F. Martitsch (13.11.2014)

          ! $1/sqrt(g)$
          fac = curr_pol + iota * curr_tor  ! (J + iota I)
          isqrg  = bmod*bmod / fac

          ! Winny for Klaus
          !s_sqrtg00_m = fac / av_b2_m
          s_sqrtg00 = fac / av_b2_m
          !PRINT *, 's_sqrtg00, s_sqrtg00_m ',s_sqrtg00, s_sqrtg00_m
          !PRINT *, 'fac, av_b2_m ',fac, av_b2_m
          !PAUSE
          ! Winny for Klaus - Ende

          !PRINT *, ' '
          !PRINT *, 'curr_pol = ',curr_pol
          !PRINT *, 'curr_tor = ',curr_tor
          !PRINT *, 'iota     = ',iota
          !PRINT *, 'fac      = ',fac
          !PRINT *, 'bmod     = ',bmod
          !PRINT *, 'isqrg    = ',isqrg
          !PRINT *, 'sqrg     = ',1.d0 / isqrg

          !PRINT *, 'sqrg11_n = ',sqrg11
          sqrg11 = sqrg11 * ABS(isqrg)
          !PRINT *, 'sqrg11   = ',sqrg11

          !! Modifications by Andreas F. Martitsch (11.03.2014)
          ! Compute the values of the additionally needed
          ! B-field components
          IF (isw_eval_bcovars .EQ. 1) THEN
             bcovar_s = isqrg*(gstb*iota+gspb)
             dbcovar_s_dtheta = (2.0d0*bmod*bb_tb/fac)*(gstb*iota+gspb) + &
                  isqrg*(gstb_tb*iota+gspb_tb)
             !PRINT *,'dbcovar_s_dtheta: ', dbcovar_s_dtheta
             dbcovar_s_dphi = (2.0d0*bmod*bb_pb/fac)*(gstb*iota+gspb) + &
                  isqrg*(gstb_pb*iota+gspb_pb)
             !PRINT *, 'dbcovar_s_dphi: ', dbcovar_s_dphi
             !STOP
          END IF
          !! End Modifications by Andreas F. Martitsch (11.03.2014)
       ELSE
          PRINT *, 'neo_magfie: s not detected!'
          STOP
       END IF
    END IF

    IF (magfie_spline .EQ. 0 .OR. s_detected .EQ. 0) THEN
       IF (magfie_spline .EQ. 1 .AND. s_detected .EQ. 0) THEN
          PRINT *, 'WARNING from neo_magfie - s out of range: ',s
          PRINT *, ' Using Fourier Summation directly'
       END IF

       PRINT *, 'magfie_spline .EQ. 0 not implemented'
       STOP

       !****************************************************************
       ! Direct summation of Fourier components
       !****************************************************************
       bmod   = 0.0_dp
       bb_s   = 0.0_dp
       bb_tb  = 0.0_dp
       bb_pb  = 0.0_dp

       !stop "Compute magnetic field components via direct Fourier summation"

       DO i = 1, mnmax
          swd = 1
          CALL splint_horner3(es,                                    &
               a_bmnc(:,i), b_bmnc(:,i), c_bmnc(:,i), d_bmnc(:,i),   &
               swd, r_mhalf(i),                                      &
               x(1), tf, tfp, tfpp, tfppp,                           &
               bmnc, bmnc_s, ypp, yppp)

          m = ixm(i)
          n = ixn(i)
          sinv = SIN(m*x(3) - n*x(2))
          cosv = COS(m*x(3) - n*x(2))

          bmod   = bmod   +     bmnc   * cosv
          bb_s   = bb_s   +     bmnc_s * cosv
          bb_tb  = bb_tb  - m * bmnc   * sinv
          bb_pb  = bb_pb  + n * bmnc   * sinv
       END DO

       swd = 1
       CALL splint_horner3(es,                                       &
            a_curr_tor, b_curr_tor, c_curr_tor, d_curr_tor,          &
            swd, m0,                                                 &
            x(1), tfone, tfzero, tfzero, tfzero,                     &
            curr_tor, curr_tor_s, ypp, yppp)
       CALL splint_horner3(es,                                       &
            a_curr_pol, b_curr_pol, c_curr_pol, d_curr_pol,          &
            swd, m0,                                                 &
            x(1), tfone, tfzero, tfzero, tfzero,                     &
            curr_pol, curr_pol_s ,ypp, yppp)
       swd = 0 ! no derivative
       CALL splint_horner3(es,                                       &
            a_iota, b_iota, c_iota, d_iota, swd, m0,                 &
            x(1), tfone, tfzero, tfzero, tfzero,                     &
            iota, yp, ypp, yppp)
    END IF

    IF (magfie_result .EQ. 1) THEN
       ! This was the original version:
       ! derived quantities in (s,theta_b,phi_b)-system
       fac   = (curr_pol + iota * curr_tor) * psi_pr
       fac1  = fac  / bmod                 ! sqrtg*bmod
       sqrtg = fac1 / bmod

       bder(1) = bb_s
       bder(2) = bb_tb
       bder(3) = bb_pb

       hcovar(1) = 0.0_dp
       hcovar(2) = curr_tor / bmod
       hcovar(3) = curr_pol / bmod

       hctrvr(1) = 0.0_dp
       hctrvr(2) = iota / fac1
       hctrvr(3) = 1.0_dp / fac1

       hcurl(1)  = (curr_pol * bb_pb      - curr_tor * bb_tb     ) / fac
       hcurl(2)  = (curr_pol * bb_s       - bmod     * curr_pol_s) / fac
       hcurl(3)  = (bmod     * curr_tor_s - curr_tor * bb_s      ) / fac
       ! Remark by Winny:
       ! The consisteny check for curr_pol shows a problem in all
       ! Greifswald (standard) input files
       ! According to the consistency check,
       ! curr_pol has to be changed to -curr_pol

    ELSEIF ( magfie_result .EQ. 0 ) THEN
       ! Modifications made by Sergie for use in SMT
       ! derived quantities in (s,theta_b,phi_b)-system
       !fac   = (curr_pol + iota * curr_tor) * psi_pr
       ! This is used in NEO2
       fac   =  curr_pol + iota * curr_tor                       !!!
       fac1  = fac  / bmod                 ! sqrtg*bmod
       fac = fac * psi_pr                                        !!!
       !    sqrtg = fac1 / bmod
       sqrtg = - fac1 / bmod * psi_pr * 1d6                      !!!
       !---------------------------------------------------------------------------
       !  iota_m = iota
       ! fac_m  =  (curr_pol + iota * curr_tor) * 1d6 * psi_pr
       !  fac_c  =  (curr_pol + iota * curr_tor) * 1d6
       !---------------------------------------------------------------------------

       bder(1) = bb_s
       bder(3) = bb_tb
       bder(2) = bb_pb
       bder=bder / bmod                                          !!!

       !! Modifications by Andreas F. Martitsch (07.03.2014)
       ! Radial covariant B-field component is now available
       IF (isw_eval_bcovars .EQ. 1) THEN
          hcovar(1) = bcovar_s / bmod
       ELSE
          hcovar(1) = 0.0_dp
       END IF
       !! End Modifications by Andreas F. Martitsch (07.03.2014)
       hcovar(3) = curr_tor / bmod
       hcovar(2) = curr_pol / bmod
       hcovar=hcovar * 1.d2                                      !!!

       hctrvr(1) = 0.0_dp
       hctrvr(3) = iota / fac1
       hctrvr(2) = 1.0_dp / fac1
       hctrvr=hctrvr * 1d-2                                      !!!

       !    hcurl(1)  = (curr_pol * bb_pb      - curr_tor * bb_tb     ) / fac
       hcurl(1)  = (curr_tor * bb_pb      - curr_pol * bb_tb     ) / fac  !!!
       hcurl(3)  = (curr_pol * bb_s       - bmod     * curr_pol_s) / fac
       hcurl(2)  = (bmod     * curr_tor_s - curr_tor * bb_s      ) / fac
       hcurl=hcurl * 1d-4                                                 !!!

       !! Modifications by Andreas F. Martitsch (12.03.2014)
       ! boozer_curr_tor, boozer_curr_pol, boozer_psi_pr,
       ! boozer_sqrtg11 and boozer_isqrg are now converted
       ! to cgs-units.
       ! This step requires changes within rhs_kin.f90 and
       ! ripple_solver.f90!
       IF (bmod0 .EQ. 0.0d0) THEN
          !PRINT *,bmod0
          boozer_curr_tor_hat=0.0d0
          boozer_curr_pol_hat=0.0d0
          boozer_curr_tor_hat_s=0.0d0
          boozer_curr_pol_hat_s=0.0d0
          boozer_psi_pr_hat=0.0d0
       ELSE
          !PRINT *,bmod0
          boozer_curr_tor_hat = (curr_tor/bmod0)*1.0d2
          boozer_curr_pol_hat = (curr_pol/bmod0)*1.0d2
          !! Modifications by Andreas F. Martitsch (24.04.2015)
          ! The following copy-paste error has been detected by
          ! Christopher and affected the versions before 24.04.2015
          ! (it is used within ripple_solver for the computation
          ! of the magnetic drift frequency)
          !boozer_curr_pol_hat_s = (curr_tor_s/bmod0)*1.0d2
          !boozer_curr_tor_hat_s = (curr_pol_s/bmod0)*1.0d2
          ! These are now the correct quantities:
          boozer_curr_tor_hat_s = (curr_tor_s/bmod0)*1.0d2
          boozer_curr_pol_hat_s = (curr_pol_s/bmod0)*1.0d2
          !! End Modifications by Andreas F. Martitsch (24.04.2015)
          boozer_psi_pr_hat = (psi_pr/bmod0)*1.0d4
       END IF
       boozer_sqrtg11 = (1.0d0/psi_pr)*sqrg11*1.0d-2
       boozer_isqrg = (1.0d0/psi_pr)*isqrg*1.0d-6
       !! End Modifications by Andreas F. Martitsch (12.03.2014)

    END IF

    boozer_iota = iota
    ! CAUTION: This quantity is only used by Klaus.
    ! Conversion from SI- to cgs-units has not yet been
    ! checked for this quantity
    boozer_sqrtg00 = s_sqrtg00
!!$    !! Modifications by Andreas F. Martitsch (12.03.2014)
!!$    ! Note! Every quantity given below is written in
!!$    ! terms of SI-units. This part has been moved to the
!!$    ! block with "magfie_result eq. 0), where the quantities are
!!$    ! converted to cgs-units. In the previous version
!!$    ! conversions were done partly within rhs_kin.f90 and
!!$    ! partly within ripple_solver.f90.
!!$    ! For this reason, please note also the changes in
!!$    ! "rhs_kin" and "ripple_solver"!
!!$    boozer_curr_tor = curr_tor
!!$    boozer_curr_pol = curr_pol
!!$    ! Compute the radial derivatives of toroidal / poloidal currents
!!$    ! (In fact these currents are already the respective
!!$    ! covariant B-field components; conversion done within
!!$    ! neo_read)
!!$    ! Caution! boozer_curr_tor/ boozer_curr_pol are given in SI-units [Tm].
!!$    ! They are not converted to cgs-units - as done within (magfie_result .EQ. 0).
!!$    boozer_curr_pol_s = curr_tor_s
!!$    boozer_curr_tor_s = curr_pol_s
!!$    boozer_psi_pr = psi_pr
!!$    boozer_sqrtg11 = sqrg11
!!$    boozer_isqrg = isqrg
!!$    !! End Modifications by Andreas F. Martitsch (12.03.2014)
  END SUBROUTINE neo_magfie_a

!!$  !! Modifications by Andreas F. Martitsch (11.03.2014)
!!$  ! Optional output (necessary for modeling the magnetic rotation)
!!$  SUBROUTINE neo_magfie_b( x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl, bcovar_s_hat_der )
!!$    ! input / output
!!$    REAL(dp), DIMENSION(:),       INTENT(in)  :: x
!!$    REAL(dp),                     INTENT(out) :: bmod
!!$    REAL(dp),                     INTENT(out) :: sqrtg
!!$    REAL(dp), DIMENSION(SIZE(x)), INTENT(out) :: bder
!!$    REAL(dp), DIMENSION(SIZE(x)), INTENT(out) :: hcovar
!!$    REAL(dp), DIMENSION(SIZE(x)), INTENT(out) :: hctrvr
!!$    REAL(dp), DIMENSION(SIZE(x)), INTENT(out) :: hcurl
!!$    REAL(dp), DIMENSION(SIZE(x)), INTENT(out) :: bcovar_s_hat_der
!!$    !
!!$    CALL neo_magfie_a( x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl )
!!$    !
!!$    ! Compute the derivatives of the radial covariant
!!$    ! B-field component (Note: cgs-units used)
!!$    bcovar_s_hat_der(1) = 0.0_dp ! not implemented at the moment
!!$    bcovar_s_hat_der(3) = dbcovar_s_dtheta / bmod0
!!$    bcovar_s_hat_der(2) = dbcovar_s_dphi / bmod0
!!$    bcovar_s_hat_der=bcovar_s_hat_der * 1.d2 ! conversion to cgs-units
!!$    !
!!$  END SUBROUTINE neo_magfie_b
!!$  !! End Modifications by Andreas F. Martitsch (11.03.2014)

!!$  !! Modifications by Andreas F. Martitsch (13.11.2014)
!!$  ! Optional output for NTV output
!!$  SUBROUTINE neo_magfie_c( x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl, bcovar_s_hat_der, R, Z )
!!$    ! input / output
!!$    REAL(dp), DIMENSION(:),       INTENT(in)  :: x
!!$    REAL(dp),                     INTENT(out) :: bmod
!!$    REAL(dp),                     INTENT(out) :: sqrtg
!!$    REAL(dp), DIMENSION(SIZE(x)), INTENT(out) :: bder
!!$    REAL(dp), DIMENSION(SIZE(x)), INTENT(out) :: hcovar
!!$    REAL(dp), DIMENSION(SIZE(x)), INTENT(out) :: hctrvr
!!$    REAL(dp), DIMENSION(SIZE(x)), INTENT(out) :: hcurl
!!$    REAL(dp), DIMENSION(SIZE(x)), INTENT(out) :: bcovar_s_hat_der
!!$    REAL(dp),                     INTENT(out) :: R
!!$    REAL(dp),                     INTENT(out) :: Z
!!$    !
!!$    CALL neo_magfie_b( x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl, bcovar_s_hat_der )
!!$    !
!!$    R = r_val * 1.d2 ! conversion to cgs-units
!!$    Z = z_val * 1.d2 ! conversion to cgs-units
!!$    !
!!$  END SUBROUTINE neo_magfie_c
!!$  !! End Modifications by Andreas F. Martitsch (13.11.2014)

  SUBROUTINE cyl_coord_a( x, x_cyl, jacobian )
    ! input / output
    REAL(dp), DIMENSION(:),       INTENT(in)            :: x
    REAL(dp), DIMENSION(SIZE(x)), INTENT(out)           :: x_cyl
    REAL(dp), DIMENSION(SIZE(x),SIZE(x)), INTENT(out)   :: jacobian
    ! local definitions
    INTEGER(i4b)                                     :: swd = 1
    INTEGER                                          :: i, m, n
    INTEGER                                          :: npsi
    REAL(dp)                                         :: m0  = 0.0_dp
    REAL(dp)                                         :: yp, ypp, yppp

    REAL(dp)                                         :: sinv, cosv

    INTEGER                                          :: k_es
    INTEGER                                          :: s_detected
    INTEGER                                          :: imn
    INTEGER                                          :: it, ip, im, in
    INTEGER                                          :: mt = 1
    INTEGER                                          :: mp = 1
    INTEGER                                          :: theta_ind, phi_ind
    INTEGER                                          :: ierr
    REAL(dp)                                         :: s
    REAL(dp)                                         :: magfie_epsi = 1.e-9
    REAL(dp)                                         :: ri, zi, li
    !! Modifications by Andreas F. Martitsch (07.03.2014)
    ! Auxiliary variables for the Fourier summation
    REAL(dp)                                         :: ri_s, zi_s, li_s
    REAL(dp)                                         :: ris, zis, lis
    REAL(dp)                                         :: ris_s, zis_s, lis_s
    !! End Modifications by Andreas F. Martitsch (07.03.2014)
    REAL(dp)                                         :: theta_d, phi_d

    REAL(dp), DIMENSION(:), ALLOCATABLE              :: s_rmnc, s_zmnc, s_lmnc
    !! Modifications by Andreas F. Martitsch (06.03.2014)
    ! Radial derivatives of (R,Z,phi)-components obtained from the 1d spline
    REAL(dp), DIMENSION(:), ALLOCATABLE              :: s_rmnc_s, s_zmnc_s, s_lmnc_s
    !! End Modifications by Andreas F. Martitsch (06.03.2014)
    !! Modifications by Andreas F. Martitsch (06.08.2014)
    ! Additional data from Boozer files without Stellarator symmetry
    REAL(kind=dp),    DIMENSION(:),     ALLOCATABLE :: s_rmns, s_zmns, s_lmns
    REAL(kind=dp),    DIMENSION(:),     ALLOCATABLE :: s_rmns_s, s_zmns_s, s_lmns_s
    !! End Modifications by Andreas F. Martitsch (06.08.2014)
    REAL(dp), DIMENSION(:,:), ALLOCATABLE            :: r,z,l
    REAL(dp), DIMENSION(:,:), ALLOCATABLE            :: r_tb,z_tb,p_tb
    REAL(dp), DIMENSION(:,:), ALLOCATABLE            :: r_pb,z_pb,p_pb
    !! Modifications by Andreas F. Martitsch (11.03.2014)
    ! Temporary storage arrays for the Fourier summations related to
    ! the radial derivatives of (R,Z,phi)-components and the
    ! additionally needed metric tensor elements
    REAL(dp), DIMENSION(:,:), ALLOCATABLE            :: r_s,z_s,p_s
    !! End Modifications by Andreas F. Martitsch (11.03.2014)

    REAL(dp), DIMENSION(:,:,:,:), POINTER            :: p_spl

    INTEGER :: ind1
    INTEGER, DIMENSION(mp_lag) :: indu
    REAL(dp), DIMENSION(mp_lag) :: xp, fp, fp1
    REAL(dp) :: fun, der, dxm1
    REAL(dp), DIMENSION(2) :: der_thph

    !*******************************************************************
    ! Initialisation if necessary
    !*******************************************************************
    IF ( .NOT. ALLOCATED(es) ) THEN
       CALL neo_read_control()
       fluxs_interp = 1
       CALL neo_init(npsi)
       PRINT *, 'theta_start,theta_end,phi_start,phi_end'
       PRINT *, theta_start,theta_end,phi_start,phi_end
       !! Modifications by Andreas F. Martitsch (27.01.2016)
       ! -> switch on/off use of splined Fourier coefficients within neo_magfie
       IF (isw_spl_fourier_cof .NE. 1) THEN
          ALLOCATE(magfie_sarray(SIZE(es)/flux_surf_dist))
          magfie_sarray = es(1:SIZE(es):flux_surf_dist)
       END IF
       !! End Modifications by Andreas F. Martitsch (27.01.2016)
    END IF
    !*******************************************************************
    ! Spline of surfaces in magfie_sarray
    !*******************************************************************
    IF (magfie_spline .EQ. 1 .AND. magfie_newspline_cyl .EQ. 1) THEN
       magfie_sarray_len =  SIZE(magfie_sarray)
       !****************************************************************
       ! Allocation
       !****************************************************************
       !! Modifications by Andreas F. Martitsch (13.11.2014)
       ! Allocate storage arrays for the 2d periodic splines
       ! of the additionally needed quantities for NTV output
       ALLOCATE( R_spl(4,4,theta_n,phi_n,magfie_sarray_len) )
       ALLOCATE( Z_spl(4,4,theta_n,phi_n,magfie_sarray_len) )
       ALLOCATE( L_spl(4,4,theta_n,phi_n,magfie_sarray_len) )
       !! End Modifications by Andreas F. Martitsch (13.11.2014)
       !****************************************************************
       ! Loop over predefined s-values
       !****************************************************************
       DO k_es = 1, magfie_sarray_len
          s = magfie_sarray(k_es)
          !*************************************************************
          ! Surface
          !*************************************************************
          IF (write_progress .EQ. 1) THEN
             PRINT *, 'Initialize Surface, k_es = ',k_es
          END IF
          ALLOCATE ( s_rmnc(mnmax) )
          ALLOCATE ( s_zmnc(mnmax) )
          ALLOCATE ( s_lmnc(mnmax) )
          !! Modifications by Andreas F. Martitsch (06.03.2014)
          ! Compute the necessary radial derivatives for the
          ! (R,Z,phi)-components obtained from the 1d spline
          ALLOCATE ( s_rmnc_s(mnmax) ) ! Allocate arrays for additional
          ALLOCATE ( s_zmnc_s(mnmax) ) ! radial derivatives
          ALLOCATE ( s_lmnc_s(mnmax) )
          !! End Modifications by Andreas F. Martitsch (06.03.2014)
          !
          !! Modifications by Andreas F. Martitsch (06.08.2014)
          ! Additional data from Boozer files without Stellarator symmetry
          IF (inp_swi .EQ. 9) THEN        ! ASDEX-U (E. Strumberger)
             ALLOCATE ( s_rmns(mnmax) )
             ALLOCATE ( s_zmns(mnmax) )
             ALLOCATE ( s_lmns(mnmax) )
             ALLOCATE ( s_rmns_s(mnmax) )
             ALLOCATE ( s_zmns_s(mnmax) )
             ALLOCATE ( s_lmns_s(mnmax) )
          END IF
          !! End Modifications by Andreas F. Martitsch (06.08.2014)
          !
          IF (isw_spl_fourier_cof .EQ. 1) THEN
             DO imn = 1, mnmax
                ! Switch swd turns on (1) / off (0) the computation of the
                ! radial derivatives within splint_horner3
                !swd = 0 ! Now we want to compute the radial derivatives
                swd = 1
                CALL splint_horner3(es,                                   &
                     a_rmnc(:,imn), b_rmnc(:,imn),                        &
                     c_rmnc(:,imn), d_rmnc(:,imn),                        &
                     swd, r_mhalf(imn),                                   &
                     s, tf, tfp, tfpp, tfppp,                             &
                     s_rmnc(imn), s_rmnc_s(imn), ypp, yppp)
                swd = 1
                CALL splint_horner3(es,                                   &
                     a_zmnc(:,imn), b_zmnc(:,imn),                        &
                     c_zmnc(:,imn), d_zmnc(:,imn),                        &
                     swd, r_mhalf(imn),                                   &
                     s, tf, tfp, tfpp, tfppp,                             &
                     s_zmnc(imn), s_zmnc_s(imn), ypp, yppp)
                swd = 1
                CALL splint_horner3(es,                                   &
                     a_lmnc(:,imn), b_lmnc(:,imn),                        &
                     c_lmnc(:,imn), d_lmnc(:,imn),                        &
                     swd, r_mhalf(imn),                                   &
                     s, tf, tfp, tfpp, tfppp,                             &
                     s_lmnc(imn), s_lmnc_s(imn), ypp, yppp)
                !
                !! Modifications by Andreas F. Martitsch (06.08.2014)
                ! Additional data from Boozer files without Stellarator symmetry
                IF (inp_swi .EQ. 9) THEN        ! ASDEX-U (E. Strumberger)
                   swd = 1
                   CALL splint_horner3(es,                                   &
                        a_rmns(:,imn), b_rmns(:,imn),                        &
                        c_rmns(:,imn), d_rmns(:,imn),                        &
                        swd, r_mhalf(imn),                                   &
                        s, tf, tfp, tfpp, tfppp,                             &
                        s_rmns(imn), s_rmns_s(imn), ypp, yppp)
                   swd = 1
                   CALL splint_horner3(es,                                   &
                        a_zmns(:,imn), b_zmns(:,imn),                        &
                        c_zmns(:,imn), d_zmns(:,imn),                        &
                        swd, r_mhalf(imn),                                   &
                        s, tf, tfp, tfpp, tfppp,                             &
                        s_zmns(imn), s_zmns_s(imn), ypp, yppp)
                   swd = 1
                   CALL splint_horner3(es,                                   &
                        a_lmns(:,imn), b_lmns(:,imn),                        &
                        c_lmns(:,imn), d_lmns(:,imn),                        &
                        swd, r_mhalf(imn),                                   &
                        s, tf, tfp, tfpp, tfppp,                             &
                        s_lmns(imn), s_lmns_s(imn), ypp, yppp)
                END IF
                !! End Modifications by Andreas F. Martitsch (06.08.2014)
                !
             END DO
          ELSE
             s_rmnc = a_rmnc(1+(k_es-1)*flux_surf_dist,:)
             s_rmnc_s = 0.0d0
             s_zmnc = a_zmnc(1+(k_es-1)*flux_surf_dist,:)
             s_zmnc_s = 0.0d0
             s_lmnc = a_lmnc(1+(k_es-1)*flux_surf_dist,:)
             s_lmnc_s = 0.0d0
             IF (inp_swi .EQ. 9) THEN        ! ASDEX-U (E. Strumberger)
                s_rmns = a_rmns(1+(k_es-1)*flux_surf_dist,:)
                s_rmns_s = 0.0d0
                s_zmns = a_zmns(1+(k_es-1)*flux_surf_dist,:)
                s_zmns_s = 0.0d0
                s_lmns = a_lmns(1+(k_es-1)*flux_surf_dist,:)
                s_lmns_s = 0.0d0
             END IF
          END IF
          !*************************************************************
          ! Fourier summation for the full theta-phi array
          !*************************************************************
          IF (write_progress .EQ. 1) THEN
             PRINT *, 'Do Fourier'
          END IF

          ALLOCATE( r(theta_n,phi_n) )  ! NEW
          ALLOCATE( z(theta_n,phi_n) )
          ALLOCATE( l(theta_n,phi_n) )
          ALLOCATE( r_tb(theta_n,phi_n) )
          ALLOCATE( z_tb(theta_n,phi_n) )
          ALLOCATE( p_tb(theta_n,phi_n) )
          ALLOCATE( r_pb(theta_n,phi_n) )
          ALLOCATE( z_pb(theta_n,phi_n) )
          ALLOCATE( p_pb(theta_n,phi_n) )
          r = 0.0d0
          z = 0.0d0
          l = 0.0d0
          r_tb = 0.0d0
          z_tb = 0.0d0
          p_tb = 0.0d0
          r_pb = 0.0d0
          z_pb = 0.0d0
          p_pb = 0.0d0
          !! Modifications by Andreas F. Martitsch (11.03.2014)
          ! Allocate temporary storage arrays for the Fourier summations
          ! related to the radial derivatives of (R,Z,phi)-components
          ALLOCATE( r_s(theta_n,phi_n) )
          ALLOCATE( z_s(theta_n,phi_n) )
          ALLOCATE( p_s(theta_n,phi_n) )
          r_s = 0.0d0
          z_s = 0.0d0
          p_s = 0.0d0
          !! End Modifications by Andreas F. Martitsch (11.03.2014)

          DO imn=1,mnmax
             ri = s_rmnc(imn) ! NEW
             zi = s_zmnc(imn)
             li = s_lmnc(imn)
             !! Modifications by Andreas F. Martitsch (07.03.2014)
             ! Auxiliary variables for the Fourier summation
             ri_s = s_rmnc_s(imn)
             zi_s = s_zmnc_s(imn)
             li_s = s_lmnc_s(imn)
             !! End Modifications by Andreas F. Martitsch (07.03.2014)
             !! Modifications by Andreas F. Martitsch (06.08.2014)
             ! Additional data from Boozer files without Stellarator symmetry
             IF (inp_swi .EQ. 9) THEN        ! ASDEX-U (E. Strumberger)
                !
                ris = s_rmns(imn)
                zis = s_zmns(imn)
                lis = s_lmns(imn)
                !
                ris_s = s_rmns_s(imn)
                zis_s = s_zmns_s(imn)
                lis_s = s_lmns_s(imn)
                !
             END IF
             !! End Modifications by Andreas F. Martitsch (06.08.2014)
             m = ixm(imn)
             n = ixn(imn)
             im = pixm(imn)
             in = pixn(imn)
             DO ip=1,phi_n
                DO it=1,theta_n
                   !! Modifications by Andreas F. Martitsch (06.08.2014)
                   ! Additional data from Boozer files without Stellarator symmetry
                   IF (inp_swi .EQ. 9) THEN        ! ASDEX-U (E. Strumberger)
                      cosv = cosmth(it,im) * cosnph(ip,in) - sinmth(it,im) * sinnph(ip,in)
                      sinv = sinmth(it,im) * cosnph(ip,in) + cosmth(it,im) * sinnph(ip,in)

                      r(it,ip) = r(it,ip) + ri*cosv + ris*sinv
                      z(it,ip) = z(it,ip) + zi*cosv + zis*sinv
                      l(it,ip) = l(it,ip) + li*cosv + lis*sinv

                      r_tb(it,ip) = r_tb(it,ip) - m*ri*sinv + m*ris*cosv
                      r_pb(it,ip) = r_pb(it,ip) - n*ri*sinv + n*ris*cosv
                      z_tb(it,ip) = z_tb(it,ip) - m*zi*sinv + m*zis*cosv
                      z_pb(it,ip) = z_pb(it,ip) - n*zi*sinv + n*zis*cosv
                      !! Modifications by Andreas F. Martitsch (12.11.2015)
                      !According to Erika Strumberger (Email 11.10.2015)
                      !the conversion from phi_b to phi is given by
                      !"\phi-phi_b = 2\pi/N_p \sum ( c \cos(2\pi (m u + n v) ) + s \sin(2\pi (m u+n v) ) )"
                      !where  \phi=2\pi/N_p v.
                      !This expression differs by a minus sign from the
                      !expression used by J. Geiger ( phi_b-\phi = ... )!
                      !-> previous versions used this definition:
                      !p_tb(it,ip) = p_tb(it,ip) + m*li*sinv - m*lis*cosv ! -l_tb
                      !p_pb(it,ip) = p_pb(it,ip) + n*li*sinv - n*lis*cosv ! -l_pb
                      !-> corrected formulas:
                      p_tb(it,ip) = p_tb(it,ip) - m*li*sinv + m*lis*cosv ! +l_tb
                      p_pb(it,ip) = p_pb(it,ip) - n*li*sinv + n*lis*cosv ! +l_pb
                      !! End Modifications by Andreas F. Martitsch (12.11.2015)

                      r_s(it,ip) = r_s(it,ip) + ri_s*cosv + ris_s*sinv
                      z_s(it,ip) = z_s(it,ip) + zi_s*cosv + zis_s*sinv
                      !! Modifications by Andreas F. Martitsch (12.11.2015)
                      !According to Erika Strumberger (Email 11.10.2015)
                      !the conversion from phi_b to phi is given by
                      !"\phi-phi_b = 2\pi/N_p \sum ( c \cos(2\pi (m u + n v) ) + s \sin(2\pi (m u+n v) ) )"
                      !where  \phi=2\pi/N_p v.
                      !This expression differs by a minus sign from the
                      !expression used by J. Geiger ( phi_b-\phi = ... )!
                      !-> previous versions used this definition:
                      !p_s(it,ip) = p_s(it,ip) - li_s*cosv - lis_s*sinv ! -l_s
                      !-> corrected formulas:
                      p_s(it,ip) = p_s(it,ip) + li_s*cosv + lis_s*sinv ! +l_s
                      !! End Modifications by Andreas F. Martitsch (12.11.2015)
                   ELSE
                      cosv = cosmth(it,im) * cosnph(ip,in) + sinmth(it,im) * sinnph(ip,in)
                      sinv = sinmth(it,im) * cosnph(ip,in) - cosmth(it,im) * sinnph(ip,in)

                      r(it,ip) = r(it,ip) + ri*cosv
                      z(it,ip) = z(it,ip) + zi*sinv
                      l(it,ip) = l(it,ip) + li*sinv

                      r_tb(it,ip) = r_tb(it,ip) - m*ri*sinv
                      r_pb(it,ip) = r_pb(it,ip) + n*ri*sinv
                      z_tb(it,ip) = z_tb(it,ip) + m*zi*cosv
                      z_pb(it,ip) = z_pb(it,ip) - n*zi*cosv
                      p_tb(it,ip) = p_tb(it,ip) - m*li*cosv ! -l_tb
                      p_pb(it,ip) = p_pb(it,ip) + n*li*cosv ! -l_pb
                      !! Modifications by Andreas F. Martitsch (07.03.2014)
                      ! Temporary storage arrays for the Fourier summations
                      ! related to the radial derivatives of (R,Z,phi)-components
                      r_s(it,ip) = r_s(it,ip) + ri_s*cosv
                      z_s(it,ip) = z_s(it,ip) + zi_s*sinv
                      p_s(it,ip) = p_s(it,ip) - li_s*sinv ! -l_s
                      !! End Modifications by Andreas F. Martitsch (07.03.2014)
                   END IF
                   !! End Modifications by Andreas F. Martitsch (06.08.2014)
                END DO
             END DO
          END DO
          DEALLOCATE( s_rmnc )
          DEALLOCATE( s_zmnc )
          DEALLOCATE( s_lmnc )
          !! Modifications by Andreas F. Martitsch (07.03.2014)
          ! Deallocate arrays for the additional radial derivatives
          DEALLOCATE( s_rmnc_s )
          DEALLOCATE( s_zmnc_s )
          DEALLOCATE( s_lmnc_s )
          !! End Modifications by Andreas F. Martitsch (07.03.2014)
          !
          !! Modifications by Andreas F. Martitsch (06.08.2014)
          ! Additional data from Boozer files without Stellarator symmetry
          IF (inp_swi .EQ. 9) THEN        ! ASDEX-U (E. Strumberger)
             DEALLOCATE ( s_rmns )
             DEALLOCATE ( s_zmns )
             DEALLOCATE ( s_lmns )
             DEALLOCATE ( s_rmns_s )
             DEALLOCATE ( s_zmns_s )
             DEALLOCATE ( s_lmns_s )
          END IF
          !! End Modifications by Andreas F. Martitsch (06.08.2014)
          !
          IF (lab_swi .EQ. 5 .OR. lab_swi .EQ. 3) THEN ! CHS, LHD
             p_tb = - p_tb
             p_pb = 1 - p_pb
             !! Modifications by Andreas F. Martitsch (07.03.2014)
             ! ToDo: Implement conversion for p_s
             PRINT *,'WARNING FROM NEO_MAGFIE: CONVERSION FOR RADIAL DERIVATIVE OF BOOZER-PHI NOT IMPLEMENTED!'
             !! End Modifications by Andreas F. Martitsch (07.03.2014)
          ELSE
             p_tb = p_tb * twopi / nfp
             p_pb = 1.0_dp + p_pb * twopi / nfp
             !! Modifications by Andreas F. Martitsch (11.03.2014)
             ! Conversion factor between Boozer-phi and L (see Boozer-file)
             p_s = p_s * twopi / nfp
             !! End Modifications by Andreas F. Martitsch (11.03.2014)
          END IF

          ! **********************************************************************
          ! Ensure periodicity boundaries to be the same
          ! **********************************************************************
          r(theta_n,:) = r(1,:)
          r(:,phi_n)   = r(:,1)
          z(theta_n,:) = z(1,:)
          z(:,phi_n)   = z(:,1)
          l(theta_n,:) = l(1,:)
          l(:,phi_n)   = l(:,1)
          r_tb(theta_n,:) = r_tb(1,:)
          r_tb(:,phi_n)   = r_tb(:,1)
          r_pb(theta_n,:) = r_pb(1,:)
          r_pb(:,phi_n)   = r_pb(:,1)
          z_tb(theta_n,:) = z_tb(1,:)
          z_tb(:,phi_n)   = z_tb(:,1)
          z_pb(theta_n,:) = z_pb(1,:)
          z_pb(:,phi_n)   = z_pb(:,1)
          p_tb(theta_n,:) = p_tb(1,:)
          p_tb(:,phi_n)   = p_tb(:,1)
          p_pb(theta_n,:) = p_pb(1,:)
          p_pb(:,phi_n)   = p_pb(:,1)
          !! Modifications by Andreas F. Martitsch (07.03.2014)
          ! Temporary storage arrays for the Fourier summations
          ! related to the radial derivatives of (R,Z,phi)-components
          r_s(theta_n,:) = r_s(1,:)
          r_s(:,phi_n)   = r_s(:,1)
          z_s(theta_n,:) = z_s(1,:)
          z_s(:,phi_n)   = z_s(:,1)
          p_s(theta_n,:) = p_s(1,:)
          p_s(:,phi_n)   = p_s(:,1)
          !! End Modifications by Andreas F. Martitsch (07.03.2014)

          !*************************************************************
          ! Do the 2-D periodic spline
          !*************************************************************
          !! Modifications by Andreas F. Martitsch (13.11.2014)
          ! Compute the 2d periodic splines (over the flux-surface)
          ! of the additionally needed quantities for NTV output
          p_spl => R_spl(:,:,:,:,k_es)
          CALL spl2d(theta_n,phi_n,theta_int,phi_int,mt,mp,            &
               r,p_spl)
          p_spl => Z_spl(:,:,:,:,k_es)
          CALL spl2d(theta_n,phi_n,theta_int,phi_int,mt,mp,            &
               z,p_spl)
          p_spl => L_spl(:,:,:,:,k_es)
          CALL spl2d(theta_n,phi_n,theta_int,phi_int,mt,mp,            &
               l,p_spl)
          !! End Modifications by Andreas F. Martitsch (13.11.2014)

          !stop "Compute magnetic field components via a spline (and not via direct Fourier summation)"

          DEALLOCATE( r )  ! NEW
          DEALLOCATE( z )
          DEALLOCATE( l )
          DEALLOCATE( r_tb )
          DEALLOCATE( z_tb )
          DEALLOCATE( p_tb )
          DEALLOCATE( r_pb )
          DEALLOCATE( z_pb )
          DEALLOCATE( p_pb )

          !! Modifications by Andreas F. Martitsch (11.03.2014)
          ! Deallocate temporary storage arrays for the Fourier summations
          ! related to the radial derivatives of (R,Z,phi)-components
          ! and the additionally needed metric tensor elements
          DEALLOCATE( r_s )
          DEALLOCATE( z_s )
          DEALLOCATE( p_s )
          !! End Modifications by Andreas F. Martitsch (11.03.2014)
       END DO
       magfie_newspline_cyl = 0
    END IF

    s_detected = 0
    IF (magfie_spline .EQ. 1) THEN
       s = x(1)
       !****************************************************************
       ! Detection of index
       !****************************************************************
       dxm1 = 1.0d0/(magfie_sarray(2)-magfie_sarray(1))
       CALL indef(s,magfie_sarray(1),dxm1,magfie_sarray_len,indu)
       xp = magfie_sarray(indu)
       s_detected = 1
       k_es = indu(1)
!!$       DO k_es = 1, magfie_sarray_len
!!$          !IF ( ABS(s-magfie_sarray(k_es)) .LT. magfie_epsi) THEN
!!$          IF ( ABS(s-magfie_sarray(k_es)) .LT. (magfie_sarray(2)-magfie_sarray(1))) THEN
!!$             s_detected = 1
!!$             EXIT
!!$          END IF
!!$       END DO
       IF (s_detected .EQ. 1) THEN
          ! ************************************************************
          ! Evaluation of 2d-splines
          ! ************************************************************
          CALL poi2d(theta_int,phi_int,mt,mp,                          &
               theta_start,theta_end,phi_start,phi_end,                &
               x(3),x(2),theta_ind,phi_ind,theta_d,phi_d,ierr)
          !
          !! Modifications by Andreas F. Martitsch (13.11.2014)
          ! Evaluate the 2d periodic splines (over the flux-surface)
          ! of the additionally needed quantities for NTV output
          DO ind1=1,4
             p_spl => R_spl(:,:,:,:,indu(ind1))
             CALL eva2d(theta_n,phi_n,theta_ind,phi_ind,theta_d,phi_d,    &
                  p_spl,r_val)
             fp(ind1) = r_val
          END DO
          CALL plag1d(s,fp,dxm1,xp,fun,der)
          r_val = fun
          r_s_val = der
          !
          DO ind1=1,4
             p_spl => Z_spl(:,:,:,:,indu(ind1))
             CALL eva2d(theta_n,phi_n,theta_ind,phi_ind,theta_d,phi_d,    &
                  p_spl,z_val)
             fp(ind1) = z_val
          END DO
          CALL plag1d(s,fp,dxm1,xp,fun,der)
          z_val = fun
          z_s_val = der
          !
          DO ind1=1,4
             p_spl => L_spl(:,:,:,:,indu(ind1))
             CALL eva2d(theta_n,phi_n,theta_ind,phi_ind,theta_d,phi_d,    &
                  p_spl,p_val)
             fp(ind1) = p_val
          END DO
          CALL plag1d(s,fp,dxm1,xp,fun,der)
          p_val = fun
          p_s_val = der
          !
          IF (lab_swi .EQ. 5 .OR. lab_swi .EQ. 3) THEN ! CHS, LHD
             p_val = x(2) - p_val
             p_s_val = - p_s_val
             !! Modifications by Andreas F. Martitsch (28.01.2016)
             ! ToDo: Implement conversion for p_s
             PRINT *,'WARNING FROM NEO_MAGFIE: CONVERSION FOR RADIAL DERIVATIVE OF BOOZER-PHI NOT CHECKED!'
             !! End Modifications by Andreas F. Martitsch (28.01.2016)
          ELSE
             p_val = x(2) + p_val * twopi / nfp
             p_s_val = p_s_val * twopi / nfp
          END IF
          !
          IF (isw_eval_spl2d_der .NE. 1) THEN
             PRINT *,"cyl_coord: Evaluation of derivatives from splined &
                  &Fourier coefficients not implemented!"
             STOP
          ELSE
             !
             DO ind1=1,4
                p_spl => R_spl(:,:,:,:,indu(ind1))
                CALL eva2d_fd(theta_n,phi_n,theta_ind,phi_ind,theta_d,phi_d,    &
                     p_spl,der_thph)
                fp(ind1) = der_thph(1)
                fp1(ind1) = der_thph(2)
             END DO
             CALL plag1d(s,fp,dxm1,xp,fun,der)
             r_tb_val = fun
             CALL plag1d(s,fp1,dxm1,xp,fun,der)
             r_pb_val = fun
             !
             DO ind1=1,4
                p_spl => Z_spl(:,:,:,:,indu(ind1))
                CALL eva2d_fd(theta_n,phi_n,theta_ind,phi_ind,theta_d,phi_d,    &
                     p_spl,der_thph)
                fp(ind1) = der_thph(1)
                fp1(ind1) = der_thph(2)
             END DO
             CALL plag1d(s,fp,dxm1,xp,fun,der)
             z_tb_val = fun
             CALL plag1d(s,fp1,dxm1,xp,fun,der)
             z_pb_val = fun
             !
             DO ind1=1,4
                p_spl => L_spl(:,:,:,:,indu(ind1))
                CALL eva2d_fd(theta_n,phi_n,theta_ind,phi_ind,theta_d,phi_d,    &
                     p_spl,der_thph)
                fp(ind1) = der_thph(1)
                fp1(ind1) = der_thph(2)
             END DO
             CALL plag1d(s,fp,dxm1,xp,fun,der)
             p_tb_val = fun
             CALL plag1d(s,fp1,dxm1,xp,fun,der)
             p_pb_val = fun
             !
             IF (lab_swi .EQ. 5 .OR. lab_swi .EQ. 3) THEN ! CHS, LHD
                p_tb_val = - p_tb_val
                p_pb_val = 1.0_dp - p_pb_val
                !! Modifications by Andreas F. Martitsch (28.01.2016)
                ! ToDo: Implement conversion for p_s
                PRINT *,'WARNING FROM NEO_MAGFIE: CONVERSION FOR RADIAL DERIVATIVE OF BOOZER-PHI NOT CHECKED!'
                !! End Modifications by Andreas F. Martitsch (28.01.2016)
             ELSE
                p_tb_val = p_tb_val * twopi / nfp
                p_pb_val = 1.0_dp + p_pb_val * twopi / nfp
             END IF
             !
          END IF
          !! End Modifications by Andreas F. Martitsch (13.11.2014)
          x_cyl = (/ r_val*1e2_dp, p_val, z_val*1e2_dp /)
          jacobian(1,:) = (/ r_s_val*1e2_dp, p_s_val, z_s_val*1e2_dp /)
          jacobian(2,:) = (/ r_tb_val*1e2_dp, p_tb_val, z_tb_val*1e2_dp /)
          jacobian(3,:) = (/ r_pb_val*1e2_dp, p_pb_val, z_pb_val*1e2_dp /)
       ELSE
          PRINT *, 'cyl_coord: s not detected!'
          STOP
       END IF
    END IF
  END SUBROUTINE cyl_coord_a

  SUBROUTINE calc_bN( x, bN, Bvec, nabla_s0_norm )
    ! input / output
    REAL(dp), DIMENSION(:),       INTENT(in)            :: x
    REAL(dp),                     INTENT(out)           :: bN
    REAL(dp), DIMENSION(SIZE(x)), INTENT(out)           :: Bvec
    REAL(dp), DIMENSION(SIZE(x)), INTENT(out)           :: nabla_s0_norm
    ! local definitions
    INTEGER(i4b)                                     :: swd = 1
    INTEGER                                          :: m, n
    INTEGER                                          :: npsi
    REAL(dp)                                         :: m0  = 0.0_dp
    REAL(dp)                                         :: yp, ypp, yppp

    REAL(dp)                                         :: sinv, cosv

    INTEGER                                          :: k_es
    INTEGER                                          :: s_detected
    INTEGER                                          :: imn
    INTEGER                                          :: it, ip, im, in
    INTEGER                                          :: mt = 1
    INTEGER                                          :: mp = 1
    INTEGER                                          :: theta_ind, phi_ind
    INTEGER                                          :: ierr
    REAL(dp)                                         :: s
    REAL(dp)                                         :: magfie_epsi = 1.e-9
    REAL(dp)                                         :: raxc, laxc, zaxc
    REAL(dp)                                         :: raxc_s, laxc_s, zaxc_s
    !! Modifications by Andreas F. Martitsch (07.03.2014)
    ! Auxiliary variables for the Fourier summation
    REAL(dp)                                         :: raxs, laxs, zaxs
    REAL(dp)                                         :: raxs_s, laxs_s, zaxs_s
    !! End Modifications by Andreas F. Martitsch (07.03.2014)
    REAL(dp)                                         :: theta_d, phi_d

    REAL(dp), DIMENSION(:), ALLOCATABLE              :: s_raxc, s_laxc, s_zaxc
    !! Modifications by Andreas F. Martitsch (06.03.2014)
    ! Radial derivatives of (R,Z,phi)-components obtained from the 1d spline
    REAL(dp), DIMENSION(:), ALLOCATABLE              :: s_raxc_s, s_laxc_s, s_zaxc_s
    !! End Modifications by Andreas F. Martitsch (06.03.2014)
    !! Modifications by Andreas F. Martitsch (06.08.2014)
    ! Additional data from Boozer files without Stellarator symmetry
    REAL(kind=dp),    DIMENSION(:),     ALLOCATABLE :: s_raxs, s_laxs, s_zaxs
    REAL(kind=dp),    DIMENSION(:),     ALLOCATABLE :: s_raxs_s, s_laxs_s, s_zaxs_s
    !! End Modifications by Andreas F. Martitsch (06.08.2014)
    REAL(dp), DIMENSION(:,:), ALLOCATABLE            :: rax, lax, zax
    REAL(dp)                                         :: rax_val, pax_val, zax_val
    REAL(dp)                                         :: rax_s_val, pax_s_val, zax_s_val
    REAL(dp)                                         :: rax_tb_val, pax_tb_val, zax_tb_val
    REAL(dp)                                         :: rax_pb_val, pax_pb_val, zax_pb_val
    REAL(dp), DIMENSION(3)                           :: nabla_s0
    REAL(dp), DIMENSION(3,3)                         :: jax, jax_inv

    REAL(dp), DIMENSION(:,:,:,:), POINTER            :: p_spl

    REAL(dp)                 :: bmod, sqrtg
    REAL(dp), DIMENSION(3)   :: bder, hcovar, hctrvr, hcurl
    REAL(dp), DIMENSION(3)   :: x_cyl
    REAL(dp), DIMENSION(3,3) :: jacobian, jacobian_cart

    INTEGER :: ind1
    INTEGER, DIMENSION(mp_lag) :: indu
    REAL(dp), DIMENSION(mp_lag) :: xp, fp, fp1
    REAL(dp) :: fun, der, dxm1
    REAL(dp), DIMENSION(2) :: der_thph

    !*******************************************************************
    ! Initialisation if necessary
    !*******************************************************************
    IF ( .NOT. ALLOCATED(es) ) THEN
       CALL neo_read_control()
       fluxs_interp = 1
       CALL neo_init(npsi)
       PRINT *, 'theta_start,theta_end,phi_start,phi_end'
       PRINT *, theta_start,theta_end,phi_start,phi_end
       !! Modifications by Andreas F. Martitsch (27.01.2016)
       ! -> switch on/off use of splined Fourier coefficients within neo_magfie
       IF (isw_spl_fourier_cof .NE. 1) THEN
          ALLOCATE(magfie_sarray(SIZE(es)/flux_surf_dist))
          magfie_sarray = es(1:SIZE(es):flux_surf_dist)
       END IF
       !! End Modifications by Andreas F. Martitsch (27.01.2016)
    END IF
    !*******************************************************************
    ! Spline of surfaces in magfie_sarray
    !*******************************************************************
    IF (magfie_spline .EQ. 1 .AND. magfie_newspline_bN .EQ. 1) THEN
       magfie_sarray_len =  SIZE(magfie_sarray)
       !****************************************************************
       ! Allocation
       !****************************************************************
       !! Modifications by Andreas F. Martitsch (13.11.2014)
       ! Allocate storage arrays for the 2d periodic splines
       ! of the additionally needed quantities for NTV output
       ALLOCATE( R_axi_spl(4,4,theta_n,phi_n,magfie_sarray_len) )
       ALLOCATE( Z_axi_spl(4,4,theta_n,phi_n,magfie_sarray_len) )
       ALLOCATE( L_axi_spl(4,4,theta_n,phi_n,magfie_sarray_len) )
       !! End Modifications by Andreas F. Martitsch (13.11.2014)
       !****************************************************************
       ! Loop over predefined s-values
       !****************************************************************
       DO k_es = 1, magfie_sarray_len
          s = magfie_sarray(k_es)
          !*************************************************************
          ! Surface
          !*************************************************************
          IF (write_progress .EQ. 1) THEN
             PRINT *, 'Initialize Surface, k_es = ',k_es
          END IF
          ALLOCATE ( s_raxc(mnmax) )
          ALLOCATE ( s_raxc_s(mnmax) )
          ALLOCATE ( s_zaxc(mnmax) )
          ALLOCATE ( s_zaxc_s(mnmax) )
          ALLOCATE ( s_laxc(mnmax) )
          ALLOCATE ( s_laxc_s(mnmax) )
          !
          !! Modifications by Andreas F. Martitsch (06.08.2014)
          ! Additional data from Boozer files without Stellarator symmetry
          IF (inp_swi .EQ. 9) THEN        ! ASDEX-U (E. Strumberger)
             ALLOCATE ( s_raxs(mnmax) )
             ALLOCATE ( s_raxs_s(mnmax) )
             ALLOCATE ( s_zaxs(mnmax) )
             ALLOCATE ( s_zaxs_s(mnmax) )
             ALLOCATE ( s_laxs(mnmax) )
             ALLOCATE ( s_laxs_s(mnmax) )
          END IF
          !! End Modifications by Andreas F. Martitsch (06.08.2014)
          !
          IF (isw_spl_fourier_cof .EQ. 1) THEN
             DO imn = 1, mnmax
                ! Switch swd turns on (1) / off (0) the computation of the
                ! radial derivatives within splint_horner3
                !swd = 0 ! Now we want to compute the radial derivatives
                swd = 1
                CALL splint_horner3(es,                                   &
                     a_rmnc(:,imn), b_rmnc(:,imn),                        &
                     c_rmnc(:,imn), d_rmnc(:,imn),                        &
                     swd, r_mhalf(imn),                                   &
                     s, tf, tfp, tfpp, tfppp,                             &
                     s_raxc(imn), s_raxc_s(imn), ypp, yppp)
                swd = 1
                CALL splint_horner3(es,                                   &
                     a_zmnc(:,imn), b_zmnc(:,imn),                        &
                     c_zmnc(:,imn), d_zmnc(:,imn),                        &
                     swd, r_mhalf(imn),                                   &
                     s, tf, tfp, tfpp, tfppp,                             &
                     s_zaxc(imn), s_zaxc_s(imn), ypp, yppp)
                swd = 1
                CALL splint_horner3(es,                                   &
                     a_lmnc(:,imn), b_lmnc(:,imn),                        &
                     c_lmnc(:,imn), d_lmnc(:,imn),                        &
                     swd, r_mhalf(imn),                                   &
                     s, tf, tfp, tfpp, tfppp,                             &
                     s_laxc(imn), s_laxc_s(imn), ypp, yppp)
                !
                !! Modifications by Andreas F. Martitsch (06.08.2014)
                ! Additional data from Boozer files without Stellarator symmetry
                IF (inp_swi .EQ. 9) THEN        ! ASDEX-U (E. Strumberger)
                   swd = 1
                   CALL splint_horner3(es,                                   &
                        a_rmns(:,imn), b_rmns(:,imn),                        &
                        c_rmns(:,imn), d_rmns(:,imn),                        &
                        swd, r_mhalf(imn),                                   &
                        s, tf, tfp, tfpp, tfppp,                             &
                        s_raxs(imn), s_raxs_s(imn), ypp, yppp)
                   swd = 1
                   CALL splint_horner3(es,                                   &
                        a_zmns(:,imn), b_zmns(:,imn),                        &
                        c_zmns(:,imn), d_zmns(:,imn),                        &
                        swd, r_mhalf(imn),                                   &
                        s, tf, tfp, tfpp, tfppp,                             &
                        s_zaxs(imn), s_zaxs_s(imn), ypp, yppp)
                   swd = 1
                   CALL splint_horner3(es,                                   &
                        a_lmns(:,imn), b_lmns(:,imn),                        &
                        c_lmns(:,imn), d_lmns(:,imn),                        &
                        swd, r_mhalf(imn),                                   &
                        s, tf, tfp, tfpp, tfppp,                             &
                        s_laxs(imn), s_laxs_s(imn), ypp, yppp)
                END IF
                !! End Modifications by Andreas F. Martitsch (06.08.2014)
                !
             END DO
          ELSE
             s_raxc = a_rmnc(1+(k_es-1)*flux_surf_dist,:)
             s_raxc_s = 0.0d0
             s_zaxc = a_zmnc(1+(k_es-1)*flux_surf_dist,:)
             s_zaxc_s = 0.0d0
             s_laxc = a_lmnc(1+(k_es-1)*flux_surf_dist,:)
             s_laxc_s = 0.0d0
             IF (inp_swi .EQ. 9) THEN        ! ASDEX-U (E. Strumberger)
                s_raxs = a_rmns(1+(k_es-1)*flux_surf_dist,:)
                s_raxs_s = 0.0d0
                s_zaxs = a_zmns(1+(k_es-1)*flux_surf_dist,:)
                s_zaxs_s = 0.0d0
                s_laxs = a_lmns(1+(k_es-1)*flux_surf_dist,:)
                s_laxs_s = 0.0d0
             END IF
          END IF
          !*************************************************************
          ! Fourier summation for the full theta-phi array
          !*************************************************************
          IF (write_progress .EQ. 1) THEN
             PRINT *, 'Do Fourier'
          END IF

          ALLOCATE( rax(theta_n,phi_n) )
          rax = 0.0d0
          ALLOCATE( zax(theta_n,phi_n) )
          zax = 0.0d0
          ALLOCATE( lax(theta_n,phi_n) )
          lax = 0.0d0
          DO imn=1,mnmax
             raxc = s_raxc(imn)
             zaxc = s_zaxc(imn)
             laxc = s_laxc(imn)
             ! Additional data from Boozer files without Stellarator symmetry
             IF (inp_swi .EQ. 9) THEN        ! ASDEX-U (E. Strumberger)
                raxs = s_raxs(imn)
                zaxs = s_zaxs(imn)
                laxs = s_laxs(imn)
             END IF
             !! End Modifications by Andreas F. Martitsch (06.08.2014)
             m = ixm(imn)
             n = ixn(imn)
             IF (n .NE. 0) CYCLE ! only axisymm part
             im = pixm(imn)
             in = pixn(imn)
             DO ip=1,phi_n
                DO it=1,theta_n
                   !! Modifications by Andreas F. Martitsch (06.08.2014)
                   ! Additional data from Boozer files without Stellarator symmetry
                   IF (inp_swi .EQ. 9) THEN        ! ASDEX-U (E. Strumberger)
                      cosv = cosmth(it,im) * cosnph(ip,in) - sinmth(it,im) * sinnph(ip,in)
                      sinv = sinmth(it,im) * cosnph(ip,in) + cosmth(it,im) * sinnph(ip,in)

                      rax(it,ip) = rax(it,ip) + raxc*cosv + raxs*sinv
                      zax(it,ip) = zax(it,ip) + zaxc*cosv + zaxs*sinv
                      lax(it,ip) = lax(it,ip) + laxc*cosv + laxs*sinv
                   ELSE
                      cosv = cosmth(it,im) * cosnph(ip,in) + sinmth(it,im) * sinnph(ip,in)
                      sinv = sinmth(it,im) * cosnph(ip,in) - cosmth(it,im) * sinnph(ip,in)

                      rax(it,ip) = rax(it,ip) + raxc*cosv
                      zax(it,ip) = zax(it,ip) + zaxc*sinv
                      lax(it,ip) = lax(it,ip) + laxc*sinv
                   END IF
                   !! End Modifications by Andreas F. Martitsch (06.08.2014)
                END DO
             END DO
          END DO
          DEALLOCATE( s_raxc )
          DEALLOCATE( s_zaxc )
          DEALLOCATE( s_laxc )
          !! Modifications by Andreas F. Martitsch (07.03.2014)
          ! Deallocate arrays for the additional radial derivatives
          DEALLOCATE( s_raxc_s )
          DEALLOCATE( s_zaxc_s )
          DEALLOCATE( s_laxc_s )
          !! End Modifications by Andreas F. Martitsch (07.03.2014)
          !
          !! Modifications by Andreas F. Martitsch (06.08.2014)
          ! Additional data from Boozer files without Stellarator symmetry
          IF (inp_swi .EQ. 9) THEN        ! ASDEX-U (E. Strumberger)
             DEALLOCATE ( s_raxs )
             DEALLOCATE ( s_raxs_s )
             DEALLOCATE ( s_zaxs )
             DEALLOCATE ( s_zaxs_s )
             DEALLOCATE ( s_laxs )
             DEALLOCATE ( s_laxs_s )
          END IF
          !! End Modifications by Andreas F. Martitsch (06.08.2014)

          ! **********************************************************************
          ! Ensure periodicity boundaries to be the same
          ! **********************************************************************
          rax(theta_n,:) = rax(1,:)
          rax(:,phi_n)   = rax(:,1)
          zax(theta_n,:) = zax(1,:)
          zax(:,phi_n)   = zax(:,1)
          lax(theta_n,:) = lax(1,:)
          lax(:,phi_n)   = lax(:,1)

          !*************************************************************
          ! Do the 2-D periodic spline
          !*************************************************************
          p_spl => R_axi_spl(:,:,:,:,k_es)
          CALL spl2d(theta_n,phi_n,theta_int,phi_int,mt,mp,            &
               rax,p_spl)
          p_spl => Z_axi_spl(:,:,:,:,k_es)
          CALL spl2d(theta_n,phi_n,theta_int,phi_int,mt,mp,            &
               zax,p_spl)
          p_spl => L_axi_spl(:,:,:,:,k_es)
          CALL spl2d(theta_n,phi_n,theta_int,phi_int,mt,mp,            &
               lax,p_spl)

          DEALLOCATE( rax )
          DEALLOCATE( zax )
          DEALLOCATE( lax )
       END DO
       magfie_newspline_bN = 0
    END IF

    s_detected = 0
    IF (magfie_spline .EQ. 1) THEN
       s = x(1)
       !****************************************************************
       ! Detection of index
       !****************************************************************
       dxm1 = 1.0d0/(magfie_sarray(2)-magfie_sarray(1))
       CALL indef(s,magfie_sarray(1),dxm1,magfie_sarray_len,indu)
       xp = magfie_sarray(indu)
       s_detected = 1
       k_es = indu(1)
!!$       DO k_es = 1, magfie_sarray_len
!!$          !IF ( ABS(s-magfie_sarray(k_es)) .LT. magfie_epsi) THEN
!!$          IF ( ABS(s-magfie_sarray(k_es)) .LT. (magfie_sarray(2)-magfie_sarray(1))) THEN
!!$             s_detected = 1
!!$             EXIT
!!$          END IF
!!$       END DO
       IF (s_detected .EQ. 1) THEN
          ! ************************************************************
          ! Evaluation of 2d-splines
          ! ************************************************************
          CALL poi2d(theta_int,phi_int,mt,mp,                          &
               theta_start,theta_end,phi_start,phi_end,                &
               x(3),x(2),theta_ind,phi_ind,theta_d,phi_d,ierr)
          !
          DO ind1=1,4
             p_spl => R_axi_spl(:,:,:,:,indu(ind1))
             CALL eva2d(theta_n,phi_n,theta_ind,phi_ind,theta_d,phi_d,    &
                  p_spl,rax_val)
             fp(ind1) = rax_val
          END DO
          CALL plag1d(s,fp,dxm1,xp,fun,der)
          rax_val = fun
          rax_s_val = der
          !
          DO ind1=1,4
             p_spl => L_axi_spl(:,:,:,:,indu(ind1))
             CALL eva2d(theta_n,phi_n,theta_ind,phi_ind,theta_d,phi_d,    &
                  p_spl,pax_val)
             fp(ind1) = pax_val
          END DO
          CALL plag1d(s,fp,dxm1,xp,fun,der)
          pax_val = fun
          pax_s_val = der
          !
          DO ind1=1,4
             p_spl => Z_axi_spl(:,:,:,:,indu(ind1))
             CALL eva2d(theta_n,phi_n,theta_ind,phi_ind,theta_d,phi_d,    &
                  p_spl,zax_val)
             fp(ind1) = zax_val
          END DO
          CALL plag1d(s,fp,dxm1,xp,fun,der)
          zax_val = fun
          zax_s_val = der
          !
          IF (lab_swi .EQ. 5 .OR. lab_swi .EQ. 3) THEN ! CHS, LHD
             pax_val = x(2) - pax_val
             pax_s_val = - pax_s_val
             !! Modifications by Andreas F. Martitsch (28.01.2016)
             ! ToDo: Implement conversion for p_s
             PRINT *,'WARNING FROM NEO_MAGFIE: CONVERSION FOR RADIAL DERIVATIVE OF BOOZER-PHI NOT CHECKED!'
             !! End Modifications by Andreas F. Martitsch (28.01.2016)
          ELSE
             pax_val = x(2) + pax_val * twopi / nfp
             pax_s_val = pax_s_val * twopi / nfp
          END IF
          !
          IF (isw_eval_spl2d_der .NE. 1) THEN
             PRINT *,"cyl_coord: Evaluation of derivatives from splined &
                  &Fourier coefficients not implemented!"
             STOP
          ELSE
             !
             DO ind1=1,4
                p_spl => R_axi_spl(:,:,:,:,indu(ind1))
                CALL eva2d_fd(theta_n,phi_n,theta_ind,phi_ind,theta_d,phi_d,    &
                     p_spl,der_thph)
                fp(ind1) = der_thph(1)
                fp1(ind1) = der_thph(2)
             END DO
             CALL plag1d(s,fp,dxm1,xp,fun,der)
             rax_tb_val = fun
             CALL plag1d(s,fp1,dxm1,xp,fun,der)
             rax_pb_val = fun
             !
             DO ind1=1,4
                p_spl => Z_axi_spl(:,:,:,:,indu(ind1))
                CALL eva2d_fd(theta_n,phi_n,theta_ind,phi_ind,theta_d,phi_d,    &
                     p_spl,der_thph)
                fp(ind1) = der_thph(1)
                fp1(ind1) = der_thph(2)
             END DO
             CALL plag1d(s,fp,dxm1,xp,fun,der)
             zax_tb_val = fun
             CALL plag1d(s,fp1,dxm1,xp,fun,der)
             zax_pb_val = fun
             !
             DO ind1=1,4
                p_spl => L_axi_spl(:,:,:,:,indu(ind1))
                CALL eva2d_fd(theta_n,phi_n,theta_ind,phi_ind,theta_d,phi_d,    &
                     p_spl,der_thph)
                fp(ind1) = der_thph(1)
                fp1(ind1) = der_thph(2)
             END DO
             CALL plag1d(s,fp,dxm1,xp,fun,der)
             pax_tb_val = fun
             CALL plag1d(s,fp1,dxm1,xp,fun,der)
             pax_pb_val = fun
             !
             IF (lab_swi .EQ. 5 .OR. lab_swi .EQ. 3) THEN ! CHS, LHD
                pax_tb_val = - pax_tb_val
                pax_pb_val = 1.0_dp - pax_pb_val
                !! Modifications by Andreas F. Martitsch (28.01.2016)
                ! ToDo: Implement conversion for p_s
                PRINT *,'WARNING FROM NEO_MAGFIE: CONVERSION FOR RADIAL DERIVATIVE OF BOOZER-PHI NOT CHECKED!'
                !! End Modifications by Andreas F. Martitsch (28.01.2016)
             ELSE
                pax_tb_val = pax_tb_val * twopi / nfp
                pax_pb_val = 1.0_dp + pax_pb_val * twopi / nfp
             END IF
             !
          END IF
          !
          jax(1,1) = (COS(pax_val)*rax_s_val-rax_val*SIN(pax_val)*pax_s_val)*1e2_dp
          jax(2,1) = (SIN(pax_val)*rax_s_val+rax_val*COS(pax_val)*pax_s_val)*1e2_dp
          jax(3,1) = zax_s_val*1e2_dp
          jax(1,2) = (COS(pax_val)*rax_tb_val-rax_val*SIN(pax_val)*pax_tb_val)*1e2_dp
          jax(2,2) = (SIN(pax_val)*rax_tb_val+rax_val*COS(pax_val)*pax_tb_val)*1e2_dp
          jax(3,2) = zax_tb_val*1e2_dp
          jax(1,3) = (COS(pax_val)*rax_pb_val-rax_val*SIN(pax_val)*pax_pb_val)*1e2_dp
          jax(2,3) = (SIN(pax_val)*rax_pb_val+rax_val*COS(pax_val)*pax_pb_val)*1e2_dp
          jax(3,3) = zax_pb_val*1e2_dp
          !
          jax_inv = jax
          CALL inv(jax_inv)
          nabla_s0 = jax_inv(1,:)
          nabla_s0_norm = nabla_s0 / SQRT(SUM(nabla_s0**2))
          !
          ! compute normal B-field component bN (w.r.t. axisymm B-field)
          CALL neo_magfie( x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl )
          CALL cyl_coord( x, x_cyl, jacobian )
          !
          jacobian_cart(1,1) = &
               COS(x_cyl(2))*jacobian(1,1)-x_cyl(1)*SIN(x_cyl(2))*jacobian(1,2)
          jacobian_cart(2,1) = &
               SIN(x_cyl(2))*jacobian(1,1)+x_cyl(1)*COS(x_cyl(2))*jacobian(1,2)
          jacobian_cart(3,1) = jacobian(1,3)
          !
          jacobian_cart(1,2) = &
               COS(x_cyl(2))*jacobian(2,1)-x_cyl(1)*SIN(x_cyl(2))*jacobian(2,2)
          jacobian_cart(2,2) = &
               SIN(x_cyl(2))*jacobian(2,1)+x_cyl(1)*COS(x_cyl(2))*jacobian(2,2)
          jacobian_cart(3,2) = jacobian(2,3)
          !
          jacobian_cart(1,3) = &
               COS(x_cyl(2))*jacobian(3,1)-x_cyl(1)*SIN(x_cyl(2))*jacobian(3,2)
          jacobian_cart(2,3) = &
               SIN(x_cyl(2))*jacobian(3,1)+x_cyl(1)*COS(x_cyl(2))*jacobian(3,2)
          jacobian_cart(3,3) = jacobian(3,3)
          !
          Bvec = (jacobian_cart(:,2)*boozer_iota+jacobian_cart(:,3))*&
               hctrvr(2)*(bmod*1.0d4)
          !
          ! conversion to cgs units
          bN = (SUM(jacobian_cart(:,2)*nabla_s0_norm)*boozer_iota + &
               SUM(jacobian_cart(:,3)*nabla_s0_norm))*hctrvr(2)*(bmod*1.0d4)
          !
       ELSE
          PRINT *, 'calc_bN: s not detected!'
          STOP
       END IF
    END IF
  END SUBROUTINE calc_bN

  SUBROUTINE inv(A)
    REAL(dp), DIMENSION(:,:) :: A
    INTEGER,  DIMENSION(SIZE(A,1)) :: ipiv
    REAL(dp), DIMENSION(SIZE(A,1)) :: work
    INTEGER :: n, info

    n = SIZE(A,1)

    !**********************************************************
    ! LU Factorization (LAPACK)
    !**********************************************************
    CALL DGETRF(n, n, A, n, ipiv, info)

    !**********************************************************
    ! Check state
    !**********************************************************
    IF (info /= 0) THEN
       STOP 'Error: Matrix is numerically singular!'
    END IF

    !**********************************************************
    ! Compute inverse matrix (LAPACK)
    !**********************************************************
    CALL DGETRI(n, A, n, ipiv, work, n, info)

    !**********************************************************
    ! Check state
    !**********************************************************
    IF (info /= 0) THEN
       STOP 'Error: Matrix inversion failed!'
    END IF
  END SUBROUTINE inv

  SUBROUTINE plot_Bfield_a( es_val, n_th, n_ph  )
    !
    REAL(dp), INTENT(in) :: es_val
    INTEGER, INTENT(in) :: n_th, n_ph
    !
    INTEGER :: i_row, i_col
    REAL(dp) :: th_val, ph_val
    REAL(dp) :: bmod, sqrtg
    REAL(dp), DIMENSION(3) :: bder, hcovar, hctrvr, hcurl
    REAL(dp), DIMENSION(3) :: x, x_cyl
    REAL(dp), DIMENSION(3,3) :: jacobian
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: R_arr, P_arr, Z_arr, B_arr
    REAL(dp) :: bN
    REAL(dp), DIMENSION(3) :: Bvec, nabla_s0_norm
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: bN_arr
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: Bx_arr, By_arr, Bz_Arr
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: nabla_s0_x_arr, nabla_s0_y_arr, nabla_s0_z_arr

    ALLOCATE(R_arr(n_th,n_ph), P_arr(n_th,n_ph), Z_arr(n_th,n_ph))
    ALLOCATE(B_arr(n_th,n_ph))
    ALLOCATE(bN_arr(n_th,n_ph))
    ALLOCATE(Bx_arr(n_th,n_ph),By_arr(n_th,n_ph),Bz_arr(n_th,n_ph))
    ALLOCATE(nabla_s0_x_arr(n_th,n_ph),nabla_s0_y_arr(n_th,n_ph),nabla_s0_z_arr(n_th,n_ph))

    DO i_row = 1,n_th
       th_val = (TWOPI/DBLE(n_th-1))*DBLE(i_row-1)
       DO i_col = 1,n_ph
          ph_val = (TWOPI/DBLE(n_ph-1))*DBLE(i_col-1)
          x = (/ es_val, ph_val, th_val /)
          CALL neo_magfie( x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl )
          B_arr(i_row,i_col) = bmod
          CALL cyl_coord( x, x_cyl, jacobian )
          R_arr(i_row,i_col) = x_cyl(1)
          P_arr(i_row,i_col) = x_cyl(2)
          Z_arr(i_row,i_col) = x_cyl(3)
          CALL calc_bN( x, bN, Bvec, nabla_s0_norm )
          bN_arr(i_row,i_col) = bN
          Bx_arr(i_row,i_col) = Bvec(1)
          By_arr(i_row,i_col) = Bvec(2)
          Bz_arr(i_row,i_col) = Bvec(3)
          nabla_s0_x_arr(i_row,i_col) = nabla_s0_norm(1)
          nabla_s0_y_arr(i_row,i_col) = nabla_s0_norm(2)
          nabla_s0_z_arr(i_row,i_col) = nabla_s0_norm(3)
       END DO
    END DO

    OPEN(unit=1111,file='B_mat_plot.dat')
    OPEN(unit=2222,file='R_mat_plot.dat')
    OPEN(unit=3333,file='P_mat_plot.dat')
    OPEN(unit=4444,file='Z_mat_plot.dat')
    OPEN(unit=5555,file='bN_mat_plot.dat')
    OPEN(unit=6666,file='Bx_mat_plot.dat')
    OPEN(unit=7777,file='By_mat_plot.dat')
    OPEN(unit=8888,file='Bz_mat_plot.dat')
    OPEN(unit=9999,file='nabla_s0_x_mat_plot.dat')
    OPEN(unit=1122,file='nabla_s0_y_mat_plot.dat')
    OPEN(unit=1133,file='nabla_s0_z_mat_plot.dat')
    !
    DO i_row = 1,n_th
       WRITE(1111,*) (B_arr(i_row,i_col),i_col=1,n_ph)
       WRITE(2222,*) (R_arr(i_row,i_col),i_col=1,n_ph)
       WRITE(3333,*) (P_arr(i_row,i_col),i_col=1,n_ph)
       WRITE(4444,*) (Z_arr(i_row,i_col),i_col=1,n_ph)
       WRITE(5555,*) (bN_arr(i_row,i_col),i_col=1,n_ph)
       WRITE(6666,*) (Bx_arr(i_row,i_col),i_col=1,n_ph)
       WRITE(7777,*) (By_arr(i_row,i_col),i_col=1,n_ph)
       WRITE(8888,*) (Bz_arr(i_row,i_col),i_col=1,n_ph)
       WRITE(9999,*) (nabla_s0_x_arr(i_row,i_col),i_col=1,n_ph)
       WRITE(1122,*) (nabla_s0_y_arr(i_row,i_col),i_col=1,n_ph)
       WRITE(1133,*) (nabla_s0_z_arr(i_row,i_col),i_col=1,n_ph)
    END DO
    !
    CLOSE(unit=1111)
    CLOSE(unit=2222)
    CLOSE(unit=3333)
    CLOSE(unit=4444)
    CLOSE(unit=5555)
    CLOSE(unit=6666)
    CLOSE(unit=7777)
    CLOSE(unit=8888)
    CLOSE(unit=9999)
    CLOSE(unit=1122)
    CLOSE(unit=1133)

    DEALLOCATE(R_arr, P_arr, Z_arr)
    DEALLOCATE(B_arr)
    DEALLOCATE(bN_arr)
    DEALLOCATE(Bx_arr, By_arr, Bz_arr)
    DEALLOCATE(nabla_s0_x_arr,nabla_s0_y_arr,nabla_s0_z_arr)

  END SUBROUTINE plot_Bfield_a

END MODULE neo_magfie_mod
