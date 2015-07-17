MODULE sbcblk_core
   !!======================================================================
   !!                       ***  MODULE  sbcblk_core  ***
   !! Ocean forcing:  momentum, heat and freshwater flux formulation
   !!=====================================================================
   !! History :  1.0  !  2004-08  (U. Schweckendiek)  Original code
   !!            2.0  !  2005-04  (L. Brodeau, A.M. Treguier) additions: 
   !!                           -  new bulk routine for efficiency
   !!                           -  WINDS ARE NOW ASSUMED TO BE AT T POINTS in input files !!!!
   !!                           -  file names and file characteristics in namelist 
   !!                           -  Implement reading of 6-hourly fields   
   !!            3.0  !  2006-06  (G. Madec) sbc rewritting   
   !!             -   !  2006-12  (L. Brodeau) Original code for TURB_CORE_2Z
   !!            3.2  !  2009-04  (B. Lemaire)  Introduce iom_put
   !!            3.3  !  2010-10  (S. Masson)  add diurnal cycle
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   sbc_blk_core  : bulk formulation as ocean surface boundary condition
   !!                   (forced mode, CORE bulk formulea)
   !!   blk_oce_core  : ocean: computes momentum, heat and freshwater fluxes
   !!   blk_ice_core  : ice  : computes momentum, heat and freshwater fluxes
   !!   turb_core     : computes the CORE turbulent transfer coefficients 
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE fldread         ! read input fields
   USE sbc_oce         ! Surface boundary condition: ocean fields
   USE sbcdcy          ! surface boundary condition: diurnal cycle
   USE iom             ! I/O manager library
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! distribued memory computing library
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE prtctl          ! Print control
#if defined key_lim3
   USE sbc_ice         ! Surface boundary condition: ice fields
#endif

   IMPLICIT NONE
   PRIVATE

   PUBLIC   sbc_blk_core         ! routine called in sbcmod module
   PUBLIC   blk_ice_core         ! routine called in sbc_ice_lim module

   INTEGER , PARAMETER ::   jp_wndi = 1           ! index of 10m wind velocity (i-component) (m/s)    at T-point
   INTEGER , PARAMETER ::   jp_wndj = 2           ! index of 10m wind velocity (j-component) (m/s)    at T-point
   INTEGER , PARAMETER ::   jp_humi = 3           ! index of specific humidity               ( - )
   INTEGER , PARAMETER ::   jp_qsr  = 4           ! index of solar heat                      (W/m2)
   INTEGER , PARAMETER ::   jp_qlw  = 5           ! index of Long wave                       (W/m2)
   INTEGER , PARAMETER ::   jp_tair = 6           ! index of 10m air temperature             (Kelvin)
   INTEGER , PARAMETER ::   jp_prec = 7           ! index of total precipitation (rain+snow) (Kg/m2/s)
   INTEGER , PARAMETER ::   jp_snow = 8           ! index of snow (solid prcipitation)       (kg/m2/s)
   INTEGER , PARAMETER ::   jp_tdif = 9           ! index of tau diff associated to HF tau   (N/m2)   at T-point
#if defined key_orca_r025
   INTEGER , PARAMETER ::   jp_swc  = 10          ! index of GEWEX correction for SW radiation  at T-point
   INTEGER , PARAMETER ::   jp_lwc  = 11          ! index of GEWEX correction for LW radiation  at T-point
   INTEGER , PARAMETER ::   jp_prc  = 12          ! index of PMWC correction forat T-point
   INTEGER , PARAMETER ::   jpfld   = 12          ! maximum number of files to read
#else
   INTEGER , PARAMETER ::   jpfld   = 9           ! maximum number of files to read
#endif
   
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf   ! structure of input fields (file informations, fields read)
         
   !                                             !!! CORE bulk parameters
   REAL(wp), PARAMETER ::   rhoa =    1.22        ! air density
   REAL(wp), PARAMETER ::   cpa  = 1000.5         ! specific heat of air
   REAL(wp), PARAMETER ::   Lv   =    2.5e6       ! latent heat of vaporization
   REAL(wp), PARAMETER ::   Ls   =    2.839e6     ! latent heat of sublimation
   REAL(wp), PARAMETER ::   Stef =    5.67e-8     ! Stefan Boltzmann constant
   REAL(wp), PARAMETER ::   Cice =    1.63e-3     ! transfer coefficient over ice
   REAL(wp), PARAMETER ::   albo =    0.066       ! ocean albedo assumed to be contant

   !                                  !!* Namelist namsbc_core : CORE bulk parameters
   LOGICAL  ::   ln_2m     = .FALSE.   ! logical flag for height of air temp. and hum
   LOGICAL  ::   ln_taudif = .FALSE.   ! logical flag to use the "mean of stress module - module of mean stress" data
   REAL(wp) ::   rn_pfac   = 1.        ! multiplication factor for precipitation
#if defined key_orca_r025
   LOGICAL  ::   ln_printdia= .TRUE.     ! logical flag for height of air temp. and hum
   LOGICAL  ::   ln_netsw   = .TRUE.     ! logical flag for height of air temp. and hum
   LOGICAL  ::   ln_core_graceopt=.FALSE., ln_core_spinup=.FALSE.
   LOGICAL  ::   ln_gwxc = .TRUE.
   LOGICAL  ::   ln_corad_antar =.FALSE., ln_corad_arc =.FALSE. , ln_cotair_arc = .FALSE.
   LOGICAL  ::   ln_coprecip =.FALSE.
   REAL(wp) ::   rn_qns_bias = 0._wp     ! heat flux bias

#endif

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO-consortium (2010) 
   !! $Id: sbcblk_core.F90 2801 2011-07-13 17:45:53Z vichi $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS
   
   SUBROUTINE sbc_blk_core( kt )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE sbc_blk_core  ***
      !!                   
      !! ** Purpose :   provide at each time step the surface ocean fluxes
      !!      (momentum, heat, freshwater and runoff) 
      !!
      !! ** Method  : (1) READ each fluxes in NetCDF files:
      !!      the 10m wind velocity (i-component) (m/s)    at T-point
      !!      the 10m wind velocity (j-component) (m/s)    at T-point
      !!      the specific humidity               ( - )
      !!      the solar heat                      (W/m2)
      !!      the Long wave                       (W/m2)
      !!      the 10m air temperature             (Kelvin)
      !!      the total precipitation (rain+snow) (Kg/m2/s)
      !!      the snow (solid prcipitation)       (kg/m2/s)
      !!   OPTIONAL parameter (see ln_taudif namelist flag):
      !!      the tau diff associated to HF tau   (N/m2)   at T-point 
      !!              (2) CALL blk_oce_core
      !!
      !!      C A U T I O N : never mask the surface stress fields
      !!                      the stress is assumed to be in the mesh referential
      !!                      i.e. the (i,j) referential
      !!
      !! ** Action  :   defined at each time-step at the air-sea interface
      !!              - utau, vtau  i- and j-component of the wind stress
      !!              - taum        wind stress module at T-point
      !!              - wndm        10m wind module at T-point
      !!              - qns, qsr    non-slor and solar heat flux
      !!              - emp, emps   evaporation minus precipitation
      !!----------------------------------------------------------------------
#if defined key_orca_r025 && key_lim2
      USE ice_2
#endif
      INTEGER, INTENT(in) ::   kt   ! ocean time step
      !!
      INTEGER  ::   ierror   ! return error code
      INTEGER  ::   ifpr     ! dummy loop indice
      INTEGER  ::   jfld     ! dummy loop arguments
      INTEGER  ::   ji, jj
      !!
      CHARACTER(len=100) ::  cn_dir   !   Root directory for location of core files
      TYPE(FLD_N), DIMENSION(jpfld) ::   slf_i     ! array of namelist informations on the fields to read
      TYPE(FLD_N) ::   sn_wndi, sn_wndj, sn_humi, sn_qsr       ! informations about the fields to be read
      TYPE(FLD_N) ::   sn_qlw , sn_tair, sn_prec, sn_snow      !   "                                 "
      TYPE(FLD_N) ::   sn_tdif                                 !   "                                 "
#if defined key_orca_r025
      TYPE(FLD_N) ::   sn_swc, sn_lwc                          !   "                                 "
      TYPE(FLD_N) ::   sn_prc
      INTEGER  ::   iter_shapiro = 250
      REAL :: zzlat, zzlat1, zzlat2, zfrld, ztmp
      REAL(wp), DIMENSION(jpi,jpj):: xyt,z_qsr,z_qlw,z_qsr1,z_qlw1,z_tair
      REAL(wp), DIMENSION(jpi,jpj):: zqsr_lr, zqsr_hr, zqlw_lr, zqlw_hr, zprec_hr, zprec_lr
      CHARACTER(len=20)  ::  c_kind='ORCA_GLOB'
      NAMELIST/namsbc_core/ cn_dir , ln_2m  , ln_taudif, rn_pfac,           &
         &                  sn_wndi, sn_wndj, sn_humi  , sn_qsr ,           &
         &                  sn_qlw , sn_tair, sn_prec  , sn_snow, sn_tdif,  &
         &                  sn_swc , sn_lwc , sn_prc   , ln_gwxc,           &
         &                  ln_corad_antar, ln_corad_arc, ln_cotair_arc, ln_coprecip ,  &
         &                  rn_qns_bias, ln_printdia, ln_netsw, ln_core_graceopt,ln_core_spinup
      !!---------------------------------------------------------------------
#else
      NAMELIST/namsbc_core/ cn_dir , ln_2m  , ln_taudif, rn_pfac,           &
         &                  sn_wndi, sn_wndj, sn_humi  , sn_qsr ,           &
         &                  sn_qlw , sn_tair, sn_prec  , sn_snow, sn_tdif
      !!---------------------------------------------------------------------
#endif

      !                                         ! ====================== !
      IF( kt == nit000 ) THEN                   !  First call kt=nit000  !
         !                                      ! ====================== !

         ! set file information (default values)
         cn_dir = './'       ! directory in which the model is executed
         !
         ! (NB: frequency positive => hours, negative => months)
         !            !    file    ! frequency ! variable ! time intep !  clim   ! 'yearly' or ! weights  ! rotation !
         !            !    name    !  (hours)  !  name    !   (T/F)    !  (T/F)  !  'monthly'  ! filename ! pairs    !
         sn_wndi = FLD_N( 'uwnd10m',    24     , 'u_10'   ,  .false.   , .false. ,   'yearly'  , ''       , ''       )
         sn_wndj = FLD_N( 'vwnd10m',    24     , 'v_10'   ,  .false.   , .false. ,   'yearly'  , ''       , ''       )
         sn_qsr  = FLD_N( 'qsw'    ,    24     , 'qsw'    ,  .false.   , .false. ,   'yearly'  , ''       , ''       )
         sn_qlw  = FLD_N( 'qlw'    ,    24     , 'qlw'    ,  .false.   , .false. ,   'yearly'  , ''       , ''       )
         sn_tair = FLD_N( 'tair10m',    24     , 't_10'   ,  .false.   , .false. ,   'yearly'  , ''       , ''       )
         sn_humi = FLD_N( 'humi10m',    24     , 'q_10'   ,  .false.   , .false. ,   'yearly'  , ''       , ''       )
         sn_prec = FLD_N( 'precip' ,    -1     , 'precip' ,  .true.    , .false. ,   'yearly'  , ''       , ''       )
         sn_snow = FLD_N( 'snow'   ,    -1     , 'snow'   ,  .true.    , .false. ,   'yearly'  , ''       , ''       )
         sn_tdif = FLD_N( 'taudif' ,    24     , 'taudif' ,  .true.    , .false. ,   'yearly'  , ''       , ''       )
#if defined key_orca_r025
         sn_swc  = FLD_N( 'swc'    ,    24     ,  'swc'   ,  .true.    , .false. ,   'yearly'  , ''       , ''       )
         sn_lwc  = FLD_N( 'lwc'    ,    24     ,  'lwc'   ,  .true.    , .false. ,   'yearly'  , ''       , ''       )
         sn_prc  = FLD_N( 'prc'    ,    24     ,  'prc'   ,  .true.    , .false. ,   'yearly'  , ''       , ''       )
#endif
         !
         REWIND( numnam )                          ! read in namlist namsbc_core
         READ  ( numnam, namsbc_core )
         !                                         ! check: do we plan to use ln_dm2dc with non-daily forcing?
         IF( ln_dm2dc .AND. sn_qsr%nfreqh /= 24 )   & 
            &   CALL ctl_stop( 'sbc_blk_core: ln_dm2dc can be activated only with daily short-wave forcing' ) 
         IF( ln_dm2dc .AND. sn_qsr%ln_tint ) THEN
            CALL ctl_warn( 'sbc_blk_core: ln_dm2dc is taking care of the temporal interpolation of daily qsr',   &
                 &         '              ==> We force time interpolation = .false. for qsr' )
            sn_qsr%ln_tint = .false.
         ENDIF
         !                                         ! store namelist information in an array
         slf_i(jp_wndi) = sn_wndi   ;   slf_i(jp_wndj) = sn_wndj
         slf_i(jp_qsr ) = sn_qsr    ;   slf_i(jp_qlw ) = sn_qlw
         slf_i(jp_tair) = sn_tair   ;   slf_i(jp_humi) = sn_humi
         slf_i(jp_prec) = sn_prec   ;   slf_i(jp_snow) = sn_snow
         slf_i(jp_tdif) = sn_tdif
         !                 
         lhftau = ln_taudif                        ! do we use HF tau information?      
         jfld = jpfld - COUNT( (/.NOT. lhftau/) )
#if defined key_orca_r025
         slf_i(jp_swc ) = sn_swc
         slf_i(jp_lwc ) = sn_lwc
         slf_i(jp_prc ) = sn_prc
         IF( .NOT. ln_gwxc )     jfld = jfld - 2
         IF( .NOT. ln_coprecip ) jfld = jfld - 1
#endif
         !
         ALLOCATE( sf(jfld), STAT=ierror )         ! set sf structure
         IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_blk_core: unable to allocate sf structure' )
         DO ifpr= 1, jfld
            ALLOCATE( sf(ifpr)%fnow(jpi,jpj,1) )
            IF( slf_i(ifpr)%ln_tint )   ALLOCATE( sf(ifpr)%fdta(jpi,jpj,1,2) )
         END DO
         !                                         ! fill sf with slf_i and control print
         CALL fld_fill( sf, slf_i, cn_dir, 'sbc_blk_core', 'flux formulation for ocean surface boundary condition', 'namsbc_core' )
         !
      ENDIF

      CALL fld_read( kt, nn_fsbc, sf )        ! input fields provided at the current time-step

      IF( MOD( kt - 1, nn_fsbc ) == 0 ) THEN

#if defined key_orca_r025
      ! Introduce ERA-Interim filtering and correction

         IF( ln_gwxc ) THEN

           call Shapiro_1D(sf(jp_qsr)%fnow(:,:,1),iter_shapiro, c_kind, zqsr_lr)
           zqsr_hr(:,:)=sf(jp_qsr)%fnow(:,:,1)-zqsr_lr(:,:)          ! We get large scale and small scale

           call Shapiro_1D(sf(jp_qlw)%fnow(:,:,1),iter_shapiro, c_kind, zqlw_lr)
           zqlw_hr(:,:)=sf(jp_qlw)%fnow(:,:,1)-zqlw_lr(:,:)          ! We get large scale and small scale

           z_qsr1(:,:)=zqsr_lr(:,:)*sf(jp_swc)%fnow(:,:,1) + zqsr_hr(:,:)
           z_qlw1(:,:)=zqlw_lr(:,:)*sf(jp_lwc)%fnow(:,:,1) + zqlw_hr(:,:)

           DO jj=1,jpj
             DO ji=1,jpi
               z_qsr1(ji,jj)=max(z_qsr1(ji,jj),0.0)
               z_qlw1(ji,jj)=max(z_qlw1(ji,jj),0.0)
             END DO
           END DO

         ENDIF

         IF( ln_coprecip ) THEN

           call Shapiro_1D(sf(jp_prec)%fnow(:,:,1),iter_shapiro,c_kind,zprec_lr)
           zprec_hr(:,:)=sf(jp_prec)%fnow(:,:,1)-zprec_lr(:,:)       ! We get large scale and small scale

           DO jj=1,jpj
             DO ji=1,jpi
               IF( zprec_lr(ji,jj) .GT. 0._wp ) THEN
                  ztmp = LOG( ( 1000._wp + sf(jp_prc)%fnow(ji,jj,1) ) * EXP( zprec_lr(ji,jj) ) / 1000._wp )
                  sf(jp_prec)%fnow(ji,jj,1) = max(ztmp+zprec_hr(ji,jj),0.0)
               ENDIF
             END DO
           END DO

         ENDIF

         IF ( ln_corad_antar ) THEN           ! correction of SW and LW in the Southern Ocean
 
           z_qsr(:,:)=0.8*z_qsr1(:,:)
           z_qlw(:,:)=1.1*z_qlw1(:,:)
           xyt(:,:) = 0.e0
           zzlat1 = -65.
           zzlat2 = -60.
           DO jj = 1, jpj
             DO ji = 1, jpi
               zzlat = gphit(ji,jj)
               IF ( zzlat >= zzlat1 .AND. zzlat <= zzlat2 ) THEN
                  xyt(ji,jj) = (zzlat2-zzlat)/(zzlat2-zzlat1)
               ELSE IF ( zzlat < zzlat1 ) THEN
                  xyt(ji,jj) = 1
               ENDIF
             END DO
           END DO
           z_qsr1(:,:)=z_qsr(:,:)*xyt(:,:)+(1.0-xyt(:,:))*z_qsr1(:,:)
           z_qlw1(:,:)=z_qlw(:,:)*xyt(:,:)+(1.0-xyt(:,:))*z_qlw1(:,:)

         ENDIF

         IF ( ln_corad_arc ) THEN         ! correction of SW in the Arctic Ocean

           z_qsr(:,:)=0.7*z_qsr1(:,:)
           xyt(:,:) = 0.e0
           zzlat1 = 78.
           zzlat2 = 82.
           DO jj = 1, jpj
             DO ji = 1, jpi
               zzlat = gphit(ji,jj)
               IF ( zzlat >= zzlat1 .AND. zzlat <= zzlat2 ) THEN
                  xyt(ji,jj) = (zzlat-zzlat1)/(zzlat2-zzlat1)
               ELSE IF ( zzlat > zzlat2 ) THEN
                  xyt(ji,jj) = 1
               ENDIF
             END DO
           END DO
           z_qsr1(:,:)=z_qsr(:,:)*xyt(:,:)+(1.0-xyt(:,:))*z_qsr1(:,:)

         ENDIF

         sf(jp_qsr)%fnow(:,:,1)=z_qsr1(:,:)
         sf(jp_qlw)%fnow(:,:,1)=z_qlw1(:,:)

#if defined key_lim2
         IF ( ln_cotair_arc ) THEN           ! correction of Air Temperature in the Arctic Ocean

           z_tair(:,:)=sf(jp_tair)%fnow(:,:,1) - 2.0
           xyt(:,:) = 0.e0 ; zzlat1 = 78. ; zzlat2 = 82.
           DO jj = 1, jpj
             DO ji = 1, jpi
               zzlat = gphit(ji,jj) ; zfrld=frld(ji,jj)
               IF ( zzlat >= zzlat1 .AND. zzlat <= zzlat2 .AND. zfrld < 0.85 ) THEN
                  xyt(ji,jj) = (zzlat-zzlat1)/(zzlat2-zzlat1)
               ELSE IF ( zzlat > zzlat2 .AND. zfrld < 0.85 ) THEN
                  xyt(ji,jj) = 1
               ENDIF
             END DO
           END DO
           sf(jp_tair)%fnow(:,:,1)=z_tair(:,:)*xyt(:,:)+(1.0-xyt(:,:))*sf(jp_tair)%fnow(:,:,1)

         ENDIF
#endif

#endif
         CALL blk_oce_core( sf, sst_m, ssu_m, ssv_m )   ! compute the surface ocean fluxes using CLIO bulk formulea

      ENDIF
      !                                               ! using CORE bulk formulea
   END SUBROUTINE sbc_blk_core
   
   
   SUBROUTINE blk_oce_core( sf, pst, pu, pv )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE blk_core  ***
      !!
      !! ** Purpose :   provide the momentum, heat and freshwater fluxes at
      !!      the ocean surface at each time step
      !!
      !! ** Method  :   CORE bulk formulea for the ocean using atmospheric
      !!      fields read in sbc_read
      !! 
      !! ** Outputs : - utau    : i-component of the stress at U-point  (N/m2)
      !!              - vtau    : j-component of the stress at V-point  (N/m2)
      !!              - taum    : Wind stress module at T-point         (N/m2)
      !!              - wndm    : Wind speed module at T-point          (m/s)
      !!              - qsr     : Solar heat flux over the ocean        (W/m2)
      !!              - qns     : Non Solar heat flux over the ocean    (W/m2)
      !!              - evap    : Evaporation over the ocean            (kg/m2/s)
      !!              - emp(s)  : evaporation minus precipitation       (kg/m2/s)
      !!
      !!  ** Nota  :   sf has to be a dummy argument for AGRIF on NEC
      !!---------------------------------------------------------------------
      USE wrk_nemo, ONLY:   wrk_in_use, wrk_not_released
      USE wrk_nemo, ONLY:   zwnd_i => wrk_2d_1  , zwnd_j => wrk_2d_2      ! wind speed components at T-point
      USE wrk_nemo, ONLY:   zqsatw => wrk_2d_3           ! specific humidity at pst
      USE wrk_nemo, ONLY:   zqlw   => wrk_2d_4  , zqsb   => wrk_2d_5      ! long wave and sensible heat fluxes
      USE wrk_nemo, ONLY:   zqla   => wrk_2d_6  , zevap  => wrk_2d_7      ! latent heat fluxes and evaporation
      USE wrk_nemo, ONLY:   Cd     => wrk_2d_8           ! transfer coefficient for momentum      (tau)
      USE wrk_nemo, ONLY:   Ch     => wrk_2d_9           ! transfer coefficient for sensible heat (Q_sens)
      USE wrk_nemo, ONLY:   Ce     => wrk_2d_10          ! transfer coefficient for evaporation   (Q_lat)
      USE wrk_nemo, ONLY:   zst    => wrk_2d_11          ! surface temperature in Kelvin
      USE wrk_nemo, ONLY:   zt_zu  => wrk_2d_12          ! air temperature at wind speed height
      USE wrk_nemo, ONLY:   zq_zu  => wrk_2d_13          ! air spec. hum.  at wind speed height
      !
      TYPE(fld), INTENT(in), DIMENSION(:)   ::   sf    ! input data
      REAL(wp) , INTENT(in), DIMENSION(:,:) ::   pst   ! surface temperature                      [Celcius]
      REAL(wp) , INTENT(in), DIMENSION(:,:) ::   pu    ! surface current at U-point (i-component) [m/s]
      REAL(wp) , INTENT(in), DIMENSION(:,:) ::   pv    ! surface current at V-point (j-component) [m/s]
      !
      INTEGER  ::   ji, jj               ! dummy loop indices
      REAL(wp) ::   zcoef_qsatw, zztmp   ! local variable
      !!---------------------------------------------------------------------

      IF( wrk_in_use(2, 1,2,3,4,5,6,7,8,9,10,11,12,13) ) THEN
         CALL ctl_stop('blk_oce_core: requested workspace arrays unavailable')   ;   RETURN
      ENDIF
      !
      ! local scalars ( place there for vector optimisation purposes)
      zcoef_qsatw = 0.98 * 640380. / rhoa
      
      zst(:,:) = pst(:,:) + rt0      ! converte Celcius to Kelvin (and set minimum value far above 0 K)

      ! ----------------------------------------------------------------------------- !
      !      0   Wind components and module at T-point relative to the moving ocean   !
      ! ----------------------------------------------------------------------------- !

      ! ... components ( U10m - U_oce ) at T-point (unmasked)
      zwnd_i(:,:) = 0.e0  
      zwnd_j(:,:) = 0.e0
#if defined key_vectopt_loop
!CDIR COLLAPSE
#endif
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1   ! vect. opt.
            zwnd_i(ji,jj) = (  sf(jp_wndi)%fnow(ji,jj,1) - 0.5 * ( pu(ji-1,jj  ) + pu(ji,jj) )  )
            zwnd_j(ji,jj) = (  sf(jp_wndj)%fnow(ji,jj,1) - 0.5 * ( pv(ji  ,jj-1) + pv(ji,jj) )  )
         END DO
      END DO
      CALL lbc_lnk( zwnd_i(:,:) , 'T', -1. )
      CALL lbc_lnk( zwnd_j(:,:) , 'T', -1. )
      ! ... scalar wind ( = | U10m - U_oce | ) at T-point (masked)
!CDIR NOVERRCHK
!CDIR COLLAPSE
      wndm(:,:) = SQRT(  zwnd_i(:,:) * zwnd_i(:,:)   &
         &             + zwnd_j(:,:) * zwnd_j(:,:)  ) * tmask(:,:,1)

      ! ----------------------------------------------------------------------------- !
      !      I   Radiative FLUXES                                                     !
      ! ----------------------------------------------------------------------------- !
    
      ! ocean albedo assumed to be constant + modify now Qsr to include the diurnal cycle                    ! Short Wave
      zztmp = 1. - albo
      IF( ln_dm2dc ) THEN   ;   qsr(:,:) = zztmp * sbc_dcy( sf(jp_qsr)%fnow(:,:,1) ) * tmask(:,:,1)
      ELSE                  ;   qsr(:,:) = zztmp *          sf(jp_qsr)%fnow(:,:,1)   * tmask(:,:,1)
      ENDIF
!CDIR COLLAPSE
      zqlw(:,:) = (  sf(jp_qlw)%fnow(:,:,1) - Stef * zst(:,:)*zst(:,:)*zst(:,:)*zst(:,:)  ) * tmask(:,:,1)   ! Long  Wave
      ! ----------------------------------------------------------------------------- !
      !     II    Turbulent FLUXES                                                    !
      ! ----------------------------------------------------------------------------- !

      ! ... specific humidity at SST and IST
!CDIR NOVERRCHK
!CDIR COLLAPSE
      zqsatw(:,:) = zcoef_qsatw * EXP( -5107.4 / zst(:,:) ) 

      ! ... NCAR Bulk formulae, computation of Cd, Ch, Ce at T-point :
      IF( ln_2m ) THEN
         !! If air temp. and spec. hum. are given at different height (2m) than wind (10m) :
         CALL TURB_CORE_2Z(2.,10., zst   , sf(jp_tair)%fnow,         &
            &                      zqsatw, sf(jp_humi)%fnow, wndm,   &
            &                      Cd    , Ch              , Ce  ,   &
            &                      zt_zu , zq_zu                   )
      ELSE
         !! If air temp. and spec. hum. are given at same height than wind (10m) :
!gm bug?  at the compiling phase, add a copy in temporary arrays...  ==> check perf
!         CALL TURB_CORE_1Z( 10., zst   (:,:), sf(jp_tair)%fnow(:,:),              &
!            &                    zqsatw(:,:), sf(jp_humi)%fnow(:,:), wndm(:,:),   &
!            &                    Cd    (:,:),             Ch  (:,:), Ce  (:,:)  )
!gm bug
! ARPDBG - this won't compile with gfortran. Fix but check performance
! as per comment above.
         CALL TURB_CORE_1Z( 10., zst   , sf(jp_tair)%fnow(:,:,1),       &
            &                    zqsatw, sf(jp_humi)%fnow(:,:,1), wndm, &
            &                    Cd    , Ch              , Ce    )
      ENDIF

      ! ... tau module, i and j component
      DO jj = 1, jpj
         DO ji = 1, jpi
            zztmp = rhoa * wndm(ji,jj) * Cd(ji,jj)
            taum  (ji,jj) = zztmp * wndm  (ji,jj)
            zwnd_i(ji,jj) = zztmp * zwnd_i(ji,jj)
            zwnd_j(ji,jj) = zztmp * zwnd_j(ji,jj)
         END DO
      END DO

      ! ... add the HF tau contribution to the wind stress module?
      IF( lhftau ) THEN 
!CDIR COLLAPSE
#if defined key_orca_r025
         ! Changed!!! Multiply by QSCAT correction
         zwnd_i(:,:) = zwnd_i(:,:) * sf(jp_tdif)%fnow(:,:,1)
         zwnd_j(:,:) = zwnd_j(:,:) * sf(jp_tdif)%fnow(:,:,1)
#endif
         taum(:,:) = taum(:,:) + sf(jp_tdif)%fnow(:,:,1)
      ENDIF
      CALL iom_put( "taum_oce", taum )   ! output wind stress module

      ! ... utau, vtau at U- and V_points, resp.
      !     Note the use of 0.5*(2-umask) in order to unmask the stress along coastlines
      DO jj = 1, jpjm1
         DO ji = 1, fs_jpim1
            utau(ji,jj) = 0.5 * ( 2. - umask(ji,jj,1) ) * ( zwnd_i(ji,jj) + zwnd_i(ji+1,jj  ) )
            vtau(ji,jj) = 0.5 * ( 2. - vmask(ji,jj,1) ) * ( zwnd_j(ji,jj) + zwnd_j(ji  ,jj+1) )
         END DO
      END DO
      CALL lbc_lnk( utau(:,:), 'U', -1. )
      CALL lbc_lnk( vtau(:,:), 'V', -1. )

      !  Turbulent fluxes over ocean
      ! -----------------------------
      IF( ln_2m ) THEN
         ! Values of temp. and hum. adjusted to 10m must be used instead of 2m values
         zevap(:,:) = MAX( 0.e0, rhoa    *Ce(:,:)*( zqsatw(:,:) - zq_zu(:,:) ) * wndm(:,:) )   ! Evaporation
         zqsb (:,:) =            rhoa*cpa*Ch(:,:)*( zst   (:,:) - zt_zu(:,:) ) * wndm(:,:)     ! Sensible Heat
      ELSE
!CDIR COLLAPSE
         zevap(:,:) = MAX( 0.e0, rhoa    *Ce(:,:)*( zqsatw(:,:) - sf(jp_humi)%fnow(:,:,1) ) * wndm(:,:) )   ! Evaporation
!CDIR COLLAPSE
         zqsb (:,:) =            rhoa*cpa*Ch(:,:)*( zst   (:,:) - sf(jp_tair)%fnow(:,:,1) ) * wndm(:,:)     ! Sensible Heat
      ENDIF
!CDIR COLLAPSE
      zqla (:,:) = Lv * zevap(:,:)                                                              ! Latent Heat

      IF(ln_ctl) THEN
         CALL prt_ctl( tab2d_1=zqla  , clinfo1=' blk_oce_core: zqla   : ', tab2d_2=Ce , clinfo2=' Ce  : ' )
         CALL prt_ctl( tab2d_1=zqsb  , clinfo1=' blk_oce_core: zqsb   : ', tab2d_2=Ch , clinfo2=' Ch  : ' )
         CALL prt_ctl( tab2d_1=zqlw  , clinfo1=' blk_oce_core: zqlw   : ', tab2d_2=qsr, clinfo2=' qsr : ' )
         CALL prt_ctl( tab2d_1=zqsatw, clinfo1=' blk_oce_core: zqsatw : ', tab2d_2=zst, clinfo2=' zst : ' )
         CALL prt_ctl( tab2d_1=utau  , clinfo1=' blk_oce_core: utau   : ', mask1=umask,   &
            &          tab2d_2=vtau  , clinfo2=              ' vtau : '  , mask2=vmask )
         CALL prt_ctl( tab2d_1=wndm  , clinfo1=' blk_oce_core: wndm   : ')
         CALL prt_ctl( tab2d_1=zst   , clinfo1=' blk_oce_core: zst    : ')
      ENDIF
       
      ! ----------------------------------------------------------------------------- !
      !     III    Total FLUXES                                                       !
      ! ----------------------------------------------------------------------------- !
     
!CDIR COLLAPSE
      qns(:,:) = zqlw(:,:) - zqsb(:,:) - zqla(:,:)      ! Downward Non Solar flux
!CDIR COLLAPSE
      emp(:,:) = zevap(:,:) - sf(jp_prec)%fnow(:,:,1) * rn_pfac * tmask(:,:,1)
!CDIR COLLAPSE
      emps(:,:) = emp(:,:)
      !
      CALL iom_put( "qlw_oce",   zqlw )                 ! output downward longwave heat over the ocean
      CALL iom_put( "qsb_oce", - zqsb )                 ! output downward sensible heat over the ocean
      CALL iom_put( "qla_oce", - zqla )                 ! output downward latent   heat over the ocean
      CALL iom_put( "qns_oce",   qns  )                 ! output downward non solar heat over the ocean
      !
      IF(ln_ctl) THEN
         CALL prt_ctl(tab2d_1=zqsb , clinfo1=' blk_oce_core: zqsb   : ', tab2d_2=zqlw , clinfo2=' zqlw  : ')
         CALL prt_ctl(tab2d_1=zqla , clinfo1=' blk_oce_core: zqla   : ', tab2d_2=qsr  , clinfo2=' qsr   : ')
         CALL prt_ctl(tab2d_1=pst  , clinfo1=' blk_oce_core: pst    : ', tab2d_2=emp  , clinfo2=' emp   : ')
         CALL prt_ctl(tab2d_1=utau , clinfo1=' blk_oce_core: utau   : ', mask1=umask,   &
            &         tab2d_2=vtau , clinfo2=              ' vtau  : ' , mask2=vmask )
      ENDIF
      !
      IF( wrk_not_released(2, 1,2,3,4,5,6,7,8,9,10,11,12,13) )   &
          CALL ctl_stop('blk_oce_core: failed to release workspace arrays')
      !
   END SUBROUTINE blk_oce_core
   
   
   SUBROUTINE blk_ice_core(  pst   , pui   , pvi   , palb ,   &
      &                      p_taui, p_tauj, p_qns , p_qsr,   &
      &                      p_qla , p_dqns, p_dqla,          &
      &                      p_tpr , p_spr ,                  &
      &                      p_fr1 , p_fr2 , cd_grid, pdim  ) 
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE blk_ice_core  ***
      !!
      !! ** Purpose :   provide the surface boundary condition over sea-ice
      !!
      !! ** Method  :   compute momentum, heat and freshwater exchanged
      !!                between atmosphere and sea-ice using CORE bulk
      !!                formulea, ice variables and read atmmospheric fields.
      !!                NB: ice drag coefficient is assumed to be a constant
      !! 
      !! caution : the net upward water flux has with mm/day unit
      !!---------------------------------------------------------------------
      USE wrk_nemo, ONLY:   wrk_in_use, wrk_not_released
      USE wrk_nemo, ONLY:   z_wnds_t => wrk_2d_1                ! wind speed ( = | U10m - U_ice | ) at T-point
      USE wrk_nemo, ONLY:   wrk_3d_4 , wrk_3d_5 , wrk_3d_6 , wrk_3d_7
      !!
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   pst      ! ice surface temperature (>0, =rt0 over land) [Kelvin]
      REAL(wp), DIMENSION(:,:)  , INTENT(in   ) ::   pui      ! ice surface velocity (i- and i- components      [m/s]
      REAL(wp), DIMENSION(:,:)  , INTENT(in   ) ::   pvi      !    at I-point (B-grid) or U & V-point (C-grid)
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   palb     ! ice albedo (clear sky) (alb_ice_cs)               [%]
      REAL(wp), DIMENSION(:,:)  , INTENT(  out) ::   p_taui   ! i- & j-components of surface ice stress        [N/m2]
      REAL(wp), DIMENSION(:,:)  , INTENT(  out) ::   p_tauj   !   at I-point (B-grid) or U & V-point (C-grid)
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   p_qns    ! non solar heat flux over ice (T-point)         [W/m2]
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   p_qsr    !     solar heat flux over ice (T-point)         [W/m2]
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   p_qla    ! latent    heat flux over ice (T-point)         [W/m2]
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   p_dqns   ! non solar heat sensistivity  (T-point)         [W/m2]
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   p_dqla   ! latent    heat sensistivity  (T-point)         [W/m2]
      REAL(wp), DIMENSION(:,:)  , INTENT(  out) ::   p_tpr    ! total precipitation          (T-point)      [Kg/m2/s]
      REAL(wp), DIMENSION(:,:)  , INTENT(  out) ::   p_spr    ! solid precipitation          (T-point)      [Kg/m2/s]
      REAL(wp), DIMENSION(:,:)  , INTENT(  out) ::   p_fr1    ! 1sr fraction of qsr penetration in ice (T-point)  [%]
      REAL(wp), DIMENSION(:,:)  , INTENT(  out) ::   p_fr2    ! 2nd fraction of qsr penetration in ice (T-point)  [%]
      CHARACTER(len=1)          , INTENT(in   ) ::   cd_grid  ! ice grid ( C or B-grid)
      INTEGER                   , INTENT(in   ) ::   pdim     ! number of ice categories
      !!
      INTEGER  ::   ji, jj, jl    ! dummy loop indices
      INTEGER  ::   ijpl          ! number of ice categories (size of 3rd dim of input arrays)
      REAL(wp) ::   zst2, zst3
      REAL(wp) ::   zcoef_wnorm, zcoef_wnorm2, zcoef_dqlw, zcoef_dqla, zcoef_dqsb
      REAL(wp) ::   zztmp                                        ! temporary variable
      REAL(wp) ::   zcoef_frca                                   ! fractional cloud amount
      REAL(wp) ::   zwnorm_f, zwndi_f , zwndj_f                  ! relative wind module and components at F-point
      REAL(wp) ::             zwndi_t , zwndj_t                  ! relative wind components at T-point
      !!
      REAL(wp), DIMENSION(:,:,:), POINTER ::   z_qlw             ! long wave heat flux over ice
      REAL(wp), DIMENSION(:,:,:), POINTER ::   z_qsb             ! sensible  heat flux over ice
      REAL(wp), DIMENSION(:,:,:), POINTER ::   z_dqlw            ! long wave heat sensitivity over ice
      REAL(wp), DIMENSION(:,:,:), POINTER ::   z_dqsb            ! sensible  heat sensitivity over ice
      !!---------------------------------------------------------------------

      ijpl  = pdim                            ! number of ice categories

      ! Set-up access to workspace arrays
      IF( wrk_in_use(2, 1) .OR. wrk_in_use(3, 4,5,6,7) ) THEN
         CALL ctl_stop('blk_ice_core: requested workspace arrays unavailable')   ;   RETURN
      ELSE IF(ijpl > jpk) THEN
         CALL ctl_stop('blk_ice_core: no. of ice categories > jpk so wrk_nemo 3D workspaces cannot be used.')
         RETURN
      END IF
      ! Set-up pointers to sub-arrays of workspaces
      z_qlw  => wrk_3d_4(:,:,1:ijpl)
      z_qsb  => wrk_3d_5(:,:,1:ijpl)
      z_dqlw => wrk_3d_6(:,:,1:ijpl)
      z_dqsb => wrk_3d_7(:,:,1:ijpl)

      ! local scalars ( place there for vector optimisation purposes)
      zcoef_wnorm  = rhoa * Cice
      zcoef_wnorm2 = rhoa * Cice * 0.5
      zcoef_dqlw   = 4.0 * 0.95 * Stef
      zcoef_dqla   = -Ls * Cice * 11637800. * (-5897.8)
      zcoef_dqsb   = rhoa * cpa * Cice
      zcoef_frca   = 1.0  - 0.3

!!gm brutal....
      z_wnds_t(:,:) = 0.e0
      p_taui  (:,:) = 0.e0
      p_tauj  (:,:) = 0.e0
!!gm end

#if defined key_lim3
      tatm_ice(:,:) = sf(jp_tair)%fnow(:,:,1)   ! LIM3: make Tair available in sea-ice. WARNING allocated after call to ice_init
#endif
      ! ----------------------------------------------------------------------------- !
      !    Wind components and module relative to the moving ocean ( U10m - U_ice )   !
      ! ----------------------------------------------------------------------------- !
      SELECT CASE( cd_grid )
      CASE( 'I' )                  ! B-grid ice dynamics :   I-point (i.e. F-point with sea-ice indexation)
         !                           and scalar wind at T-point ( = | U10m - U_ice | ) (masked)
!CDIR NOVERRCHK
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! B grid : NO vector opt
               ! ... scalar wind at I-point (fld being at T-point)
               zwndi_f = 0.25 * (  sf(jp_wndi)%fnow(ji-1,jj  ,1) + sf(jp_wndi)%fnow(ji  ,jj  ,1)   &
                  &              + sf(jp_wndi)%fnow(ji-1,jj-1,1) + sf(jp_wndi)%fnow(ji  ,jj-1,1)  ) - pui(ji,jj)
               zwndj_f = 0.25 * (  sf(jp_wndj)%fnow(ji-1,jj  ,1) + sf(jp_wndj)%fnow(ji  ,jj  ,1)   &
                  &              + sf(jp_wndj)%fnow(ji-1,jj-1,1) + sf(jp_wndj)%fnow(ji  ,jj-1,1)  ) - pvi(ji,jj)
               zwnorm_f = zcoef_wnorm * SQRT( zwndi_f * zwndi_f + zwndj_f * zwndj_f )
               ! ... ice stress at I-point
               p_taui(ji,jj) = zwnorm_f * zwndi_f
               p_tauj(ji,jj) = zwnorm_f * zwndj_f
               ! ... scalar wind at T-point (fld being at T-point)
               zwndi_t = sf(jp_wndi)%fnow(ji,jj,1) - 0.25 * (  pui(ji,jj+1) + pui(ji+1,jj+1)   &
                  &                                          + pui(ji,jj  ) + pui(ji+1,jj  )  )
               zwndj_t = sf(jp_wndj)%fnow(ji,jj,1) - 0.25 * (  pvi(ji,jj+1) + pvi(ji+1,jj+1)   &
                  &                                          + pvi(ji,jj  ) + pvi(ji+1,jj  )  )
               z_wnds_t(ji,jj)  = SQRT( zwndi_t * zwndi_t + zwndj_t * zwndj_t ) * tmask(ji,jj,1)
            END DO
         END DO
         CALL lbc_lnk( p_taui  , 'I', -1. )
         CALL lbc_lnk( p_tauj  , 'I', -1. )
         CALL lbc_lnk( z_wnds_t, 'T',  1. )
         !
      CASE( 'C' )                  ! C-grid ice dynamics :   U & V-points (same as ocean)
#if defined key_vectopt_loop
!CDIR COLLAPSE
#endif
         DO jj = 2, jpj
            DO ji = fs_2, jpi   ! vect. opt.
               zwndi_t = (  sf(jp_wndi)%fnow(ji,jj,1) - 0.5 * ( pui(ji-1,jj  ) + pui(ji,jj) )  )
               zwndj_t = (  sf(jp_wndj)%fnow(ji,jj,1) - 0.5 * ( pvi(ji  ,jj-1) + pvi(ji,jj) )  )
               z_wnds_t(ji,jj)  = SQRT( zwndi_t * zwndi_t + zwndj_t * zwndj_t ) * tmask(ji,jj,1)
            END DO
         END DO
#if defined key_vectopt_loop
!CDIR COLLAPSE
#endif
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vect. opt.
               p_taui(ji,jj) = zcoef_wnorm2 * ( z_wnds_t(ji+1,jj  ) + z_wnds_t(ji,jj) )                          &
                  &          * ( 0.5 * (sf(jp_wndi)%fnow(ji+1,jj,1) + sf(jp_wndi)%fnow(ji,jj,1) ) - pui(ji,jj) )
               p_tauj(ji,jj) = zcoef_wnorm2 * ( z_wnds_t(ji,jj+1  ) + z_wnds_t(ji,jj) )                          &
                  &          * ( 0.5 * (sf(jp_wndj)%fnow(ji,jj+1,1) + sf(jp_wndj)%fnow(ji,jj,1) ) - pvi(ji,jj) )
            END DO
         END DO
         CALL lbc_lnk( p_taui  , 'U', -1. )
         CALL lbc_lnk( p_tauj  , 'V', -1. )
         CALL lbc_lnk( z_wnds_t, 'T',  1. )
         !
      END SELECT

      zztmp = 1. / ( 1. - albo )
      !                                     ! ========================== !
      DO jl = 1, ijpl                       !  Loop over ice categories  !
         !                                  ! ========================== !
!CDIR NOVERRCHK
!CDIR COLLAPSE
         DO jj = 1 , jpj
!CDIR NOVERRCHK
            DO ji = 1, jpi
               ! ----------------------------!
               !      I   Radiative FLUXES   !
               ! ----------------------------!
               zst2 = pst(ji,jj,jl) * pst(ji,jj,jl)
               zst3 = pst(ji,jj,jl) * zst2
               ! Short Wave (sw)
               p_qsr(ji,jj,jl) = zztmp * ( 1. - palb(ji,jj,jl) ) * qsr(ji,jj)
               ! Long  Wave (lw)
               z_qlw(ji,jj,jl) = 0.95 * (  sf(jp_qlw)%fnow(ji,jj,1) - Stef * pst(ji,jj,jl) * zst3  ) * tmask(ji,jj,1)
               ! lw sensitivity
               z_dqlw(ji,jj,jl) = zcoef_dqlw * zst3                                               

               ! ----------------------------!
               !     II    Turbulent FLUXES  !
               ! ----------------------------!

               ! ... turbulent heat fluxes
               ! Sensible Heat
               z_qsb(ji,jj,jl) = rhoa * cpa * Cice * z_wnds_t(ji,jj) * ( pst(ji,jj,jl) - sf(jp_tair)%fnow(ji,jj,1) )
               ! Latent Heat
               p_qla(ji,jj,jl) = MAX( 0.e0, rhoa * Ls  * Cice * z_wnds_t(ji,jj)   &                           
                  &                    * (  11637800. * EXP( -5897.8 / pst(ji,jj,jl) ) / rhoa - sf(jp_humi)%fnow(ji,jj,1)  ) )
               ! Latent heat sensitivity for ice (Dqla/Dt)
               p_dqla(ji,jj,jl) = zcoef_dqla * z_wnds_t(ji,jj) / ( zst2 ) * EXP( -5897.8 / pst(ji,jj,jl) )
               ! Sensible heat sensitivity (Dqsb_ice/Dtn_ice)
               z_dqsb(ji,jj,jl) = zcoef_dqsb * z_wnds_t(ji,jj)

               ! ----------------------------!
               !     III    Total FLUXES     !
               ! ----------------------------!
               ! Downward Non Solar flux
               p_qns (ji,jj,jl) =     z_qlw (ji,jj,jl) - z_qsb (ji,jj,jl) - p_qla (ji,jj,jl)      
               ! Total non solar heat flux sensitivity for ice
               p_dqns(ji,jj,jl) = - ( z_dqlw(ji,jj,jl) + z_dqsb(ji,jj,jl) + p_dqla(ji,jj,jl) )    
            END DO
            !
         END DO
         !
      END DO
      !
      !--------------------------------------------------------------------
      ! FRACTIONs of net shortwave radiation which is not absorbed in the
      ! thin surface layer and penetrates inside the ice cover
      ! ( Maykut and Untersteiner, 1971 ; Ebert and Curry, 1993 )
    
!CDIR COLLAPSE
      p_fr1(:,:) = ( 0.18 * ( 1.0 - zcoef_frca ) + 0.35 * zcoef_frca )
!CDIR COLLAPSE
      p_fr2(:,:) = ( 0.82 * ( 1.0 - zcoef_frca ) + 0.65 * zcoef_frca )
       
!CDIR COLLAPSE
      p_tpr(:,:) = sf(jp_prec)%fnow(:,:,1) * rn_pfac      ! total precipitation [kg/m2/s]
!CDIR COLLAPSE
      p_spr(:,:) = sf(jp_snow)%fnow(:,:,1) * rn_pfac      ! solid precipitation [kg/m2/s]
      CALL iom_put( 'snowpre', p_spr )                  ! Snow precipitation 
      !
      IF(ln_ctl) THEN
         CALL prt_ctl(tab3d_1=p_qla   , clinfo1=' blk_ice_core: p_qla  : ', tab3d_2=z_qsb   , clinfo2=' z_qsb  : ', kdim=ijpl)
         CALL prt_ctl(tab3d_1=z_qlw   , clinfo1=' blk_ice_core: z_qlw  : ', tab3d_2=p_dqla  , clinfo2=' p_dqla : ', kdim=ijpl)
         CALL prt_ctl(tab3d_1=z_dqsb  , clinfo1=' blk_ice_core: z_dqsb : ', tab3d_2=z_dqlw  , clinfo2=' z_dqlw : ', kdim=ijpl)
         CALL prt_ctl(tab3d_1=p_dqns  , clinfo1=' blk_ice_core: p_dqns : ', tab3d_2=p_qsr   , clinfo2=' p_qsr  : ', kdim=ijpl)
         CALL prt_ctl(tab3d_1=pst     , clinfo1=' blk_ice_core: pst    : ', tab3d_2=p_qns   , clinfo2=' p_qns  : ', kdim=ijpl)
         CALL prt_ctl(tab2d_1=p_tpr   , clinfo1=' blk_ice_core: p_tpr  : ', tab2d_2=p_spr   , clinfo2=' p_spr  : ')
         CALL prt_ctl(tab2d_1=p_taui  , clinfo1=' blk_ice_core: p_taui : ', tab2d_2=p_tauj  , clinfo2=' p_tauj : ')
         CALL prt_ctl(tab2d_1=z_wnds_t, clinfo1=' blk_ice_core: z_wnds_t : ')
      ENDIF

      IF( wrk_not_released(2, 1)       .OR.   &
          wrk_not_released(3, 4,5,6,7) )   CALL ctl_stop('blk_ice_core: failed to release workspace arrays')
      !
   END SUBROUTINE blk_ice_core
  

   SUBROUTINE TURB_CORE_1Z(zu, sst, T_a, q_sat, q_a,   &
      &                        dU , Cd , Ch   , Ce   )
      !!----------------------------------------------------------------------
      !!                      ***  ROUTINE  turb_core  ***
      !!
      !! ** Purpose :   Computes turbulent transfert coefficients of surface
      !!                fluxes according to Large & Yeager (2004)
      !!
      !! ** Method  :   I N E R T I A L   D I S S I P A T I O N   M E T H O D
      !!      Momentum, Latent and sensible heat exchange coefficients
      !!      Caution: this procedure should only be used in cases when air
      !!      temperature (T_air), air specific humidity (q_air) and wind (dU)
      !!      are provided at the same height 'zzu'!
      !!
      !! References :   Large & Yeager, 2004 : ???
      !!----------------------------------------------------------------------
      USE wrk_nemo, ONLY: wrk_in_use, wrk_not_released, iwrk_in_use, iwrk_not_released
      USE wrk_nemo, ONLY: dU10 => wrk_2d_14        ! dU                             [m/s]
      USE wrk_nemo, ONLY: dT => wrk_2d_15          ! air/sea temperature difference   [K]
      USE wrk_nemo, ONLY: dq => wrk_2d_16          ! air/sea humidity difference      [K]
      USE wrk_nemo, ONLY: Cd_n10 => wrk_2d_17      ! 10m neutral drag coefficient
      USE wrk_nemo, ONLY: Ce_n10 => wrk_2d_18      ! 10m neutral latent coefficient
      USE wrk_nemo, ONLY: Ch_n10 => wrk_2d_19      ! 10m neutral sensible coefficient
      USE wrk_nemo, ONLY: sqrt_Cd_n10 => wrk_2d_20 ! root square of Cd_n10
      USE wrk_nemo, ONLY: sqrt_Cd => wrk_2d_21     ! root square of Cd
      USE wrk_nemo, ONLY: T_vpot => wrk_2d_22      ! virtual potential temperature    [K]
      USE wrk_nemo, ONLY: T_star => wrk_2d_23      ! turbulent scale of tem. fluct.
      USE wrk_nemo, ONLY: q_star => wrk_2d_24      ! turbulent humidity of temp. fluct.
      USE wrk_nemo, ONLY: U_star => wrk_2d_25      ! turb. scale of velocity fluct.
      USE wrk_nemo, ONLY: L => wrk_2d_26           ! Monin-Obukov length              [m]
      USE wrk_nemo, ONLY: zeta => wrk_2d_27        ! stability parameter at height zu
      USE wrk_nemo, ONLY: U_n10 => wrk_2d_28       ! neutral wind velocity at 10m     [m]
      USE wrk_nemo, ONLY: xlogt  => wrk_2d_29,    xct => wrk_2d_30,   &
                          zpsi_h => wrk_2d_31, zpsi_m => wrk_2d_32
      USE wrk_nemo, ONLY: stab => iwrk_2d_1      ! 1st guess stability test integer
      !
      REAL(wp)                , INTENT(in   ) ::   zu      ! altitude of wind measurement       [m]
      REAL(wp), DIMENSION(:,:), INTENT(in   ) ::   sst     ! sea surface temperature         [Kelvin]
      REAL(wp), DIMENSION(:,:), INTENT(in   ) ::   T_a     ! potential air temperature       [Kelvin]
      REAL(wp), DIMENSION(:,:), INTENT(in   ) ::   q_sat   ! sea surface specific humidity   [kg/kg]
      REAL(wp), DIMENSION(:,:), INTENT(in   ) ::   q_a     ! specific air humidity           [kg/kg]
      REAL(wp), DIMENSION(:,:), INTENT(in   ) ::   dU      ! wind module |U(zu)-U(0)|        [m/s]
      REAL(wp), DIMENSION(:,:), INTENT(  out) ::   Cd      ! transfert coefficient for momentum       (tau)
      REAL(wp), DIMENSION(:,:), INTENT(  out) ::   Ch      ! transfert coefficient for temperature (Q_sens)
      REAL(wp), DIMENSION(:,:), INTENT(  out) ::   Ce      ! transfert coefficient for evaporation  (Q_lat)
      !!
      INTEGER :: j_itt
      INTEGER , PARAMETER ::   nb_itt = 3
      REAL(wp), PARAMETER ::   grav   = 9.8   ! gravity                       
      REAL(wp), PARAMETER ::   kappa  = 0.4   ! von Karman s constant
      !!----------------------------------------------------------------------

      IF(  wrk_in_use(2,             14,15,16,17,18,19,        &
                         20,21,22,23,24,25,26,27,28,29,        &
                         30,31,32)                      .OR.   &
          iwrk_in_use(2, 1)                               ) THEN
         CALL ctl_stop('TURB_CORE_1Z: requested workspace arrays unavailable')   ;   RETURN
      ENDIF

      !! * Start
      !! Air/sea differences
      dU10 = max(0.5, dU)     ! we don't want to fall under 0.5 m/s
      dT = T_a - sst          ! assuming that T_a is allready the potential temp. at zzu
      dq = q_a - q_sat
      !!    
      !! Virtual potential temperature
      T_vpot = T_a*(1. + 0.608*q_a)
      !!
      !! Neutral Drag Coefficient
      stab    = 0.5 + sign(0.5,dT)    ! stable : stab = 1 ; unstable : stab = 0 
      Cd_n10  = 1E-3 * ( 2.7/dU10 + 0.142 + dU10/13.09 )    !   L & Y eq. (6a)
      sqrt_Cd_n10 = sqrt(Cd_n10)
      Ce_n10  = 1E-3 * ( 34.6 * sqrt_Cd_n10 )               !   L & Y eq. (6b)
      Ch_n10  = 1E-3*sqrt_Cd_n10*(18*stab + 32.7*(1-stab)) !   L & Y eq. (6c), (6d)
      !!
      !! Initializing transfert coefficients with their first guess neutral equivalents :
      Cd = Cd_n10 ;  Ce = Ce_n10 ;  Ch = Ch_n10 ;  sqrt_Cd = sqrt(Cd)

      !! * Now starting iteration loop
      DO j_itt=1, nb_itt
         !! Turbulent scales :
         U_star  = sqrt_Cd*dU10                !   L & Y eq. (7a)
         T_star  = Ch/sqrt_Cd*dT               !   L & Y eq. (7b)
         q_star  = Ce/sqrt_Cd*dq               !   L & Y eq. (7c)

         !! Estimate the Monin-Obukov length :
         L  = (U_star**2)/( kappa*grav*(T_star/T_vpot + q_star/(q_a + 1./0.608)) )

         !! Stability parameters :
         zeta  = zu/L ;  zeta   = sign( min(abs(zeta),10.0), zeta )
         zpsi_h  = psi_h(zeta)
         zpsi_m  = psi_m(zeta)

         !! Shifting the wind speed to 10m and neutral stability :
         U_n10 = dU10*1./(1. + sqrt_Cd_n10/kappa*(log(zu/10.) - zpsi_m)) !  L & Y eq. (9a)

         !! Updating the neutral 10m transfer coefficients :
         Cd_n10  = 1E-3 * (2.7/U_n10 + 0.142 + U_n10/13.09)              !  L & Y eq. (6a)
         sqrt_Cd_n10 = sqrt(Cd_n10)
         Ce_n10  = 1E-3 * (34.6 * sqrt_Cd_n10)                           !  L & Y eq. (6b)
         stab    = 0.5 + sign(0.5,zeta)
         Ch_n10  = 1E-3*sqrt_Cd_n10*(18.*stab + 32.7*(1-stab))           !  L & Y eq. (6c), (6d)

         !! Shifting the neutral  10m transfer coefficients to ( zu , zeta ) :
         !!
         xct = 1. + sqrt_Cd_n10/kappa*(log(zu/10) - zpsi_m)
         Cd  = Cd_n10/(xct*xct) ;  sqrt_Cd = sqrt(Cd)
         !!
         xlogt = log(zu/10.) - zpsi_h
         !!
         xct = 1. + Ch_n10*xlogt/kappa/sqrt_Cd_n10
         Ch  = Ch_n10*sqrt_Cd/sqrt_Cd_n10/xct
         !!
         xct = 1. + Ce_n10*xlogt/kappa/sqrt_Cd_n10
         Ce  = Ce_n10*sqrt_Cd/sqrt_Cd_n10/xct
         !!
      END DO
      !!
      IF( wrk_not_released(2,             14,15,16,17,18,19,          &
         &                    20,21,22,23,24,25,26,27,28,29,          &
         &                    30,31,32                      )   .OR.  &      
         iwrk_not_released(2, 1)                                  )   &
         CALL ctl_stop('TURB_CORE_1Z: failed to release workspace arrays')
      !
    END SUBROUTINE TURB_CORE_1Z


    SUBROUTINE TURB_CORE_2Z(zt, zu, sst, T_zt, q_sat, q_zt, dU, Cd, Ch, Ce, T_zu, q_zu)
      !!----------------------------------------------------------------------
      !!                      ***  ROUTINE  turb_core  ***
      !!
      !! ** Purpose :   Computes turbulent transfert coefficients of surface 
      !!                fluxes according to Large & Yeager (2004).
      !!
      !! ** Method  :   I N E R T I A L   D I S S I P A T I O N   M E T H O D
      !!      Momentum, Latent and sensible heat exchange coefficients
      !!      Caution: this procedure should only be used in cases when air
      !!      temperature (T_air) and air specific humidity (q_air) are at 2m
      !!      whereas wind (dU) is at 10m.
      !!
      !! References :   Large & Yeager, 2004 : ???
      !!----------------------------------------------------------------------
      USE wrk_nemo, ONLY: wrk_in_use, wrk_not_released, iwrk_in_use, iwrk_not_released
      USE wrk_nemo, ONLY: dU10 => wrk_2d_14        ! dU                             [m/s]
      USE wrk_nemo, ONLY: dT => wrk_2d_15          ! air/sea temperature difference   [K]
      USE wrk_nemo, ONLY: dq => wrk_2d_16          ! air/sea humidity difference      [K]
      USE wrk_nemo, ONLY: Cd_n10 => wrk_2d_17      ! 10m neutral drag coefficient
      USE wrk_nemo, ONLY: Ce_n10 => wrk_2d_18      ! 10m neutral latent coefficient
      USE wrk_nemo, ONLY: Ch_n10 => wrk_2d_19      ! 10m neutral sensible coefficient
      USE wrk_nemo, ONLY: sqrt_Cd_n10 => wrk_2d_20 ! root square of Cd_n10
      USE wrk_nemo, ONLY: sqrt_Cd => wrk_2d_21     ! root square of Cd
      USE wrk_nemo, ONLY: T_vpot => wrk_2d_22      ! virtual potential temperature    [K]
      USE wrk_nemo, ONLY: T_star => wrk_2d_23     ! turbulent scale of tem. fluct.
      USE wrk_nemo, ONLY: q_star => wrk_2d_24     ! turbulent humidity of temp. fluct.
      USE wrk_nemo, ONLY: U_star => wrk_2d_25     ! turb. scale of velocity fluct.
      USE wrk_nemo, ONLY: L => wrk_2d_26          ! Monin-Obukov length              [m]
      USE wrk_nemo, ONLY: zeta_u => wrk_2d_27     ! stability parameter at height zu
      USE wrk_nemo, ONLY: zeta_t => wrk_2d_28     ! stability parameter at height zt
      USE wrk_nemo, ONLY: U_n10 => wrk_2d_29      ! neutral wind velocity at 10m     [m]
      USE wrk_nemo, ONLY: xlogt => wrk_2d_30, xct => wrk_2d_31, zpsi_hu => wrk_2d_32, zpsi_ht => wrk_2d_33, zpsi_m => wrk_2d_34
      USE wrk_nemo, ONLY: stab => iwrk_2d_1      ! 1st guess stability test integer
      !!
      REAL(wp), INTENT(in)   :: &
         zt,      &     ! height for T_zt and q_zt                   [m]
         zu             ! height for dU                              [m]
      REAL(wp), INTENT(in), DIMENSION(jpi,jpj) ::  &
         sst,      &     ! sea surface temperature              [Kelvin]
         T_zt,     &     ! potential air temperature            [Kelvin]
         q_sat,    &     ! sea surface specific humidity         [kg/kg]
         q_zt,     &     ! specific air humidity                 [kg/kg]
         dU              ! relative wind module |U(zu)-U(0)|       [m/s]
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj)  ::  &
         Cd,       &     ! transfer coefficient for momentum         (tau)
         Ch,       &     ! transfer coefficient for sensible heat (Q_sens)
         Ce,       &     ! transfert coefficient for evaporation   (Q_lat)
         T_zu,     &     ! air temp. shifted at zu                     [K]
         q_zu            ! spec. hum.  shifted at zu               [kg/kg]

      INTEGER :: j_itt
      INTEGER, PARAMETER :: nb_itt = 3   ! number of itterations
      REAL(wp), PARAMETER ::                        &
         grav   = 9.8,      &  ! gravity                       
         kappa  = 0.4          ! von Karman's constant
      !!----------------------------------------------------------------------
      !!  * Start

      IF(  wrk_in_use(2,             14,15,16,17,18,19,        &
                         20,21,22,23,24,25,26,27,28,29,        &         
                         30,31,32,33,34)                .OR.   &
          iwrk_in_use(2, 1)                               ) THEN
         CALL ctl_stop('TURB_CORE_2Z: requested workspace arrays unavailable')   ;   RETURN
      ENDIF

      !! Initial air/sea differences
      dU10 = max(0.5, dU)      !  we don't want to fall under 0.5 m/s
      dT = T_zt - sst 
      dq = q_zt - q_sat

      !! Neutral Drag Coefficient :
      stab = 0.5 + sign(0.5,dT)                 ! stab = 1  if dT > 0  -> STABLE
      Cd_n10  = 1E-3*( 2.7/dU10 + 0.142 + dU10/13.09 ) 
      sqrt_Cd_n10 = sqrt(Cd_n10)
      Ce_n10  = 1E-3*( 34.6 * sqrt_Cd_n10 )
      Ch_n10  = 1E-3*sqrt_Cd_n10*(18*stab + 32.7*(1 - stab))

      !! Initializing transf. coeff. with their first guess neutral equivalents :
      Cd = Cd_n10 ;  Ce = Ce_n10 ;  Ch = Ch_n10 ;  sqrt_Cd = sqrt(Cd)

      !! Initializing z_u values with z_t values :
      T_zu = T_zt ;  q_zu = q_zt

      !!  * Now starting iteration loop
      DO j_itt=1, nb_itt
         dT = T_zu - sst ;  dq = q_zu - q_sat ! Updating air/sea differences
         T_vpot = T_zu*(1. + 0.608*q_zu)    ! Updating virtual potential temperature at zu
         U_star = sqrt_Cd*dU10                ! Updating turbulent scales :   (L & Y eq. (7))
         T_star  = Ch/sqrt_Cd*dT              !
         q_star  = Ce/sqrt_Cd*dq              !
         !!
         L = (U_star*U_star) &                ! Estimate the Monin-Obukov length at height zu
              & / (kappa*grav/T_vpot*(T_star*(1.+0.608*q_zu) + 0.608*T_zu*q_star))
         !! Stability parameters :
         zeta_u  = zu/L  ;  zeta_u = sign( min(abs(zeta_u),10.0), zeta_u )
         zeta_t  = zt/L  ;  zeta_t = sign( min(abs(zeta_t),10.0), zeta_t )
         zpsi_hu = psi_h(zeta_u)
         zpsi_ht = psi_h(zeta_t)
         zpsi_m  = psi_m(zeta_u)
         !!
         !! Shifting the wind speed to 10m and neutral stability : (L & Y eq.(9a))
!        U_n10 = dU10/(1. + sqrt_Cd_n10/kappa*(log(zu/10.) - psi_m(zeta_u)))
         U_n10 = dU10/(1. + sqrt_Cd_n10/kappa*(log(zu/10.) - zpsi_m))
         !!
         !! Shifting temperature and humidity at zu :          (L & Y eq. (9b-9c))
!        T_zu = T_zt - T_star/kappa*(log(zt/zu) + psi_h(zeta_u) - psi_h(zeta_t))
         T_zu = T_zt - T_star/kappa*(log(zt/zu) + zpsi_hu - zpsi_ht)
!        q_zu = q_zt - q_star/kappa*(log(zt/zu) + psi_h(zeta_u) - psi_h(zeta_t))
         q_zu = q_zt - q_star/kappa*(log(zt/zu) + zpsi_hu - zpsi_ht)
         !!
         !! q_zu cannot have a negative value : forcing 0
         stab = 0.5 + sign(0.5,q_zu) ;  q_zu = stab*q_zu
         !!
         !! Updating the neutral 10m transfer coefficients :
         Cd_n10  = 1E-3 * (2.7/U_n10 + 0.142 + U_n10/13.09)    ! L & Y eq. (6a)
         sqrt_Cd_n10 = sqrt(Cd_n10)
         Ce_n10  = 1E-3 * (34.6 * sqrt_Cd_n10)                 ! L & Y eq. (6b)
         stab    = 0.5 + sign(0.5,zeta_u)
         Ch_n10  = 1E-3*sqrt_Cd_n10*(18.*stab + 32.7*(1-stab)) ! L & Y eq. (6c-6d)
         !!
         !!
         !! Shifting the neutral 10m transfer coefficients to (zu,zeta_u) :
!        xct = 1. + sqrt_Cd_n10/kappa*(log(zu/10.) - psi_m(zeta_u))
         xct = 1. + sqrt_Cd_n10/kappa*(log(zu/10.) - zpsi_m)
         Cd = Cd_n10/(xct*xct) ; sqrt_Cd = sqrt(Cd)
         !!
!        xlogt = log(zu/10.) - psi_h(zeta_u)
         xlogt = log(zu/10.) - zpsi_hu
         !!
         xct = 1. + Ch_n10*xlogt/kappa/sqrt_Cd_n10
         Ch  = Ch_n10*sqrt_Cd/sqrt_Cd_n10/xct
         !!
         xct = 1. + Ce_n10*xlogt/kappa/sqrt_Cd_n10
         Ce  = Ce_n10*sqrt_Cd/sqrt_Cd_n10/xct
         !!
         !!
      END DO
      !!
     IF( wrk_not_released(2,              14,15,16,17,18,19,          &
         &                    20,21,22,23,24,25,26,27,28,29,          &
         &                    30,31,32,33,34                )   .OR.  &  
         iwrk_not_released(2, 1)                                  )   &
         CALL ctl_stop('TURB_CORE_2Z: failed to release workspace arrays')
      !
    END SUBROUTINE TURB_CORE_2Z


    FUNCTION psi_m(zta)   !! Psis, L & Y eq. (8c), (8d), (8e)
      !-------------------------------------------------------------------------------
      USE wrk_nemo, ONLY: wrk_in_use, wrk_not_released
      USE wrk_nemo, ONLY:     X2 => wrk_2d_35
      USE wrk_nemo, ONLY:     X  => wrk_2d_36
      USE wrk_nemo, ONLY: stabit => wrk_2d_37
      !!
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: zta

      REAL(wp), PARAMETER :: pi = 3.141592653589793_wp
      REAL(wp), DIMENSION(jpi,jpj)             :: psi_m
      !-------------------------------------------------------------------------------

      IF( wrk_in_use(2, 35,36,37) ) THEN
         CALL ctl_stop('psi_m: requested workspace arrays unavailable')   ;   RETURN
      ENDIF

      X2 = sqrt(abs(1. - 16.*zta))  ;  X2 = max(X2 , 1.0) ;  X  = sqrt(X2)
      stabit    = 0.5 + sign(0.5,zta)
      psi_m = -5.*zta*stabit  &                                                          ! Stable
         &    + (1. - stabit)*(2*log((1. + X)/2) + log((1. + X2)/2) - 2*atan(X) + pi/2)  ! Unstable 

      IF( wrk_not_released(2, 35,36,37) )   CALL ctl_stop('psi_m: failed to release workspace arrays')
      !
    END FUNCTION psi_m


    FUNCTION psi_h( zta )    !! Psis, L & Y eq. (8c), (8d), (8e)
      !-------------------------------------------------------------------------------
      USE wrk_nemo, ONLY: wrk_in_use, wrk_not_released
      USE wrk_nemo, ONLY:     X2 => wrk_2d_35
      USE wrk_nemo, ONLY:     X  => wrk_2d_36
      USE wrk_nemo, ONLY: stabit => wrk_2d_37
      !
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   zta
      !
      REAL(wp), DIMENSION(jpi,jpj)             ::   psi_h
      !-------------------------------------------------------------------------------

      IF( wrk_in_use(2, 35,36,37) ) THEN
         CALL ctl_stop('psi_h: requested workspace arrays unavailable')   ;   RETURN
      ENDIF

      X2 = sqrt(abs(1. - 16.*zta))  ;  X2 = max(X2 , 1.) ;  X  = sqrt(X2)
      stabit    = 0.5 + sign(0.5,zta)
      psi_h = -5.*zta*stabit  &                                       ! Stable
         &    + (1. - stabit)*(2.*log( (1. + X2)/2. ))                 ! Unstable

      IF( wrk_not_released(2, 35,36,37) )   CALL ctl_stop('psi_h: failed to release workspace arrays')
      !
    END FUNCTION psi_h

    SUBROUTINE Shapiro_1D(rla_varin,id_np, cd_overlap, rlpa_varout) !GIG
      !!=====================================================================
      !!
      !! Description: This function applies a 1D Shapiro filter
      !!              (3 points filter) horizontally to a 2D field
      !!              in regular grid
      !! Arguments :
      !!            rla_varin   : Input variable to filter
      !!            zla_mask    : Input mask variable
      !!            id_np       : Number of Shapiro filter iterations
      !!            cd_overlap  : Logical argument for periodical condition
      !!                          (global ocean case)
      !!            rlpa_varout : Output filtered variable
      !!
      !! History : 08/2009  S. CAILLEAU : from 1st version of N. FERRY
      !!           09/2009  C. REGNIER  : Corrections
      !!
      !!=====================================================================
      IMPLICIT NONE
      INTEGER, INTENT(IN)                       :: id_np
      REAL(wp), DIMENSION(jpi,jpj), INTENT(IN)  :: rla_varin !GIG
      CHARACTER(len=20), INTENT(IN)             :: cd_overlap !GIG
      REAL(wp), DIMENSION(jpi,jpj), INTENT(OUT) :: rlpa_varout !GIG

      REAL(wp), DIMENSION(jpi,jpj)              :: rlpa_varout_tmp
      REAL, PARAMETER                           :: rl_alpha = 1./2.    ! fixed stability coefficient (isotrope case)
      REAL, parameter                           :: rap_aniso_diff_XY=2.25 ! anisotrope case
      REAL                                      :: alphax,alphay, znum, zden,test
      INTEGER                                   :: ji, jj, jn, nn
!
!! rap_aniso_diff_XY=2.25 : valeur trouvée empiriquement pour 140 itération pour le filtre de Shapiro et
!! pour un rapport d'anisotopie de 1.5 : on filtre de plus rapidement en x qu'eny.
!------------------------------------------------------------------------------
!
! Loop on several filter iterations

!     Global ocean case
      IF (( cd_overlap == 'MERCA_GLOB' )   .OR.   &
          ( cd_overlap == 'REGULAR_GLOB' ) .OR.   &
          ( cd_overlap == 'ORCA_GLOB' )) THEN
             rlpa_varout(:,:) = rla_varin(:,:)
             rlpa_varout_tmp(:,:) = rlpa_varout(:,:)
!

       alphax=1./2.
       alphay=1./2.
!  Dx/Dy=rap_aniso_diff_XY  , D_ = vitesse de diffusion
!  140 passes du fitre, Lx/Ly=1.5, le rap_aniso_diff_XY correspondant est:
       IF ( rap_aniso_diff_XY .GE. 1. ) alphay=alphay/rap_aniso_diff_XY
       IF ( rap_aniso_diff_XY .LT. 1. ) alphax=alphax*rap_aniso_diff_XY

        DO jn = 1,id_np   ! number of passes of the filter
            DO ji = 2,jpim1
               DO jj = 2,jpjm1
                  ! We crop on the coast
                   znum = rlpa_varout_tmp(ji,jj)   &
                          + 0.25*alphax*(rlpa_varout_tmp(ji-1,jj  )-rlpa_varout_tmp(ji,jj))*tmask(ji-1,jj  ,1)  &
                          + 0.25*alphax*(rlpa_varout_tmp(ji+1,jj  )-rlpa_varout_tmp(ji,jj))*tmask(ji+1,jj  ,1)  &
                          + 0.25*alphay*(rlpa_varout_tmp(ji  ,jj-1)-rlpa_varout_tmp(ji,jj))*tmask(ji  ,jj-1,1)  &
                          + 0.25*alphay*(rlpa_varout_tmp(ji  ,jj+1)-rlpa_varout_tmp(ji,jj))*tmask(ji  ,jj+1,1)
                   rlpa_varout(ji,jj)=znum*tmask(ji,jj,1)+rla_varin(ji,jj)*(1.-tmask(ji,jj,1))
                ENDDO  ! end loop ji
            ENDDO  ! end loop jj
!
!
!           Periodical condition in case of cd_overlap (global ocean)
!           - on a mercator projection grid we consider that singular point at poles
!             are a mean of the values at points of the previous latitude
!           - on ORCA and regular grid we copy the values at points of the previous latitude
            IF ( cd_overlap == 'MERCAT_GLOB' ) THEN
!GIG case unchecked
               rlpa_varout(1,1) = SUM(rlpa_varout(:,2)) / jpi
               rlpa_varout(jpi,jpj) = SUM(rlpa_varout(:,jpj-1)) / jpi
            ELSE
               call lbc_lnk(rlpa_varout, 'T', 1.) ! Boundary condition
            ENDIF
            rlpa_varout_tmp(:,:) = rlpa_varout(:,:)
         ENDDO  ! end loop jn
      ENDIF

!
    END SUBROUTINE Shapiro_1D

  
   !!======================================================================
END MODULE sbcblk_core
