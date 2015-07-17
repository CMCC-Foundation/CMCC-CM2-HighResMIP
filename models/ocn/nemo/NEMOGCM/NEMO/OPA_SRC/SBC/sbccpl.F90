MODULE sbccpl
   !!======================================================================
   !!                       ***  MODULE  sbccpl  ***
   !! Surface Boundary Condition :  momentum, heat and freshwater fluxes in coupled mode
   !!======================================================================
   !! History :  2.0  ! 2007-06  (R. Redler, N. Keenlyside, W. Park) Original code split into flxmod & taumod
   !!            3.0  ! 2008-02  (G. Madec, C Talandier)  surface module
   !!            3.1  ! 2009_02  (G. Madec, S. Masson, E. Maisonave, A. Caubel) generic coupled interface
   !!----------------------------------------------------------------------
#if defined key_oasis3 || defined key_oasis4
   !!----------------------------------------------------------------------
   !!   'key_oasis3' or 'key_oasis4'   Coupled Ocean/Atmosphere formulation
   !!----------------------------------------------------------------------
   !!   namsbc_cpl      : coupled formulation namlist
   !!   sbc_cpl_init    : initialisation of the coupled exchanges
   !!   sbc_cpl_rcv     : receive fields from the atmosphere over the ocean (ocean only)
   !!                     receive stress from the atmosphere over the ocean (ocean-ice case)
   !!   sbc_cpl_ice_tau : receive stress from the atmosphere over ice
   !!   sbc_cpl_ice_flx : receive fluxes from the atmosphere over ice
   !!   sbc_cpl_snd     : send     fields to the atmosphere
   !!----------------------------------------------------------------------
   USE dom_oce         ! ocean space and time domain
   USE sbc_oce         ! Surface boundary condition: ocean fields
   USE sbc_ice         ! Surface boundary condition: ice fields
   USE sbcdcy          ! surface boundary condition: diurnal cycle
   USE phycst          ! physical constants
#if defined key_lim3
   USE par_ice         ! ice parameters
   USE ice             ! ice variables
#endif
#if defined key_lim2
   USE par_ice_2       ! ice parameters
   USE ice_2           ! ice variables
#endif
#if defined key_oasis3
   USE cpl_oasis3      ! OASIS3 coupling
#endif
#if defined key_oasis4
   USE cpl_oasis4      ! OASIS4 coupling
#endif
   USE geo2ocean       ! 
   USE restart         !
   USE oce   , ONLY : tn, un, vn
   USE albedo          !
   USE in_out_manager  ! I/O manager
   USE iom             ! NetCDF library
   USE lib_mpp         ! distribued memory computing library
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
#if defined key_cpl_carbon_cycle
   USE p4zflx, ONLY : oce_co2
#endif
   USE diaar5, ONLY :   lk_diaar5
   IMPLICIT NONE
   PRIVATE

   PUBLIC   sbc_cpl_rcv        ! routine called by sbc_ice_lim(_2).F90
   PUBLIC   sbc_cpl_snd        ! routine called by step.F90
   PUBLIC   sbc_cpl_ice_tau    ! routine called by sbc_ice_lim(_2).F90
   PUBLIC   sbc_cpl_ice_flx    ! routine called by sbc_ice_lim(_2).F90

   INTEGER, PARAMETER ::   jpr_otx1   =  1            ! 3 atmosphere-ocean stress components on grid 1
   INTEGER, PARAMETER ::   jpr_oty1   =  2            ! 
   INTEGER, PARAMETER ::   jpr_otz1   =  3            ! 
   INTEGER, PARAMETER ::   jpr_otx2   =  4            ! 3 atmosphere-ocean stress components on grid 2
   INTEGER, PARAMETER ::   jpr_oty2   =  5            ! 
   INTEGER, PARAMETER ::   jpr_otz2   =  6            ! 
   INTEGER, PARAMETER ::   jpr_itx1   =  7            ! 3 atmosphere-ice   stress components on grid 1
   INTEGER, PARAMETER ::   jpr_ity1   =  8            ! 
   INTEGER, PARAMETER ::   jpr_itz1   =  9            ! 
   INTEGER, PARAMETER ::   jpr_itx2   = 10            ! 3 atmosphere-ice   stress components on grid 2
   INTEGER, PARAMETER ::   jpr_ity2   = 11            ! 
   INTEGER, PARAMETER ::   jpr_itz2   = 12            ! 
   INTEGER, PARAMETER ::   jpr_qsroce = 13            ! Qsr above the ocean
   INTEGER, PARAMETER ::   jpr_qsrice = 14            ! Qsr above the ice
   INTEGER, PARAMETER ::   jpr_qsrmix = 15 
   INTEGER, PARAMETER ::   jpr_qnsoce = 16            ! Qns above the ocean
   INTEGER, PARAMETER ::   jpr_qnsice = 17            ! Qns above the ice
   INTEGER, PARAMETER ::   jpr_qnsmix = 18
   INTEGER, PARAMETER ::   jpr_rain   = 19            ! total liquid precipitation (rain)
   INTEGER, PARAMETER ::   jpr_snow   = 20            ! solid precipitation over the ocean (snow)
   INTEGER, PARAMETER ::   jpr_tevp   = 21            ! total evaporation
   INTEGER, PARAMETER ::   jpr_ievp   = 22            ! solid evaporation (sublimation)
   INTEGER, PARAMETER ::   jpr_sbpr   = 23            ! sublimation - liquid precipitation - solid precipitation
   INTEGER, PARAMETER ::   jpr_semp   = 24            ! solid freshwater budget (sublimation - snow)
   INTEGER, PARAMETER ::   jpr_oemp   = 25            ! ocean freshwater budget (evap - precip)
   INTEGER, PARAMETER ::   jpr_w10m   = 26            ! 10m wind
   INTEGER, PARAMETER ::   jpr_dqnsdt = 27            ! d(Q non solar)/d(temperature)
   INTEGER, PARAMETER ::   jpr_rnf    = 28            ! runoffs
   INTEGER, PARAMETER ::   jpr_cal    = 29            ! calving
   INTEGER, PARAMETER ::   jpr_taum   = 30            ! wind stress module
#if ! defined key_cpl_carbon_cycle
   INTEGER, PARAMETER ::   jprcv      = 30            ! total number of fields received
#else
   INTEGER, PARAMETER ::   jpr_co2    = 31
   INTEGER, PARAMETER ::   jprcv      = 31            ! total number of fields received
#endif   
   INTEGER, PARAMETER ::   jps_fice   =  1            ! ice fraction 
   INTEGER, PARAMETER ::   jps_toce   =  2            ! ocean temperature
   INTEGER, PARAMETER ::   jps_tice   =  3            ! ice   temperature
   INTEGER, PARAMETER ::   jps_tmix   =  4            ! mixed temperature (ocean+ice)
   INTEGER, PARAMETER ::   jps_albice =  5            ! ice   albedo
   INTEGER, PARAMETER ::   jps_albmix =  6            ! mixed albedo
   INTEGER, PARAMETER ::   jps_hice   =  7            ! ice  thickness
   INTEGER, PARAMETER ::   jps_hsnw   =  8            ! snow thickness
   INTEGER, PARAMETER ::   jps_ocx1   =  9            ! ocean current on grid 1
   INTEGER, PARAMETER ::   jps_ocy1   = 10            !
   INTEGER, PARAMETER ::   jps_ocz1   = 11            !
   INTEGER, PARAMETER ::   jps_ivx1   = 12            ! ice   current on grid 1
   INTEGER, PARAMETER ::   jps_ivy1   = 13            !
   INTEGER, PARAMETER ::   jps_ivz1   = 14            !
#if ! defined key_cpl_carbon_cycle
   INTEGER, PARAMETER ::   jpsnd      = 14            ! total number of fields sended
#else
   INTEGER, PARAMETER ::   jps_co2    = 15
   INTEGER, PARAMETER ::   jpsnd      = 15            ! total number of fields sended
#endif   
   !                                                         !!** namelist namsbc_cpl **
   ! Send to the atmosphere                                   !
   CHARACTER(len=100) ::   cn_snd_temperature = 'oce only'    ! 'oce only' 'weighted oce and ice' or 'mixed oce-ice'
   CHARACTER(len=100) ::   cn_snd_albedo      = 'none'        ! 'none' 'weighted ice' or 'mixed oce-ice'
   CHARACTER(len=100) ::   cn_snd_thickness   = 'none'        ! 'none' or 'weighted ice and snow'
   CHARACTER(len=100) ::   cn_snd_crt_nature  = 'none'        ! 'none' 'oce only' 'weighted oce and ice' or 'mixed oce-ice'   
   CHARACTER(len=100) ::   cn_snd_crt_refere  = 'spherical'   ! 'spherical' or 'cartesian'
   CHARACTER(len=100) ::   cn_snd_crt_orient  = 'local grid'  ! 'eastward-northward' or 'local grid'
   CHARACTER(len=100) ::   cn_snd_crt_grid    = 'T'           ! always at 'T' point
#if defined key_cpl_carbon_cycle 
   CHARACTER(len=100) ::   cn_snd_co2         = 'none'        ! 'none' or 'coupled'
#endif
   ! Received from the atmosphere                             !
   CHARACTER(len=100) ::   cn_rcv_tau_nature  = 'oce only'    ! 'oce only' 'oce and ice' or 'mixed oce-ice'
   CHARACTER(len=100) ::   cn_rcv_tau_refere  = 'spherical'   ! 'spherical' or 'cartesian'
   CHARACTER(len=100) ::   cn_rcv_tau_orient  = 'local grid'  ! 'eastward-northward' or 'local grid'
   CHARACTER(len=100) ::   cn_rcv_tau_grid    = 'T'           ! 'T', 'U,V', 'U,V,I', 'T,I', or 'T,U,V'
   CHARACTER(len=100) ::   cn_rcv_w10m        = 'none'        ! 'none' or 'coupled'
   CHARACTER(len=100) ::   cn_rcv_dqnsdt      = 'none'        ! 'none' or 'coupled'
   CHARACTER(len=100) ::   cn_rcv_qsr         = 'oce only'    ! 'oce only' 'conservative' 'oce and ice' or 'mixed oce-ice'
   CHARACTER(len=100) ::   cn_rcv_qns         = 'oce only'    ! 'oce only' 'conservative' 'oce and ice' or 'mixed oce-ice'
   CHARACTER(len=100) ::   cn_rcv_emp         = 'oce only'    ! 'oce only' 'conservative' or 'oce and ice'
   CHARACTER(len=100) ::   cn_rcv_rnf         = 'coupled'     ! 'coupled' 'climato' or 'mixed'
   CHARACTER(len=100) ::   cn_rcv_cal         = 'none'        ! 'none' or 'coupled'
   CHARACTER(len=100) ::   cn_rcv_taumod      = 'none'        ! 'none' or 'coupled'
#if defined key_cpl_carbon_cycle 
   CHARACTER(len=100) ::   cn_rcv_co2         = 'none'        ! 'none' or 'coupled'
#endif

!!   CHARACTER(len=100), PUBLIC ::   cn_rcv_rnf   !: ???             ==>>  !!gm   treat this case in a different maner
   
   CHARACTER(len=100), DIMENSION(4) ::   cn_snd_crt           ! array combining cn_snd_crt_*
   CHARACTER(len=100), DIMENSION(4) ::   cn_rcv_tau           ! array combining cn_rcv_tau_*

   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   albedo_oce_mix     ! ocean albedo sent to atmosphere (mix clear/overcast sky)

   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   frcv               ! all fields recieved from the atmosphere
   INTEGER , ALLOCATABLE, SAVE, DIMENSION(    :) ::   nrcvinfo           ! OASIS info argument

#if ! defined key_lim2   &&   ! defined key_lim3
   ! quick patch to be able to run the coupled model without sea-ice...
   INTEGER, PARAMETER ::   jpl = 1 
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   hicif, hsnif, u_ice, v_ice,fr1_i0,fr2_i0 ! jpi, jpj
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   tn_ice, alb_ice ! (jpi,jpj,jpl)
   REAL(wp) ::  lfus
#endif

   !! Substitution
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: sbccpl.F90 2715 2011-03-30 15:58:35Z rblod $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS
  
   INTEGER FUNCTION sbc_cpl_alloc()
      !!----------------------------------------------------------------------
      !!             ***  FUNCTION sbc_cpl_alloc  ***
      !!----------------------------------------------------------------------
      INTEGER :: ierr(2)
      !!----------------------------------------------------------------------
      ierr(:) = 0
      !
      ALLOCATE( albedo_oce_mix(jpi,jpj), frcv(jpi,jpj,jprcv), nrcvinfo(jprcv),  STAT=ierr(1) )
      !
#if ! defined key_lim2 && ! defined key_lim3
      ! quick patch to be able to run the coupled model without sea-ice...
      ALLOCATE( hicif(jpi,jpj) , u_ice(jpi,jpj) , fr1_i0(jpi,jpj) , tn_ice (jpi,jpj,jpl) ,     &
                hsnif(jpi,jpj) , v_ice(jpi,jpj) , fr2_i0(jpi,jpj) , alb_ice(jpi,jpj,jpl) , STAT=ierr(2) )
#endif
      sbc_cpl_alloc = MAXVAL( ierr )
      IF( lk_mpp            )   CALL mpp_sum ( sbc_cpl_alloc )
      IF( sbc_cpl_alloc > 0 )   CALL ctl_warn('sbc_cpl_alloc: allocation of arrays failed')
      !
   END FUNCTION sbc_cpl_alloc


   SUBROUTINE sbc_cpl_init( k_ice )     
      !!----------------------------------------------------------------------
      !!             ***  ROUTINE sbc_cpl_init  ***
      !!
      !! ** Purpose :   Initialisation of send and recieved information from
      !!                the atmospheric component
      !!
      !! ** Method  : * Read namsbc_cpl namelist 
      !!              * define the receive interface
      !!              * define the send    interface
      !!              * initialise the OASIS coupler
      !!----------------------------------------------------------------------
      USE wrk_nemo, ONLY:   wrk_in_use, wrk_not_released
      USE wrk_nemo, ONLY:   zacs => wrk_2d_3 , zaos => wrk_2d_4   ! clear & overcast sky albedos
      !!
      INTEGER, INTENT(in) ::   k_ice    ! ice management in the sbc (=0/1/2/3)
      !!
      INTEGER ::   jn   ! dummy loop index
      !!
      NAMELIST/namsbc_cpl/  cn_snd_temperature, cn_snd_albedo    , cn_snd_thickness,                 &          
         cn_snd_crt_nature, cn_snd_crt_refere , cn_snd_crt_orient, cn_snd_crt_grid ,                 &
         cn_rcv_w10m      , cn_rcv_taumod     ,                                                      &
         cn_rcv_tau_nature, cn_rcv_tau_refere , cn_rcv_tau_orient, cn_rcv_tau_grid ,                 &
         cn_rcv_dqnsdt    , cn_rcv_qsr        , cn_rcv_qns       , cn_rcv_emp      , cn_rcv_rnf , cn_rcv_cal
#if defined key_cpl_carbon_cycle 
      NAMELIST/namsbc_cpl_co2/  cn_snd_co2, cn_rcv_co2
#endif
      !!---------------------------------------------------------------------

      IF( wrk_in_use(2, 3,4) ) THEN
         CALL ctl_stop('sbc_cpl_init: requested workspace arrays unavailable')   ;   RETURN
      ENDIF

      ! ================================ !
      !      Namelist informations       !
      ! ================================ !

      REWIND( numnam )                    ! ... read namlist namsbc_cpl
      READ  ( numnam, namsbc_cpl )

      IF(lwp) THEN                        ! control print
         WRITE(numout,*)
         WRITE(numout,*)'sbc_cpl_init : namsbc_cpl namelist '
         WRITE(numout,*)'~~~~~~~~~~~~'
         WRITE(numout,*)'   received fields'
         WRITE(numout,*)'       10m wind module                    cn_rcv_w10m        = ', cn_rcv_w10m 
         WRITE(numout,*)'       surface stress - nature            cn_rcv_tau_nature  = ', cn_rcv_tau_nature
         WRITE(numout,*)'                      - referential       cn_rcv_tau_refere  = ', cn_rcv_tau_refere
         WRITE(numout,*)'                      - orientation       cn_rcv_tau_orient  = ', cn_rcv_tau_orient
         WRITE(numout,*)'                      - mesh              cn_rcv_tau_grid    = ', cn_rcv_tau_grid
         WRITE(numout,*)'       non-solar heat flux sensitivity    cn_rcv_dqnsdt      = ', cn_rcv_dqnsdt
         WRITE(numout,*)'       solar heat flux                    cn_rcv_qsr         = ', cn_rcv_qsr  
         WRITE(numout,*)'       non-solar heat flux                cn_rcv_qns         = ', cn_rcv_qns
         WRITE(numout,*)'       freshwater budget                  cn_rcv_emp         = ', cn_rcv_emp
         WRITE(numout,*)'       runoffs                            cn_rcv_rnf         = ', cn_rcv_rnf
         WRITE(numout,*)'       calving                            cn_rcv_cal         = ', cn_rcv_cal 
         WRITE(numout,*)'       stress module                      cn_rcv_taumod      = ', cn_rcv_taumod
         WRITE(numout,*)'   sent fields'
         WRITE(numout,*)'       surface temperature                cn_snd_temperature = ', cn_snd_temperature
         WRITE(numout,*)'       albedo                             cn_snd_albedo      = ', cn_snd_albedo
         WRITE(numout,*)'       ice/snow thickness                 cn_snd_thickness   = ', cn_snd_thickness  
         WRITE(numout,*)'       surface current - nature           cn_snd_crt_nature  = ', cn_snd_crt_nature 
         WRITE(numout,*)'                       - referential      cn_snd_crt_refere  = ', cn_snd_crt_refere 
         WRITE(numout,*)'                       - orientation      cn_snd_crt_orient  = ', cn_snd_crt_orient
         WRITE(numout,*)'                       - mesh             cn_snd_crt_grid    = ', cn_snd_crt_grid 
      ENDIF

#if defined key_cpl_carbon_cycle 
      REWIND( numnam )                    ! read namlist namsbc_cpl_co2
      READ  ( numnam, namsbc_cpl_co2 )
      IF(lwp) THEN                        ! control print
         WRITE(numout,*)
         WRITE(numout,*)'sbc_cpl_init : namsbc_cpl_co2 namelist '
         WRITE(numout,*)'~~~~~~~~~~~~'
         WRITE(numout,*)'   received fields'
         WRITE(numout,*)'       atm co2                            cn_rcv_co2         = ', cn_rcv_co2
         WRITE(numout,*)'   sent fields'
         WRITE(numout,*)'      oce co2 flux                        cn_snd_co2         = ', cn_snd_co2
          WRITE(numout,*)
      ENDIF
#endif
      ! save current & stress in an array and suppress possible blank in the name
      cn_snd_crt(1) = TRIM( cn_snd_crt_nature )   ;   cn_snd_crt(2) = TRIM( cn_snd_crt_refere )
      cn_snd_crt(3) = TRIM( cn_snd_crt_orient )   ;   cn_snd_crt(4) = TRIM( cn_snd_crt_grid   )
      cn_rcv_tau(1) = TRIM( cn_rcv_tau_nature )   ;   cn_rcv_tau(2) = TRIM( cn_rcv_tau_refere )
      cn_rcv_tau(3) = TRIM( cn_rcv_tau_orient )   ;   cn_rcv_tau(4) = TRIM( cn_rcv_tau_grid   )

      !                                   ! allocate zdfric arrays
      IF( sbc_cpl_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'sbc_cpl_alloc : unable to allocate arrays' )
     
      ! ================================ !
      !   Define the receive interface   !
      ! ================================ !
      nrcvinfo(:) = OASIS_idle   ! needed by nrcvinfo(jpr_otx1) if we do not receive ocean stress 

      ! for each field: define the OASIS name                              (srcv(:)%clname)
      !                 define receive or not from the namelist parameters (srcv(:)%laction)
      !                 define the north fold type of lbc                  (srcv(:)%nsgn)

      ! default definitions of srcv
      srcv(:)%laction = .FALSE.   ;   srcv(:)%clgrid = 'T'   ;   srcv(:)%nsgn = 1.

      !                                                      ! ------------------------- !
      !                                                      ! ice and ocean wind stress !   
      !                                                      ! ------------------------- !
      !                                                           ! Name 
      srcv(jpr_otx1)%clname = 'O_OTaux1'      ! 1st ocean component on grid ONE (T or U)
      srcv(jpr_oty1)%clname = 'O_OTauy1'      ! 2nd   -      -         -     - 
      srcv(jpr_otz1)%clname = 'O_OTauz1'      ! 3rd   -      -         -     - 
      srcv(jpr_otx2)%clname = 'O_OTaux2'      ! 1st ocean component on grid TWO (V)
      srcv(jpr_oty2)%clname = 'O_OTauy2'      ! 2nd   -      -         -     - 
      srcv(jpr_otz2)%clname = 'O_OTauz2'      ! 3rd   -      -         -     - 
      !
      srcv(jpr_itx1)%clname = 'O_ITaux1'      ! 1st  ice  component on grid ONE (T, F, I or U)
      srcv(jpr_ity1)%clname = 'O_ITauy1'      ! 2nd   -      -         -     - 
      srcv(jpr_itz1)%clname = 'O_ITauz1'      ! 3rd   -      -         -     - 
      srcv(jpr_itx2)%clname = 'O_ITaux2'      ! 1st  ice  component on grid TWO (V)
      srcv(jpr_ity2)%clname = 'O_ITauy2'      ! 2nd   -      -         -     - 
      srcv(jpr_itz2)%clname = 'O_ITauz2'      ! 3rd   -      -         -     - 
      ! 
      ! Vectors: change of sign at north fold ONLY if on the local grid
      IF( TRIM( cn_rcv_tau(3) ) == 'local grid' )   srcv(jpr_otx1:jpr_itz2)%nsgn = -1.
      
      !                                                           ! Set grid and action
      SELECT CASE( TRIM( cn_rcv_tau(4) ) )      !  'T', 'U,V', 'U,V,I', 'U,V,F', 'T,I', 'T,F', or 'T,U,V'
      CASE( 'T' ) 
         srcv(jpr_otx1:jpr_itz2)%clgrid  = 'T'        ! oce and ice components given at T-point
         srcv(jpr_otx1:jpr_otz1)%laction = .TRUE.     ! receive oce components on grid 1 
         srcv(jpr_itx1:jpr_itz1)%laction = .TRUE.     ! receive ice components on grid 1 
      CASE( 'U,V' ) 
         srcv(jpr_otx1:jpr_otz1)%clgrid  = 'U'        ! oce components given at U-point
         srcv(jpr_otx2:jpr_otz2)%clgrid  = 'V'        !           and           V-point
         srcv(jpr_itx1:jpr_itz1)%clgrid  = 'U'        ! ice components given at U-point
         srcv(jpr_itx2:jpr_itz2)%clgrid  = 'V'        !           and           V-point
         srcv(jpr_otx1:jpr_itz2)%laction = .TRUE.     ! receive oce and ice components on both grid 1 & 2
      CASE( 'U,V,T' )
         srcv(jpr_otx1:jpr_otz1)%clgrid  = 'U'        ! oce components given at U-point
         srcv(jpr_otx2:jpr_otz2)%clgrid  = 'V'        !           and           V-point
         srcv(jpr_itx1:jpr_itz1)%clgrid  = 'T'        ! ice components given at T-point
         srcv(jpr_otx1:jpr_otz2)%laction = .TRUE.     ! receive oce components on grid 1 & 2
         srcv(jpr_itx1:jpr_itz1)%laction = .TRUE.     ! receive ice components on grid 1 only
      CASE( 'U,V,I' )
         srcv(jpr_otx1:jpr_otz1)%clgrid  = 'U'        ! oce components given at U-point
         srcv(jpr_otx2:jpr_otz2)%clgrid  = 'V'        !           and           V-point
         srcv(jpr_itx1:jpr_itz1)%clgrid  = 'I'        ! ice components given at I-point
         srcv(jpr_otx1:jpr_otz2)%laction = .TRUE.     ! receive oce components on grid 1 & 2
         srcv(jpr_itx1:jpr_itz1)%laction = .TRUE.     ! receive ice components on grid 1 only
      CASE( 'U,V,F' )
         srcv(jpr_otx1:jpr_otz1)%clgrid  = 'U'        ! oce components given at U-point
         srcv(jpr_otx2:jpr_otz2)%clgrid  = 'V'        !           and           V-point
         srcv(jpr_itx1:jpr_itz1)%clgrid  = 'F'        ! ice components given at F-point
         srcv(jpr_otx1:jpr_otz2)%laction = .TRUE.     ! receive oce components on grid 1 & 2
         srcv(jpr_itx1:jpr_itz1)%laction = .TRUE.     ! receive ice components on grid 1 only
      CASE( 'T,I' ) 
         srcv(jpr_otx1:jpr_itz2)%clgrid  = 'T'        ! oce and ice components given at T-point
         srcv(jpr_itx1:jpr_itz1)%clgrid  = 'I'        ! ice components given at I-point
         srcv(jpr_otx1:jpr_otz1)%laction = .TRUE.     ! receive oce components on grid 1 
         srcv(jpr_itx1:jpr_itz1)%laction = .TRUE.     ! receive ice components on grid 1 
      CASE( 'T,F' ) 
         srcv(jpr_otx1:jpr_itz2)%clgrid  = 'T'        ! oce and ice components given at T-point
         srcv(jpr_itx1:jpr_itz1)%clgrid  = 'F'        ! ice components given at F-point
         srcv(jpr_otx1:jpr_otz1)%laction = .TRUE.     ! receive oce components on grid 1 
         srcv(jpr_itx1:jpr_itz1)%laction = .TRUE.     ! receive ice components on grid 1 
      CASE( 'T,U,V' )
         srcv(jpr_otx1:jpr_otz1)%clgrid  = 'T'        ! oce components given at T-point
         srcv(jpr_itx1:jpr_itz1)%clgrid  = 'U'        ! ice components given at U-point
         srcv(jpr_itx2:jpr_itz2)%clgrid  = 'V'        !           and           V-point
         srcv(jpr_otx1:jpr_otz1)%laction = .TRUE.     ! receive oce components on grid 1 only
         srcv(jpr_itx1:jpr_itz2)%laction = .TRUE.     ! receive ice components on grid 1 & 2
      CASE default   
         CALL ctl_stop( 'sbc_cpl_init: wrong definition of cn_rcv_tau(4)' )
      END SELECT
      !
      IF( TRIM( cn_rcv_tau(2) ) == 'spherical' )   &           ! spherical: 3rd component not received
         &     srcv( (/jpr_otz1, jpr_otz2, jpr_itz1, jpr_itz2/) )%laction = .FALSE. 
      !
      IF( TRIM( cn_rcv_tau(1) ) /= 'oce and ice' ) THEN        ! 'oce and ice' case ocean stress on ocean mesh used
         srcv(jpr_itx1:jpr_itz2)%laction = .FALSE.    ! ice components not received
         srcv(jpr_itx1)%clgrid = 'U'                  ! ocean stress used after its transformation
         srcv(jpr_ity1)%clgrid = 'V'                  ! i.e. it is always at U- & V-points for i- & j-comp. resp.
      ENDIF
       
      !                                                      ! ------------------------- !
      !                                                      !    freshwater budget      !   E-P
      !                                                      ! ------------------------- !
      ! we suppose that atmosphere modele do not make the difference between precipiration (liquide or solid)
      ! over ice of free ocean within the same atmospheric cell.cd 
      srcv(jpr_rain)%clname = 'OTotRain'      ! Rain = liquid precipitation
      srcv(jpr_snow)%clname = 'OTotSnow'      ! Snow = solid precipitation
      srcv(jpr_tevp)%clname = 'OTotEvap'      ! total evaporation (over oce + ice sublimation)
      srcv(jpr_ievp)%clname = 'OIceEvap'      ! evaporation over ice = sublimation
      srcv(jpr_sbpr)%clname = 'OSubMPre'      ! sublimation - liquid precipitation - solid precipitation 
      srcv(jpr_semp)%clname = 'OISubMSn'      ! ice solid water budget = sublimation - solid precipitation
      srcv(jpr_oemp)%clname = 'OOEvaMPr'      ! ocean water budget = ocean Evap - ocean precip
      SELECT CASE( TRIM( cn_rcv_emp ) )
      CASE( 'oce only'      )   ;   srcv(                                 jpr_oemp   )%laction = .TRUE. 
      CASE( 'conservative'  )   ;   srcv( (/jpr_rain, jpr_snow, jpr_ievp, jpr_tevp/) )%laction = .TRUE.
      CASE( 'oce and ice'   )   ;   srcv( (/jpr_ievp, jpr_sbpr, jpr_semp, jpr_oemp/) )%laction = .TRUE.
      CASE default              ;   CALL ctl_stop( 'sbc_cpl_init: wrong definition of cn_rcv_emp' )
      END SELECT

      !                                                      ! ------------------------- !
      !                                                      !     Runoffs & Calving     !   
      !                                                      ! ------------------------- !
      srcv(jpr_rnf   )%clname = 'O_Runoff'   ;   IF( TRIM( cn_rcv_rnf ) == 'coupled' )   srcv(jpr_rnf)%laction = .TRUE.
                                                 IF( TRIM( cn_rcv_rnf ) == 'climato' )   THEN   ;   ln_rnf = .TRUE.
                                                 ELSE                                           ;   ln_rnf = .FALSE.
                                                 ENDIF
      srcv(jpr_cal   )%clname = 'OCalving'   ;   IF( TRIM( cn_rcv_cal ) == 'coupled' )   srcv(jpr_cal)%laction = .TRUE.

      !                                                      ! ------------------------- !
      !                                                      !    non solar radiation    !   Qns
      !                                                      ! ------------------------- !
      srcv(jpr_qnsoce)%clname = 'O_QnsOce'
      srcv(jpr_qnsice)%clname = 'O_QnsIce'
      srcv(jpr_qnsmix)%clname = 'O_QnsMix'
      SELECT CASE( TRIM( cn_rcv_qns ) )
      CASE( 'oce only'      )   ;   srcv(               jpr_qnsoce   )%laction = .TRUE.
      CASE( 'conservative'  )   ;   srcv( (/jpr_qnsice, jpr_qnsmix/) )%laction = .TRUE.
      CASE( 'oce and ice'   )   ;   srcv( (/jpr_qnsice, jpr_qnsoce/) )%laction = .TRUE.
      CASE( 'mixed oce-ice' )   ;   srcv(               jpr_qnsmix   )%laction = .TRUE. 
      CASE default              ;   CALL ctl_stop( 'sbc_cpl_init: wrong definition of cn_rcv_qns' )
      END SELECT

      !                                                      ! ------------------------- !
      !                                                      !    solar radiation        !   Qsr
      !                                                      ! ------------------------- !
      srcv(jpr_qsroce)%clname = 'O_QsrOce'
      srcv(jpr_qsrice)%clname = 'O_QsrIce'
      srcv(jpr_qsrmix)%clname = 'O_QsrMix'
      SELECT CASE( TRIM( cn_rcv_qsr ) )
      CASE( 'oce only'      )   ;   srcv(               jpr_qsroce   )%laction = .TRUE.
      CASE( 'conservative'  )   ;   srcv( (/jpr_qsrice, jpr_qsrmix/) )%laction = .TRUE.
      CASE( 'oce and ice'   )   ;   srcv( (/jpr_qsrice, jpr_qsroce/) )%laction = .TRUE.
      CASE( 'mixed oce-ice' )   ;   srcv(               jpr_qsrmix   )%laction = .TRUE. 
      CASE default              ;   CALL ctl_stop( 'sbc_cpl_init: wrong definition of cn_rcv_qsr' )
      END SELECT

      !                                                      ! ------------------------- !
      !                                                      !   non solar sensitivity   !   d(Qns)/d(T)
      !                                                      ! ------------------------- !
      srcv(jpr_dqnsdt)%clname = 'O_dQnsdT'   
      IF( TRIM( cn_rcv_dqnsdt ) == 'coupled' )   srcv(jpr_dqnsdt)%laction = .TRUE.
      !
      ! non solar sensitivity mandatory for ice model
      IF( TRIM( cn_rcv_dqnsdt ) == 'none' .AND. k_ice /= 0 ) &
         CALL ctl_stop( 'sbc_cpl_init: cn_rcv_dqnsdt must be coupled in namsbc_cpl namelist' )
      ! non solar sensitivity mandatory for mixed oce-ice solar radiation coupling technique
      IF( TRIM( cn_rcv_dqnsdt ) == 'none' .AND. TRIM( cn_rcv_qns ) == 'mixed oce-ice' ) &
         CALL ctl_stop( 'sbc_cpl_init: namsbc_cpl namelist mismatch between cn_rcv_qns and cn_rcv_dqnsdt' )
      !                                                      ! ------------------------- !
      !                                                      !    Ice Qsr penetration    !   
      !                                                      ! ------------------------- !
      ! fraction of net shortwave radiation which is not absorbed in the thin surface layer 
      ! and penetrates inside the ice cover ( Maykut and Untersteiner, 1971 ; Elbert anbd Curry, 1993 )
      ! Coupled case: since cloud cover is not received from atmosphere 
      !               ===> defined as constant value -> definition done in sbc_cpl_init
      fr1_i0(:,:) = 0.18
      fr2_i0(:,:) = 0.82
      !                                                      ! ------------------------- !
      !                                                      !      10m wind module      !   
      !                                                      ! ------------------------- !
      srcv(jpr_w10m)%clname = 'O_Wind10'   ;   IF( TRIM(cn_rcv_w10m  ) == 'coupled' )   srcv(jpr_w10m)%laction = .TRUE. 
      !
      !                                                      ! ------------------------- !
      !                                                      !   wind stress module      !   
      !                                                      ! ------------------------- !
      srcv(jpr_taum)%clname = 'O_TauMod'   ;   IF( TRIM(cn_rcv_taumod) == 'coupled' )   srcv(jpr_taum)%laction = .TRUE.
      lhftau = srcv(jpr_taum)%laction

#if defined key_cpl_carbon_cycle
      !                                                      ! ------------------------- !
      !                                                      !      Atmospheric CO2      !
      !                                                      ! ------------------------- !
      srcv(jpr_co2 )%clname = 'O_AtmCO2'   ;   IF( TRIM(cn_rcv_co2   ) == 'coupled' )    srcv(jpr_co2 )%laction = .TRUE.
#endif
     
      ! ================================ !
      !     Define the send interface    !
      ! ================================ !
      ! for each field: define the OASIS name                           (srcv(:)%clname)
      !                 define send or not from the namelist parameters (srcv(:)%laction)
      !                 define the north fold type of lbc               (srcv(:)%nsgn)
      
      ! default definitions of nsnd
      ssnd(:)%laction = .FALSE.   ;   ssnd(:)%clgrid = 'T'   ;   ssnd(:)%nsgn = 1.
         
      !                                                      ! ------------------------- !
      !                                                      !    Surface temperature    !
      !                                                      ! ------------------------- !
      ssnd(jps_toce)%clname = 'O_SSTSST'
      ssnd(jps_tice)%clname = 'O_TepIce'
      ssnd(jps_tmix)%clname = 'O_TepMix'
      SELECT CASE( TRIM( cn_snd_temperature ) )
      CASE( 'oce only'             )   ;   ssnd(   jps_toce             )%laction = .TRUE.
      CASE( 'weighted oce and ice' )   ;   ssnd( (/jps_toce, jps_tice/) )%laction = .TRUE.
      CASE( 'mixed oce-ice'        )   ;   ssnd(   jps_tmix             )%laction = .TRUE.
      CASE default   ;   CALL ctl_stop( 'sbc_cpl_init: wrong definition of cn_snd_temperature' )
      END SELECT
     
      !                                                      ! ------------------------- !
      !                                                      !          Albedo           !
      !                                                      ! ------------------------- !
      ssnd(jps_albice)%clname = 'O_AlbIce' 
      ssnd(jps_albmix)%clname = 'O_AlbMix'
      SELECT CASE( TRIM( cn_snd_albedo ) )
      CASE( 'none'          )       ! nothing to do
      CASE( 'weighted ice'  )   ;   ssnd(jps_albice)%laction = .TRUE.
      CASE( 'mixed oce-ice' )   ;   ssnd(jps_albmix)%laction = .TRUE.
      CASE default   ;   CALL ctl_stop( 'sbc_cpl_init: wrong definition of cn_snd_albedo' )
      END SELECT
      !
      ! Need to calculate oceanic albedo if
      !     1. sending mixed oce-ice albedo or
      !     2. receiving mixed oce-ice solar radiation 
      IF ( TRIM ( cn_snd_albedo ) == 'mixed oce-ice' .OR. TRIM ( cn_rcv_qsr ) == 'mixed oce-ice' ) THEN
         CALL albedo_oce( zaos, zacs )
         ! Due to lack of information on nebulosity : mean clear/overcast sky
         albedo_oce_mix(:,:) = ( zacs(:,:) + zaos(:,:) ) * 0.5
      ENDIF

      !                                                      ! ------------------------- !
      !                                                      !  Ice fraction & Thickness ! 
      !                                                      ! ------------------------- !
      ssnd(jps_fice)%clname = 'OIceFrac'   
      ssnd(jps_hice)%clname = 'O_IceTck'
      ssnd(jps_hsnw)%clname = 'O_SnwTck'
      IF( k_ice /= 0 )   ssnd(jps_fice)%laction = .TRUE.       ! if ice treated in the ocean (even in climato case)
      IF( TRIM( cn_snd_thickness ) == 'weighted ice and snow' )   ssnd( (/jps_hice, jps_hsnw/) )%laction = .TRUE.
         
      !                                                      ! ------------------------- !
      !                                                      !      Surface current      !
      !                                                      ! ------------------------- !
      !        ocean currents              !            ice velocities
      ssnd(jps_ocx1)%clname = 'O_OCurx1'   ;   ssnd(jps_ivx1)%clname = 'O_IVelx1'
      ssnd(jps_ocy1)%clname = 'O_OCury1'   ;   ssnd(jps_ivy1)%clname = 'O_IVely1'
      ssnd(jps_ocz1)%clname = 'O_OCurz1'   ;   ssnd(jps_ivz1)%clname = 'O_IVelz1'
      !
      ssnd(jps_ocx1:jps_ivz1)%nsgn = -1.   ! vectors: change of the sign at the north fold

      IF( cn_snd_crt(4) /= 'T' )   CALL ctl_stop( 'cn_snd_crt(4) must be equal to T' )
      ssnd(jps_ocx1:jps_ivz1)%clgrid  = 'T'      ! all oce and ice components on the same unique grid

      ssnd(jps_ocx1:jps_ivz1)%laction = .TRUE.   ! default: all are send
      IF( TRIM( cn_snd_crt(2) ) == 'spherical' )   ssnd( (/jps_ocz1, jps_ivz1/) )%laction = .FALSE. 
      SELECT CASE( TRIM( cn_snd_crt(1) ) )
      CASE( 'none'                 )   ;   ssnd(jps_ocx1:jps_ivz1)%laction = .FALSE.
      CASE( 'oce only'             )   ;   ssnd(jps_ivx1:jps_ivz1)%laction = .FALSE.
      CASE( 'weighted oce and ice' )   !   nothing to do
      CASE( 'mixed oce-ice'        )   ;   ssnd(jps_ivx1:jps_ivz1)%laction = .FALSE.
      CASE default   ;   CALL ctl_stop( 'sbc_cpl_init: wrong definition of cn_snd_crt(1)' )
      END SELECT

#if defined key_cpl_carbon_cycle
      !                                                      ! ------------------------- !
      !                                                      !          CO2 flux         !
      !                                                      ! ------------------------- !
      ssnd(jps_co2)%clname = 'O_CO2FLX' ;  IF( TRIM(cn_snd_co2) == 'coupled' )    ssnd(jps_co2 )%laction = .TRUE.
#endif
      !
      ! ================================ !
      !   initialisation of the coupler  !
      ! ================================ !

      CALL cpl_prism_define(jprcv, jpsnd)            
      !
      IF( ln_dm2dc .AND. ( cpl_prism_freq( jpr_qsroce ) + cpl_prism_freq( jpr_qsrmix ) /= 86400 ) )   &
         &   CALL ctl_stop( 'sbc_cpl_init: diurnal cycle reconstruction (ln_dm2dc) needs daily couping for solar radiation' )

      IF( wrk_not_released(2, 3,4) )   CALL ctl_stop('sbc_cpl_init: failed to release workspace arrays')
      !
   END SUBROUTINE sbc_cpl_init


   SUBROUTINE sbc_cpl_rcv( kt, k_fsbc, k_ice )     
      !!----------------------------------------------------------------------
      !!             ***  ROUTINE sbc_cpl_rcv  ***
      !!
      !! ** Purpose :   provide the stress over the ocean and, if no sea-ice,
      !!                provide the ocean heat and freshwater fluxes.
      !!
      !! ** Method  : - Receive all the atmospheric fields (stored in frcv array). called at each time step.
      !!                OASIS controls if there is something do receive or not. nrcvinfo contains the info
      !!                to know if the field was really received or not
      !!
      !!              --> If ocean stress was really received:
      !!
      !!                  - transform the received ocean stress vector from the received
      !!                 referential and grid into an atmosphere-ocean stress in 
      !!                 the (i,j) ocean referencial and at the ocean velocity point. 
      !!                    The received stress are :
      !!                     - defined by 3 components (if cartesian coordinate)
      !!                            or by 2 components (if spherical)
      !!                     - oriented along geographical   coordinate (if eastward-northward)
      !!                            or  along the local grid coordinate (if local grid)
      !!                     - given at U- and V-point, resp.   if received on 2 grids
      !!                            or at T-point               if received on 1 grid
      !!                    Therefore and if necessary, they are successively 
      !!                  processed in order to obtain them 
      !!                     first  as  2 components on the sphere 
      !!                     second as  2 components oriented along the local grid
      !!                     third  as  2 components on the U,V grid 
      !!
      !!              --> 
      !!
      !!              - In 'ocean only' case, non solar and solar ocean heat fluxes 
      !!             and total ocean freshwater fluxes  
      !!
      !! ** Method  :   receive all fields from the atmosphere and transform 
      !!              them into ocean surface boundary condition fields 
      !!
      !! ** Action  :   update  utau, vtau   ocean stress at U,V grid 
      !!                        taum, wndm   wind stres and wind speed module at T-point
      !!                        qns , qsr    non solar and solar ocean heat fluxes   ('ocean only case)
      !!                        emp = emps   evap. - precip. (- runoffs) (- calving) ('ocean only case)
      !!----------------------------------------------------------------------
      USE wrk_nemo, ONLY:   wrk_in_use, wrk_not_released
      USE wrk_nemo, ONLY:   ztx => wrk_2d_1 , zty => wrk_2d_2
      !!
      INTEGER, INTENT(in) ::   kt       ! ocean model time step index
      INTEGER, INTENT(in) ::   k_fsbc   ! frequency of sbc (-> ice model) computation 
      INTEGER, INTENT(in) ::   k_ice    ! ice management in the sbc (=0/1/2/3)
      !!
      LOGICAL ::    llnewtx, llnewtau      ! update wind stress components and module??
      INTEGER  ::   ji, jj, jn             ! dummy loop indices
      INTEGER  ::   isec                   ! number of seconds since nit000 (assuming rdttra did not change since nit000)
      REAL(wp) ::   zcumulneg, zcumulpos   ! temporary scalars     
      REAL(wp) ::   zcoef                  ! temporary scalar
      REAL(wp) ::   zrhoa  = 1.22          ! Air density kg/m3
      REAL(wp) ::   zcdrag = 1.5e-3        ! drag coefficient
      REAL(wp) ::   zzx, zzy               ! temporary variables
      !!----------------------------------------------------------------------

      IF( wrk_in_use(2, 1,2) ) THEN
         CALL ctl_stop('sbc_cpl_rcv: requested workspace arrays unavailable')   ;   RETURN
      ENDIF

      IF( kt == nit000 )   CALL sbc_cpl_init( k_ice )          ! initialisation

      !                                                 ! Receive all the atmos. fields (including ice information)
      isec = ( kt - nit000 ) * NINT( rdttra(1) )             ! date of exchanges
      DO jn = 1, jprcv                                       ! received fields sent by the atmosphere
         IF( srcv(jn)%laction )   CALL cpl_prism_rcv( jn, isec, frcv(:,:,jn), nrcvinfo(jn) )
      END DO

      !                                                      ! ========================= !
      IF( srcv(jpr_otx1)%laction ) THEN                      !  ocean stress components  !
         !                                                   ! ========================= !
         ! define frcv(:,:,jpr_otx1) and frcv(:,:,jpr_oty1): stress at U/V point along model grid
         ! => need to be done only when we receive the field
         IF(  nrcvinfo(jpr_otx1) == OASIS_Rcv ) THEN
            !
            IF( TRIM( cn_rcv_tau(2) ) == 'cartesian' ) THEN            ! 2 components on the sphere
               !                                                       ! (cartesian to spherical -> 3 to 2 components)
               !
               CALL geo2oce( frcv(:,:,jpr_otx1), frcv(:,:,jpr_oty1), frcv(:,:,jpr_otz1),   &
                  &          srcv(jpr_otx1)%clgrid, ztx, zty )
               frcv(:,:,jpr_otx1) = ztx(:,:)   ! overwrite 1st comp. on the 1st grid
               frcv(:,:,jpr_oty1) = zty(:,:)   ! overwrite 2nd comp. on the 1st grid
               !
               IF( srcv(jpr_otx2)%laction ) THEN
                  CALL geo2oce( frcv(:,:,jpr_otx2), frcv(:,:,jpr_oty2), frcv(:,:,jpr_otz2),   &
                     &          srcv(jpr_otx2)%clgrid, ztx, zty )
                  frcv(:,:,jpr_otx2) = ztx(:,:)   ! overwrite 1st comp. on the 2nd grid
                  frcv(:,:,jpr_oty2) = zty(:,:)   ! overwrite 2nd comp. on the 2nd grid
               ENDIF
               !
            ENDIF
            !
            IF( TRIM( cn_rcv_tau(3) ) == 'eastward-northward' ) THEN   ! 2 components oriented along the local grid
               !                                                       ! (geographical to local grid -> rotate the components)
               CALL rot_rep( frcv(:,:,jpr_otx1), frcv(:,:,jpr_oty1), srcv(jpr_otx1)%clgrid, 'en->i', ztx )   
               frcv(:,:,jpr_otx1) = ztx(:,:)      ! overwrite 1st component on the 1st grid
               IF( srcv(jpr_otx2)%laction ) THEN
                  CALL rot_rep( frcv(:,:,jpr_otx2), frcv(:,:,jpr_oty2), srcv(jpr_otx2)%clgrid, 'en->j', zty )   
               ELSE
                  CALL rot_rep( frcv(:,:,jpr_otx1), frcv(:,:,jpr_oty1), srcv(jpr_otx1)%clgrid, 'en->j', zty )  
               ENDIF
               frcv(:,:,jpr_oty1) = zty(:,:)      ! overwrite 2nd component on the 2nd grid
            ENDIF
            !                              
            IF( srcv(jpr_otx1)%clgrid == 'T' ) THEN
               DO jj = 2, jpjm1                                          ! T ==> (U,V)
                  DO ji = fs_2, fs_jpim1   ! vector opt.
                     frcv(ji,jj,jpr_otx1) = 0.5 * ( frcv(ji+1,jj  ,jpr_otx1) + frcv(ji,jj,jpr_otx1) )
                     frcv(ji,jj,jpr_oty1) = 0.5 * ( frcv(ji  ,jj+1,jpr_oty1) + frcv(ji,jj,jpr_oty1) )
                  END DO
               END DO
               CALL lbc_lnk( frcv(:,:,jpr_otx1), 'U',  -1. )   ;   CALL lbc_lnk( frcv(:,:,jpr_oty1), 'V',  -1. )
            ENDIF
            llnewtx = .TRUE.
         ELSE
            llnewtx = .FALSE.
         ENDIF
         !                                                   ! ========================= !
      ELSE                                                   !   No dynamical coupling   !
         !                                                   ! ========================= !
         frcv(:,:,jpr_otx1) = 0.e0                               ! here simply set to zero 
         frcv(:,:,jpr_oty1) = 0.e0                               ! an external read in a file can be added instead
         llnewtx = .TRUE.
         !
      ENDIF
      
      !                                                      ! ========================= !
      !                                                      !    wind stress module     !   (taum)
      !                                                      ! ========================= !
      !
      IF( .NOT. srcv(jpr_taum)%laction ) THEN                    ! compute wind stress module from its components if not received 
         ! => need to be done only when otx1 was changed
         IF( llnewtx ) THEN
!CDIR NOVERRCHK
            DO jj = 2, jpjm1
!CDIR NOVERRCHK
               DO ji = fs_2, fs_jpim1   ! vect. opt.
                  zzx = frcv(ji-1,jj  ,jpr_otx1) + frcv(ji,jj,jpr_otx1) 
                  zzy = frcv(ji  ,jj-1,jpr_oty1) + frcv(ji,jj,jpr_oty1) 
                  frcv(ji,jj,jpr_taum) = 0.5 * SQRT( zzx * zzx + zzy * zzy )
               END DO
            END DO
            CALL lbc_lnk( frcv(:,:,jpr_taum), 'T', 1. )
            llnewtau = .TRUE.
         ELSE
            llnewtau = .FALSE.
         ENDIF
      ELSE
         llnewtau = nrcvinfo(jpr_taum) == OASIS_Rcv
         ! Stress module can be negative when received (interpolation problem)
         IF( llnewtau ) THEN 
            DO jj = 1, jpj
               DO ji = 1, jpi 
                  frcv(ji,jj,jpr_taum) = MAX( 0.0e0, frcv(ji,jj,jpr_taum) )
               END DO
            END DO
         ENDIF
      ENDIF
      
      !                                                      ! ========================= !
      !                                                      !      10 m wind speed      !   (wndm)
      !                                                      ! ========================= !
      !
      IF( .NOT. srcv(jpr_w10m)%laction ) THEN                    ! compute wind spreed from wind stress module if not received  
         ! => need to be done only when taumod was changed
         IF( llnewtau ) THEN 
            zcoef = 1. / ( zrhoa * zcdrag ) 
!CDIR NOVERRCHK
            DO jj = 1, jpj
!CDIR NOVERRCHK
               DO ji = 1, jpi 
                  frcv(ji,jj,jpr_w10m) = SQRT( frcv(ji,jj,jpr_taum) * zcoef )
               END DO
            END DO
         ENDIF
      ENDIF

      ! u(v)tau and taum will be modified by ice model (wndm will be changed by PISCES)
      ! -> need to be reset before each call of the ice/fsbc      
      IF( MOD( kt-1, k_fsbc ) == 0 ) THEN
         !
         utau(:,:) = frcv(:,:,jpr_otx1)                   
         vtau(:,:) = frcv(:,:,jpr_oty1)
         taum(:,:) = frcv(:,:,jpr_taum)
         wndm(:,:) = frcv(:,:,jpr_w10m)
         CALL iom_put( "taum_oce", taum )   ! output wind stress module
         !  
      ENDIF
      !                                                      ! ========================= !
      IF( k_ice <= 1 ) THEN                                  !  heat & freshwater fluxes ! (Ocean only case)
         !                                                   ! ========================= !
         !
         !                                                       ! non solar heat flux over the ocean (qns)
         IF( srcv(jpr_qnsoce)%laction )   qns(:,:) = frcv(:,:,jpr_qnsoce)
         IF( srcv(jpr_qnsmix)%laction )   qns(:,:) = frcv(:,:,jpr_qnsmix)  
         ! add the latent heat of solid precip. melting
         IF( srcv(jpr_snow  )%laction )   qns(:,:) = qns(:,:) - frcv(:,:,jpr_snow) * lfus              

         !                                                       ! solar flux over the ocean          (qsr)
         IF( srcv(jpr_qsroce)%laction )   qsr(:,:) = frcv(:,:,jpr_qsroce) 
         IF( srcv(jpr_qsrmix)%laction )   qsr(:,:) = frcv(:,:,jpr_qsrmix)
         IF( ln_dm2dc )   qsr(:,:) = sbc_dcy( qsr )                           ! modify qsr to include the diurnal cycle
         !
         !                                                       ! total freshwater fluxes over the ocean (emp, emps)
         SELECT CASE( TRIM( cn_rcv_emp ) )                                    ! evaporation - precipitation
         CASE( 'conservative' )
            emp(:,:) = frcv(:,:,jpr_tevp) - ( frcv(:,:,jpr_rain) + frcv(:,:,jpr_snow) )
         CASE( 'oce only', 'oce and ice' )
            emp(:,:) = frcv(:,:,jpr_oemp)
         CASE default
            CALL ctl_stop( 'sbc_cpl_rcv: wrong definition of cn_rcv_emp' )
         END SELECT
         !
         !                                                        ! runoffs and calving (added in emp)
         IF( srcv(jpr_rnf)%laction )   emp(:,:) = emp(:,:) - frcv(:,:,jpr_rnf)        
         IF( srcv(jpr_cal)%laction )   emp(:,:) = emp(:,:) - frcv(:,:,jpr_cal)
         !
!!gm :  this seems to be internal cooking, not sure to need that in a generic interface 
!!gm                                       at least should be optional...
!!         IF( TRIM( cn_rcv_rnf ) == 'coupled' ) THEN     ! add to the total freshwater budget
!!            ! remove negative runoff
!!            zcumulpos = SUM( MAX( frcv(:,:,jpr_rnf), 0.e0 ) * e1t(:,:) * e2t(:,:) * tmask_i(:,:) ) 
!!            zcumulneg = SUM( MIN( frcv(:,:,jpr_rnf), 0.e0 ) * e1t(:,:) * e2t(:,:) * tmask_i(:,:) )
!!            IF( lk_mpp )   CALL mpp_sum( zcumulpos )   ! sum over the global domain
!!            IF( lk_mpp )   CALL mpp_sum( zcumulneg ) 
!!            IF( zcumulpos /= 0. ) THEN                 ! distribute negative runoff on positive runoff grid points
!!               zcumulneg = 1.e0 + zcumulneg / zcumulpos
!!               frcv(:,:,jpr_rnf) = MAX( frcv(:,:,jpr_rnf), 0.e0 ) * zcumulneg
!!            ENDIF     
!!            ! add runoff to e-p 
!!            emp(:,:) = emp(:,:) - frcv(:,:,jpr_rnf)
!!         ENDIF
!!gm  end of internal cooking
         !
         emps(:,:) = emp(:,:)                                        ! concentration/dilution = emp
  
         !                                                           ! 10 m wind speed
         IF( srcv(jpr_w10m)%laction )   wndm(:,:) = frcv(:,:,jpr_w10m)
         !
#if defined  key_cpl_carbon_cycle
         !                                                              ! atmosph. CO2 (ppm)
         IF( srcv(jpr_co2)%laction )   atm_co2(:,:) = frcv(:,:,jpr_co2)
#endif

      ENDIF
      !
      IF( wrk_not_released(2, 1,2) )   CALL ctl_stop('sbc_cpl_rcv: failed to release workspace arrays')
      !
   END SUBROUTINE sbc_cpl_rcv
   

   SUBROUTINE sbc_cpl_ice_tau( p_taui, p_tauj )     
      !!----------------------------------------------------------------------
      !!             ***  ROUTINE sbc_cpl_ice_tau  ***
      !!
      !! ** Purpose :   provide the stress over sea-ice in coupled mode 
      !!
      !! ** Method  :   transform the received stress from the atmosphere into
      !!             an atmosphere-ice stress in the (i,j) ocean referencial
      !!             and at the velocity point of the sea-ice model (cp_ice_msh):
      !!                'C'-grid : i- (j-) components given at U- (V-) point 
      !!                'I'-grid : B-grid lower-left corner: both components given at I-point 
      !!
      !!                The received stress are :
      !!                 - defined by 3 components (if cartesian coordinate)
      !!                        or by 2 components (if spherical)
      !!                 - oriented along geographical   coordinate (if eastward-northward)
      !!                        or  along the local grid coordinate (if local grid)
      !!                 - given at U- and V-point, resp.   if received on 2 grids
      !!                        or at a same point (T or I) if received on 1 grid
      !!                Therefore and if necessary, they are successively 
      !!             processed in order to obtain them 
      !!                 first  as  2 components on the sphere 
      !!                 second as  2 components oriented along the local grid
      !!                 third  as  2 components on the cp_ice_msh point 
      !!
      !!                In 'oce and ice' case, only one vector stress field 
      !!             is received. It has already been processed in sbc_cpl_rcv
      !!             so that it is now defined as (i,j) components given at U-
      !!             and V-points, respectively. Therefore, here only the third
      !!             transformation is done and only if the ice-grid is a 'I'-grid. 
      !!
      !! ** Action  :   return ptau_i, ptau_j, the stress over the ice at cp_ice_msh point
      !!----------------------------------------------------------------------
      USE wrk_nemo, ONLY:   wrk_in_use, wrk_not_released
      USE wrk_nemo, ONLY:   ztx => wrk_2d_1 , zty => wrk_2d_2
      !!
      REAL(wp), INTENT(out), DIMENSION(:,:) ::   p_taui   ! i- & j-components of atmos-ice stress [N/m2]
      REAL(wp), INTENT(out), DIMENSION(:,:) ::   p_tauj   ! at I-point (B-grid) or U & V-point (C-grid)
      !!
      INTEGER ::   ji, jj                          ! dummy loop indices
      INTEGER ::   itx                             ! index of taux over ice
      !!----------------------------------------------------------------------

      IF( wrk_in_use(2, 1,2) ) THEN
         CALL ctl_stop('sbc_cpl_ice_tau: requested workspace arrays unavailable')   ;   RETURN
      ENDIF

      IF( srcv(jpr_itx1)%laction ) THEN   ;   itx =  jpr_itx1   
      ELSE                                ;   itx =  jpr_otx1
      ENDIF

      ! do something only if we just received the stress from atmosphere
      IF(  nrcvinfo(itx) == OASIS_Rcv ) THEN

         !                                                      ! ======================= !
         IF( srcv(jpr_itx1)%laction ) THEN                      !   ice stress received   !
            !                                                   ! ======================= !
            !  
            IF( TRIM( cn_rcv_tau(2) ) == 'cartesian' ) THEN            ! 2 components on the sphere
               !                                                       ! (cartesian to spherical -> 3 to 2 components)
               CALL geo2oce( frcv(:,:,jpr_itx1), frcv(:,:,jpr_ity1), frcv(:,:,jpr_itz1),   &
                  &          srcv(jpr_itx1)%clgrid, ztx, zty )
               frcv(:,:,jpr_itx1) = ztx(:,:)   ! overwrite 1st comp. on the 1st grid
               frcv(:,:,jpr_itx1) = zty(:,:)   ! overwrite 2nd comp. on the 1st grid
               !
               IF( srcv(jpr_itx2)%laction ) THEN
                  CALL geo2oce( frcv(:,:,jpr_itx2), frcv(:,:,jpr_ity2), frcv(:,:,jpr_itz2),   &
                     &          srcv(jpr_itx2)%clgrid, ztx, zty )
                  frcv(:,:,jpr_itx2) = ztx(:,:)   ! overwrite 1st comp. on the 2nd grid
                  frcv(:,:,jpr_ity2) = zty(:,:)   ! overwrite 2nd comp. on the 2nd grid
               ENDIF
               !
            ENDIF
            !
            IF( TRIM( cn_rcv_tau(3) ) == 'eastward-northward' ) THEN   ! 2 components oriented along the local grid
               !                                                       ! (geographical to local grid -> rotate the components)
               CALL rot_rep( frcv(:,:,jpr_itx1), frcv(:,:,jpr_ity1), srcv(jpr_itx1)%clgrid, 'en->i', ztx )   
               frcv(:,:,jpr_itx1) = ztx(:,:)      ! overwrite 1st component on the 1st grid
               IF( srcv(jpr_itx2)%laction ) THEN
                  CALL rot_rep( frcv(:,:,jpr_itx2), frcv(:,:,jpr_ity2), srcv(jpr_itx2)%clgrid, 'en->j', zty )   
               ELSE
                  CALL rot_rep( frcv(:,:,jpr_itx1), frcv(:,:,jpr_ity1), srcv(jpr_itx1)%clgrid, 'en->j', zty )  
               ENDIF
               frcv(:,:,jpr_ity1) = zty(:,:)      ! overwrite 2nd component on the 1st grid
            ENDIF
            !                                                   ! ======================= !
         ELSE                                                   !     use ocean stress    !
            !                                                   ! ======================= !
            frcv(:,:,jpr_itx1) = frcv(:,:,jpr_otx1)
            frcv(:,:,jpr_ity1) = frcv(:,:,jpr_oty1)
            !
         ENDIF

         !                                                      ! ======================= !
         !                                                      !     put on ice grid     !
         !                                                      ! ======================= !
         !    
         !                                                  j+1   j     -----V---F
         ! ice stress on ice velocity point (cp_ice_msh)                 !       |
         ! (C-grid ==>(U,V) or B-grid ==> I or F)                 j      |   T   U
         !                                                               |       |
         !                                                   j    j-1   -I-------|
         !                                               (for I)         |       |
         !                                                              i-1  i   i
         !                                                               i      i+1 (for I)
         SELECT CASE ( cp_ice_msh )
            !
         CASE( 'I' )                                         ! B-grid ==> I
            SELECT CASE ( srcv(jpr_itx1)%clgrid )
            CASE( 'U' )
               DO jj = 2, jpjm1                                   ! (U,V) ==> I
                  DO ji = 2, jpim1   ! NO vector opt.
                     p_taui(ji,jj) = 0.5 * ( frcv(ji-1,jj  ,jpr_itx1) + frcv(ji-1,jj-1,jpr_itx1) )
                     p_tauj(ji,jj) = 0.5 * ( frcv(ji  ,jj-1,jpr_ity1) + frcv(ji-1,jj-1,jpr_ity1) )
                  END DO
               END DO
            CASE( 'F' )
               DO jj = 2, jpjm1                                   ! F ==> I
                  DO ji = 2, jpim1   ! NO vector opt.
                     p_taui(ji,jj) = frcv(ji-1,jj-1,jpr_itx1) 
                     p_tauj(ji,jj) = frcv(ji-1,jj-1,jpr_ity1)  
                  END DO
               END DO
            CASE( 'T' )
               DO jj = 2, jpjm1                                   ! T ==> I
                  DO ji = 2, jpim1   ! NO vector opt.
                     p_taui(ji,jj) = 0.25 * ( frcv(ji,jj  ,jpr_itx1) + frcv(ji-1,jj  ,jpr_itx1)   &
                        &                   + frcv(ji,jj-1,jpr_itx1) + frcv(ji-1,jj-1,jpr_itx1) ) 
                     p_tauj(ji,jj) = 0.25 * ( frcv(ji,jj  ,jpr_ity1) + frcv(ji-1,jj  ,jpr_ity1)   &
                        &                   + frcv(ji,jj-1,jpr_ity1) + frcv(ji-1,jj-1,jpr_ity1) )
                  END DO
               END DO
            CASE( 'I' )
               p_taui(:,:) = frcv(:,:,jpr_itx1)                   ! I ==> I
               p_tauj(:,:) = frcv(:,:,jpr_ity1)
            END SELECT
            IF( srcv(jpr_itx1)%clgrid /= 'I' ) THEN 
               CALL lbc_lnk( p_taui, 'I',  -1. )   ;   CALL lbc_lnk( p_tauj, 'I',  -1. )
            ENDIF
            !
         CASE( 'F' )                                         ! B-grid ==> F
            SELECT CASE ( srcv(jpr_itx1)%clgrid )
            CASE( 'U' )
               DO jj = 2, jpjm1                                   ! (U,V) ==> F
                  DO ji = fs_2, fs_jpim1   ! vector opt.
                     p_taui(ji,jj) = 0.5 * ( frcv(ji,jj,jpr_itx1) + frcv(ji  ,jj+1,jpr_itx1) )
                     p_tauj(ji,jj) = 0.5 * ( frcv(ji,jj,jpr_ity1) + frcv(ji+1,jj  ,jpr_ity1) )
                  END DO
               END DO
            CASE( 'I' )
               DO jj = 2, jpjm1                                   ! I ==> F
                  DO ji = 2, jpim1   ! NO vector opt.
                     p_taui(ji,jj) = frcv(ji+1,jj+1,jpr_itx1) 
                     p_tauj(ji,jj) = frcv(ji+1,jj+1,jpr_ity1)  
                  END DO
               END DO
            CASE( 'T' )
               DO jj = 2, jpjm1                                   ! T ==> F
                  DO ji = 2, jpim1   ! NO vector opt.
                     p_taui(ji,jj) = 0.25 * ( frcv(ji,jj  ,jpr_itx1) + frcv(ji+1,jj  ,jpr_itx1)   &
                        &                   + frcv(ji,jj+1,jpr_itx1) + frcv(ji+1,jj+1,jpr_itx1) ) 
                     p_tauj(ji,jj) = 0.25 * ( frcv(ji,jj  ,jpr_ity1) + frcv(ji+1,jj  ,jpr_ity1)   &
                        &                   + frcv(ji,jj+1,jpr_ity1) + frcv(ji+1,jj+1,jpr_ity1) )
                  END DO
               END DO
            CASE( 'F' )
               p_taui(:,:) = frcv(:,:,jpr_itx1)                   ! F ==> F
               p_tauj(:,:) = frcv(:,:,jpr_ity1)
            END SELECT
            IF( srcv(jpr_itx1)%clgrid /= 'F' ) THEN 
               CALL lbc_lnk( p_taui, 'F',  -1. )   ;   CALL lbc_lnk( p_tauj, 'F',  -1. )
            ENDIF
            !
         CASE( 'C' )                                         ! C-grid ==> U,V
            SELECT CASE ( srcv(jpr_itx1)%clgrid )
            CASE( 'U' )
               p_taui(:,:) = frcv(:,:,jpr_itx1)                   ! (U,V) ==> (U,V)
               p_tauj(:,:) = frcv(:,:,jpr_ity1)
            CASE( 'F' )
               DO jj = 2, jpjm1                                   ! F ==> (U,V)
                  DO ji = fs_2, fs_jpim1   ! vector opt.
                     p_taui(ji,jj) = 0.5 * ( frcv(ji,jj,jpr_itx1) + frcv(ji  ,jj-1,jpr_itx1) )
                     p_tauj(ji,jj) = 0.5 * ( frcv(ji,jj,jpr_ity1) + frcv(ji-1,jj  ,jpr_ity1) )
                  END DO
               END DO
            CASE( 'T' )
               DO jj = 2, jpjm1                                   ! T ==> (U,V)
                  DO ji = fs_2, fs_jpim1   ! vector opt.
                     p_taui(ji,jj) = 0.5 * ( frcv(ji+1,jj  ,jpr_itx1) + frcv(ji,jj,jpr_itx1) )
                     p_tauj(ji,jj) = 0.5 * ( frcv(ji  ,jj+1,jpr_ity1) + frcv(ji,jj,jpr_ity1) )
                  END DO
               END DO
            CASE( 'I' )
               DO jj = 2, jpjm1                                   ! I ==> (U,V)
                  DO ji = 2, jpim1   ! NO vector opt.
                     p_taui(ji,jj) = 0.5 * ( frcv(ji+1,jj+1,jpr_itx1) + frcv(ji+1,jj  ,jpr_itx1) )
                     p_tauj(ji,jj) = 0.5 * ( frcv(ji+1,jj+1,jpr_ity1) + frcv(ji  ,jj+1,jpr_ity1) )
                  END DO
               END DO
            END SELECT
            IF( srcv(jpr_itx1)%clgrid /= 'U' ) THEN 
               CALL lbc_lnk( p_taui, 'U',  -1. )   ;   CALL lbc_lnk( p_tauj, 'V',  -1. )
            ENDIF
         END SELECT

         !!gm Should be useless as sbc_cpl_ice_tau only called at coupled frequency
         ! The receive stress are transformed such that in all case frcv(:,:,jpr_itx1), frcv(:,:,jpr_ity1)
         ! become the i-component and j-component of the stress at the right grid point 
         !!gm  frcv(:,:,jpr_itx1) = p_taui(:,:)
         !!gm  frcv(:,:,jpr_ity1) = p_tauj(:,:)
         !!gm
      ENDIF
      !   
      IF( wrk_not_released(2, 1,2) )   CALL ctl_stop('sbc_cpl_ice_tau: failed to release workspace arrays')
      !
   END SUBROUTINE sbc_cpl_ice_tau
   

   SUBROUTINE sbc_cpl_ice_flx( p_frld  ,                                  &
      &                        pqns_tot, pqns_ice, pqsr_tot , pqsr_ice,   &
      &                        pemp_tot, pemp_ice, pdqns_ice, psprecip,   &
      &                        palbi   , psst    , pist                 )
      !!----------------------------------------------------------------------
      !!             ***  ROUTINE sbc_cpl_ice_flx_rcv  ***
      !!
      !! ** Purpose :   provide the heat and freshwater fluxes of the 
      !!              ocean-ice system.
      !!
      !! ** Method  :   transform the fields received from the atmosphere into
      !!             surface heat and fresh water boundary condition for the 
      !!             ice-ocean system. The following fields are provided:
      !!              * total non solar, solar and freshwater fluxes (qns_tot, 
      !!             qsr_tot and emp_tot) (total means weighted ice-ocean flux)
      !!             NB: emp_tot include runoffs and calving.
      !!              * fluxes over ice (qns_ice, qsr_ice, emp_ice) where
      !!             emp_ice = sublimation - solid precipitation as liquid
      !!             precipitation are re-routed directly to the ocean and 
      !!             runoffs and calving directly enter the ocean.
      !!              * solid precipitation (sprecip), used to add to qns_tot 
      !!             the heat lost associated to melting solid precipitation
      !!             over the ocean fraction.
      !!       ===>> CAUTION here this changes the net heat flux received from
      !!             the atmosphere
      !!
      !!             N.B. - fields over sea-ice are passed in argument so that
      !!                 the module can be compile without sea-ice.
      !!                  - the fluxes have been separated from the stress as
      !!                 (a) they are updated at each ice time step compare to
      !!                 an update at each coupled time step for the stress, and
      !!                 (b) the conservative computation of the fluxes over the
      !!                 sea-ice area requires the knowledge of the ice fraction
      !!                 after the ice advection and before the ice thermodynamics,
      !!                 so that the stress is updated before the ice dynamics
      !!                 while the fluxes are updated after it.
      !!
      !! ** Action  :   update at each nf_ice time step:
      !!                   pqns_tot, pqsr_tot  non-solar and solar total heat fluxes
      !!                   pqns_ice, pqsr_ice  non-solar and solar heat fluxes over the ice
      !!                   pemp_tot            total evaporation - precipitation(liquid and solid) (-runoff)(-calving)
      !!                   pemp_ice            ice sublimation - solid precipitation over the ice
      !!                   pdqns_ice           d(non-solar heat flux)/d(Temperature) over the ice
      !!                   sprecip             solid precipitation over the ocean  
      !!----------------------------------------------------------------------
      USE wrk_nemo, ONLY:   wrk_in_use, wrk_not_released
      USE wrk_nemo, ONLY:   zcptn  => wrk_2d_1   ! rcp * tn(:,:,1)
      USE wrk_nemo, ONLY:   ztmp   => wrk_2d_2   ! temporary array
      USE wrk_nemo, ONLY:   zsnow  => wrk_2d_3   ! snow precipitation 
      USE wrk_nemo, ONLY:   zicefr => wrk_3d_4   ! ice fraction 
      !!
      REAL(wp), INTENT(in   ), DIMENSION(:,:,:) ::   p_frld     ! lead fraction                [0 to 1]
      REAL(wp), INTENT(  out), DIMENSION(:,:  ) ::   pqns_tot   ! total non solar heat flux    [W/m2]
      REAL(wp), INTENT(  out), DIMENSION(:,:,:) ::   pqns_ice   ! ice   non solar heat flux    [W/m2]
      REAL(wp), INTENT(  out), DIMENSION(:,:  ) ::   pqsr_tot   ! total     solar heat flux    [W/m2]
      REAL(wp), INTENT(  out), DIMENSION(:,:,:) ::   pqsr_ice   ! ice       solar heat flux    [W/m2]
      REAL(wp), INTENT(  out), DIMENSION(:,:  ) ::   pemp_tot   ! total     freshwater budget        [Kg/m2/s]
      REAL(wp), INTENT(  out), DIMENSION(:,:  ) ::   pemp_ice   ! solid freshwater budget over ice   [Kg/m2/s]
      REAL(wp), INTENT(  out), DIMENSION(:,:  ) ::   psprecip   ! Net solid precipitation (=emp_ice) [Kg/m2/s]
      REAL(wp), INTENT(  out), DIMENSION(:,:,:) ::   pdqns_ice  ! d(Q non solar)/d(Temperature) over ice
      ! optional arguments, used only in 'mixed oce-ice' case
      REAL(wp), INTENT(in   ), DIMENSION(:,:,:), OPTIONAL ::   palbi   ! ice albedo 
      REAL(wp), INTENT(in   ), DIMENSION(:,:  ), OPTIONAL ::   psst    ! sea surface temperature     [Celcius]
      REAL(wp), INTENT(in   ), DIMENSION(:,:,:), OPTIONAL ::   pist    ! ice surface temperature     [Kelvin]
      !!
      INTEGER ::   ji, jj           ! dummy loop indices
      INTEGER ::   isec, info       ! temporary integer
      REAL(wp)::   zcoef, ztsurf    ! temporary scalar
      !!----------------------------------------------------------------------

      IF( wrk_in_use(2, 1,2,3) .OR. wrk_in_use(3, 4) ) THEN
         CALL ctl_stop('sbc_cpl_ice_flx: requested workspace arrays unavailable')   ;   RETURN
      ENDIF

      zicefr(:,:,1) = 1.- p_frld(:,:,1)
      IF( lk_diaar5 )   zcptn(:,:) = rcp * tn(:,:,1)
      !
      !                                                      ! ========================= !
      !                                                      !    freshwater budget      !   (emp)
      !                                                      ! ========================= !
      !
      !                                                           ! total Precipitations - total Evaporation (emp_tot)
      !                                                           ! solid precipitation  - sublimation       (emp_ice)
      !                                                           ! solid Precipitation                      (sprecip)
      SELECT CASE( TRIM( cn_rcv_emp ) )
      CASE( 'conservative'  )   ! received fields: jpr_rain, jpr_snow, jpr_ievp, jpr_tevp
         pemp_tot(:,:) = frcv(:,:,jpr_tevp) - frcv(:,:,jpr_rain) - frcv(:,:,jpr_snow)
         pemp_ice(:,:) = frcv(:,:,jpr_ievp) - frcv(:,:,jpr_snow)
         zsnow   (:,:) = frcv(:,:,jpr_snow)
                           CALL iom_put( 'rain'         , frcv(:,:,jpr_rain)              )   ! liquid precipitation 
         IF( lk_diaar5 )   CALL iom_put( 'hflx_rain_cea', frcv(:,:,jpr_rain) * zcptn(:,:) )   ! heat flux from liq. precip. 
         ztmp(:,:) = frcv(:,:,jpr_tevp) - frcv(:,:,jpr_ievp) * zicefr(:,:,1)
                           CALL iom_put( 'evap_ao_cea'  , ztmp                            )   ! ice-free oce evap (cell average)
         IF( lk_diaar5 )   CALL iom_put( 'hflx_evap_cea', ztmp(:,:         ) * zcptn(:,:) )   ! heat flux from from evap (cell ave)
      CASE( 'oce and ice'   )   ! received fields: jpr_sbpr, jpr_semp, jpr_oemp
         pemp_tot(:,:) = p_frld(:,:,1) * frcv(:,:,jpr_oemp) + zicefr(:,:,1) * frcv(:,:,jpr_sbpr) 
         pemp_ice(:,:) = frcv(:,:,jpr_semp)
         zsnow   (:,:) = - frcv(:,:,jpr_semp) + frcv(:,:,jpr_ievp)
      END SELECT
      psprecip(:,:) = - pemp_ice(:,:)
      CALL iom_put( 'snowpre'    , zsnow                               )   ! Snow
      CALL iom_put( 'snow_ao_cea', zsnow(:,:         ) * p_frld(:,:,1) )   ! Snow        over ice-free ocean  (cell average)
      CALL iom_put( 'snow_ai_cea', zsnow(:,:         ) * zicefr(:,:,1) )   ! Snow        over sea-ice         (cell average)
      CALL iom_put( 'subl_ai_cea', frcv (:,:,jpr_ievp) * zicefr(:,:,1) )   ! Sublimation over sea-ice         (cell average)
      !   
      !                                                           ! runoffs and calving (put in emp_tot)
      IF( srcv(jpr_rnf)%laction ) THEN 
         pemp_tot(:,:) = pemp_tot(:,:) - frcv(:,:,jpr_rnf)
                           CALL iom_put( 'runoffs'      , frcv(:,:,jpr_rnf )              )   ! rivers
         IF( lk_diaar5 )   CALL iom_put( 'hflx_rnf_cea' , frcv(:,:,jpr_rnf ) * zcptn(:,:) )   ! heat flux from rivers
      ENDIF
      IF( srcv(jpr_cal)%laction ) THEN 
         pemp_tot(:,:) = pemp_tot(:,:) - frcv(:,:,jpr_cal)
         CALL iom_put( 'calving', frcv(:,:,jpr_cal) )
      ENDIF
      !
!!gm :  this seems to be internal cooking, not sure to need that in a generic interface 
!!gm                                       at least should be optional...
!!       ! remove negative runoff                            ! sum over the global domain
!!       zcumulpos = SUM( MAX( frcv(:,:,jpr_rnf), 0.e0 ) * e1t(:,:) * e2t(:,:) * tmask_i(:,:) ) 
!!       zcumulneg = SUM( MIN( frcv(:,:,jpr_rnf), 0.e0 ) * e1t(:,:) * e2t(:,:) * tmask_i(:,:) )
!!       IF( lk_mpp )   CALL mpp_sum( zcumulpos )
!!       IF( lk_mpp )   CALL mpp_sum( zcumulneg ) 
!!       IF( zcumulpos /= 0. ) THEN                          ! distribute negative runoff on positive runoff grid points
!!          zcumulneg = 1.e0 + zcumulneg / zcumulpos
!!          frcv(:,:,jpr_rnf) = MAX( frcv(:,:,jpr_rnf), 0.e0 ) * zcumulneg
!!       ENDIF     
!!       pemp_tot(:,:) = pemp_tot(:,:) - frcv(:,:,jpr_rnf)   ! add runoff to e-p 
!!
!!gm  end of internal cooking


      !                                                      ! ========================= !
      SELECT CASE( TRIM( cn_rcv_qns ) )                      !   non solar heat fluxes   !   (qns)
      !                                                      ! ========================= !
      CASE( 'conservative' )                                      ! the required fields are directly provided
         pqns_tot(:,:  ) = frcv(:,:,jpr_qnsmix)
         pqns_ice(:,:,1) = frcv(:,:,jpr_qnsice)
      CASE( 'oce and ice' )       ! the total flux is computed from ocean and ice fluxes
         pqns_tot(:,:  ) =  p_frld(:,:,1) * frcv(:,:,jpr_qnsoce) + zicefr(:,:,1) * frcv(:,:,jpr_qnsice)
         pqns_ice(:,:,1) =  frcv(:,:,jpr_qnsice)
      CASE( 'mixed oce-ice' )     ! the ice flux is cumputed from the total flux, the SST and ice informations
         pqns_tot(:,:  ) = frcv(:,:,jpr_qnsmix)
         pqns_ice(:,:,1) = frcv(:,:,jpr_qnsmix)    &
            &            + frcv(:,:,jpr_dqnsdt) * ( pist(:,:,1) - ( (rt0 + psst(:,:  ) ) * p_frld(:,:,1)   &
            &                                                   +          pist(:,:,1)   * zicefr(:,:,1) ) )
      END SELECT
      ztmp(:,:) = p_frld(:,:,1) * zsnow(:,:) * lfus               ! add the latent heat of solid precip. melting
      pqns_tot(:,:) = pqns_tot(:,:) - ztmp(:,:)                   ! over free ocean 
      IF( lk_diaar5 )   CALL iom_put( 'hflx_snow_cea', ztmp + zsnow(:,:) * zcptn(:,:) )   ! heat flux from snow (cell average)
!!gm
!!    currently it is taken into account in leads budget but not in the qns_tot, and thus not in 
!!    the flux that enter the ocean....
!!    moreover 1 - it is not diagnose anywhere.... 
!!             2 - it is unclear for me whether this heat lost is taken into account in the atmosphere or not...
!!
!! similar job should be done for snow and precipitation temperature
      !                                     
      IF( srcv(jpr_cal)%laction ) THEN                            ! Iceberg melting 
         ztmp(:,:) = frcv(:,:,jpr_cal) * lfus                     ! add the latent heat of iceberg melting 
         pqns_tot(:,:) = pqns_tot(:,:) - ztmp(:,:)
         IF( lk_diaar5 )   CALL iom_put( 'hflx_cal_cea', ztmp + frcv(:,:,jpr_cal) * zcptn(:,:) )   ! heat flux from calving
      ENDIF

      !                                                      ! ========================= !
      SELECT CASE( TRIM( cn_rcv_qsr ) )                      !      solar heat fluxes    !   (qsr)
      !                                                      ! ========================= !
      CASE( 'conservative' )
         pqsr_tot(:,:  ) = frcv(:,:,jpr_qsrmix)
         pqsr_ice(:,:,1) = frcv(:,:,jpr_qsrice)
      CASE( 'oce and ice' )
         pqsr_tot(:,:  ) =  p_frld(:,:,1) * frcv(:,:,jpr_qsroce) + zicefr(:,:,1) * frcv(:,:,jpr_qsrice)
         pqsr_ice(:,:,1) =  frcv(:,:,jpr_qsrice)
      CASE( 'mixed oce-ice' )
         pqsr_tot(:,:  ) = frcv(:,:,jpr_qsrmix)
!       Create solar heat flux over ice using incoming solar heat flux and albedos
!       ( see OASIS3 user guide, 5th edition, p39 )
         pqsr_ice(:,:,1) = frcv(:,:,jpr_qsrmix) * ( 1.- palbi(:,:,1) )   &
            &            / (  1.- ( albedo_oce_mix(:,:  ) * p_frld(:,:,1)   &
            &                     + palbi         (:,:,1) * zicefr(:,:,1) ) )
      END SELECT
      IF( ln_dm2dc ) THEN   ! modify qsr to include the diurnal cycle
         pqsr_tot(:,:  ) = sbc_dcy( pqsr_tot(:,:  ) )
         pqsr_ice(:,:,1) = sbc_dcy( pqsr_ice(:,:,1) )
      ENDIF

      SELECT CASE( TRIM( cn_rcv_dqnsdt ) )
      CASE ('coupled')
          pdqns_ice(:,:,1) = frcv(:,:,jpr_dqnsdt)
      END SELECT

      IF( wrk_not_released(2, 1,2,3)  .OR.   &
          wrk_not_released(3, 4)      )   CALL ctl_stop('sbc_cpl_ice_flx: failed to release workspace arrays')
      !
   END SUBROUTINE sbc_cpl_ice_flx
   
   
   SUBROUTINE sbc_cpl_snd( kt )
      !!----------------------------------------------------------------------
      !!             ***  ROUTINE sbc_cpl_snd  ***
      !!
      !! ** Purpose :   provide the ocean-ice informations to the atmosphere
      !!
      !! ** Method  :   send to the atmosphere through a call to cpl_prism_snd
      !!              all the needed fields (as defined in sbc_cpl_init)
      !!----------------------------------------------------------------------
      USE wrk_nemo, ONLY:   wrk_in_use, wrk_not_released
      USE wrk_nemo, ONLY:   zfr_l => wrk_2d_1   ! 1. - fr_i(:,:)
      USE wrk_nemo, ONLY:   ztmp1 => wrk_2d_2 , ztmp2 => wrk_2d_3
      USE wrk_nemo, ONLY:   zotx1 => wrk_2d_4 , zoty1 => wrk_2d_5 , zotz1 => wrk_2d_6
      USE wrk_nemo, ONLY:   zitx1 => wrk_2d_7 , zity1 => wrk_2d_8 , zitz1 => wrk_2d_9
      !
      INTEGER, INTENT(in) ::   kt
      !
      INTEGER ::   ji, jj       ! dummy loop indices
      INTEGER ::   isec, info   ! local integer
      !!----------------------------------------------------------------------

      IF( wrk_in_use(2, 1,2,3,4,5,6,7,8,9) ) THEN
         CALL ctl_stop('sbc_cpl_snd: requested workspace arrays are unavailable')   ;   RETURN
      ENDIF

      isec = ( kt - nit000 ) * NINT(rdttra(1))        ! date of exchanges

      zfr_l(:,:) = 1.- fr_i(:,:)

      !                                                      ! ------------------------- !
      !                                                      !    Surface temperature    !   in Kelvin
      !                                                      ! ------------------------- !
      SELECT CASE( cn_snd_temperature)
      CASE( 'oce only'             )   ;   ztmp1(:,:) =   tn(:,:,1) + rt0
      CASE( 'weighted oce and ice' )   ;   ztmp1(:,:) = ( tn(:,:,1) + rt0 ) * zfr_l(:,:)   
                                           ztmp2(:,:) =   tn_ice(:,:,1)     *  fr_i(:,:)
      CASE( 'mixed oce-ice'        )   ;   ztmp1(:,:) = ( tn(:,:,1) + rt0 ) * zfr_l(:,:) + tn_ice(:,:,1) * fr_i(:,:)
      CASE default                     ;   CALL ctl_stop( 'sbc_cpl_snd: wrong definition of cn_snd_temperature' )
      END SELECT
      IF( ssnd(jps_toce)%laction )   CALL cpl_prism_snd( jps_toce, isec, ztmp1, info )
      IF( ssnd(jps_tice)%laction )   CALL cpl_prism_snd( jps_tice, isec, ztmp2, info )
      IF( ssnd(jps_tmix)%laction )   CALL cpl_prism_snd( jps_tmix, isec, ztmp1, info )
      !
      !                                                      ! ------------------------- !
      !                                                      !           Albedo          !
      !                                                      ! ------------------------- !
      IF( ssnd(jps_albice)%laction ) THEN                         ! ice 
         ztmp1(:,:) = alb_ice(:,:,1) * fr_i(:,:)
         CALL cpl_prism_snd( jps_albice, isec, ztmp1, info )
      ENDIF
      IF( ssnd(jps_albmix)%laction ) THEN                         ! mixed ice-ocean
         ztmp1(:,:) = albedo_oce_mix(:,:) * zfr_l(:,:) + alb_ice(:,:,1) * fr_i(:,:)
         CALL cpl_prism_snd( jps_albmix, isec, ztmp1, info )
      ENDIF
      !                                                      ! ------------------------- !
      !                                                      !  Ice fraction & Thickness ! 
      !                                                      ! ------------------------- !
      IF( ssnd(jps_fice)%laction )   CALL cpl_prism_snd( jps_fice, isec, fr_i                  , info )
      IF( ssnd(jps_hice)%laction )   CALL cpl_prism_snd( jps_hice, isec, hicif(:,:) * fr_i(:,:), info )
      IF( ssnd(jps_hsnw)%laction )   CALL cpl_prism_snd( jps_hsnw, isec, hsnif(:,:) * fr_i(:,:), info )
      !
#if defined key_cpl_carbon_cycle
      !                                                      ! ------------------------- !
      !                                                      !  CO2 flux from PISCES     ! 
      !                                                      ! ------------------------- !
      IF( ssnd(jps_co2)%laction )   CALL cpl_prism_snd( jps_co2, isec, oce_co2 , info )
      !
#endif
      IF( ssnd(jps_ocx1)%laction ) THEN                      !      Surface current      !
         !                                                   ! ------------------------- !
         !    
         !                                                  j+1   j     -----V---F
         ! surface velocity always sent from T point                     !       |
         !                                                        j      |   T   U
         !                                                               |       |
         !                                                   j    j-1   -I-------|
         !                                               (for I)         |       |
         !                                                              i-1  i   i
         !                                                               i      i+1 (for I)
         SELECT CASE( TRIM( cn_snd_crt(1) ) )
         CASE( 'oce only'             )      ! C-grid ==> T
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zotx1(ji,jj) = 0.5 * ( un(ji,jj,1) + un(ji-1,jj  ,1) )
                  zoty1(ji,jj) = 0.5 * ( vn(ji,jj,1) + vn(ji  ,jj-1,1) ) 
               END DO
            END DO
         CASE( 'weighted oce and ice' )   
            SELECT CASE ( cp_ice_msh )
            CASE( 'C' )                      ! Ocean and Ice on C-grid ==> T
               DO jj = 2, jpjm1
                  DO ji = fs_2, fs_jpim1   ! vector opt.
                     zotx1(ji,jj) = 0.5 * ( un   (ji,jj,1) + un   (ji-1,jj  ,1) ) * zfr_l(ji,jj)  
                     zoty1(ji,jj) = 0.5 * ( vn   (ji,jj,1) + vn   (ji  ,jj-1,1) ) * zfr_l(ji,jj)
                     zitx1(ji,jj) = 0.5 * ( u_ice(ji,jj  ) + u_ice(ji-1,jj    ) ) *  fr_i(ji,jj)
                     zity1(ji,jj) = 0.5 * ( v_ice(ji,jj  ) + v_ice(ji  ,jj-1  ) ) *  fr_i(ji,jj)
                  END DO
               END DO
            CASE( 'I' )                      ! Ocean on C grid, Ice on I-point (B-grid) ==> T
               DO jj = 2, jpjm1
                  DO ji = 2, jpim1   ! NO vector opt.
                     zotx1(ji,jj) = 0.5  * ( un(ji,jj,1)      + un(ji-1,jj  ,1) ) * zfr_l(ji,jj)  
                     zoty1(ji,jj) = 0.5  * ( vn(ji,jj,1)      + vn(ji  ,jj-1,1) ) * zfr_l(ji,jj)  
                     zitx1(ji,jj) = 0.25 * ( u_ice(ji+1,jj+1) + u_ice(ji,jj+1)                     &
                        &                  + u_ice(ji+1,jj  ) + u_ice(ji,jj  )  ) *  fr_i(ji,jj)
                     zity1(ji,jj) = 0.25 * ( v_ice(ji+1,jj+1) + v_ice(ji,jj+1)                     &
                        &                  + v_ice(ji+1,jj  ) + v_ice(ji,jj  )  ) *  fr_i(ji,jj)
                  END DO
               END DO
            CASE( 'F' )                      ! Ocean on C grid, Ice on F-point (B-grid) ==> T
               DO jj = 2, jpjm1
                  DO ji = 2, jpim1   ! NO vector opt.
                     zotx1(ji,jj) = 0.5  * ( un(ji,jj,1)      + un(ji-1,jj  ,1) ) * zfr_l(ji,jj)  
                     zoty1(ji,jj) = 0.5  * ( vn(ji,jj,1)      + vn(ji  ,jj-1,1) ) * zfr_l(ji,jj)  
                     zitx1(ji,jj) = 0.25 * ( u_ice(ji-1,jj-1) + u_ice(ji,jj-1)                     &
                        &                  + u_ice(ji-1,jj  ) + u_ice(ji,jj  )  ) *  fr_i(ji,jj)
                     zity1(ji,jj) = 0.25 * ( v_ice(ji-1,jj-1) + v_ice(ji,jj-1)                     &
                        &                  + v_ice(ji-1,jj  ) + v_ice(ji,jj  )  ) *  fr_i(ji,jj)
                  END DO
               END DO
            END SELECT
            CALL lbc_lnk( zitx1, 'T', -1. )   ;   CALL lbc_lnk( zity1, 'T', -1. )
         CASE( 'mixed oce-ice'        )
            SELECT CASE ( cp_ice_msh )
            CASE( 'C' )                      ! Ocean and Ice on C-grid ==> T
               DO jj = 2, jpjm1
                  DO ji = fs_2, fs_jpim1   ! vector opt.
                     zotx1(ji,jj) = 0.5 * ( un   (ji,jj,1) + un   (ji-1,jj  ,1) ) * zfr_l(ji,jj)   &
                        &         + 0.5 * ( u_ice(ji,jj  ) + u_ice(ji-1,jj    ) ) *  fr_i(ji,jj)
                     zoty1(ji,jj) = 0.5 * ( vn   (ji,jj,1) + vn   (ji  ,jj-1,1) ) * zfr_l(ji,jj)   &
                        &         + 0.5 * ( v_ice(ji,jj  ) + v_ice(ji  ,jj-1  ) ) *  fr_i(ji,jj)
                  END DO
               END DO
            CASE( 'I' )                      ! Ocean on C grid, Ice on I-point (B-grid) ==> T
               DO jj = 2, jpjm1
                  DO ji = 2, jpim1   ! NO vector opt.
                     zotx1(ji,jj) = 0.5  * ( un(ji,jj,1)      + un(ji-1,jj  ,1) ) * zfr_l(ji,jj)   &   
                        &         + 0.25 * ( u_ice(ji+1,jj+1) + u_ice(ji,jj+1)                     &
                        &                  + u_ice(ji+1,jj  ) + u_ice(ji,jj  )  ) *  fr_i(ji,jj)
                     zoty1(ji,jj) = 0.5  * ( vn(ji,jj,1)      + vn(ji  ,jj-1,1) ) * zfr_l(ji,jj)   & 
                        &         + 0.25 * ( v_ice(ji+1,jj+1) + v_ice(ji,jj+1)                     &
                        &                  + v_ice(ji+1,jj  ) + v_ice(ji,jj  )  ) *  fr_i(ji,jj)
                  END DO
               END DO
            CASE( 'F' )                      ! Ocean on C grid, Ice on F-point (B-grid) ==> T
               DO jj = 2, jpjm1
                  DO ji = 2, jpim1   ! NO vector opt.
                     zotx1(ji,jj) = 0.5  * ( un(ji,jj,1)      + un(ji-1,jj  ,1) ) * zfr_l(ji,jj)   &   
                        &         + 0.25 * ( u_ice(ji-1,jj-1) + u_ice(ji,jj-1)                     &
                        &                  + u_ice(ji-1,jj  ) + u_ice(ji,jj  )  ) *  fr_i(ji,jj)
                     zoty1(ji,jj) = 0.5  * ( vn(ji,jj,1)      + vn(ji  ,jj-1,1) ) * zfr_l(ji,jj)   & 
                        &         + 0.25 * ( v_ice(ji-1,jj-1) + v_ice(ji,jj-1)                     &
                        &                  + v_ice(ji-1,jj  ) + v_ice(ji,jj  )  ) *  fr_i(ji,jj)
                  END DO
               END DO
            END SELECT
         END SELECT
         CALL lbc_lnk( zotx1, 'T', -1. )   ;   CALL lbc_lnk( zoty1, 'T', -1. )
         !
         !
         IF( TRIM( cn_snd_crt(3) ) == 'eastward-northward' ) THEN             ! Rotation of the components
            !                                                                     ! Ocean component
            CALL rot_rep( zotx1, zoty1, ssnd(jps_ocx1)%clgrid, 'ij->e', ztmp1 )       ! 1st component 
            CALL rot_rep( zotx1, zoty1, ssnd(jps_ocx1)%clgrid, 'ij->n', ztmp2 )       ! 2nd component 
            zotx1(:,:) = ztmp1(:,:)                                                   ! overwrite the components 
            zoty1(:,:) = ztmp2(:,:)
            IF( ssnd(jps_ivx1)%laction ) THEN                                     ! Ice component
               CALL rot_rep( zitx1, zity1, ssnd(jps_ivx1)%clgrid, 'ij->e', ztmp1 )    ! 1st component 
               CALL rot_rep( zitx1, zity1, ssnd(jps_ivx1)%clgrid, 'ij->n', ztmp2 )    ! 2nd component 
               zitx1(:,:) = ztmp1(:,:)                                                ! overwrite the components 
               zity1(:,:) = ztmp2(:,:)
            ENDIF
         ENDIF
         !
         ! spherical coordinates to cartesian -> 2 components to 3 components
         IF( TRIM( cn_snd_crt(2) ) == 'cartesian' ) THEN
            ztmp1(:,:) = zotx1(:,:)                     ! ocean currents
            ztmp2(:,:) = zoty1(:,:)
            CALL oce2geo ( ztmp1, ztmp2, 'T', zotx1, zoty1, zotz1 )
            !
            IF( ssnd(jps_ivx1)%laction ) THEN           ! ice velocities
               ztmp1(:,:) = zitx1(:,:)
               ztmp1(:,:) = zity1(:,:)
               CALL oce2geo ( ztmp1, ztmp2, 'T', zitx1, zity1, zitz1 )
            ENDIF
         ENDIF
         !
         IF( ssnd(jps_ocx1)%laction )   CALL cpl_prism_snd( jps_ocx1, isec, zotx1, info )   ! ocean x current 1st grid
         IF( ssnd(jps_ocy1)%laction )   CALL cpl_prism_snd( jps_ocy1, isec, zoty1, info )   ! ocean y current 1st grid
         IF( ssnd(jps_ocz1)%laction )   CALL cpl_prism_snd( jps_ocz1, isec, zotz1, info )   ! ocean z current 1st grid
         !
         IF( ssnd(jps_ivx1)%laction )   CALL cpl_prism_snd( jps_ivx1, isec, zitx1, info )   ! ice   x current 1st grid
         IF( ssnd(jps_ivy1)%laction )   CALL cpl_prism_snd( jps_ivy1, isec, zity1, info )   ! ice   y current 1st grid
         IF( ssnd(jps_ivz1)%laction )   CALL cpl_prism_snd( jps_ivz1, isec, zitz1, info )   ! ice   z current 1st grid
         ! 
      ENDIF
      !
      IF( wrk_not_released(2, 1,2,3,4,5,6,7,8,9) )   CALL ctl_stop('sbc_cpl_snd: failed to release workspace arrays')
      !
   END SUBROUTINE sbc_cpl_snd

#else
   !!----------------------------------------------------------------------
   !!   Dummy module                                            NO coupling
   !!----------------------------------------------------------------------
   USE par_kind        ! kind definition
CONTAINS
   SUBROUTINE sbc_cpl_snd( kt )
      WRITE(*,*) 'sbc_cpl_snd: You should not have seen this print! error?', kt
   END SUBROUTINE sbc_cpl_snd
   !
   SUBROUTINE sbc_cpl_rcv( kt, k_fsbc, k_ice )     
      WRITE(*,*) 'sbc_cpl_snd: You should not have seen this print! error?', kt, k_fsbc, k_ice
   END SUBROUTINE sbc_cpl_rcv
   !
   SUBROUTINE sbc_cpl_ice_tau( p_taui, p_tauj )     
      REAL(wp), INTENT(out), DIMENSION(:,:) ::   p_taui   ! i- & j-components of atmos-ice stress [N/m2]
      REAL(wp), INTENT(out), DIMENSION(:,:) ::   p_tauj   ! at I-point (B-grid) or U & V-point (C-grid)
      p_taui(:,:) = 0.   ;   p_tauj(:,:) = 0. ! stupid definition to avoid warning message when compiling...
      WRITE(*,*) 'sbc_cpl_snd: You should not have seen this print! error?'
   END SUBROUTINE sbc_cpl_ice_tau
   !
   SUBROUTINE sbc_cpl_ice_flx( p_frld  ,                                  &
      &                        pqns_tot, pqns_ice, pqsr_tot , pqsr_ice,   &
      &                        pemp_tot, pemp_ice, pdqns_ice, psprecip,   &
      &                        palbi   , psst    , pist                )
      REAL(wp), INTENT(in   ), DIMENSION(:,:,:) ::   p_frld     ! lead fraction                [0 to 1]
      REAL(wp), INTENT(  out), DIMENSION(:,:  ) ::   pqns_tot   ! total non solar heat flux    [W/m2]
      REAL(wp), INTENT(  out), DIMENSION(:,:,:) ::   pqns_ice   ! ice   non solar heat flux    [W/m2]
      REAL(wp), INTENT(  out), DIMENSION(:,:  ) ::   pqsr_tot   ! total     solar heat flux    [W/m2]
      REAL(wp), INTENT(  out), DIMENSION(:,:,:) ::   pqsr_ice   ! ice       solar heat flux    [W/m2]
      REAL(wp), INTENT(  out), DIMENSION(:,:  ) ::   pemp_tot   ! total     freshwater budget  [Kg/m2/s]
      REAL(wp), INTENT(  out), DIMENSION(:,:  ) ::   pemp_ice   ! ice solid freshwater budget  [Kg/m2/s]
      REAL(wp), INTENT(  out), DIMENSION(:,:,:) ::   pdqns_ice  ! d(Q non solar)/d(Temperature) over ice
      REAL(wp), INTENT(  out), DIMENSION(:,:  ) ::   psprecip   ! solid precipitation          [Kg/m2/s]
      REAL(wp), INTENT(in   ), DIMENSION(:,:,:), OPTIONAL ::   palbi   ! ice albedo
      REAL(wp), INTENT(in   ), DIMENSION(:,:  ), OPTIONAL ::   psst    ! sea surface temperature      [Celcius]
      REAL(wp), INTENT(in   ), DIMENSION(:,:,:), OPTIONAL ::   pist    ! ice surface temperature      [Kelvin]
      WRITE(*,*) 'sbc_cpl_snd: You should not have seen this print! error?', p_frld(1,1,1), palbi(1,1,1), psst(1,1), pist(1,1,1) 
      ! stupid definition to avoid warning message when compiling...
      pqns_tot(:,:) = 0. ; pqns_ice(:,:,:) = 0. ; pdqns_ice(:,:,:) = 0.
      pqsr_tot(:,:) = 0. ; pqsr_ice(:,:,:) = 0. 
      pemp_tot(:,:) = 0. ; pemp_ice(:,:)   = 0. ; psprecip(:,:) = 0.
   END SUBROUTINE sbc_cpl_ice_flx
   
#endif

   !!======================================================================
END MODULE sbccpl
