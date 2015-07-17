MODULE trcini
   !!======================================================================
   !!                         ***  MODULE trcini  ***
   !! TOP :   Manage the passive tracer initialization
   !!======================================================================
   !! History :   -   ! 1991-03 (O. Marti)  original code
   !!            1.0  ! 2005-03 (O. Aumont, A. El Moussaoui) F90
   !!            2.0  ! 2005-10 (C. Ethe, G. Madec) revised architecture
   !!            4.0  ! 2011-01 (A. R. Porter, STFC Daresbury) dynamical allocation
   !!----------------------------------------------------------------------
#if defined key_top
   !!----------------------------------------------------------------------
   !!   'key_top'                                                TOP models
   !!----------------------------------------------------------------------
   !!   trc_init  :   Initialization for passive tracer
   !!   top_alloc :   allocate the TOP arrays
   !!----------------------------------------------------------------------
   USE oce_trc
   USE trc
   USE trcrst
   USE trcnam          ! Namelist read
   USE trcini_cfc      ! CFC      initialisation
   USE trcini_lobster  ! LOBSTER  initialisation
   USE trcini_pisces   ! PISCES   initialisation
   USE trcini_c14b     ! C14 bomb initialisation
   USE trcini_my_trc   ! MY_TRC   initialisation
#if defined key_bfm
   USE par_bfm
#endif
   USE trcdta   
   USE daymod
   USE zpshde          ! partial step: hor. derivative   (zps_hde routine)
   USE prtctl_trc      ! Print control passive tracers (prt_ctl_trc_init routine)
   
   IMPLICIT NONE
   PRIVATE
   
   PUBLIC   trc_init   ! called by opa

    !! * Substitutions
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2011)
   !! $Id: trcini.F90 2715 2011-03-30 15:58:35Z rblod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS
   
   SUBROUTINE trc_init
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_init  ***
      !!
      !! ** Purpose :   Initialization of the passive tracer fields 
      !!
      !! ** Method  : - read namelist
      !!              - control the consistancy 
      !!              - compute specific initialisations
      !!              - set initial tracer fields (either read restart 
      !!                or read data or analytical formulation
      !!---------------------------------------------------------------------
      INTEGER ::   jk, jn    ! dummy loop indices
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'trc_init : initial set up of the passive tracers'
      IF(lwp) WRITE(numout,*) '~~~~~~~'

      CALL top_alloc()              ! allocate TOP arrays

      !                             ! masked grid volume
      DO jk = 1, jpk
         cvol(:,:,jk) = e1t(:,:) * e2t(:,:) * fse3t(:,:,jk) * tmask(:,:,jk) 
      END DO

      !                             ! total volume of the ocean
#if ! defined key_degrad
      areatot = glob_sum( cvol(:,:,:) )
#else
      areatot = glob_sum( cvol(:,:,:) * facvol(:,:,:) )  ! degrad option: reduction by facvol
#endif

      CALL trc_nam                  ! read passive tracers namelists

      !                             ! restart for passive tracer (input)
      IF( ln_rsttr ) THEN
         IF(lwp) WRITE(numout,*) '       read a restart file for passive tracer : ', cn_trcrst_in
         IF(lwp) WRITE(numout,*) ' '
      ELSE
         IF( lwp .AND. lk_dtatrc ) THEN
            DO jn = 1, jptra
               IF( lutini(jn) )  &                  ! open input FILE only IF lutini(jn) is true
                  &  WRITE(numout,*) ' read an initial file for passive tracer number :', jn, ' traceur : ', ctrcnm(jn) 
             END DO
          ENDIF
          IF( lwp ) WRITE(numout,*)
      ENDIF

      IF( ln_dm2dc .AND. ( lk_pisces .OR. lk_lobster .OR. lk_bfm) )    &
         &       CALL ctl_stop( ' The diurnal cycle is not compatible with PISCES or LOBSTER or BFM ' )

      IF( nn_cla == 1 )   &
         &       CALL ctl_stop( ' Cross Land Advection not yet implemented with passive tracer ; nn_cla must be 0' )

      IF( lk_lobster ) THEN   ;   CALL trc_ini_lobster      ! LOBSTER bio-model
      ELSE                    ;   IF(lwp) WRITE(numout,*) '          LOBSTER not used'
      ENDIF
      
      IF( lk_pisces  ) THEN   ;   CALL trc_ini_pisces       ! PISCES  bio-model
      ELSE                    ;   IF(lwp) WRITE(numout,*) '          PISCES not used'
      ENDIF
      
      IF( lk_cfc     ) THEN   ;   CALL trc_ini_cfc          ! CFC     tracers
      ELSE                    ;   IF(lwp) WRITE(numout,*) '          CFC not used'
      ENDIF

      IF( lk_c14b    ) THEN   ;   CALL trc_ini_c14b         ! C14 bomb  tracer
      ELSE                    ;   IF(lwp) WRITE(numout,*) '          C14 not used'
      ENDIF
      
      IF( lk_my_trc  ) THEN   ;   CALL trc_ini_my_trc       ! MY_TRC  tracers
      ELSE                    ;   IF(lwp) WRITE(numout,*) '          MY_TRC not used'
      ENDIF

      IF( lk_bfm     ) THEN   ;   CALL trc_ini_bfm          ! BFM  tracers
      ELSE                    ;   IF(lwp) WRITE(numout,*) '          BFM not used'
      ENDIF

      IF( ln_rsttr ) THEN
        !
        IF( lk_offline )  neuler = 1   ! Set time-step indicator at nit000 (leap-frog)
        CALL trc_rst_read              ! restart from a file
        !
      ELSE
        IF( lk_offline )  THEN
           neuler = 0                  ! Set time-step indicator at nit000 (euler)
           CALL day_init               ! set calendar
        ENDIF
#if defined key_dtatrc
        CALL trc_dta( nit000 )      ! Initialization of tracer from a file that may also be used for damping
        DO jn = 1, jptra
           IF( lutini(jn) )   trn(:,:,:,jn) = trdta(:,:,:,jn) * tmask(:,:,:)   ! initialisation from file if required
        END DO
#endif
        trb(:,:,:,:) = trn(:,:,:,:)
        ! 
      ENDIF
 
      tra(:,:,:,:) = 0._wp
      
      IF( ln_zps .AND. .NOT. lk_c1d )   &              ! Partial steps: before horizontal gradient of passive
        &    CALL zps_hde( nit000, jptra, trn, gtru, gtrv )       ! tracers at the bottom ocean level


      !           
      trai = 0._wp         ! Computation content of all tracers
      DO jn = 1, jptra
#if ! defined key_degrad
         trai = trai + glob_sum( trn(:,:,:,jn) * cvol(:,:,:) )
#else
         trai = trai + glob_sum( trn(:,:,:,jn) * cvol(:,:,:) * facvol(:,:,:) ) ! degrad option: reduction by facvol
#endif
      END DO      

      IF(lwp) THEN               ! control print
         WRITE(numout,*)
         WRITE(numout,*)
         WRITE(numout,*) '          *** Total number of passive tracer jptra = ', jptra
         WRITE(numout,*) '          *** Total volume of ocean                = ', areatot
         WRITE(numout,*) '          *** Total inital content of all tracers  = ', trai
         WRITE(numout,*)
      ENDIF

      IF(ln_ctl) THEN            ! print mean trends (used for debugging)
         CALL prt_ctl_trc_init
         WRITE(charout, FMT="('ini ')")
         CALL prt_ctl_trc_info( charout )
         CALL prt_ctl_trc( tab4d=trn, mask=tmask, clinfo=ctrcnm )
      ENDIF
      !
   END SUBROUTINE trc_init


   SUBROUTINE top_alloc
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE top_alloc  ***
      !!
      !! ** Purpose :   Allocate all the dynamic arrays of the OPA modules
      !!----------------------------------------------------------------------
      USE trcadv        , ONLY:   trc_adv_alloc          ! TOP-related alloc routines...
      USE trc           , ONLY:   trc_alloc
#ifdef key_bfm
      USE trcnxtbfm     , ONLY:   trc_nxt_alloc
#else
      USE trcnxt        , ONLY:   trc_nxt_alloc
#endif
      USE trczdf        , ONLY:   trc_zdf_alloc
      USE trdmod_trc_oce, ONLY:   trd_mod_trc_oce_alloc
#if ! defined key_iomput
      USE trcdia        , ONLY:   trc_dia_alloc
#endif
#if defined key_trcdmp 
      USE trcdmp        , ONLY:   trc_dmp_alloc
#endif
#if defined key_dtatrc
      USE trcdta        , ONLY:   trc_dta_alloc
#endif
#if defined key_trdmld_trc   ||   defined key_esopa
      USE trdmld_trc    , ONLY:   trd_mld_trc_alloc
#endif
      !
      INTEGER :: ierr
      !!----------------------------------------------------------------------
      !
      ierr =        trc_adv_alloc()          ! Start of TOP-related alloc routines...
      ierr = ierr + trc_alloc    ()
      ierr = ierr + trc_nxt_alloc()
      ierr = ierr + trc_zdf_alloc()
      ierr = ierr + trd_mod_trc_oce_alloc()
#if ! defined key_iomput
      ierr = ierr + trc_dia_alloc()
#endif
#if defined key_trcdmp 
      ierr = ierr + trc_dmp_alloc()
#endif
#if defined key_dtatrc
      ierr = ierr + trc_dta_alloc()
#endif
#if defined key_trdmld_trc   ||   defined key_esopa
      ierr = ierr + trd_mld_trc_alloc()
#endif
      !
      IF( lk_mpp    )   CALL mpp_sum( ierr )
      IF( ierr /= 0 )   CALL ctl_stop( 'STOP', 'top_alloc : unable to allocate standard ocean arrays' )
      !
   END SUBROUTINE top_alloc

#else
   !!----------------------------------------------------------------------
   !!  Empty module :                                     No passive tracer
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_init                      ! Dummy routine   
   END SUBROUTINE trc_init
#endif

   !!======================================================================
END MODULE trcini
