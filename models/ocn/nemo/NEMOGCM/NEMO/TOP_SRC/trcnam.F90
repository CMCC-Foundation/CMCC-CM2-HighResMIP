MODULE trcnam
   !!======================================================================
   !!                       ***  MODULE trcnam  ***
   !! TOP :   Read and print options for the passive tracer run (namelist)
   !!======================================================================
   !! History :    -   !  1996-11  (M.A. Foujols, M. Levy)  original code
   !!              -   !  1998-04  (M.A Foujols, L. Bopp) ahtrb0 for isopycnal mixing
   !!              -   !  1999-10  (M.A. Foujols, M. Levy) separation of sms
   !!              -   !  2000-07  (A. Estublier) add TVD and MUSCL : Tests on ndttrc
   !!              -   !  2000-11  (M.A Foujols, E Kestenare) trcrat, ahtrc0 and aeivtr0
   !!              -   !  2001-01 (E Kestenare) suppress ndttrc=1 for CEN2 and TVD schemes
   !!             1.0  !  2005-03 (O. Aumont, A. El Moussaoui) F90
   !!----------------------------------------------------------------------
#if defined key_top
   !!----------------------------------------------------------------------
   !!   'key_top'                                                TOP models
   !!----------------------------------------------------------------------
   !!   trc_nam    :  Read and print options for the passive tracer run (namelist)
   !!----------------------------------------------------------------------
   USE oce_trc
   USE trc
   USE trcnam_trp        ! Transport namelist
   USE trcnam_lobster    ! LOBSTER namelist
   USE trcnam_pisces     ! PISCES namelist
   USE trcnam_cfc        ! CFC SMS namelist
   USE trcnam_c14b       ! C14 SMS namelist
   USE trcnam_my_trc     ! MY_TRC SMS namelist
   USE trcnam_bfm        ! BFM namelist 
   USE trdmod_trc_oce

   IMPLICIT NONE
   PRIVATE 

   PUBLIC trc_nam      ! called in trcini

   !! * Substitutions
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcnam.F90 2715 2011-03-30 15:58:35Z rblod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trc_nam
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_nam  ***
      !!
      !! ** Purpose :   READ and PRINT options for the passive tracer run (namelist) 
      !!
      !! ** Method  : - read passive tracer namelist 
      !!              - read namelist of each defined SMS model
      !!                ( (LOBSTER, PISCES, CFC, MY_TRC )
      !!---------------------------------------------------------------------
      INTEGER ::  jn

      ! Definition of a tracer as a structure
      TYPE PTRACER
         CHARACTER(len = 20)  :: clsname  !: short name
         CHARACTER(len = 80 ) :: cllname  !: long name
         CHARACTER(len = 20 ) :: clunit   !: unit
         LOGICAL              :: llinit   !: read in a file or not
         LOGICAL              :: llsave   !: save the tracer or not
      END TYPE PTRACER

      TYPE(PTRACER) , DIMENSION(jptra) :: sn_tracer

      !!
      NAMELIST/namtrc/    nn_dttrc, nn_writetrc, ln_rsttr, nn_rsttr, &
                          cn_trcrst_in, cn_trcrst_out, sn_tracer
#if defined key_trdmld_trc  || defined key_trdtrc
      NAMELIST/namtrc_trd/ nn_trd_trc, nn_ctls_trc, rn_ucf_trc, &
                         ln_trdmld_trc_restart, ln_trdmld_trc_instant, &
                         cn_trdrst_trc_in, cn_trdrst_trc_out, ln_trdtrc
#endif

      !!---------------------------------------------------------------------

      IF(lwp) WRITE(numout,*) 'trc_nam : read the passive tracer namelists'
      IF(lwp) WRITE(numout,*) '~~~~~~~'

      CALL ctl_opn( numnat, 'namelist_top', 'OLD', 'FORMATTED', 'SEQUENTIAL', 1, numout, .FALSE. )

      ! Namelist nattrc (files)
      ! ----------------------------------------------
      nn_dttrc    = 1                 ! default values
      nn_writetrc = 10      
      ln_rsttr    = .FALSE.
      nn_rsttr    =  0
      cn_trcrst_in  = 'restart_trc'
      cn_trcrst_out = 'restart_trc'
      DO jn = 1, jptra
         WRITE(ctrcnm(jn),'("TR_",I1)'           ) jn
         WRITE(ctrcnl(jn),'("TRACER NUMBER ",I1)') jn
         ctrcun(jn) = 'mmole/m3'
         lutini(jn) = .FALSE. 
         lutsav(jn) = .TRUE. 
      END DO

      REWIND( numnat )               ! read nattrc
      READ  ( numnat, namtrc )

      DO jn = 1, jptra
         ctrcnm(jn) = TRIM( sn_tracer(jn)%clsname )
         ctrcnl(jn) = TRIM( sn_tracer(jn)%cllname )
         ctrcun(jn) = TRIM( sn_tracer(jn)%clunit  )
         lutini(jn) =       sn_tracer(jn)%llinit 
         lutsav(jn) =       sn_tracer(jn)%llsave
      END DO


      IF(lwp) THEN                   ! control print
         WRITE(numout,*)
         WRITE(numout,*) ' Namelist : namtrc'
         WRITE(numout,*) '    time step freq. for pass. trac. nn_dttrc             = ', nn_dttrc
         WRITE(numout,*) '    frequency of outputs for passive tracers nn_writetrc = ', nn_writetrc  
         WRITE(numout,*) '    restart LOGICAL for passive tr. ln_rsttr             = ', ln_rsttr
         WRITE(numout,*) '    control of time step for p. tr. nn_rsttr             = ', nn_rsttr
         WRITE(numout,*) ' '
         DO jn = 1, jptra
            WRITE(numout,*) '   tracer nb             : ', jn 
            WRITE(numout,*) '   short name            : ', ctrcnm(jn)
            WRITE(numout,*) '   long name             : ', ctrcnl(jn)
            WRITE(numout,*) '   unit                  : ', ctrcun(jn)
            WRITE(numout,*) '   initial value in FILE : ', lutini(jn) 
            WRITE(numout,*) ' '
         END DO
      ENDIF

      rdttrc(:) = rdttra(:) * FLOAT( nn_dttrc )   ! vertical profile of passive tracer time-step
  
      IF(lwp) WRITE(numout,*) 
      IF(lwp) WRITE(numout,*) '    Passive Tracer  time step    rdttrc  = ', rdttrc(1)
      IF(lwp) WRITE(numout,*) 

#if defined key_trdmld_trc || defined key_trdtrc
      nn_trd_trc  = 20
      nn_ctls_trc =  9
      rn_ucf_trc   =  1.
      ln_trdmld_trc_instant = .TRUE.
      ln_trdmld_trc_restart =.FALSE.
      cn_trdrst_trc_in  = "restart_mld_trc"
      cn_trdrst_trc_out = "restart_mld_trc"
      ln_trdtrc(:) = .FALSE.

      REWIND( numnat )               !  namelist namtoptrd : passive tracer trends diagnostic
      READ  ( numnat, namtrc_trd )

     IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) ' trd_mld_trc_init : read namelist namtrc_trd                    '
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~                                               '
         WRITE(numout,*) '   * frequency of trends diagnostics   nn_trd_trc             = ', nn_trd_trc
         WRITE(numout,*) '   * control surface type              nn_ctls_trc            = ', nn_ctls_trc
         WRITE(numout,*) '   * restart for ML diagnostics        ln_trdmld_trc_restart  = ', ln_trdmld_trc_restart
         WRITE(numout,*) '   * flag to diagnose trends of                                 '
         WRITE(numout,*) '     instantantaneous or mean ML T/S   ln_trdmld_trc_instant  = ', ln_trdmld_trc_instant
         WRITE(numout,*) '   * unit conversion factor            rn_ucf_trc             = ', rn_ucf_trc
         DO jn = 1, jptra
            IF( ln_trdtrc(jn) ) WRITE(numout,*) '    compute ML trends for tracer number :', jn
         END DO
      ENDIF
#endif

      ! namelist of transport
      ! ---------------------
      CALL trc_nam_trp


      ! namelist of SMS
      ! ---------------      
      IF( lk_lobster ) THEN   ;   CALL trc_nam_lobster      ! LOBSTER bio-model
      ELSE                    ;   IF(lwp) WRITE(numout,*) '          LOBSTER not used'
      ENDIF

      IF( lk_pisces  ) THEN   ;   CALL trc_nam_pisces      ! PISCES  bio-model
      ELSE                    ;   IF(lwp) WRITE(numout,*) '          PISCES not used'
      ENDIF

      IF( lk_cfc     ) THEN   ;   CALL trc_nam_cfc         ! CFC     tracers
      ELSE                    ;   IF(lwp) WRITE(numout,*) '          CFC not used'
      ENDIF

      IF( lk_c14b     ) THEN   ;   CALL trc_nam_c14b         ! C14 bomb     tracers
      ELSE                    ;   IF(lwp) WRITE(numout,*) '          C14 not used'
      ENDIF

      IF( lk_my_trc  ) THEN   ;   CALL trc_nam_my_trc      ! MY_TRC  tracers
      ELSE                    ;   IF(lwp) WRITE(numout,*) '          MY_TRC not used'
      ENDIF
      
      IF( lk_bfm  ) THEN      ;   CALL trc_nam_bfm         ! BFM  tracers
      ELSE                    ;   IF(lwp) WRITE(numout,*) '          BFM not used'
      ENDIF
      !
   END SUBROUTINE trc_nam

#else
   !!----------------------------------------------------------------------
   !!  Dummy module :                                     No passive tracer
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_nam                      ! Empty routine   
   END SUBROUTINE trc_nam
#endif

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcnam.F90 2715 2011-03-30 15:58:35Z rblod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!======================================================================
END MODULE  trcnam
