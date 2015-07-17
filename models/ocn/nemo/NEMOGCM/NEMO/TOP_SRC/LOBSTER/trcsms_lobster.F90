MODULE trcsms_lobster
   !!======================================================================
   !!                         ***  MODULE trcsms_lobster  ***
   !! TOP :   Time loop of LOBSTER model
   !!======================================================================
   !! History :   1.0  !            M. Levy
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  revised architecture
   !!----------------------------------------------------------------------
#if defined key_lobster
   !!----------------------------------------------------------------------
   !!   'key_lobster'                                       LOBSTER bio-model
   !!----------------------------------------------------------------------
   !!   trcsms_lobster        :  Time loop of passive tracers sms
   !!----------------------------------------------------------------------
   USE oce_trc          !
   USE trc
   USE trcbio
   USE trcopt
   USE trcsed
   USE trcexp
   USE trdmod_oce
   USE trdmod_trc_oce
   USE trdmod_trc
   USE trdmld_trc

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_sms_lobster    ! called in trcsms.F90

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcsms_lobster.F90 2715 2011-03-30 15:58:35Z rblod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_sms_lobster( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_sms_lobster  ***
      !!
      !! ** Purpose :  Managment of the call to Biological sources and sinks 
      !!               routines of LOBSTER bio-model 
      !!
      !! ** Method  : - ???
      !! --------------------------------------------------------------------
      USE wrk_nemo, ONLY: wrk_in_use, wrk_not_released
      USE wrk_nemo, ONLY: ztrlob => wrk_3d_1   ! used for lobster sms trends
      !!
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index      
      INTEGER :: jn
      !! --------------------------------------------------------------------

      IF( wrk_in_use(3, 1) ) THEN
         CALL ctl_stop('trc_sms_lobster : requested workspace array unavailable')   ;   RETURN
      ENDIF

      CALL trc_opt( kt )      ! optical model
      CALL trc_bio( kt )      ! biological model
      CALL trc_sed( kt )      ! sedimentation model
      CALL trc_exp( kt )      ! export

      IF( l_trdtrc ) THEN
          DO jn = jp_lob0, jp_lob1
            ztrlob(:,:,:) = tra(:,:,:,jn)
            CALL trd_mod_trc( ztrlob, jn, jptra_trd_sms, kt )   ! save trends
          END DO
      END IF

      IF( lk_trdmld_trc )  CALL trd_mld_bio( kt )   ! trends: Mixed-layer

      IF( wrk_not_released(3, 1) )   CALL ctl_stop('trc_sms_lobster : failed to release workspace array.')
      !
   END SUBROUTINE trc_sms_lobster

#else
   !!======================================================================
   !!  Dummy module :                                     No passive tracer
   !!======================================================================
CONTAINS
   SUBROUTINE trc_sms_lobster( kt )                   ! Empty routine
      INTEGER, INTENT( in ) ::   kt
      WRITE(*,*) 'trc_sms_lobster: You should not have seen this print! error?', kt
   END SUBROUTINE trc_sms_lobster
#endif 

   !!======================================================================
END MODULE  trcsms_lobster
