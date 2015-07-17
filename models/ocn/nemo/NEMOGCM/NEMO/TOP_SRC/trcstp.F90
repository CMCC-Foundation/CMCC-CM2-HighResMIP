MODULE trcstp
   !!======================================================================
   !!                       ***  MODULE trcstp  ***
   !! Time-stepping    : time loop of opa for passive tracer
   !!======================================================================
   !! History :  1.0  !  2004-03  (C. Ethe)  Original
   !!----------------------------------------------------------------------
#if defined key_top
   !!----------------------------------------------------------------------
   !!   trc_stp      : passive tracer system time-stepping
   !!----------------------------------------------------------------------
   USE oce_trc          ! ocean dynamics and active tracers variables
#if ! defined key_bfm
   USE trc
   USE trctrp           ! passive tracers transport
   USE trcsms           ! passive tracers sources and sinks
   USE prtctl_trc       ! Print control for debbuging
   USE trcdia
   USE trcwri
   USE trcrst
   USE trdmod_trc_oce
   USE trdmld_trc
#else
   USE api_bfm, ONLY: bio_calc
   USE par_bfm
#ifdef key_obcbfm
   USE trcobcdta_bfm
#endif
#endif
   USE iom
   USE in_out_manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_stp    ! called by step
   
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcstp.F90 2528 2010-12-27 17:33:53Z rblod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_stp( kt )
      !!-------------------------------------------------------------------
      !!                     ***  ROUTINE trc_stp  ***
      !!                      
      !! ** Purpose : Time loop of opa for passive tracer
      !! 
      !! ** Method  : 
      !!              Compute the passive tracers trends 
      !!              Update the passive tracers
      !!-------------------------------------------------------------------
      INTEGER, INTENT( in ) ::  kt  ! ocean time-step index
      CHARACTER (len=25)    ::  charout
      !!-------------------------------------------------------------------

#if defined key_bfm
   !---------------------------------------------
   ! Check the main BFM flag
   !---------------------------------------------
      IF (bio_calc) THEN

                             CALL trc_bfm( kt )           ! main call to BFM

                             CALL trc_trp_bfm( kt )       ! transport of BFM tracers

                             CALL trc_dia_bfm( kt )       ! diagnostic output for BFM
      END IF
#else

      IF( MOD( kt - 1 , nn_dttrc ) == 0 ) THEN      ! only every nn_dttrc time step
         !
         IF(ln_ctl) THEN
            WRITE(charout,FMT="('kt =', I4,'  d/m/y =',I2,I2,I4)") kt, nday, nmonth, nyear
            CALL prt_ctl_trc_info(charout)
         ENDIF
         !
         tra(:,:,:,:) = 0.e0
         !
         IF( kt == nit000 .AND. lk_trdmld_trc  )  &
            &                      CALL trd_mld_trc_init        ! trends: Mixed-layer
                                   CALL trc_rst_opn( kt )       ! Open tracer restart file 
         IF( lk_iomput ) THEN  ;   CALL trc_wri( kt )           ! output of passive tracers
         ELSE                  ;   CALL trc_dia( kt )
         ENDIF
                                   CALL trc_sms( kt )           ! tracers: sink and source
                                   CALL trc_trp( kt )           ! transport of passive tracers
         IF( kt == nit000 )     CALL iom_close( numrtr )     ! close input  passive tracers restart file
         IF( lrst_trc )            CALL trc_rst_wri( kt )       ! write tracer restart file
         IF( lk_trdmld_trc  )      CALL trd_mld_trc( kt )       ! trends: Mixed-layer
         !
      ENDIF
#endif

   END SUBROUTINE trc_stp

#else
   !!----------------------------------------------------------------------
   !!   Default key                                     NO passive tracers
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_stp( kt )        ! Empty routine
      WRITE(*,*) 'trc_stp: You should not have seen this print! error?', kt
   END SUBROUTINE trc_stp
#endif

   !!======================================================================
END MODULE trcstp
