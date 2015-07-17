MODULE trcsms_pisces
   !!======================================================================
   !!                         ***  MODULE trcsms_pisces  ***
   !! TOP :   PISCES Source Minus Sink manager
   !!======================================================================
   !! History :   1.0  !  2004-03 (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   'key_pisces'                                       PISCES bio-model
   !!----------------------------------------------------------------------
   !!   trcsms_pisces        :  Time loop of passive tracers sms
   !!----------------------------------------------------------------------
   USE oce_trc         !
   USE trc
   USE sms_pisces
   
   USE p4zint          ! 
   USE p4zche          ! 
   USE p4zbio          ! 
   USE p4zsink         ! 
   USE p4zopt          ! 
   USE p4zlim          ! 
   USE p4zprod         !
   USE p4zmort         !
   USE p4zmicro        ! 
   USE p4zmeso         ! 
   USE p4zrem          ! 
   USE p4zsed          ! 
   USE p4zlys          ! 
   USE p4zflx          ! 

   USE prtctl_trc

   USE trdmod_oce
   USE trdmod_trc

   USE sedmodel

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_sms_pisces    ! called in trcsms.F90

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcsms_pisces.F90 2715 2011-03-30 15:58:35Z rblod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trc_sms_pisces( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_sms_pisces  ***
      !!
      !! ** Purpose :   Managment of the call to Biological sources and sinks 
      !!              routines of PISCES bio-model
      !!
      !! ** Method  : - at each new day ...
      !!              - several calls of bio and sed ???
      !!              - ...
      !!---------------------------------------------------------------------
      USE wrk_nemo, ONLY: wrk_in_use, wrk_not_released
      USE wrk_nemo, ONLY: ztrpis => wrk_3d_1   ! used for pisces sms trends
      !
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index      
      !!
      INTEGER ::   jnt, jn
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------

      IF( kt == nit000 )   CALL trc_sms_pisces_init    ! Initialization (first time-step only)

      IF( wrk_in_use(3,1) )  THEN
        CALL ctl_stop('trc_sms_pisces : requested workspace array unavailable.')  ;  RETURN
      ENDIF

      IF( ndayflxtr /= nday_year ) THEN      ! New days
         !
         ndayflxtr = nday_year

         IF(lwp) write(numout,*)
         IF(lwp) write(numout,*) ' New chemical constants and various rates for biogeochemistry at new day : ', nday_year
         IF(lwp) write(numout,*) '~~~~~~'

         CALL p4z_che          ! computation of chemical constants
         CALL p4z_int          ! computation of various rates for biogeochemistry
         !
      ENDIF

      DO jnt = 1, nrdttrc          ! Potential time splitting if requested
         !
         CALL p4z_bio (kt, jnt)    ! Compute soft tissue production (POC)
         CALL p4z_sed (kt, jnt)    ! compute soft tissue remineralisation
         !
         trb(:,:,:,:) = trn(:,:,:,:)
         !
      END DO

      CALL p4z_lys( kt )             ! Compute CaCO3 saturation
      CALL p4z_flx( kt )             ! Compute surface fluxes

      DO jn = jp_pcs0, jp_pcs1
        CALL lbc_lnk( trn(:,:,:,jn), 'T', 1. )
        CALL lbc_lnk( trb(:,:,:,jn), 'T', 1. )
        CALL lbc_lnk( tra(:,:,:,jn), 'T', 1. )
      END DO


      IF( l_trdtrc ) THEN
          DO jn = jp_pcs0, jp_pcs1
            ztrpis(:,:,:) = tra(:,:,:,jn)
            CALL trd_mod_trc( ztrpis, jn, jptra_trd_sms, kt )   ! save trends
          END DO
          DEALLOCATE( ztrpis )
      END IF

      IF( lk_sed ) THEN 
         !
         CALL sed_model( kt )     !  Main program of Sediment model
         !
         DO jn = jp_pcs0, jp_pcs1
           CALL lbc_lnk( trn(:,:,:,jn), 'T', 1. )
         END DO
         !
      ENDIF

      IF( wrk_not_released(3,1) ) CALL ctl_stop('trc_sms_pisces : failed to release workspace array.') 

   END SUBROUTINE trc_sms_pisces

   SUBROUTINE trc_sms_pisces_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trc_sms_pisces_init  ***
      !!
      !! ** Purpose :   Initialization of PH variable
      !!
      !!----------------------------------------------------------------------
      INTEGER  ::  ji, jj, jk
      REAL(wp) ::  zcaralk, zbicarb, zco3
      REAL(wp) ::  ztmas, ztmas1

      IF( .NOT. ln_rsttr ) THEN
         ! Initialization of chemical variables of the carbon cycle
         ! --------------------------------------------------------
         DO jk = 1, jpk
            DO jj = 1, jpj
               DO ji = 1, jpi
                  ztmas   = tmask(ji,jj,jk)
                  ztmas1  = 1. - tmask(ji,jj,jk)
                  zcaralk = trn(ji,jj,jk,jptal) - borat(ji,jj,jk) / (  1. + 1.E-8 / ( rtrn + akb3(ji,jj,jk) )  )
                  zco3    = ( zcaralk - trn(ji,jj,jk,jpdic) ) * ztmas + 0.5e-3 * ztmas1
                  zbicarb = ( 2. * trn(ji,jj,jk,jpdic) - zcaralk )
                  hi(ji,jj,jk) = ( ak23(ji,jj,jk) * zbicarb / zco3 ) * ztmas + 1.e-9 * ztmas1
               END DO
            END DO
         END DO
         !
      END IF

      ! Time step duration for biology
      xstep = rfact2 / rday

      CALL p4z_sink_init      ! vertical flux of particulate organic matter
      CALL p4z_opt_init       ! Optic: PAR in the water column
      CALL p4z_lim_init       ! co-limitations by the various nutrients
      CALL p4z_prod_init      ! phytoplankton growth rate over the global ocean. 
      CALL p4z_rem_init       ! remineralisation
      CALL p4z_mort_init      ! phytoplankton mortality
      CALL p4z_micro_init     ! microzooplankton
      CALL p4z_meso_init      ! mesozooplankton
      CALL p4z_sed_init       ! sedimentation
      CALL p4z_lys_init       ! calcite saturation
      CALL p4z_flx_init       ! gas exchange

      ndayflxtr = 0

   END SUBROUTINE trc_sms_pisces_init

#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE trc_sms_pisces( kt )                   ! Empty routine
      INTEGER, INTENT( in ) ::   kt
      WRITE(*,*) 'trc_sms_pisces: You should not have seen this print! error?', kt
   END SUBROUTINE trc_sms_pisces
#endif 

   !!======================================================================
END MODULE trcsms_pisces 
