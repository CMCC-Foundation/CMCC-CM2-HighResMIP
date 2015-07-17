MODULE trcsed
   !!======================================================================
   !!                         ***  MODULE p4sed  ***
   !! TOP :   PISCES Compute loss of organic matter in the sediments
   !!======================================================================
   !! History :    -   !  1995-06 (M. Levy)  original code
   !!              -   !  2000-12 (E. Kestenare)  clean up
   !!             2.0  !  2007-12  (C. Deltel, G. Madec)  F90 + simplifications
   !!----------------------------------------------------------------------
#if defined key_lobster
   !!----------------------------------------------------------------------
   !!   'key_lobster'                                     LOBSTER bio-model
   !!----------------------------------------------------------------------
   !!   trc_sed        :  Compute loss of organic matter in the sediments
   !!----------------------------------------------------------------------
   USE oce_trc         !
   USE trc
   USE sms_lobster
   USE lbclnk
   USE trdmod_oce
   USE trdmod_trc
   USE iom
   USE prtctl_trc      ! Print control for debbuging

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_sed    ! called in ???

   !!* Substitution
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcsed.F90 2715 2011-03-30 15:58:35Z rblod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trc_sed( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_sed  ***
      !!
      !! ** Purpose :   compute the now trend due to the vertical sedimentation of
      !!              detritus and add it to the general trend of detritus equations
      !!
      !! ** Method  :   this ROUTINE compute not exactly the advection but the
      !!              transport term, i.e.  dz(wt) and dz(ws)., dz(wtr)
      !!              using an upstream scheme
      !!              the now vertical advection of tracers is given by:
      !!                      dz(trn wn) = 1/bt dk+1( e1t e2t vsed (trn) )
      !!              add this trend now to the general trend of tracer (ta,sa,tra):
      !!                             tra = tra + dz(trn wn)
      !!        
      !!              IF 'key_diabio' is defined, the now vertical advection
      !!              trend of passive tracers is saved for futher diagnostics.
      !!---------------------------------------------------------------------
      USE wrk_nemo, ONLY: wrk_in_use, wrk_not_released
      USE wrk_nemo, ONLY: zwork => wrk_3d_2
      USE wrk_nemo, ONLY: zw2d  => wrk_2d_1 ! only used (if defined 
                                            ! key_diatrc && defined key_iomput)
      !!
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index      
      !!
      INTEGER  ::   ji, jj, jk, jl
      REAL(wp) ::   ztra
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::   ztrbio
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------

      IF( ( wrk_in_use(3,2)) .OR. ( wrk_in_use(2,1)) ) THEN
         CALL ctl_stop('trc_sed : requested workspace arrays unavailable.')
         RETURN
      END IF

      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) ' trc_sed: LOBSTER sedimentation'
         IF(lwp) WRITE(numout,*) ' ~~~~~~~'
      ENDIF

      ! sedimentation of detritus  : upstream scheme
      ! --------------------------------------------

      ! for detritus sedimentation only - jp_lob_det
      zwork(:,:,1  ) = 0.e0      ! surface value set to zero
      zwork(:,:,jpk) = 0.e0      ! bottom value  set to zero

#if defined key_diatrc && defined key_iomput
      zw2d(:,:) = 0.
# endif

      IF( l_trdtrc )THEN
         ALLOCATE( ztrbio(jpi,jpj,jpk) )
         ztrbio(:,:,:) = tra(:,:,:,jp_lob_det)
      ENDIF

      ! tracer flux at w-point: we use -vsed (downward flux)  with simplification : no e1*e2
      DO jk = 2, jpkm1
         zwork(:,:,jk) = -vsed * trn(:,:,jk-1,jp_lob_det)
      END DO

      ! tracer flux divergence at t-point added to the general trend
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1,jpi
               ztra  = - ( zwork(ji,jj,jk) - zwork(ji,jj,jk+1) ) / fse3t(ji,jj,jk)
               tra(ji,jj,jk,jp_lob_det) = tra(ji,jj,jk,jp_lob_det) + ztra
#if defined key_diabio
               trbio(ji,jj,jk,jp_lob0_trd + 7) = ztra
#endif
#if defined key_diatrc
# if ! defined key_iomput
               trc2d(ji,jj,jp_lob0_2d + 7) = trc2d(ji,jj,jp_lob0_2d + 7) + ztra * fse3t(ji,jj,jk) * 86400.
# else
               zw2d(ji,jj) = zw2d(ji,jj) + ztra * fse3t(ji,jj,jk) * 86400.
# endif
#endif
            END DO
         END DO
      END DO

#if defined key_diabio
      jl = jp_lob0_trd + 7
      CALL lbc_lnk (trbio(:,:,1,jl), 'T', 1. )    ! Lateral boundary conditions on trcbio
#endif
#if defined key_diatrc
# if ! defined key_iomput
      jl = jp_lob0_2d + 7
      CALL lbc_lnk( trc2d(:,:,jl), 'T', 1. )      ! Lateral boundary conditions on trc2d
# else
      CALL lbc_lnk( zw2d(:,:), 'T', 1. )      ! Lateral boundary conditions on zw2d
      CALL iom_put( "TDETSED", zw2d )
# endif
#endif
      !

      IF( l_trdtrc ) THEN
         ztrbio(:,:,:) = tra(:,:,:,jp_lob_det) - ztrbio(:,:,:)
         jl = jp_lob0_trd + 7
         CALL trd_mod_trc( ztrbio, jl, kt )   ! handle the trend
      ENDIF

      IF( l_trdtrc ) DEALLOCATE( ztrbio )

      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('sed')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF

      IF( ( wrk_not_released(3, 2) ) .OR. ( wrk_not_released(2, 1) ) )  &
       &         CALL ctl_stop('trc_sed : failed to release workspace arrays.')

   END SUBROUTINE trc_sed

#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE trc_sed( kt )                   ! Empty routine
      INTEGER, INTENT( in ) ::   kt
      WRITE(*,*) 'trc_sed: You should not have seen this print! error?', kt
   END SUBROUTINE trc_sed
#endif 

   !!======================================================================
END MODULE  trcsed
