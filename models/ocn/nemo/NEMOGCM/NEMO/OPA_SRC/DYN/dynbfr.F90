MODULE dynbfr
   !!==============================================================================
   !!                 ***  MODULE  dynbfr  ***
   !! Ocean dynamics :  bottom friction component of the momentum mixing trend
   !!==============================================================================
   !! History :  3.2  ! 2008-11  (A. C. Coward)  Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dyn_bfr      : Update the momentum trend with the bottom friction contribution
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain variables 
   USE zdf_oce         ! ocean vertical physics variables
   USE trdmod          ! ocean active dynamics and tracers trends 
   USE trdmod_oce      ! ocean variables trends
   USE in_out_manager  ! I/O manager
   USE prtctl          ! Print control

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dyn_bfr    !  routine called by step.F90

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "zdfddm_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: dynbfr.F90 2715 2011-03-30 15:58:35Z rblod $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS
   
   SUBROUTINE dyn_bfr( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_bfr  ***
      !!
      !! ** Purpose :   compute the bottom friction ocean dynamics physics.
      !!
      !! ** Action  :   (ua,va)   momentum trend increased by bottom friction trend
      !!---------------------------------------------------------------------
      USE oce, ONLY:   ztrduv => tsa   ! tsa used as 4D workspace
      !!
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      !! 
      INTEGER  ::   ji, jj       ! dummy loop indexes
      INTEGER  ::   ikbu, ikbv   ! local integers
      REAL(wp) ::   zm1_2dt      ! local scalar
      !!---------------------------------------------------------------------
      !
      zm1_2dt = - 1._wp / ( 2._wp * rdt )

      IF( l_trddyn )   THEN                      ! temporary save of ua and va trends
         ztrduv(:,:,:,1) = ua(:,:,:)
         ztrduv(:,:,:,2) = va(:,:,:)
      ENDIF

# if defined key_vectopt_loop
      DO jj = 1, 1
         DO ji = jpi+2, jpij-jpi-1   ! vector opt. (forced unrolling)
# else
      DO jj = 2, jpjm1
         DO ji = 2, jpim1
# endif
            ikbu = mbku(ji,jj)          ! deepest ocean u- & v-levels
            ikbv = mbkv(ji,jj)
            !
            ! Apply stability criteria on absolute value  : abs(bfr/e3) < 1/(2dt) => bfr/e3 > -1/(2dt)
            ua(ji,jj,ikbu) = ua(ji,jj,ikbu) + MAX(  bfrua(ji,jj) / fse3u(ji,jj,ikbu) , zm1_2dt  ) * ub(ji,jj,ikbu)
            va(ji,jj,ikbv) = va(ji,jj,ikbv) + MAX(  bfrva(ji,jj) / fse3v(ji,jj,ikbv) , zm1_2dt  ) * vb(ji,jj,ikbv)
         END DO
      END DO

      !
      IF( l_trddyn )   THEN                      ! save the vertical diffusive trends for further diagnostics
         ztrduv(:,:,:,1) = ua(:,:,:) - ztrduv(:,:,:,1)
         ztrduv(:,:,:,2) = va(:,:,:) - ztrduv(:,:,:,2)
         CALL trd_mod( ztrduv(:,:,:,1), ztrduv(:,:,:,2), jpdyn_trd_bfr, 'DYN', kt )
      ENDIF
      !                                          ! print mean trends (used for debugging)
      IF(ln_ctl)   CALL prt_ctl( tab3d_1=ua, clinfo1=' bfr  - Ua: ', mask1=umask,               &
         &                       tab3d_2=va, clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
      !
   END SUBROUTINE dyn_bfr

   !!==============================================================================
END MODULE dynbfr
