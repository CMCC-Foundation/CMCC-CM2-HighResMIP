MODULE dynkeg
   !!======================================================================
   !!                       ***  MODULE  dynkeg  ***
   !! Ocean dynamics:  kinetic energy gradient trend
   !!======================================================================
   !! History :  1.0  !  87-09  (P. Andrich, m.-a. Foujols)  Original code
   !!            7.0  !  97-05  (G. Madec)  Split dynber into dynkeg and dynhpg
   !!            9.0  !  02-07  (G. Madec)  F90: Free form and module
   !!----------------------------------------------------------------------
   
   !!----------------------------------------------------------------------
   !!   dyn_keg      : update the momentum trend with the horizontal tke
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE trdmod          ! ocean dynamics trends 
   USE trdmod_oce      ! ocean variables trends
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! MPP library
   USE prtctl          ! Print control

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dyn_keg    ! routine called by step module
   
   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: dynkeg.F90 2777 2011-06-07 09:55:02Z smasson $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dyn_keg( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_keg  ***
      !!
      !! ** Purpose :   Compute the now momentum trend due to the horizontal
      !!      gradient of the horizontal kinetic energy and add it to the 
      !!      general momentum trend.
      !!
      !! ** Method  :   Compute the now horizontal kinetic energy 
      !!         zhke = 1/2 [ mi-1( un^2 ) + mj-1( vn^2 ) ]
      !!      Take its horizontal gradient and add it to the general momentum
      !!      trend (ua,va).
      !!         ua = ua - 1/e1u di[ zhke ]
      !!         va = va - 1/e2v dj[ zhke ]
      !!
      !! ** Action : - Update the (ua, va) with the hor. ke gradient trend
      !!             - save this trends (l_trddyn=T) for post-processing
      !!----------------------------------------------------------------------
      USE wrk_nemo, ONLY:   wrk_in_use, wrk_not_released
      USE oce     , ONLY:   ztrdu => ta       , ztrdv => sa   ! (ta,sa) used as 3D workspace   
      USE wrk_nemo, ONLY:   zhke  => wrk_3d_1                 ! 3D workspace
      !!
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step index
      !!
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   zu, zv       ! temporary scalars
      !!----------------------------------------------------------------------

      IF( wrk_in_use(3,1) ) THEN
         CALL ctl_stop('dyn_keg: requested workspace array is unavailable')   ;   RETURN
      ENDIF

      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_keg : kinetic energy gradient trend'
         IF(lwp) WRITE(numout,*) '~~~~~~~'
      ENDIF

      IF( l_trddyn ) THEN           ! Save ua and va trends
         ztrdu(:,:,:) = ua(:,:,:) 
         ztrdv(:,:,:) = va(:,:,:) 
      ENDIF
      
      !                                                ! ===============
      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============
         DO jj = 2, jpj         ! Horizontal kinetic energy at T-point
            DO ji = fs_2, jpi   ! vector opt.
               zu = 0.25 * (  un(ji-1,jj  ,jk) * un(ji-1,jj  ,jk)   &
                  &         + un(ji  ,jj  ,jk) * un(ji  ,jj  ,jk)  )
               zv = 0.25 * (  vn(ji  ,jj-1,jk) * vn(ji  ,jj-1,jk)   &
                  &         + vn(ji  ,jj  ,jk) * vn(ji  ,jj  ,jk)  )
               zhke(ji,jj,jk) = zv + zu
!!gm simplier coding  ==>>   ~ faster
!    don't forget to suppress local zu zv scalars
!               zhke(ji,jj,jk) = 0.25 * (   un(ji-1,jj  ,jk) * un(ji-1,jj  ,jk)   &
!                  &                      + un(ji  ,jj  ,jk) * un(ji  ,jj  ,jk)   &
!                  &                      + vn(ji  ,jj-1,jk) * vn(ji  ,jj-1,jk)   &
!                  &                      + vn(ji  ,jj  ,jk) * vn(ji  ,jj  ,jk)   )
!!gm end <<==
            END DO  
         END DO  
         DO jj = 2, jpjm1       ! add the gradient of kinetic energy to the general momentum trends
            DO ji = fs_2, fs_jpim1   ! vector opt.
               ua(ji,jj,jk) = ua(ji,jj,jk) - ( zhke(ji+1,jj  ,jk) - zhke(ji,jj,jk) ) / e1u(ji,jj)
               va(ji,jj,jk) = va(ji,jj,jk) - ( zhke(ji  ,jj+1,jk) - zhke(ji,jj,jk) ) / e2v(ji,jj)
            END DO 
         END DO
!!gm idea to be tested  ==>>   is it faster on scalar computers ?
!         DO jj = 2, jpjm1       ! add the gradient of kinetic energy to the general momentum trends
!            DO ji = fs_2, fs_jpim1   ! vector opt.
!               ua(ji,jj,jk) = ua(ji,jj,jk) - 0.25 * ( + un(ji+1,jj  ,jk) * un(ji+1,jj  ,jk)   &
!                  &                                   + vn(ji+1,jj-1,jk) * vn(ji+1,jj-1,jk)   &
!                  &                                   + vn(ji+1,jj  ,jk) * vn(ji+1,jj  ,jk)   &
!                  !
!                  &                                   - un(ji-1,jj  ,jk) * un(ji-1,jj  ,jk)   &
!                  &                                   - vn(ji  ,jj-1,jk) * vn(ji  ,jj-1,jk)   &
!                  &                                   - vn(ji  ,jj  ,jk) * vn(ji  ,jj  ,jk)   ) / e1u(ji,jj)
!                  !
!               va(ji,jj,jk) = va(ji,jj,jk) - 0.25 * (   un(ji-1,jj+1,jk) * un(ji-1,jj+1,jk)   &
!                  &                                   + un(ji  ,jj+1,jk) * un(ji  ,jj+1,jk)   &
!                  &                                   + vn(ji  ,jj+1,jk) * vn(ji  ,jj+1,jk)   &
!                  !
!                  &                                   - un(ji-1,jj  ,jk) * un(ji-1,jj  ,jk)   &
!                  &                                   - un(ji  ,jj  ,jk) * un(ji  ,jj  ,jk)   &
!                  &                                   - vn(ji  ,jj  ,jk) * vn(ji  ,jj  ,jk)   ) / e2v(ji,jj)
!            END DO 
!         END DO
!!gm en idea            <<==
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============

      IF( l_trddyn ) THEN      ! save the Kinetic Energy trends for diagnostic
         ztrdu(:,:,:) = ua(:,:,:) - ztrdu(:,:,:)
         ztrdv(:,:,:) = va(:,:,:) - ztrdv(:,:,:)
         CALL trd_mod( ztrdu, ztrdv, jpdyn_trd_keg, 'DYN', kt )
      ENDIF
      !
      IF(ln_ctl)   CALL prt_ctl( tab3d_1=ua, clinfo1=' keg  - Ua: ', mask1=umask,   &
         &                       tab3d_2=va, clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
      !
      IF( wrk_not_released(3, 1) )   CALL ctl_stop('dyn_keg: failed to release workspace array')
      !
   END SUBROUTINE dyn_keg

   !!======================================================================
END MODULE dynkeg
