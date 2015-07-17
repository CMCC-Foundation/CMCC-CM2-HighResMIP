MODULE domvvl
   !!======================================================================
   !!                       ***  MODULE domvvl   ***
   !! Ocean : 
   !!======================================================================
   !! History :  2.0  !  2006-06  (B. Levier, L. Marie)  original code
   !!            3.1  !  2009-02  (G. Madec, M. Leclair, R. Benshila)  pure z* coordinate
   !!----------------------------------------------------------------------
#if defined key_vvl
   !!----------------------------------------------------------------------
   !!   'key_vvl'                              variable volume
   !!----------------------------------------------------------------------
   !!   dom_vvl     : defined coefficients to distribute ssh on each layers
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE sbc_oce         ! surface boundary condition: ocean
   USE phycst          ! physical constants
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! distributed memory computing library
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dom_vvl       ! called by domain.F90
   PUBLIC   dom_vvl_alloc ! called by nemogcm.F90

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   ee_t, ee_u, ee_v, ee_f   !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   mut , muu , muv , muf    !: ??? 

   REAL(wp),         ALLOCATABLE, SAVE, DIMENSION(:)     ::   r2dt   ! vertical profile time-step, = 2 rdttra 
      !                                                              ! except at nit000 (=rdttra) if neuler=0

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 4.0 , NEMO Consortium (2011)
   !! $Id: domvvl.F90 2779 2011-06-07 17:09:43Z acc $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS       

   INTEGER FUNCTION dom_vvl_alloc()
      !!----------------------------------------------------------------------
      !!                ***  ROUTINE dom_vvl_alloc  ***
      !!----------------------------------------------------------------------
      !
      ALLOCATE( mut (jpi,jpj,jpk) , muu (jpi,jpj,jpk) , muv (jpi,jpj,jpk) , muf (jpi,jpj,jpk) ,     &
         &      ee_t(jpi,jpj)     , ee_u(jpi,jpj)     , ee_v(jpi,jpj)     , ee_f(jpi,jpj)     ,     &
         &      r2dt        (jpk)                                                             , STAT=dom_vvl_alloc )
         !
      IF( lk_mpp             )   CALL mpp_sum ( dom_vvl_alloc )
      IF( dom_vvl_alloc /= 0 )   CALL ctl_warn('dom_vvl_alloc: failed to allocate arrays')
      !
   END FUNCTION dom_vvl_alloc


   SUBROUTINE dom_vvl
      !!----------------------------------------------------------------------
      !!                ***  ROUTINE dom_vvl  ***
      !!                   
      !! ** Purpose :  compute coefficients muX at T-U-V-F points to spread
      !!               ssh over the whole water column (scale factors)
      !!----------------------------------------------------------------------
      USE wrk_nemo, ONLY:   wrk_in_use, wrk_not_released
      USE wrk_nemo, ONLY:   zs_t   => wrk_2d_1 , zs_u_1 => wrk_2d_2 , zs_v_1 => wrk_2d_3     ! 2D workspace
      !
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   zcoefu , zcoefv   , zcoeff                   ! local scalars
      REAL(wp) ::   zv_t_ij, zv_t_ip1j, zv_t_ijp1, zv_t_ip1jp1   !   -      -
      !!----------------------------------------------------------------------

      IF( wrk_in_use(2, 1,2,3) ) THEN
         CALL ctl_stop('dom_vvl: requested workspace arrays unavailable')   ;   RETURN
      ENDIF

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'dom_vvl : Variable volume initialization'
         WRITE(numout,*) '~~~~~~~~  compute coef. used to spread ssh over each layers'
      ENDIF
      
      IF( dom_vvl_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'dom_vvl : unable to allocate arrays' )

      fsdept(:,:,:) = gdept (:,:,:)
      fsdepw(:,:,:) = gdepw (:,:,:)
      fsde3w(:,:,:) = gdep3w(:,:,:)
      fse3t (:,:,:) = e3t   (:,:,:)
      fse3u (:,:,:) = e3u   (:,:,:)
      fse3v (:,:,:) = e3v   (:,:,:)
      fse3f (:,:,:) = e3f   (:,:,:)
      fse3w (:,:,:) = e3w   (:,:,:)
      fse3uw(:,:,:) = e3uw  (:,:,:)
      fse3vw(:,:,:) = e3vw  (:,:,:)

      !                                 !==  mu computation  ==!
      ee_t(:,:) = fse3t_0(:,:,1)                ! Lower bound : thickness of the first model level
      ee_u(:,:) = fse3u_0(:,:,1)
      ee_v(:,:) = fse3v_0(:,:,1)
      ee_f(:,:) = fse3f_0(:,:,1)
      DO jk = 2, jpkm1                          ! Sum of the masked vertical scale factors
         ee_t(:,:) = ee_t(:,:) + fse3t_0(:,:,jk) * tmask(:,:,jk)
         ee_u(:,:) = ee_u(:,:) + fse3u_0(:,:,jk) * umask(:,:,jk)
         ee_v(:,:) = ee_v(:,:) + fse3v_0(:,:,jk) * vmask(:,:,jk)
         DO jj = 1, jpjm1                      ! f-point : fmask=shlat at coasts, use the product of umask
            ee_f(:,jj) = ee_f(:,jj) + fse3f_0(:,jj,jk) *  umask(:,jj,jk) * umask(:,jj+1,jk)
         END DO
      END DO  
      !                                         ! Compute and mask the inverse of the local depth at T, U, V and F points
      ee_t(:,:) = 1. / ee_t(:,:) * tmask(:,:,1)
      ee_u(:,:) = 1. / ee_u(:,:) * umask(:,:,1)
      ee_v(:,:) = 1. / ee_v(:,:) * vmask(:,:,1)
      DO jj = 1, jpjm1                               ! f-point case fmask cannot be used 
         ee_f(:,jj) = 1. / ee_f(:,jj) * umask(:,jj,1) * umask(:,jj+1,1)
      END DO
      CALL lbc_lnk( ee_f, 'F', 1. )                  ! lateral boundary condition on ee_f
      !
      DO jk = 1, jpk                            ! mu coefficients
         mut(:,:,jk) = ee_t(:,:) * tmask(:,:,jk)     ! T-point at T levels
         muu(:,:,jk) = ee_u(:,:) * umask(:,:,jk)     ! U-point at T levels
         muv(:,:,jk) = ee_v(:,:) * vmask(:,:,jk)     ! V-point at T levels
      END DO
      DO jk = 1, jpk                                 ! F-point : fmask=shlat at coasts, use the product of umask
         DO jj = 1, jpjm1
               muf(:,jj,jk) = ee_f(:,jj) * umask(:,jj,jk) * umask(:,jj+1,jk)   ! at T levels
         END DO
         muf(:,jpj,jk) = 0.e0
      END DO
      CALL lbc_lnk( muf, 'F', 1. )                   ! lateral boundary condition


      hu_0(:,:) = 0.e0                          ! Reference ocean depth at U- and V-points
      hv_0(:,:) = 0.e0
      DO jk = 1, jpk
         hu_0(:,:) = hu_0(:,:) + fse3u_0(:,:,jk) * umask(:,:,jk)
         hv_0(:,:) = hv_0(:,:) + fse3v_0(:,:,jk) * vmask(:,:,jk)
      END DO
      
      ! surface at t-points and inverse surface at (u/v)-points used in surface averaging computations
      ! for ssh and scale factors
      zs_t  (:,:) =         e1t(:,:) * e2t(:,:)
      zs_u_1(:,:) = 0.5 / ( e1u(:,:) * e2u(:,:) )
      zs_v_1(:,:) = 0.5 / ( e1v(:,:) * e2v(:,:) )

      DO jj = 1, jpjm1                          ! initialise before and now Sea Surface Height at u-, v-, f-points
         DO ji = 1, jpim1   ! NO vector opt.
            zcoefu = umask(ji,jj,1) * zs_u_1(ji,jj)
            zcoefv = vmask(ji,jj,1) * zs_v_1(ji,jj)
            zcoeff = 0.5 * umask(ji,jj,1) * umask(ji,jj+1,1) / ( e1f(ji,jj) * e2f(ji,jj) )
            ! before fields
            zv_t_ij       = zs_t(ji  ,jj  ) * sshb(ji  ,jj  )
            zv_t_ip1j     = zs_t(ji+1,jj  ) * sshb(ji+1,jj  )
            zv_t_ijp1     = zs_t(ji  ,jj+1) * sshb(ji  ,jj+1)
            sshu_b(ji,jj) = zcoefu * ( zv_t_ij + zv_t_ip1j )
            sshv_b(ji,jj) = zcoefv * ( zv_t_ij + zv_t_ijp1 )
            ! now fields
            zv_t_ij       = zs_t(ji  ,jj  ) * sshn(ji  ,jj  )
            zv_t_ip1j     = zs_t(ji+1,jj  ) * sshn(ji+1,jj  )
            zv_t_ijp1     = zs_t(ji  ,jj+1) * sshn(ji  ,jj+1)
            zv_t_ip1jp1   = zs_t(ji  ,jj+1) * sshn(ji  ,jj+1)
            sshu_n(ji,jj) = zcoefu * ( zv_t_ij + zv_t_ip1j )
            sshv_n(ji,jj) = zcoefv * ( zv_t_ij + zv_t_ijp1 )
            sshf_n(ji,jj) = zcoeff * ( zv_t_ij + zv_t_ip1j + zv_t_ijp1 + zv_t_ip1jp1 )
         END DO
      END DO
      CALL lbc_lnk( sshu_n, 'U', 1. )   ;   CALL lbc_lnk( sshu_b, 'U', 1. )      ! lateral boundary conditions
      CALL lbc_lnk( sshv_n, 'V', 1. )   ;   CALL lbc_lnk( sshv_b, 'V', 1. )
      CALL lbc_lnk( sshf_n, 'F', 1. )

                                                ! initialise before scale factors at (u/v)-points
      ! Scale factor anomaly at (u/v)-points: surface averaging of scale factor at t-points
      DO jk = 1, jpkm1
         DO jj = 1, jpjm1
            DO ji = 1, jpim1
               zv_t_ij           = zs_t(ji  ,jj  ) * fse3t_b(ji  ,jj  ,jk)
               zv_t_ip1j         = zs_t(ji+1,jj  ) * fse3t_b(ji+1,jj  ,jk)
               zv_t_ijp1         = zs_t(ji  ,jj+1) * fse3t_b(ji  ,jj+1,jk)
               fse3u_b(ji,jj,jk) = umask(ji,jj,jk) * ( zs_u_1(ji,jj) * ( zv_t_ij + zv_t_ip1j ) - fse3u_0(ji,jj,jk) )
               fse3v_b(ji,jj,jk) = vmask(ji,jj,jk) * ( zs_v_1(ji,jj) * ( zv_t_ij + zv_t_ijp1 ) - fse3v_0(ji,jj,jk) )
            END DO
         END DO
      END DO
      CALL lbc_lnk( fse3u_b(:,:,:), 'U', 1. )               ! lateral boundary conditions
      CALL lbc_lnk( fse3v_b(:,:,:), 'V', 1. )
      ! Add initial scale factor to scale factor anomaly
      fse3u_b(:,:,:) = fse3u_b(:,:,:) + fse3u_0(:,:,:)
      fse3v_b(:,:,:) = fse3v_b(:,:,:) + fse3v_0(:,:,:)
      !
      IF( wrk_not_released(2, 1,2,3) )   CALL ctl_stop('dom_vvl: failed to release workspace arrays')
      !
   END SUBROUTINE dom_vvl

#else
   !!----------------------------------------------------------------------
   !!   Default option :                                      Empty routine
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE dom_vvl
   END SUBROUTINE dom_vvl
#endif

   !!======================================================================
END MODULE domvvl
