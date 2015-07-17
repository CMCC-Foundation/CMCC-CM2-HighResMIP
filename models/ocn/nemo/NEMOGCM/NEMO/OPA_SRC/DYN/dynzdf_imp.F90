MODULE dynzdf_imp
   !!======================================================================
   !!                    ***  MODULE  dynzdf_imp  ***
   !! Ocean dynamics:  vertical component(s) of the momentum mixing trend
   !!======================================================================
   !! History :  OPA  !  1990-10  (B. Blanke)  Original code
   !!            8.0  !  1997-05  (G. Madec)  vertical component of isopycnal
   !!   NEMO     0.5  !  2002-08  (G. Madec)  F90: Free form and module
   !!            3.3  !  2010-04  (M. Leclair, G. Madec)  Forcing averaged over 2 time steps
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dyn_zdf_imp  : update the momentum trend with the vertical diffusion using a implicit time-stepping
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE sbc_oce         ! surface boundary condition: ocean
   USE zdf_oce         ! ocean vertical physics
   USE phycst          ! physical constants
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! MPP library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dyn_zdf_imp   ! called by step.F90

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: dynzdf_imp.F90 2715 2011-03-30 15:58:35Z rblod $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dyn_zdf_imp( kt, p2dt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_zdf_imp  ***
      !!                   
      !! ** Purpose :   Compute the trend due to the vert. momentum diffusion
      !!      and the surface forcing, and add it to the general trend of 
      !!      the momentum equations.
      !!
      !! ** Method  :   The vertical momentum mixing trend is given by :
      !!             dz( avmu dz(u) ) = 1/e3u dk+1( avmu/e3uw dk(ua) )
      !!      backward time stepping
      !!      Surface boundary conditions: wind stress input (averaged over kt-1/2 & kt+1/2)
      !!      Bottom boundary conditions : bottom stress (cf zdfbfr.F)
      !!      Add this trend to the general trend ua :
      !!         ua = ua + dz( avmu dz(u) )
      !!
      !! ** Action : - Update (ua,va) arrays with the after vertical diffusive mixing trend.
      !!---------------------------------------------------------------------
      USE wrk_nemo, ONLY:   wrk_in_use, wrk_not_released
      USE oce     , ONLY:  zwd  => ta       , zws   => sa   ! (ta,sa) used as 3D workspace
      USE wrk_nemo, ONLY:   zwi => wrk_3d_3                 ! 3D workspace
      !!
      INTEGER , INTENT(in) ::   kt    ! ocean time-step index
      REAL(wp), INTENT(in) ::  p2dt   ! vertical profile of tracer time-step
      !!
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   z1_p2dt, zcoef, zzwi, zzws, zrhs   ! local scalars
      !!----------------------------------------------------------------------

      IF( wrk_in_use(3, 3) ) THEN
         CALL ctl_stop('dyn_zdf_imp: requested workspace array unavailable')   ;   RETURN
      END IF

      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_zdf_imp : vertical momentum diffusion implicit operator'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~ '
      ENDIF

      ! 0. Local constant initialization
      ! --------------------------------
      z1_p2dt = 1._wp / p2dt      ! inverse of the timestep

      ! 1. Vertical diffusion on u
      ! ---------------------------
      ! Matrix and second member construction
      ! bottom boundary condition: both zwi and zws must be masked as avmu can take
      ! non zero value at the ocean bottom depending on the bottom friction
      ! used but the bottom velocities have already been updated with the bottom
      ! friction velocity in dyn_bfr using values from the previous timestep. There
      ! is no need to include these in the implicit calculation.
      !
      DO jk = 1, jpkm1        ! Matrix
         DO jj = 2, jpjm1 
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zcoef = - p2dt / fse3u(ji,jj,jk)
               zzwi          = zcoef * avmu (ji,jj,jk  ) / fse3uw(ji,jj,jk  )
               zwi(ji,jj,jk) = zzwi  * umask(ji,jj,jk)
               zzws          = zcoef * avmu (ji,jj,jk+1) / fse3uw(ji,jj,jk+1)
               zws(ji,jj,jk) = zzws  * umask(ji,jj,jk+1)
               zwd(ji,jj,jk) = 1._wp - zwi(ji,jj,jk) - zzws
            END DO
         END DO
      END DO
      DO jj = 2, jpjm1        ! Surface boudary conditions
         DO ji = fs_2, fs_jpim1   ! vector opt.
            zwi(ji,jj,1) = 0._wp
            zwd(ji,jj,1) = 1._wp - zws(ji,jj,1)
         END DO
      END DO

      ! Matrix inversion starting from the first level
      !-----------------------------------------------------------------------
      !   solve m.x = y  where m is a tri diagonal matrix ( jpk*jpk )
      !
      !        ( zwd1 zws1   0    0    0  )( zwx1 ) ( zwy1 )
      !        ( zwi2 zwd2 zws2   0    0  )( zwx2 ) ( zwy2 )
      !        (  0   zwi3 zwd3 zws3   0  )( zwx3 )=( zwy3 )
      !        (        ...               )( ...  ) ( ...  )
      !        (  0    0    0   zwik zwdk )( zwxk ) ( zwyk )
      !
      !   m is decomposed in the product of an upper and a lower triangular matrix
      !   The 3 diagonal terms are in 2d arrays: zwd, zws, zwi
      !   The solution (the after velocity) is in ua
      !-----------------------------------------------------------------------
      !
      DO jk = 2, jpkm1        !==  First recurrence : Dk = Dk - Lk * Uk-1 / Dk-1   (increasing k)  ==
         DO jj = 2, jpjm1   
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zwd(ji,jj,jk) = zwd(ji,jj,jk) - zwi(ji,jj,jk) * zws(ji,jj,jk-1) / zwd(ji,jj,jk-1)
            END DO
         END DO
      END DO
      !
      DO jj = 2, jpjm1        !==  second recurrence:    SOLk = RHSk - Lk / Dk-1  Lk-1  ==
         DO ji = fs_2, fs_jpim1   ! vector opt.
            ua(ji,jj,1) = ub(ji,jj,1) + p2dt * (  ua(ji,jj,1) + 0.5_wp * ( utau_b(ji,jj) + utau(ji,jj) )   &
               &                                                       / ( fse3u(ji,jj,1) * rau0       )  )
         END DO
      END DO
      DO jk = 2, jpkm1
         DO jj = 2, jpjm1   
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zrhs = ub(ji,jj,jk) + p2dt * ua(ji,jj,jk)   ! zrhs=right hand side
               ua(ji,jj,jk) = zrhs - zwi(ji,jj,jk) / zwd(ji,jj,jk-1) * ua(ji,jj,jk-1)
            END DO
         END DO
      END DO
      !
      DO jj = 2, jpjm1        !==  thrid recurrence : SOLk = ( Lk - Uk * Ek+1 ) / Dk  ==
         DO ji = fs_2, fs_jpim1   ! vector opt.
            ua(ji,jj,jpkm1) = ua(ji,jj,jpkm1) / zwd(ji,jj,jpkm1)
         END DO
      END DO
      DO jk = jpk-2, 1, -1
         DO jj = 2, jpjm1   
            DO ji = fs_2, fs_jpim1   ! vector opt.
               ua(ji,jj,jk) = ( ua(ji,jj,jk) - zws(ji,jj,jk) * ua(ji,jj,jk+1) ) / zwd(ji,jj,jk)
            END DO
         END DO
      END DO

      ! Normalization to obtain the general momentum trend ua
      DO jk = 1, jpkm1
         DO jj = 2, jpjm1   
            DO ji = fs_2, fs_jpim1   ! vector opt.
               ua(ji,jj,jk) = ( ua(ji,jj,jk) - ub(ji,jj,jk) ) * z1_p2dt
            END DO
         END DO
      END DO


      ! 2. Vertical diffusion on v
      ! ---------------------------
      ! Matrix and second member construction
      ! bottom boundary condition: both zwi and zws must be masked as avmv can take
      ! non zero value at the ocean bottom depending on the bottom friction
      ! used but the bottom velocities have already been updated with the bottom
      ! friction velocity in dyn_bfr using values from the previous timestep. There
      ! is no need to include these in the implicit calculation.
      !
      DO jk = 1, jpkm1        ! Matrix
         DO jj = 2, jpjm1   
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zcoef = -p2dt / fse3v(ji,jj,jk)
               zzwi          = zcoef * avmv (ji,jj,jk  ) / fse3vw(ji,jj,jk  )
               zwi(ji,jj,jk) =  zzwi * vmask(ji,jj,jk)
               zzws          = zcoef * avmv (ji,jj,jk+1) / fse3vw(ji,jj,jk+1)
               zws(ji,jj,jk) =  zzws * vmask(ji,jj,jk+1)
               zwd(ji,jj,jk) = 1._wp - zwi(ji,jj,jk) - zzws
            END DO
         END DO
      END DO
      DO jj = 2, jpjm1        ! Surface boudary conditions
         DO ji = fs_2, fs_jpim1   ! vector opt.
            zwi(ji,jj,1) = 0._wp
            zwd(ji,jj,1) = 1._wp - zws(ji,jj,1)
         END DO
      END DO

      ! Matrix inversion
      !-----------------------------------------------------------------------
      !   solve m.x = y  where m is a tri diagonal matrix ( jpk*jpk )
      !
      !        ( zwd1 zws1   0    0    0  )( zwx1 ) ( zwy1 )
      !        ( zwi2 zwd2 zws2   0    0  )( zwx2 ) ( zwy2 )
      !        (  0   zwi3 zwd3 zws3   0  )( zwx3 )=( zwy3 )
      !        (        ...               )( ...  ) ( ...  )
      !        (  0    0    0   zwik zwdk )( zwxk ) ( zwyk )
      !
      !   m is decomposed in the product of an upper and lower triangular matrix
      !   The 3 diagonal terms are in 2d arrays: zwd, zws, zwi
      !   The solution (after velocity) is in 2d array va
      !-----------------------------------------------------------------------
      !
      DO jk = 2, jpkm1        !==  First recurrence : Dk = Dk - Lk * Uk-1 / Dk-1   (increasing k)  ==
         DO jj = 2, jpjm1   
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zwd(ji,jj,jk) = zwd(ji,jj,jk) - zwi(ji,jj,jk) * zws(ji,jj,jk-1) / zwd(ji,jj,jk-1)
            END DO
         END DO
      END DO
      !
      DO jj = 2, jpjm1        !==  second recurrence:    SOLk = RHSk - Lk / Dk-1  Lk-1  ==
         DO ji = fs_2, fs_jpim1   ! vector opt.
            va(ji,jj,1) = vb(ji,jj,1) + p2dt * (  va(ji,jj,1) + 0.5_wp * ( vtau_b(ji,jj) + vtau(ji,jj) )   &
               &                                                       / ( fse3v(ji,jj,1) * rau0       )  )
         END DO
      END DO
      DO jk = 2, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zrhs = vb(ji,jj,jk) + p2dt * va(ji,jj,jk)   ! zrhs=right hand side
               va(ji,jj,jk) = zrhs - zwi(ji,jj,jk) / zwd(ji,jj,jk-1) * va(ji,jj,jk-1)
            END DO
         END DO
      END DO
      !
      DO jj = 2, jpjm1        !==  thrid recurrence : SOLk = ( Lk - Uk * SOLk+1 ) / Dk  ==
         DO ji = fs_2, fs_jpim1   ! vector opt.
            va(ji,jj,jpkm1) = va(ji,jj,jpkm1) / zwd(ji,jj,jpkm1)
         END DO
      END DO
      DO jk = jpk-2, 1, -1
         DO jj = 2, jpjm1   
            DO ji = fs_2, fs_jpim1   ! vector opt.
               va(ji,jj,jk) = ( va(ji,jj,jk) - zws(ji,jj,jk) * va(ji,jj,jk+1) ) / zwd(ji,jj,jk)
            END DO
         END DO
      END DO

      ! Normalization to obtain the general momentum trend va
      DO jk = 1, jpkm1
         DO jj = 2, jpjm1   
            DO ji = fs_2, fs_jpim1   ! vector opt.
               va(ji,jj,jk) = ( va(ji,jj,jk) - vb(ji,jj,jk) ) * z1_p2dt
            END DO
         END DO
      END DO
      !
      IF( wrk_not_released(3, 3) )   CALL ctl_stop('dyn_zdf_imp: failed to release workspace array')
      !
   END SUBROUTINE dyn_zdf_imp

   !!==============================================================================
END MODULE dynzdf_imp
