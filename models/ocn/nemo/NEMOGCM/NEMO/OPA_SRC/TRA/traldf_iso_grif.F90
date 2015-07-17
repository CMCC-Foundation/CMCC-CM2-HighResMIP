MODULE traldf_iso_grif
   !!======================================================================
   !!                   ***  MODULE  traldf_iso_grif  ***
   !! Ocean  tracers:  horizontal component of the lateral tracer mixing trend
   !!======================================================================
   !! History : 3.3  ! 2010-10  (G. Nurser, C. Harris, G. Madec)  
   !!                !          Griffies operator version adapted from traldf_iso.F90
   !!----------------------------------------------------------------------
#if   defined key_ldfslp   ||   defined key_esopa
   !!----------------------------------------------------------------------
   !!   'key_ldfslp'               slope of the lateral diffusive direction
   !!----------------------------------------------------------------------
   !!   tra_ldf_iso_grif  : update the tracer trend with the horizontal component  
   !!                       of the Griffies iso-neutral laplacian operator 
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and active tracers
   USE dom_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE trc_oce         ! share passive tracers/Ocean variables
   USE zdf_oce         ! ocean vertical physics
   USE ldftra_oce      ! ocean active tracers: lateral physics
   USE ldfslp          ! iso-neutral slopes
   USE diaptr          ! poleward transport diagnostics
   USE in_out_manager  ! I/O manager
   USE iom             ! I/O library
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp         ! MPP library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   tra_ldf_iso_grif   ! routine called by traldf.F90

   REAL(wp), PUBLIC, DIMENSION(:,:,:), ALLOCATABLE, SAVE ::   psix_eiv, psiy_eiv   !: eiv stream function (diag only)
   REAL(wp), PUBLIC, DIMENSION(:,:,:), ALLOCATABLE, SAVE ::   ah_wslp2             !: aeiv*w-slope^2
   REAL(wp),         DIMENSION(:,:,:), ALLOCATABLE, SAVE ::   zdkt                 !  atypic workspace

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "ldftra_substitute.h90"
#  include "vectopt_loop_substitute.h90"
#  include "ldfeiv_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: traldf_iso_grif.F90 2715 2011-03-30 15:58:35Z rblod $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

  SUBROUTINE tra_ldf_iso_grif( kt, cdtype, pgu, pgv,              &
       &                                   ptb, pta, kjpt, pahtb0 )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_ldf_iso_grif  ***
      !!
      !! ** Purpose :   Compute the before horizontal tracer (t & s) diffusive 
      !!      trend for a laplacian tensor (ezxcept the dz[ dz[.] ] term) and 
      !!      add it to the general trend of tracer equation.
      !!
      !! ** Method  :   The horizontal component of the lateral diffusive trends 
      !!      is provided by a 2nd order operator rotated along neural or geopo-
      !!      tential surfaces to which an eddy induced advection can be added
      !!      It is computed using before fields (forward in time) and isopyc-
      !!      nal or geopotential slopes computed in routine ldfslp.
      !!
      !!      1st part :  masked horizontal derivative of T  ( di[ t ] )
      !!      ========    with partial cell update if ln_zps=T.
      !!
      !!      2nd part :  horizontal fluxes of the lateral mixing operator
      !!      ========    
      !!         zftu = (aht+ahtb0) e2u*e3u/e1u di[ tb ]
      !!               - aht       e2u*uslp    dk[ mi(mk(tb)) ]
      !!         zftv = (aht+ahtb0) e1v*e3v/e2v dj[ tb ]
      !!               - aht       e2u*vslp    dk[ mj(mk(tb)) ]
      !!      take the horizontal divergence of the fluxes:
      !!         difft = 1/(e1t*e2t*e3t) {  di-1[ zftu ] +  dj-1[ zftv ]  }
      !!      Add this trend to the general trend (ta,sa):
      !!         ta = ta + difft
      !!
      !!      3rd part: vertical trends of the lateral mixing operator
      !!      ========  (excluding the vertical flux proportional to dk[t] )
      !!      vertical fluxes associated with the rotated lateral mixing:
      !!         zftw =-aht {  e2t*wslpi di[ mi(mk(tb)) ]
      !!                     + e1t*wslpj dj[ mj(mk(tb)) ]  }
      !!      take the horizontal divergence of the fluxes:
      !!         difft = 1/(e1t*e2t*e3t) dk[ zftw ]
      !!      Add this trend to the general trend (ta,sa):
      !!         pta = pta + difft
      !!
      !! ** Action :   Update pta arrays with the before rotated diffusion
      !!----------------------------------------------------------------------
      USE wrk_nemo, ONLY:   wrk_in_use, wrk_not_released
      USE oce     , ONLY:   zftu => ua       , zftv => va            ! (ua,va) used as 3D workspace
      USE wrk_nemo, ONLY:   zdit => wrk_3d_6 , zdjt => wrk_3d_7 , ztfw => wrk_3d_8   ! 3D workspace
      USE wrk_nemo, ONLY:   z2d  => wrk_2d_1                                         ! 2D workspace
      !
      INTEGER                              , INTENT(in   ) ::   kt         ! ocean time-step index
      CHARACTER(len=3)                     , INTENT(in   ) ::   cdtype     ! =TRA or TRC (tracer indicator)
      INTEGER                              , INTENT(in   ) ::   kjpt       ! number of tracers
      REAL(wp), DIMENSION(jpi,jpj    ,kjpt), INTENT(in   ) ::   pgu, pgv   ! tracer gradient at pstep levels
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(in   ) ::   ptb        ! before and now tracer fields
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(inout) ::   pta        ! tracer trend 
      REAL(wp)                             , INTENT(in   ) ::   pahtb0     ! background diffusion coef
      !
      INTEGER  ::  ji, jj, jk,jn   ! dummy loop indices
      INTEGER  ::  ip,jp,kp        ! dummy loop indices
      INTEGER  ::  ierr            ! temporary integer
      REAL(wp) ::  zmsku, zabe1, zcof1, zcoef3   ! local scalars
      REAL(wp) ::  zmskv, zabe2, zcof2, zcoef4   !   -      -
      REAL(wp) ::  zcoef0, zbtr                  !   -      -
      !REAL(wp), POINTER, DIMENSION(:,:,:) ::   zdkt           ! 2D+1 workspace
      !
      REAL(wp) ::   zslope_skew, zslope_iso, zslope2, zbu, zbv
      REAL(wp) ::   ze1ur, zdxt, ze2vr, ze3wr, zdyt, zdzt
      REAL(wp) ::   zah, zah_slp, zaei_slp
#if defined key_diaar5
      REAL(wp) ::   zztmp              ! local scalar
#endif
      !!----------------------------------------------------------------------

      IF( wrk_in_use(3, 6,7,8) .OR. wrk_in_use(2, 1) ) THEN
         CALL ctl_stop('tra_ldf_iso_grif: requested workspace arrays unavailable.')   ;   RETURN
      ENDIF
      ! ARP - line below uses 'bounds re-mapping' which is only defined in
      ! Fortran 2003 and up. We would be OK if code was written to use
      ! zdkt(:,:,1:2) instead as then wouldn't need to re-map bounds.
      ! As it is, we make zdkt a module array and allocate it in _alloc().
      !zdkt(1:jpi,1:jpj,0:1) => wrk_3d_9(:,:,1:2)

      IF( kt == nit000 )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tra_ldf_iso_grif : rotated laplacian diffusion operator on ', cdtype
         IF(lwp) WRITE(numout,*) '                   WARNING: STILL UNDER TEST, NOT RECOMMENDED. USE AT YOUR OWN PERIL'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
         ALLOCATE( ah_wslp2(jpi,jpj,jpk) , zdkt(jpi,jpj,0:1), STAT=ierr )
         IF( lk_mpp   )   CALL mpp_sum ( ierr )
         IF( ierr > 0 )   CALL ctl_stop('STOP', 'tra_ldf_iso_grif: unable to allocate arrays')
         IF( ln_traldf_gdia ) THEN
            ALLOCATE( psix_eiv(jpi,jpj,jpk) , psiy_eiv(jpi,jpj,jpk) , STAT=ierr )
            IF( lk_mpp   )   CALL mpp_sum ( ierr )
            IF( ierr > 0 )   CALL ctl_stop('STOP', 'tra_ldf_iso_grif: unable to allocate diagnostics')
         ENDIF
      ENDIF

      !!----------------------------------------------------------------------
      !!   0 - calculate  ah_wslp2, psix_eiv, psiy_eiv 
      !!----------------------------------------------------------------------

!!gm Future development: consider using Ah defined at T-points and attached to the 4 t-point triads

      ah_wslp2(:,:,:) = 0._wp
      IF( ln_traldf_gdia ) THEN
         psix_eiv(:,:,:) = 0._wp
         psiy_eiv(:,:,:) = 0._wp
      ENDIF

      DO ip = 0, 1
         DO kp = 0, 1
            DO jk = 1, jpkm1
               DO jj = 1, jpjm1
                  DO ji = 1, fs_jpim1
                     ze3wr = 1._wp / fse3w(ji+ip,jj,jk+kp)
                     zbu   = 0.25_wp * e1u(ji,jj) * e2u(ji,jj) * fse3u(ji,jj,jk)
                     zah   = fsahtu(ji,jj,jk)                                  !  fsaht(ji+ip,jj,jk)
                     zslope_skew = triadi_g(ji+ip,jj,jk,1-ip,kp)
                     zslope2 = zslope_skew - ( fsdept(ji+1,jj,jk) - fsdept(ji ,jj ,jk) ) * ze1ur * umask(ji,jj,jk+kp)
                     zslope2 = zslope2 *zslope2
                     ah_wslp2(ji+ip,jj,jk+kp) = ah_wslp2(ji+ip,jj,jk+kp)    &
                        &                     + zah * ( zbu * ze3wr / ( e1t(ji+ip,jj) * e2t(ji+ip,jj) ) ) * zslope2
                     IF( ln_traldf_gdia ) THEN
                        zaei_slp = fsaeiw(ji+ip,jj,jk) * zslope_skew        !fsaeit(ji+ip,jj,jk)*zslope_skew
                        psix_eiv(ji,jj,jk+kp) = psix_eiv(ji,jj,jk+kp) + 0.25_wp * zaei_slp
                     ENDIF
                  END DO
               END DO
            END DO
         END DO
      END DO
      !
      DO jp = 0, 1
         DO kp = 0, 1
            DO jk = 1, jpkm1
               DO jj = 1, jpjm1
                  DO ji=1,fs_jpim1
                     ze3wr = 1.0_wp / fse3w(ji,jj+jp,jk+kp)
                     zbv   = 0.25_wp * e1v(ji,jj) * e2v(ji,jj) * fse3v(ji,jj,jk)
                     zah   = fsahtu(ji,jj,jk)                                       !fsaht(ji,jj+jp,jk)
                     zslope_skew = triadj_g(ji,jj+jp,jk,1-jp,kp)
                     zslope2 = zslope_skew - ( fsdept(ji,jj+1,jk) - fsdept(ji,jj,jk) ) * ze2vr * vmask(ji,jj,jk+kp)
                     zslope2 = zslope2 * zslope2
                     ah_wslp2(ji,jj+jp,jk+kp) = ah_wslp2(ji,jj+jp,jk+kp)   &
                        &                     + zah * ( zbv * ze3wr / ( e1t(ji,jj+jp) * e2t(ji,jj+jp) ) ) * zslope2
                     IF( ln_traldf_gdia ) THEN
                        zaei_slp = fsaeiw(ji,jj+jp,jk) * zslope_skew     !fsaeit(ji,jj+jp,jk)*zslope_skew
                        psiy_eiv(ji,jj,jk+kp) = psiy_eiv(ji,jj,jk+kp) + 0.25_wp * zaei_slp
                     ENDIF
                  END DO
               END DO
            END DO
         END DO
      END DO
      !
      !                                                          ! ===========
      DO jn = 1, kjpt                                            ! tracer loop
         !                                                       ! ===========
         ! Zero fluxes for each tracer
         ztfw(:,:,:) = 0._wp
         zftu(:,:,:) = 0._wp
         zftv(:,:,:) = 0._wp
         !                                               
         DO jk = 1, jpkm1                          !==  before lateral T & S gradients at T-level jk  ==!
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1   ! vector opt.
                  zdit(ji,jj,jk) = ( ptb(ji+1,jj  ,jk,jn) - ptb(ji,jj,jk,jn) ) * umask(ji,jj,jk)
                  zdjt(ji,jj,jk) = ( ptb(ji  ,jj+1,jk,jn) - ptb(ji,jj,jk,jn) ) * vmask(ji,jj,jk)
               END DO
            END DO
         END DO
         IF( ln_zps ) THEN                               ! partial steps: correction at the last level
# if defined key_vectopt_loop
            DO jj = 1, 1
               DO ji = 1, jpij-jpi   ! vector opt. (forced unrolling)
# else
            DO jj = 1, jpjm1
               DO ji = 1, jpim1
# endif
                  zdit(ji,jj,mbku(ji,jj)) = pgu(ji,jj,jn)          
                  zdjt(ji,jj,mbkv(ji,jj)) = pgv(ji,jj,jn)      
               END DO
            END DO
         ENDIF

         !!----------------------------------------------------------------------
         !!   II - horizontal trend  (full)
         !!----------------------------------------------------------------------
         !
         DO jk = 1, jpkm1
            !
            !                    !==  Vertical tracer gradient at level jk and jk+1
            zdkt(:,:,1) = ( ptb(:,:,jk,jn) - ptb(:,:,jk+1,jn) ) * tmask(:,:,jk+1)
            !
            !                          ! surface boundary condition: zdkt(jk=1)=zdkt(jk=2)
            IF( jk == 1 ) THEN   ;   zdkt(:,:,0) = zdkt(:,:,1)
            ELSE                 ;   zdkt(:,:,0) = ( ptb(:,:,jk-1,jn) - ptb(:,:,jk,jn) ) * tmask(:,:,jk)
            ENDIF

            IF( .NOT. l_triad_iso ) THEN
               triadi = triadi_g
               triadj = triadj_g
            ENDIF

            DO ip = 0, 1              !==  Horizontal & vertical fluxes
               DO kp = 0, 1
                  DO jj = 1, jpjm1
                     DO ji = 1, fs_jpim1
                        ze1ur = 1._wp / e1u(ji,jj)
                        zdxt = zdit(ji,jj,jk) * ze1ur
                        ze3wr = 1._wp / fse3w(ji+ip,jj,jk+kp)
                        zdzt  = zdkt(ji+ip,jj,kp) * ze3wr
                        zslope_skew = triadi_g(ji+ip,jj,jk,1-ip,kp)
                        zslope_iso  = triadi(ji+ip,jj,jk,1-ip,kp)

                        zbu = 0.25_wp * e1u(ji,jj) * e2u(ji,jj) * fse3u(ji,jj,jk)
                        zah = fsahtu(ji,jj,jk)   !*umask(ji,jj,jk+kp)         !fsaht(ji+ip,jj,jk)           ===>>  ????
                        zah_slp  = zah * zslope_iso
                        zaei_slp = fsaeiw(ji+ip,jj,jk) * zslope_skew    !fsaeit(ji+ip,jj,jk)*zslope_skew
                        zftu(ji,jj,jk) = zftu(ji,jj,jk) - ( zah * zdxt + (zah_slp - zaei_slp) * zdzt ) * zbu * ze1ur
                        ztfw(ji+ip,jj,jk+kp) = ztfw(ji+ip,jj,jk+kp) - (zah_slp + zaei_slp) * zdxt * zbu * ze3wr
                     END DO
                  END DO
               END DO
            END DO

            DO jp = 0, 1
               DO kp = 0, 1
                  DO jj = 1, jpjm1
                     DO ji = 1, fs_jpim1
                        ze2vr = 1._wp / e2v(ji,jj)
                        zdyt = zdjt(ji,jj,jk) * ze2vr
                        ze3wr = 1._wp / fse3w(ji,jj+jp,jk+kp)
                        zdzt = zdkt(ji,jj+jp,kp) * ze3wr
                        zslope_skew = triadj_g(ji,jj+jp,jk,1-jp,kp)
                        zslope_iso  = triadj(ji,jj+jp,jk,1-jp,kp)
                        zbv = 0.25_wp * e1v(ji,jj) * e2v(ji,jj) * fse3v(ji,jj,jk)
                        zah = fsahtv(ji,jj,jk)        !*vmask(ji,jj,jk+kp)         !fsaht(ji,jj+jp,jk)
                        zah_slp = zah * zslope_iso
                        zaei_slp = fsaeiw(ji,jj+jp,jk) * zslope_skew    !fsaeit(ji,jj+jp,jk)*zslope_skew
                        zftv(ji,jj,jk) = zftv(ji,jj,jk) - ( zah * zdyt + (zah_slp - zaei_slp) * zdzt ) * zbv * ze2vr
                        ztfw(ji,jj+jp,jk+kp) = ztfw(ji,jj+jp,jk+kp) - (zah_slp + zaei_slp) * zdyt * zbv * ze3wr
                     END DO
                  END DO
               END DO
            END DO

            !                        !==  divergence and add to the general trend  ==!
            DO jj = 2 , jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zbtr = 1._wp / ( e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,jk) )
                  pta(ji,jj,jk,jn) = pta(ji,jj,jk,jn) + zbtr * (   zftu(ji-1,jj,jk) - zftu(ji,jj,jk)   &
                     &                                           + zftv(ji,jj-1,jk) - zftv(ji,jj,jk)   )
               END DO
            END DO
            !
         END DO
         !
         DO jk = 1, jpkm1            !== Divergence of vertical fluxes added to the general tracer trend
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  pta(ji,jj,jk,jn) = pta(ji,jj,jk,jn) + (  ztfw(ji,jj,jk+1) - ztfw(ji,jj,jk)  )   &
                     &                                / ( e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,jk) )
               END DO
            END DO
         END DO
         !
         !                            ! "Poleward" diffusive heat or salt transports (T-S case only)
         IF( cdtype == 'TRA' .AND. ln_diaptr .AND. ( MOD( kt, nn_fptr ) == 0 ) ) THEN
            IF( jn == jp_tem)   htr_ldf(:) = ptr_vj( zftv(:,:,:) )        ! 3.3  names
            IF( jn == jp_sal)   str_ldf(:) = ptr_vj( zftv(:,:,:) )
         ENDIF

#if defined key_diaar5
         IF( cdtype == 'TRA' .AND. jn == jp_tem  ) THEN
            z2d(:,:) = 0._wp 
            zztmp = rau0 * rcp 
            DO jk = 1, jpkm1
               DO jj = 2, jpjm1
                  DO ji = fs_2, fs_jpim1   ! vector opt.
                     z2d(ji,jj) = z2d(ji,jj) + zftu(ji,jj,jk) 
                  END DO
               END DO
            END DO
            z2d(:,:) = zztmp * z2d(:,:)
            CALL lbc_lnk( z2d, 'U', -1. )
            CALL iom_put( "udiff_heattr", z2d )                  ! heat transport in i-direction
            z2d(:,:) = 0._wp 
            DO jk = 1, jpkm1
               DO jj = 2, jpjm1
                  DO ji = fs_2, fs_jpim1   ! vector opt.
                     z2d(ji,jj) = z2d(ji,jj) + zftv(ji,jj,jk) 
                  END DO
               END DO
            END DO
            z2d(:,:) = zztmp * z2d(:,:)
            CALL lbc_lnk( z2d, 'V', -1. )
            CALL iom_put( "vdiff_heattr", z2d )                  !  heat transport in i-direction
         END IF
#endif
         !
      END DO
      !
      IF( wrk_not_released(3, 6,7,8) .OR.   &
          wrk_not_released(2, 1)       )   CALL ctl_stop('tra_ldf_iso_grif: failed to release workspace arrays')
      !
  END SUBROUTINE tra_ldf_iso_grif

#else
   !!----------------------------------------------------------------------
   !!   default option :   Dummy code   NO rotation of the diffusive tensor
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE tra_ldf_iso_grif( kt, cdtype, pgu, pgv, ptb, pta, kjpt, pahtb0 )      ! Empty routine
      CHARACTER(len=3) ::   cdtype
      REAL, DIMENSION(:,:,:) ::   pgu, pgv   ! tracer gradient at pstep levels
      REAL, DIMENSION(:,:,:,:) ::   ptb, pta
      WRITE(*,*) 'tra_ldf_iso_grif: You should not have seen this print! error?', kt, cdtype,    &
         &                  pgu(1,1,1), pgv(1,1,1), ptb(1,1,1,1), pta(1,1,1,1), kjpt, pahtb0
   END SUBROUTINE tra_ldf_iso_grif
#endif

   !!==============================================================================
END MODULE traldf_iso_grif
