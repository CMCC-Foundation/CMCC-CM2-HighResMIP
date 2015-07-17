MODULE dynhpg
   !!======================================================================
   !!                       ***  MODULE  dynhpg  ***
   !! Ocean dynamics:  hydrostatic pressure gradient trend
   !!======================================================================
   !! History :  OPA  !  1987-09  (P. Andrich, M.-A. Foujols)  hpg_zco: Original code
   !!            5.0  !  1991-11  (G. Madec)
   !!            7.0  !  1996-01  (G. Madec)  hpg_sco: Original code for s-coordinates
   !!            8.0  !  1997-05  (G. Madec)  split dynber into dynkeg and dynhpg
   !!            8.5  !  2002-07  (G. Madec)  F90: Free form and module
   !!            8.5  !  2002-08  (A. Bozec)  hpg_zps: Original code
   !!   NEMO     1.0  !  2005-10  (A. Beckmann, B.W. An)  various s-coordinate options
   !!                 !         Original code for hpg_ctl, hpg_hel hpg_wdj, hpg_djc, hpg_rot 
   !!             -   !  2005-11  (G. Madec) style & small optimisation
   !!            3.3  !  2010-10  (C. Ethe, G. Madec) reorganisation of initialisation phase
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dyn_hpg      : update the momentum trend with the now horizontal
   !!                  gradient of the hydrostatic pressure
   !!   dyn_hpg_init : initialisation and control of options
   !!       hpg_zco  : z-coordinate scheme
   !!       hpg_zps  : z-coordinate plus partial steps (interpolation)
   !!       hpg_sco  : s-coordinate (standard jacobian formulation)
   !!       hpg_hel  : s-coordinate (helsinki modification)
   !!       hpg_wdj  : s-coordinate (weighted density jacobian)
   !!       hpg_djc  : s-coordinate (Density Jacobian with Cubic polynomial)
   !!       hpg_rot  : s-coordinate (ROTated axes scheme)
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE trdmod          ! ocean dynamics trends 
   USE trdmod_oce      ! ocean variables trends
   USE in_out_manager  ! I/O manager
   USE prtctl          ! Print control
   USE lbclnk          ! lateral boundary condition 
   USE lib_mpp         ! MPP library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dyn_hpg        ! routine called by step module
   PUBLIC   dyn_hpg_init   ! routine called by opa module

   !                                              !!* Namelist namdyn_hpg : hydrostatic pressure gradient 
   LOGICAL , PUBLIC ::   ln_hpg_zco    = .TRUE.    !: z-coordinate - full steps
   LOGICAL , PUBLIC ::   ln_hpg_zps    = .FALSE.   !: z-coordinate - partial steps (interpolation)
   LOGICAL , PUBLIC ::   ln_hpg_sco    = .FALSE.   !: s-coordinate (standard jacobian formulation)
   LOGICAL , PUBLIC ::   ln_hpg_hel    = .FALSE.   !: s-coordinate (helsinki modification)
   LOGICAL , PUBLIC ::   ln_hpg_wdj    = .FALSE.   !: s-coordinate (weighted density jacobian)
   LOGICAL , PUBLIC ::   ln_hpg_djc    = .FALSE.   !: s-coordinate (Density Jacobian with Cubic polynomial)
   LOGICAL , PUBLIC ::   ln_hpg_rot    = .FALSE.   !: s-coordinate (ROTated axes scheme)
   REAL(wp), PUBLIC ::   rn_gamma      = 0._wp     !: weighting coefficient
   LOGICAL , PUBLIC ::   ln_dynhpg_imp = .FALSE.   !: semi-implicite hpg flag

   INTEGER  ::   nhpg  =  0   ! = 0 to 6, type of pressure gradient scheme used ! (deduced from ln_hpg_... flags)

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: dynhpg.F90 2715 2011-03-30 15:58:35Z rblod $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dyn_hpg( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_hpg  ***
      !!
      !! ** Method  :   Call the hydrostatic pressure gradient routine 
      !!              using the scheme defined in the namelist
      !!   
      !! ** Action : - Update (ua,va) with the now hydrastatic pressure trend
      !!             - Save the trend (l_trddyn=T)
      !!----------------------------------------------------------------------
      USE wrk_nemo, ONLY:   wrk_in_use, wrk_not_released
      USE wrk_nemo, ONLY:   ztrdu => wrk_3d_1 , ztrdv => wrk_3d_2   ! 3D workspace
      !!
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      !!----------------------------------------------------------------------
      !
      IF( wrk_in_use(3, 1,2) ) THEN
         CALL ctl_stop('dyn_hpg: requested workspace arrays are unavailable')   ;   RETURN
      ENDIF
      !
      IF( l_trddyn ) THEN                    ! Temporary saving of ua and va trends (l_trddyn)
         ztrdu(:,:,:) = ua(:,:,:)  
         ztrdv(:,:,:) = va(:,:,:) 
      ENDIF      
      !
      SELECT CASE ( nhpg )      ! Hydrastatic pressure gradient computation
      CASE (  0 )   ;   CALL hpg_zco    ( kt )      ! z-coordinate
      CASE (  1 )   ;   CALL hpg_zps    ( kt )      ! z-coordinate plus partial steps (interpolation)
      CASE (  2 )   ;   CALL hpg_sco    ( kt )      ! s-coordinate (standard jacobian formulation)
      CASE (  3 )   ;   CALL hpg_hel    ( kt )      ! s-coordinate (helsinki modification)
      CASE (  4 )   ;   CALL hpg_wdj    ( kt )      ! s-coordinate (weighted density jacobian)
      CASE (  5 )   ;   CALL hpg_djc    ( kt )      ! s-coordinate (Density Jacobian with Cubic polynomial)
      CASE (  6 )   ;   CALL hpg_rot    ( kt )      ! s-coordinate (ROTated axes scheme)
      END SELECT
      !
      IF( l_trddyn ) THEN      ! save the hydrostatic pressure gradient trends for momentum trend diagnostics
         ztrdu(:,:,:) = ua(:,:,:) - ztrdu(:,:,:)
         ztrdv(:,:,:) = va(:,:,:) - ztrdv(:,:,:)
         CALL trd_mod( ztrdu, ztrdv, jpdyn_trd_hpg, 'DYN', kt )
      ENDIF          
      !
      IF(ln_ctl)   CALL prt_ctl( tab3d_1=ua, clinfo1=' hpg  - Ua: ', mask1=umask,   &
         &                       tab3d_2=va, clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
      !
      IF( wrk_not_released(3, 1,2) )   CALL ctl_stop('dyn_hpg: failed to release workspace arrays')
      !
   END SUBROUTINE dyn_hpg


   SUBROUTINE dyn_hpg_init
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE dyn_hpg_init  ***
      !!
      !! ** Purpose :   initializations for the hydrostatic pressure gradient
      !!              computation and consistency control
      !!
      !! ** Action  :   Read the namelist namdyn_hpg and check the consistency
      !!      with the type of vertical coordinate used (zco, zps, sco)
      !!----------------------------------------------------------------------
      INTEGER ::   ioptio = 0      ! temporary integer
      !!
      NAMELIST/namdyn_hpg/ ln_hpg_zco, ln_hpg_zps, ln_hpg_sco, ln_hpg_hel,    &
         &                 ln_hpg_wdj, ln_hpg_djc, ln_hpg_rot, rn_gamma  , ln_dynhpg_imp
      !!----------------------------------------------------------------------
      !
      REWIND( numnam )               ! Read Namelist namdyn_hpg
      READ  ( numnam, namdyn_hpg )
      !
      IF(lwp) THEN                   ! Control print
         WRITE(numout,*)
         WRITE(numout,*) 'dyn_hpg_init : hydrostatic pressure gradient initialisation'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namdyn_hpg : choice of hpg scheme'
         WRITE(numout,*) '      z-coord. - full steps                             ln_hpg_zco    = ', ln_hpg_zco
         WRITE(numout,*) '      z-coord. - partial steps (interpolation)          ln_hpg_zps    = ', ln_hpg_zps
         WRITE(numout,*) '      s-coord. (standard jacobian formulation)          ln_hpg_sco    = ', ln_hpg_sco
         WRITE(numout,*) '      s-coord. (helsinki modification)                  ln_hpg_hel    = ', ln_hpg_hel
         WRITE(numout,*) '      s-coord. (weighted density jacobian)              ln_hpg_wdj    = ', ln_hpg_wdj
         WRITE(numout,*) '      s-coord. (Density Jacobian: Cubic polynomial)     ln_hpg_djc    = ', ln_hpg_djc
         WRITE(numout,*) '      s-coord. (ROTated axes scheme)                    ln_hpg_rot    = ', ln_hpg_rot
         WRITE(numout,*) '      weighting coeff. (wdj scheme)                     rn_gamma      = ', rn_gamma
         WRITE(numout,*) '      time stepping: centered (F) or semi-implicit (T)  ln_dynhpg_imp = ', ln_dynhpg_imp
      ENDIF
      !
      IF( lk_vvl .AND. .NOT. ln_hpg_sco )   &
         &   CALL ctl_stop('dyn_hpg_init : variable volume key_vvl require the standard jacobian formulation hpg_sco')
      !
      !                               ! Set nhpg from ln_hpg_... flags
      IF( ln_hpg_zco )   nhpg = 0
      IF( ln_hpg_zps )   nhpg = 1
      IF( ln_hpg_sco )   nhpg = 2
      IF( ln_hpg_hel )   nhpg = 3
      IF( ln_hpg_wdj )   nhpg = 4
      IF( ln_hpg_djc )   nhpg = 5
      IF( ln_hpg_rot )   nhpg = 6
      !
      !                               ! Consitency check
      ioptio = 0 
      IF( ln_hpg_zco )   ioptio = ioptio + 1
      IF( ln_hpg_zps )   ioptio = ioptio + 1
      IF( ln_hpg_sco )   ioptio = ioptio + 1
      IF( ln_hpg_hel )   ioptio = ioptio + 1
      IF( ln_hpg_wdj )   ioptio = ioptio + 1
      IF( ln_hpg_djc )   ioptio = ioptio + 1
      IF( ln_hpg_rot )   ioptio = ioptio + 1
      IF( ioptio /= 1 )   CALL ctl_stop( 'NO or several hydrostatic pressure gradient options used' )
      !
   END SUBROUTINE dyn_hpg_init


   SUBROUTINE hpg_zco( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE hpg_zco  ***
      !!
      !! ** Method  :   z-coordinate case, levels are horizontal surfaces.
      !!      The now hydrostatic pressure gradient at a given level, jk,
      !!      is computed by taking the vertical integral of the in-situ
      !!      density gradient along the model level from the suface to that
      !!      level:    zhpi = grav .....
      !!                zhpj = grav .....
      !!      add it to the general momentum trend (ua,va).
      !!            ua = ua - 1/e1u * zhpi
      !!            va = va - 1/e2v * zhpj
      !! 
      !! ** Action : - Update (ua,va) with the now hydrastatic pressure trend
      !!----------------------------------------------------------------------
      USE oce, ONLY:   zhpi => ta , zhpj => sa   ! (ta,sa) used as 3D workspace
      !!
      INTEGER, INTENT(in) ::   kt    ! ocean time-step index
      !!
      INTEGER  ::   ji, jj, jk       ! dummy loop indices
      REAL(wp) ::   zcoef0, zcoef1   ! temporary scalars
      !!----------------------------------------------------------------------
      
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn:hpg_zco : hydrostatic pressure gradient trend'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~   z-coordinate case '
      ENDIF
      
      zcoef0 = - grav * 0.5_wp      ! Local constant initialization 

      ! Surface value
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1   ! vector opt.
            zcoef1 = zcoef0 * fse3w(ji,jj,1)
            ! hydrostatic pressure gradient
            zhpi(ji,jj,1) = zcoef1 * ( rhd(ji+1,jj,1) - rhd(ji,jj,1) ) / e1u(ji,jj)
            zhpj(ji,jj,1) = zcoef1 * ( rhd(ji,jj+1,1) - rhd(ji,jj,1) ) / e2v(ji,jj)
            ! add to the general momentum trend
            ua(ji,jj,1) = ua(ji,jj,1) + zhpi(ji,jj,1)
            va(ji,jj,1) = va(ji,jj,1) + zhpj(ji,jj,1)
         END DO
      END DO
      !
      ! interior value (2=<jk=<jpkm1)
      DO jk = 2, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zcoef1 = zcoef0 * fse3w(ji,jj,jk)
               ! hydrostatic pressure gradient
               zhpi(ji,jj,jk) = zhpi(ji,jj,jk-1)   &
                  &           + zcoef1 * (  ( rhd(ji+1,jj,jk)+rhd(ji+1,jj,jk-1) )   &
                  &                       - ( rhd(ji  ,jj,jk)+rhd(ji  ,jj,jk-1) )  ) / e1u(ji,jj)

               zhpj(ji,jj,jk) = zhpj(ji,jj,jk-1)   &
                  &           + zcoef1 * (  ( rhd(ji,jj+1,jk)+rhd(ji,jj+1,jk-1) )   &
                  &                       - ( rhd(ji,jj,  jk)+rhd(ji,jj  ,jk-1) )  ) / e2v(ji,jj)
               ! add to the general momentum trend
               ua(ji,jj,jk) = ua(ji,jj,jk) + zhpi(ji,jj,jk)
               va(ji,jj,jk) = va(ji,jj,jk) + zhpj(ji,jj,jk)
            END DO
         END DO
      END DO
      !
   END SUBROUTINE hpg_zco


   SUBROUTINE hpg_zps( kt )
      !!---------------------------------------------------------------------
      !!                 ***  ROUTINE hpg_zps  ***
      !!                    
      !! ** Method  :   z-coordinate plus partial steps case.  blahblah...
      !! 
      !! ** Action  : - Update (ua,va) with the now hydrastatic pressure trend
      !!---------------------------------------------------------------------- 
      USE oce, ONLY:   zhpi => ta , zhpj => sa   ! (ta,sa) used as 3D workspace
      !!
      INTEGER, INTENT(in) ::   kt    ! ocean time-step index
      !!
      INTEGER  ::   ji, jj, jk                       ! dummy loop indices
      INTEGER  ::   iku, ikv                         ! temporary integers
      REAL(wp) ::   zcoef0, zcoef1, zcoef2, zcoef3   ! temporary scalars
      !!----------------------------------------------------------------------

      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn:hpg_zps : hydrostatic pressure gradient trend'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~   z-coordinate with partial steps - vector optimization'
      ENDIF

      ! Local constant initialization
      zcoef0 = - grav * 0.5_wp

      !  Surface value (also valid in partial step case)
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1   ! vector opt.
            zcoef1 = zcoef0 * fse3w(ji,jj,1)
            ! hydrostatic pressure gradient
            zhpi(ji,jj,1) = zcoef1 * ( rhd(ji+1,jj  ,1) - rhd(ji,jj,1) ) / e1u(ji,jj)
            zhpj(ji,jj,1) = zcoef1 * ( rhd(ji  ,jj+1,1) - rhd(ji,jj,1) ) / e2v(ji,jj)
            ! add to the general momentum trend
            ua(ji,jj,1) = ua(ji,jj,1) + zhpi(ji,jj,1)
            va(ji,jj,1) = va(ji,jj,1) + zhpj(ji,jj,1)
         END DO
      END DO

      ! interior value (2=<jk=<jpkm1)
      DO jk = 2, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zcoef1 = zcoef0 * fse3w(ji,jj,jk)
               ! hydrostatic pressure gradient
               zhpi(ji,jj,jk) = zhpi(ji,jj,jk-1)   &
                  &           + zcoef1 * (  ( rhd(ji+1,jj,jk) + rhd(ji+1,jj,jk-1) )   &
                  &                       - ( rhd(ji  ,jj,jk) + rhd(ji  ,jj,jk-1) )  ) / e1u(ji,jj)

               zhpj(ji,jj,jk) = zhpj(ji,jj,jk-1)   &
                  &           + zcoef1 * (  ( rhd(ji,jj+1,jk) + rhd(ji,jj+1,jk-1) )   &
                  &                       - ( rhd(ji,jj,  jk) + rhd(ji,jj  ,jk-1) )  ) / e2v(ji,jj)
               ! add to the general momentum trend
               ua(ji,jj,jk) = ua(ji,jj,jk) + zhpi(ji,jj,jk)
               va(ji,jj,jk) = va(ji,jj,jk) + zhpj(ji,jj,jk)
            END DO
         END DO
      END DO

      ! partial steps correction at the last level  (use gru & grv computed in zpshde.F90)
# if defined key_vectopt_loop
         jj = 1
         DO ji = jpi+2, jpij-jpi-1   ! vector opt. (forced unrolling)
# else
      DO jj = 2, jpjm1
         DO ji = 2, jpim1
# endif
            iku = mbku(ji,jj)
            ikv = mbkv(ji,jj)
            zcoef2 = zcoef0 * MIN( fse3w(ji,jj,iku), fse3w(ji+1,jj  ,iku) )
            zcoef3 = zcoef0 * MIN( fse3w(ji,jj,ikv), fse3w(ji  ,jj+1,ikv) )
            IF( iku > 1 ) THEN            ! on i-direction (level 2 or more)
               ua  (ji,jj,iku) = ua(ji,jj,iku) - zhpi(ji,jj,iku)         ! subtract old value
               zhpi(ji,jj,iku) = zhpi(ji,jj,iku-1)                   &   ! compute the new one
                  &            + zcoef2 * ( rhd(ji+1,jj,iku-1) - rhd(ji,jj,iku-1) + gru(ji,jj) ) / e1u(ji,jj)
               ua  (ji,jj,iku) = ua(ji,jj,iku) + zhpi(ji,jj,iku)         ! add the new one to the general momentum trend
            ENDIF
            IF( ikv > 1 ) THEN            ! on j-direction (level 2 or more)
               va  (ji,jj,ikv) = va(ji,jj,ikv) - zhpj(ji,jj,ikv)         ! subtract old value
               zhpj(ji,jj,ikv) = zhpj(ji,jj,ikv-1)                   &   ! compute the new one
                  &            + zcoef3 * ( rhd(ji,jj+1,ikv-1) - rhd(ji,jj,ikv-1) + grv(ji,jj) ) / e2v(ji,jj)
               va  (ji,jj,ikv) = va(ji,jj,ikv) + zhpj(ji,jj,ikv)         ! add the new one to the general momentum trend
            ENDIF
# if ! defined key_vectopt_loop
         END DO
# endif
      END DO
      !
   END SUBROUTINE hpg_zps


   SUBROUTINE hpg_sco( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE hpg_sco  ***
      !!
      !! ** Method  :   s-coordinate case. Jacobian scheme.
      !!      The now hydrostatic pressure gradient at a given level, jk,
      !!      is computed by taking the vertical integral of the in-situ
      !!      density gradient along the model level from the suface to that
      !!      level. s-coordinates (ln_sco): a corrective term is added
      !!      to the horizontal pressure gradient :
      !!         zhpi = grav .....  + 1/e1u mi(rhd) di[ grav dep3w ]
      !!         zhpj = grav .....  + 1/e2v mj(rhd) dj[ grav dep3w ]
      !!      add it to the general momentum trend (ua,va).
      !!         ua = ua - 1/e1u * zhpi
      !!         va = va - 1/e2v * zhpj
      !!
      !! ** Action : - Update (ua,va) with the now hydrastatic pressure trend
      !!----------------------------------------------------------------------
      USE oce, ONLY:   zhpi => ta , zhpj => sa   ! (ta,sa) used as 3D workspace
      !!
      INTEGER, INTENT(in) ::   kt    ! ocean time-step index
      !!
      INTEGER  ::   ji, jj, jk                 ! dummy loop indices
      REAL(wp) ::   zcoef0, zuap, zvap, znad   ! temporary scalars
      !!----------------------------------------------------------------------

      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn:hpg_sco : hydrostatic pressure gradient trend'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~   s-coordinate case, OPA original scheme used'
      ENDIF

      ! Local constant initialization
      zcoef0 = - grav * 0.5_wp
      ! To use density and not density anomaly
      IF ( lk_vvl ) THEN   ;     znad = 1._wp          ! Variable volume
      ELSE                 ;     znad = 0._wp         ! Fixed volume
      ENDIF

      ! Surface value
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1   ! vector opt.   
            ! hydrostatic pressure gradient along s-surfaces
            zhpi(ji,jj,1) = zcoef0 / e1u(ji,jj) * ( fse3w(ji+1,jj  ,1) * ( znad + rhd(ji+1,jj  ,1) )   &
               &                                  - fse3w(ji  ,jj  ,1) * ( znad + rhd(ji  ,jj  ,1) ) )
            zhpj(ji,jj,1) = zcoef0 / e2v(ji,jj) * ( fse3w(ji  ,jj+1,1) * ( znad + rhd(ji  ,jj+1,1) )   &
               &                                  - fse3w(ji  ,jj  ,1) * ( znad + rhd(ji  ,jj  ,1) ) )
            ! s-coordinate pressure gradient correction
            zuap = -zcoef0 * ( rhd   (ji+1,jj,1) + rhd   (ji,jj,1) + 2._wp * znad )   &
               &           * ( fsde3w(ji+1,jj,1) - fsde3w(ji,jj,1) ) / e1u(ji,jj)
            zvap = -zcoef0 * ( rhd   (ji,jj+1,1) + rhd   (ji,jj,1) + 2._wp * znad )   &
               &           * ( fsde3w(ji,jj+1,1) - fsde3w(ji,jj,1) ) / e2v(ji,jj)
            ! add to the general momentum trend
            ua(ji,jj,1) = ua(ji,jj,1) + zhpi(ji,jj,1) + zuap
            va(ji,jj,1) = va(ji,jj,1) + zhpj(ji,jj,1) + zvap
         END DO  
      END DO   
            
      ! interior value (2=<jk=<jpkm1)
      DO jk = 2, jpkm1                                  
         DO jj = 2, jpjm1     
            DO ji = fs_2, fs_jpim1   ! vector opt.      
               ! hydrostatic pressure gradient along s-surfaces
               zhpi(ji,jj,jk) = zhpi(ji,jj,jk-1) + zcoef0 / e1u(ji,jj)   & 
                  &           * (  fse3w(ji+1,jj,jk) * ( rhd(ji+1,jj,jk) + rhd(ji+1,jj,jk-1) + 2*znad )   & 
                  &              - fse3w(ji  ,jj,jk) * ( rhd(ji  ,jj,jk) + rhd(ji  ,jj,jk-1) + 2*znad )  )
               zhpj(ji,jj,jk) = zhpj(ji,jj,jk-1) + zcoef0 / e2v(ji,jj)   &
                  &           * (  fse3w(ji,jj+1,jk) * ( rhd(ji,jj+1,jk) + rhd(ji,jj+1,jk-1) + 2*znad )   &
                  &              - fse3w(ji,jj  ,jk) * ( rhd(ji,jj,  jk) + rhd(ji,jj  ,jk-1) + 2*znad )  )
               ! s-coordinate pressure gradient correction
               zuap = -zcoef0 * ( rhd   (ji+1,jj  ,jk) + rhd   (ji,jj,jk) + 2._wp * znad )   &
                  &           * ( fsde3w(ji+1,jj  ,jk) - fsde3w(ji,jj,jk) ) / e1u(ji,jj)
               zvap = -zcoef0 * ( rhd   (ji  ,jj+1,jk) + rhd   (ji,jj,jk) + 2._wp * znad )   &
                  &           * ( fsde3w(ji  ,jj+1,jk) - fsde3w(ji,jj,jk) ) / e2v(ji,jj)
               ! add to the general momentum trend
               ua(ji,jj,jk) = ua(ji,jj,jk) + zhpi(ji,jj,jk) + zuap
               va(ji,jj,jk) = va(ji,jj,jk) + zhpj(ji,jj,jk) + zvap
            END DO
         END DO
      END DO
      !
   END SUBROUTINE hpg_sco


   SUBROUTINE hpg_hel( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE hpg_hel  ***
      !!
      !! ** Method  :   s-coordinate case.
      !!      The now hydrostatic pressure gradient at a given level
      !!      jk is computed by taking the vertical integral of the in-situ 
      !!      density gradient along the model level from the suface to that 
      !!      level. s-coordinates (ln_sco): a corrective term is added
      !!      to the horizontal pressure gradient :
      !!         zhpi = grav .....  + 1/e1u mi(rhd) di[ grav dep3w ]
      !!         zhpj = grav .....  + 1/e2v mj(rhd) dj[ grav dep3w ]
      !!      add it to the general momentum trend (ua,va).
      !!         ua = ua - 1/e1u * zhpi
      !!         va = va - 1/e2v * zhpj
      !!
      !! ** Action : - Update (ua,va) with the now hydrastatic pressure trend
      !!             - Save the trend (l_trddyn=T)
      !!----------------------------------------------------------------------
      USE oce, ONLY:   zhpi => ta , zhpj => sa   ! (ta,sa) used as 3D workspace
      !!
      INTEGER, INTENT(in) ::   kt    ! ocean time-step index
      !!
      INTEGER  ::   ji, jj, jk           ! dummy loop indices
      REAL(wp) ::   zcoef0, zuap, zvap   ! temporary scalars
      !!----------------------------------------------------------------------

      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn:hpg_hel : hydrostatic pressure gradient trend'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~   s-coordinate case, helsinki modified scheme'
      ENDIF

      ! Local constant initialization
      zcoef0 = - grav * 0.5_wp
 
      ! Surface value
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1   ! vector opt.
            ! hydrostatic pressure gradient along s-surfaces
            zhpi(ji,jj,1) = zcoef0 / e1u(ji,jj) * ( fse3t(ji+1,jj  ,1) * rhd(ji+1,jj  ,1)  &
               &                                  - fse3t(ji  ,jj  ,1) * rhd(ji  ,jj  ,1) )
            zhpj(ji,jj,1) = zcoef0 / e2v(ji,jj) * ( fse3t(ji  ,jj+1,1) * rhd(ji  ,jj+1,1)  &
               &                                  - fse3t(ji  ,jj  ,1) * rhd(ji  ,jj  ,1) )
            ! s-coordinate pressure gradient correction
            zuap = -zcoef0 * ( rhd   (ji+1,jj,1) + rhd   (ji,jj,1) )   &
               &           * ( fsdept(ji+1,jj,1) - fsdept(ji,jj,1) ) / e1u(ji,jj)
            zvap = -zcoef0 * ( rhd   (ji,jj+1,1) + rhd   (ji,jj,1) )   &
               &           * ( fsdept(ji,jj+1,1) - fsdept(ji,jj,1) ) / e2v(ji,jj)
            ! add to the general momentum trend
            ua(ji,jj,1) = ua(ji,jj,1) + zhpi(ji,jj,1) + zuap
            va(ji,jj,1) = va(ji,jj,1) + zhpj(ji,jj,1) + zvap
         END DO
      END DO
      !
      ! interior value (2=<jk=<jpkm1)
      DO jk = 2, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               ! hydrostatic pressure gradient along s-surfaces
               zhpi(ji,jj,jk) = zhpi(ji,jj,jk-1) &
                  &           +  zcoef0 / e1u(ji,jj) * ( fse3t(ji+1,jj,jk  ) * rhd(ji+1,jj,jk)     &
                  &                                     -fse3t(ji  ,jj,jk  ) * rhd(ji  ,jj,jk)   ) &
                  &           +  zcoef0 / e1u(ji,jj) * ( fse3t(ji+1,jj,jk-1) * rhd(ji+1,jj,jk-1)   &
                  &                                     -fse3t(ji  ,jj,jk-1) * rhd(ji  ,jj,jk-1) )
               zhpj(ji,jj,jk) = zhpj(ji,jj,jk-1) &
                  &           +  zcoef0 / e2v(ji,jj) * ( fse3t(ji,jj+1,jk  ) * rhd(ji,jj+1,jk)     &
                  &                                     -fse3t(ji,jj  ,jk  ) * rhd(ji,jj,  jk)   ) &
                  &           +  zcoef0 / e2v(ji,jj) * ( fse3t(ji,jj+1,jk-1) * rhd(ji,jj+1,jk-1)   &
                  &                                     -fse3t(ji,jj  ,jk-1) * rhd(ji,jj,  jk-1) )
               ! s-coordinate pressure gradient correction
               zuap = - zcoef0 * ( rhd   (ji+1,jj,jk) + rhd   (ji,jj,jk) )   &
                  &            * ( fsdept(ji+1,jj,jk) - fsdept(ji,jj,jk) ) / e1u(ji,jj)
               zvap = - zcoef0 * ( rhd   (ji,jj+1,jk) + rhd   (ji,jj,jk) )   &
                  &            * ( fsdept(ji,jj+1,jk) - fsdept(ji,jj,jk) ) / e2v(ji,jj)
               ! add to the general momentum trend
               ua(ji,jj,jk) = ua(ji,jj,jk) + zhpi(ji,jj,jk) + zuap
               va(ji,jj,jk) = va(ji,jj,jk) + zhpj(ji,jj,jk) + zvap
            END DO
         END DO
      END DO
      !
   END SUBROUTINE hpg_hel


   SUBROUTINE hpg_wdj( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE hpg_wdj  ***
      !!
      !! ** Method  :   Weighted Density Jacobian (wdj) scheme (song 1998)
      !!      The weighting coefficients from the namelist parameter rn_gamma
      !!      (alpha=0.5-rn_gamma ; beta=1-alpha=0.5+rn_gamma
      !!
      !! Reference : Song, Mon. Wea. Rev., 126, 3213-3230, 1998.
      !!----------------------------------------------------------------------
      USE oce, ONLY:   zhpi => ta , zhpj => sa   ! (ta,sa) used as 3D workspace
      !!
      INTEGER, INTENT(in) ::   kt    ! ocean time-step index
      !!
      INTEGER  ::   ji, jj, jk           ! dummy loop indices
      REAL(wp) ::   zcoef0, zuap, zvap   ! temporary scalars
      REAL(wp) ::   zalph , zbeta        !    "         "
      !!----------------------------------------------------------------------

      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn:hpg_wdj : hydrostatic pressure gradient trend'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~   Weighted Density Jacobian'
      ENDIF

      ! Local constant initialization
      zcoef0 = - grav * 0.5_wp
      zalph  = 0.5_wp - rn_gamma    ! weighting coefficients (alpha=0.5-rn_gamma
      zbeta  = 0.5_wp + rn_gamma    !                        (beta =1-alpha=0.5+rn_gamma

      ! Surface value (no ponderation)
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1   ! vector opt.
            ! hydrostatic pressure gradient along s-surfaces
            zhpi(ji,jj,1) = zcoef0 / e1u(ji,jj) * (  fse3w(ji+1,jj  ,1) * rhd(ji+1,jj  ,1)   &
               &                                   - fse3w(ji  ,jj  ,1) * rhd(ji  ,jj  ,1)  )
            zhpj(ji,jj,1) = zcoef0 / e2v(ji,jj) * (  fse3w(ji  ,jj+1,1) * rhd(ji  ,jj+1,1)   &
               &                                   - fse3w(ji  ,jj  ,1) * rhd(ji,  jj  ,1)  )
            ! s-coordinate pressure gradient correction
            zuap = -zcoef0 * ( rhd   (ji+1,jj,1) + rhd   (ji,jj,1) )   &
               &           * ( fsde3w(ji+1,jj,1) - fsde3w(ji,jj,1) ) / e1u(ji,jj)
            zvap = -zcoef0 * ( rhd   (ji,jj+1,1) + rhd   (ji,jj,1) )   &
               &           * ( fsde3w(ji,jj+1,1) - fsde3w(ji,jj,1) ) / e2v(ji,jj)
            ! add to the general momentum trend
            ua(ji,jj,1) = ua(ji,jj,1) + zhpi(ji,jj,1) + zuap
            va(ji,jj,1) = va(ji,jj,1) + zhpj(ji,jj,1) + zvap
         END DO
      END DO

      ! Interior value (2=<jk=<jpkm1) (weighted with zalph & zbeta)
      DO jk = 2, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zhpi(ji,jj,jk) = zhpi(ji,jj,jk-1) + zcoef0 / e1u(ji,jj)                            &
                  &           * (   (            fsde3w(ji+1,jj,jk  ) + fsde3w(ji,jj,jk  )        &
                  &                            - fsde3w(ji+1,jj,jk-1) - fsde3w(ji,jj,jk-1)    )   &
                  &               * (  zalph * ( rhd   (ji+1,jj,jk-1) - rhd   (ji,jj,jk-1) )      &
                  &                  + zbeta * ( rhd   (ji+1,jj,jk  ) - rhd   (ji,jj,jk  ) )  )   &
                  &             -   (            rhd   (ji+1,jj,jk  ) + rhd   (ji,jj,jk  )        &
                  &                           - rhd   (ji+1,jj,jk-1) - rhd   (ji,jj,jk-1)     )   &
                  &               * (  zalph * ( fsde3w(ji+1,jj,jk-1) - fsde3w(ji,jj,jk-1) )      &
                  &                  + zbeta * ( fsde3w(ji+1,jj,jk  ) - fsde3w(ji,jj,jk  ) )  )  )
               zhpj(ji,jj,jk) = zhpj(ji,jj,jk-1) + zcoef0 / e2v(ji,jj)                            &
                  &           * (   (           fsde3w(ji,jj+1,jk  ) + fsde3w(ji,jj,jk  )         &
                  &                           - fsde3w(ji,jj+1,jk-1) - fsde3w(ji,jj,jk-1)     )   &
                  &               * (  zalph * ( rhd   (ji,jj+1,jk-1) - rhd   (ji,jj,jk-1) )      &
                  &                  + zbeta * ( rhd   (ji,jj+1,jk  ) - rhd   (ji,jj,jk  ) )  )   &
                  &             -   (            rhd   (ji,jj+1,jk  ) + rhd   (ji,jj,jk  )        &
                  &                            - rhd   (ji,jj+1,jk-1) - rhd   (ji,jj,jk-1)    )   &
                  &               * (  zalph * ( fsde3w(ji,jj+1,jk-1) - fsde3w(ji,jj,jk-1) )      &
                  &                  + zbeta * ( fsde3w(ji,jj+1,jk  ) - fsde3w(ji,jj,jk  ) )  )  )
               ! add to the general momentum trend
               ua(ji,jj,jk) = ua(ji,jj,jk) + zhpi(ji,jj,jk)
               va(ji,jj,jk) = va(ji,jj,jk) + zhpj(ji,jj,jk)
            END DO
         END DO
      END DO
      !
   END SUBROUTINE hpg_wdj


   SUBROUTINE hpg_djc( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE hpg_djc  ***
      !!
      !! ** Method  :   Density Jacobian with Cubic polynomial scheme
      !! 
      !! Reference: Shchepetkin and McWilliams, J. Geophys. Res., 108(C3), 3090, 2003
      !!----------------------------------------------------------------------
      USE wrk_nemo, ONLY:   wrk_in_use, wrk_not_released
      USE oce     , ONLY:   zhpi  => ta        , zhpj => sa       ! (ta,sa) used as 3D workspace
      USE wrk_nemo, ONLY:   drhox => wrk_3d_1  , dzx  => wrk_3d_2
      USE wrk_nemo, ONLY:   drhou => wrk_3d_3  , dzu  => wrk_3d_4 , rho_i => wrk_3d_5
      USE wrk_nemo, ONLY:   drhoy => wrk_3d_6  , dzy  => wrk_3d_7
      USE wrk_nemo, ONLY:   drhov => wrk_3d_8  , dzv  => wrk_3d_9 , rho_j => wrk_3d_10
      USE wrk_nemo, ONLY:   drhoz => wrk_3d_11 , dzz  => wrk_3d_12 
      USE wrk_nemo, ONLY:   drhow => wrk_3d_13 , dzw  => wrk_3d_14
      USE wrk_nemo, ONLY:   rho_k => wrk_3d_15
      !!
      INTEGER, INTENT(in) ::   kt    ! ocean time-step index
      !!
      INTEGER  ::   ji, jj, jk          ! dummy loop indices
      REAL(wp) ::   zcoef0, zep, cffw   ! temporary scalars
      REAL(wp) ::   z1_10, cffu, cffx   !    "         "
      REAL(wp) ::   z1_12, cffv, cffy   !    "         "
      !!----------------------------------------------------------------------

      IF( wrk_in_use(3, 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15) ) THEN
         CALL ctl_stop('dyn:hpg_djc: requested workspace arrays unavailable')   ;   RETURN
      ENDIF

      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn:hpg_djc : hydrostatic pressure gradient trend'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~   s-coordinate case, density Jacobian with cubic polynomial scheme'
      ENDIF

      ! Local constant initialization
      zcoef0 = - grav * 0.5_wp
      z1_10  = 1._wp / 10._wp
      z1_12  = 1._wp / 12._wp

      !----------------------------------------------------------------------------------------
      !  compute and store in provisional arrays elementary vertical and horizontal differences
      !----------------------------------------------------------------------------------------

!!bug gm   Not a true bug, but... dzz=e3w  for dzx, dzy verify what it is really

      DO jk = 2, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               drhoz(ji,jj,jk) = rhd   (ji  ,jj  ,jk) - rhd   (ji,jj,jk-1)
               dzz  (ji,jj,jk) = fsde3w(ji  ,jj  ,jk) - fsde3w(ji,jj,jk-1)
               drhox(ji,jj,jk) = rhd   (ji+1,jj  ,jk) - rhd   (ji,jj,jk  )
               dzx  (ji,jj,jk) = fsde3w(ji+1,jj  ,jk) - fsde3w(ji,jj,jk  )
               drhoy(ji,jj,jk) = rhd   (ji  ,jj+1,jk) - rhd   (ji,jj,jk  )
               dzy  (ji,jj,jk) = fsde3w(ji  ,jj+1,jk) - fsde3w(ji,jj,jk  )
            END DO
         END DO
      END DO

      !-------------------------------------------------------------------------
      ! compute harmonic averages using eq. 5.18
      !-------------------------------------------------------------------------
      zep = 1.e-15

!!bug  gm  drhoz not defined at level 1 and used (jk-1 with jk=2)
!!bug  gm  idem for drhox, drhoy et ji=jpi and jj=jpj

      DO jk = 2, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               cffw = 2._wp * drhoz(ji  ,jj  ,jk) * drhoz(ji,jj,jk-1)

               cffu = 2._wp * drhox(ji+1,jj  ,jk) * drhox(ji,jj,jk  )
               cffx = 2._wp * dzx  (ji+1,jj  ,jk) * dzx  (ji,jj,jk  )
  
               cffv = 2._wp * drhoy(ji  ,jj+1,jk) * drhoy(ji,jj,jk  )
               cffy = 2._wp * dzy  (ji  ,jj+1,jk) * dzy  (ji,jj,jk  )

               IF( cffw > zep) THEN
                  drhow(ji,jj,jk) = 2._wp *   drhoz(ji,jj,jk) * drhoz(ji,jj,jk-1)   &
                     &                    / ( drhoz(ji,jj,jk) + drhoz(ji,jj,jk-1) )
               ELSE
                  drhow(ji,jj,jk) = 0._wp
               ENDIF

               dzw(ji,jj,jk) = 2._wp *   dzz(ji,jj,jk) * dzz(ji,jj,jk-1)   &
                  &                  / ( dzz(ji,jj,jk) + dzz(ji,jj,jk-1) )

               IF( cffu > zep ) THEN
                  drhou(ji,jj,jk) = 2._wp *   drhox(ji+1,jj,jk) * drhox(ji,jj,jk)   &
                     &                    / ( drhox(ji+1,jj,jk) + drhox(ji,jj,jk) )
               ELSE
                  drhou(ji,jj,jk ) = 0._wp
               ENDIF

               IF( cffx > zep ) THEN
                  dzu(ji,jj,jk) = 2._wp *   dzx(ji+1,jj,jk) * dzx(ji,jj,jk)   &
                     &                  / ( dzx(ji+1,jj,jk) + dzx(ji,jj,jk) )
               ELSE
                  dzu(ji,jj,jk) = 0._wp
               ENDIF

               IF( cffv > zep ) THEN
                  drhov(ji,jj,jk) = 2._wp *   drhoy(ji,jj+1,jk) * drhoy(ji,jj,jk)   &
                     &                    / ( drhoy(ji,jj+1,jk) + drhoy(ji,jj,jk) )
               ELSE
                  drhov(ji,jj,jk) = 0._wp
               ENDIF

               IF( cffy > zep ) THEN
                  dzv(ji,jj,jk) = 2._wp *   dzy(ji,jj+1,jk) * dzy(ji,jj,jk)   &
                     &                  / ( dzy(ji,jj+1,jk) + dzy(ji,jj,jk) )
               ELSE
                  dzv(ji,jj,jk) = 0._wp
               ENDIF

            END DO
         END DO
      END DO

      !----------------------------------------------------------------------------------
      ! apply boundary conditions at top and bottom using 5.36-5.37
      !----------------------------------------------------------------------------------
      drhow(:,:, 1 ) = 1.5_wp * ( drhoz(:,:, 2 ) - drhoz(:,:,  1  ) ) - 0.5_wp * drhow(:,:,  2  )
      drhou(:,:, 1 ) = 1.5_wp * ( drhox(:,:, 2 ) - drhox(:,:,  1  ) ) - 0.5_wp * drhou(:,:,  2  )
      drhov(:,:, 1 ) = 1.5_wp * ( drhoy(:,:, 2 ) - drhoy(:,:,  1  ) ) - 0.5_wp * drhov(:,:,  2  )

      drhow(:,:,jpk) = 1.5_wp * ( drhoz(:,:,jpk) - drhoz(:,:,jpkm1) ) - 0.5_wp * drhow(:,:,jpkm1)
      drhou(:,:,jpk) = 1.5_wp * ( drhox(:,:,jpk) - drhox(:,:,jpkm1) ) - 0.5_wp * drhou(:,:,jpkm1)
      drhov(:,:,jpk) = 1.5_wp * ( drhoy(:,:,jpk) - drhoy(:,:,jpkm1) ) - 0.5_wp * drhov(:,:,jpkm1)


      !--------------------------------------------------------------
      ! Upper half of top-most grid box, compute and store
      !-------------------------------------------------------------

!!bug gm   :  e3w-de3w = 0.5*e3w  ....  and de3w(2)-de3w(1)=e3w(2) ....   to be verified
!          true if de3w is really defined as the sum of the e3w scale factors as, it seems to me, it should be

      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1   ! vector opt.
            rho_k(ji,jj,1) = -grav * ( fse3w(ji,jj,1) - fsde3w(ji,jj,1) )               &
               &                   * (  rhd(ji,jj,1)                                    &
               &                     + 0.5_wp * ( rhd(ji,jj,2) - rhd(ji,jj,1) )         &
               &                              * ( fse3w (ji,jj,1) - fsde3w(ji,jj,1) )   &
               &                              / ( fsde3w(ji,jj,2) - fsde3w(ji,jj,1) )  ) 
         END DO
      END DO

!!bug gm    : here also, simplification is possible
!!bug gm    : optimisation: 1/10 and 1/12 the division should be done before the loop

      DO jk = 2, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.

               rho_k(ji,jj,jk) = zcoef0 * ( rhd   (ji,jj,jk) + rhd   (ji,jj,jk-1) )                                   &
                  &                     * ( fsde3w(ji,jj,jk) - fsde3w(ji,jj,jk-1) )                                   &
                  &            - grav * z1_10 * (                                                                     &
                  &     ( drhow (ji,jj,jk) - drhow (ji,jj,jk-1) )                                                     &
                  &   * ( fsde3w(ji,jj,jk) - fsde3w(ji,jj,jk-1) - z1_12 * ( dzw  (ji,jj,jk) + dzw  (ji,jj,jk-1) ) )   &
                  &   - ( dzw   (ji,jj,jk) - dzw   (ji,jj,jk-1) )                                                     &
                  &   * ( rhd   (ji,jj,jk) - rhd   (ji,jj,jk-1) - z1_12 * ( drhow(ji,jj,jk) + drhow(ji,jj,jk-1) ) )   &
                  &                             )

               rho_i(ji,jj,jk) = zcoef0 * ( rhd   (ji+1,jj,jk) + rhd   (ji,jj,jk) )                                   &
                  &                     * ( fsde3w(ji+1,jj,jk) - fsde3w(ji,jj,jk) )                                   &
                  &            - grav* z1_10 * (                                                                      &
                  &     ( drhou (ji+1,jj,jk) - drhou (ji,jj,jk) )                                                     &
                  &   * ( fsde3w(ji+1,jj,jk) - fsde3w(ji,jj,jk) - z1_12 * ( dzu  (ji+1,jj,jk) + dzu  (ji,jj,jk) ) )   &
                  &   - ( dzu   (ji+1,jj,jk) - dzu   (ji,jj,jk) )                                                     &
                  &   * ( rhd   (ji+1,jj,jk) - rhd   (ji,jj,jk) - z1_12 * ( drhou(ji+1,jj,jk) + drhou(ji,jj,jk) ) )   &
                  &                            )

               rho_j(ji,jj,jk) = zcoef0 * ( rhd   (ji,jj+1,jk) + rhd   (ji,jj,jk) )                                   &
                  &                     * ( fsde3w(ji,jj+1,jk) - fsde3w(ji,jj,jk) )                                   &
                  &            - grav* z1_10 * (                                                                      &
                  &     ( drhov (ji,jj+1,jk) - drhov (ji,jj,jk) )                                                     &
                  &   * ( fsde3w(ji,jj+1,jk) - fsde3w(ji,jj,jk) - z1_12 * ( dzv  (ji,jj+1,jk) + dzv  (ji,jj,jk) ) )   &
                  &   - ( dzv   (ji,jj+1,jk) - dzv   (ji,jj,jk) )                                                     &
                  &   * ( rhd   (ji,jj+1,jk) - rhd   (ji,jj,jk) - z1_12 * ( drhov(ji,jj+1,jk) + drhov(ji,jj,jk) ) )   &
                  &                            )

            END DO
         END DO
      END DO
      CALL lbc_lnk(rho_k,'W',1.)
      CALL lbc_lnk(rho_i,'U',1.)
      CALL lbc_lnk(rho_j,'V',1.)


      ! ---------------
      !  Surface value
      ! ---------------
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1   ! vector opt.
            zhpi(ji,jj,1) = ( rho_k(ji+1,jj  ,1) - rho_k(ji,jj,1) - rho_i(ji,jj,1) ) / e1u(ji,jj)
            zhpj(ji,jj,1) = ( rho_k(ji  ,jj+1,1) - rho_k(ji,jj,1) - rho_j(ji,jj,1) ) / e2v(ji,jj)
            ! add to the general momentum trend
            ua(ji,jj,1) = ua(ji,jj,1) + zhpi(ji,jj,1)
            va(ji,jj,1) = va(ji,jj,1) + zhpj(ji,jj,1)
         END DO
      END DO

      ! ----------------
      !  interior value   (2=<jk=<jpkm1)
      ! ----------------
      DO jk = 2, jpkm1
         DO jj = 2, jpjm1 
            DO ji = fs_2, fs_jpim1   ! vector opt.
               ! hydrostatic pressure gradient along s-surfaces
               zhpi(ji,jj,jk) = zhpi(ji,jj,jk-1)                                &
                  &           + (  ( rho_k(ji+1,jj,jk) - rho_k(ji,jj,jk  ) )    &
                  &              - ( rho_i(ji  ,jj,jk) - rho_i(ji,jj,jk-1) )  ) / e1u(ji,jj)
               zhpj(ji,jj,jk) = zhpj(ji,jj,jk-1)                                &
                  &           + (  ( rho_k(ji,jj+1,jk) - rho_k(ji,jj,jk  ) )    &
                  &               -( rho_j(ji,jj  ,jk) - rho_j(ji,jj,jk-1) )  ) / e2v(ji,jj)
               ! add to the general momentum trend
               ua(ji,jj,jk) = ua(ji,jj,jk) + zhpi(ji,jj,jk)
               va(ji,jj,jk) = va(ji,jj,jk) + zhpj(ji,jj,jk)
            END DO
         END DO
      END DO
      !
      IF( wrk_not_released(3, 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15) )   &
         CALL ctl_stop('dyn:hpg_djc: failed to release workspace arrays')
      !
   END SUBROUTINE hpg_djc


   SUBROUTINE hpg_rot( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE hpg_rot  ***
      !!
      !! ** Method  :   rotated axes scheme (Thiem and Berntsen 2005)
      !!
      !! Reference: Thiem & Berntsen, Ocean Modelling, In press, 2005.
      !!----------------------------------------------------------------------
      USE wrk_nemo, ONLY:   wrk_in_use, wrk_not_released
      USE oce     , ONLY:   zhpi    => ta       , zhpj    => sa       ! (ta,sa) used as 3D workspace
      USE wrk_nemo, ONLY:   zdistr  => wrk_2d_1 , zsina   => wrk_2d_2 , zcosa  => wrk_2d_3
      USE wrk_nemo, ONLY:   zhpiorg => wrk_3d_1 , zhpirot => wrk_3d_2
      USE wrk_nemo, ONLY:   zhpitra => wrk_3d_3 , zhpine  => wrk_3d_4
      USE wrk_nemo, ONLY:   zhpjorg => wrk_3d_5 , zhpjrot => wrk_3d_6
      USE wrk_nemo, ONLY:   zhpjtra => wrk_3d_7 , zhpjne  => wrk_3d_8
      !!
      INTEGER, INTENT(in) ::   kt    ! ocean time-step index
      !!
      INTEGER  ::   ji, jj, jk          ! dummy loop indices
      REAL(wp) ::   zforg, zcoef0, zuap, zmskd1, zmskd1m   ! temporary scalar
      REAL(wp) ::   zfrot        , zvap, zmskd2, zmskd2m   !    "         "
      !!----------------------------------------------------------------------

      IF( wrk_in_use(2, 1,2,3)             .OR.   &
          wrk_in_use(3, 1,2,3,4,5,6,7,8) ) THEN
         CALL ctl_stop('dyn:hpg_rot: requested workspace arrays unavailable')   ;   RETURN
      ENDIF

      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn:hpg_rot : hydrostatic pressure gradient trend'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~   s-coordinate case, rotated axes scheme used'
      ENDIF

      ! -------------------------------
      !  Local constant initialization
      ! -------------------------------
      zcoef0 = - grav * 0.5_wp
      zforg  = 0.95_wp
      zfrot  = 1._wp - zforg

      ! inverse of the distance between 2 diagonal T-points (defined at F-point) (here zcoef0/distance)
      zdistr(:,:) = zcoef0 / SQRT( e1f(:,:)*e1f(:,:) + e2f(:,:)*e1f(:,:) )

      ! sinus and cosinus of diagonal angle at F-point
      zsina(:,:) = ATAN2( e2f(:,:), e1f(:,:) )
      zcosa(:,:) = COS( zsina(:,:) )
      zsina(:,:) = SIN( zsina(:,:) )

      ! ---------------
      !  Surface value
      ! ---------------
      ! compute and add to the general trend the pressure gradients along the axes
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1   ! vector opt.
            ! hydrostatic pressure gradient along s-surfaces
            zhpiorg(ji,jj,1) = zcoef0 / e1u(ji,jj) * (  fse3t(ji+1,jj,1) * rhd(ji+1,jj,1)   &
               &                                      - fse3t(ji  ,jj,1) * rhd(ji  ,jj,1)  )
            zhpjorg(ji,jj,1) = zcoef0 / e2v(ji,jj) * (  fse3t(ji,jj+1,1) * rhd(ji,jj+1,1)   &
               &                                      - fse3t(ji,jj  ,1) * rhd(ji,jj  ,1)  )
            ! s-coordinate pressure gradient correction
            zuap = -zcoef0 * ( rhd   (ji+1,jj  ,1) + rhd   (ji,jj,1) )   &
               &           * ( fsdept(ji+1,jj  ,1) - fsdept(ji,jj,1) ) / e1u(ji,jj)
            zvap = -zcoef0 * ( rhd   (ji  ,jj+1,1) + rhd   (ji,jj,1) )   &
               &           * ( fsdept(ji  ,jj+1,1) - fsdept(ji,jj,1) ) / e2v(ji,jj)
            ! add to the general momentum trend
            ua(ji,jj,1) = ua(ji,jj,1) + zforg * ( zhpiorg(ji,jj,1) + zuap )
            va(ji,jj,1) = va(ji,jj,1) + zforg * ( zhpjorg(ji,jj,1) + zvap )
         END DO
      END DO

      ! compute the pressure gradients in the diagonal directions
      DO jj = 1, jpjm1
         DO ji = 1, fs_jpim1   ! vector opt.
            zmskd1 = tmask(ji+1,jj+1,1) * tmask(ji  ,jj,1)      ! mask in the 1st diagnonal
            zmskd2 = tmask(ji  ,jj+1,1) * tmask(ji+1,jj,1)      ! mask in the 2nd diagnonal
            ! hydrostatic pressure gradient along s-surfaces
            zhpitra(ji,jj,1) = zdistr(ji,jj) * zmskd1 * (  fse3t(ji+1,jj+1,1) * rhd(ji+1,jj+1,1)   &
               &                                         - fse3t(ji  ,jj  ,1) * rhd(ji  ,jj  ,1)  )
            zhpjtra(ji,jj,1) = zdistr(ji,jj) * zmskd2 * (  fse3t(ji  ,jj+1,1) * rhd(ji  ,jj+1,1)   &
               &                                         - fse3t(ji+1,jj  ,1) * rhd(ji+1,jj  ,1)  )
            ! s-coordinate pressure gradient correction
            zuap = -zdistr(ji,jj) * zmskd1 * ( rhd   (ji+1,jj+1,1) + rhd   (ji  ,jj,1) )   &
               &                           * ( fsdept(ji+1,jj+1,1) - fsdept(ji  ,jj,1) )
            zvap = -zdistr(ji,jj) * zmskd2 * ( rhd   (ji  ,jj+1,1) + rhd   (ji+1,jj,1) )   &
               &                           * ( fsdept(ji  ,jj+1,1) - fsdept(ji+1,jj,1) )
            ! back rotation
            zhpine(ji,jj,1) = zcosa(ji,jj) * ( zhpitra(ji,jj,1) + zuap )   &
               &            - zsina(ji,jj) * ( zhpjtra(ji,jj,1) + zvap )
            zhpjne(ji,jj,1) = zsina(ji,jj) * ( zhpitra(ji,jj,1) + zuap )   &
               &            + zcosa(ji,jj) * ( zhpjtra(ji,jj,1) + zvap )
         END DO
      END DO

      ! interpolate and add to the general trend the diagonal gradient
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1   ! vector opt.
            ! averaging
            zhpirot(ji,jj,1) = 0.5 * ( zhpine(ji,jj,1) + zhpine(ji  ,jj-1,1) )
            zhpjrot(ji,jj,1) = 0.5 * ( zhpjne(ji,jj,1) + zhpjne(ji-1,jj  ,1) )
            ! add to the general momentum trend
            ua(ji,jj,1) = ua(ji,jj,1) + zfrot * zhpirot(ji,jj,1) 
            va(ji,jj,1) = va(ji,jj,1) + zfrot * zhpjrot(ji,jj,1) 
         END DO
      END DO

      ! -----------------
      ! 2. interior value (2=<jk=<jpkm1)
      ! -----------------
      ! compute and add to the general trend the pressure gradients along the axes
      DO jk = 2, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               ! hydrostatic pressure gradient along s-surfaces
               zhpiorg(ji,jj,jk) = zhpiorg(ji,jj,jk-1)                                                 &
                  &              +  zcoef0 / e1u(ji,jj) * (  fse3t(ji+1,jj,jk  ) * rhd(ji+1,jj,jk  )   &
                  &                                        - fse3t(ji  ,jj,jk  ) * rhd(ji  ,jj,jk  )   &
                  &                                        + fse3t(ji+1,jj,jk-1) * rhd(ji+1,jj,jk-1)   &
                  &                                        - fse3t(ji  ,jj,jk-1) * rhd(ji  ,jj,jk-1)  )
               zhpjorg(ji,jj,jk) = zhpjorg(ji,jj,jk-1)                                                 &
                  &              +  zcoef0 / e2v(ji,jj) * (  fse3t(ji,jj+1,jk  ) * rhd(ji,jj+1,jk  )   &
                  &                                        - fse3t(ji,jj  ,jk  ) * rhd(ji,jj,  jk  )   &
                  &                                        + fse3t(ji,jj+1,jk-1) * rhd(ji,jj+1,jk-1)   &
                  &                                        - fse3t(ji,jj  ,jk-1) * rhd(ji,jj,  jk-1)  )
               ! s-coordinate pressure gradient correction
               zuap = - zcoef0 * ( rhd   (ji+1,jj  ,jk) + rhd   (ji,jj,jk) )   &
                  &            * ( fsdept(ji+1,jj  ,jk) - fsdept(ji,jj,jk) ) / e1u(ji,jj)
               zvap = - zcoef0 * ( rhd   (ji  ,jj+1,jk) + rhd   (ji,jj,jk) )   &
                  &            * ( fsdept(ji  ,jj+1,jk) - fsdept(ji,jj,jk) ) / e2v(ji,jj)
               ! add to the general momentum trend
               ua(ji,jj,jk) = ua(ji,jj,jk) + zforg*( zhpiorg(ji,jj,jk) + zuap )
               va(ji,jj,jk) = va(ji,jj,jk) + zforg*( zhpjorg(ji,jj,jk) + zvap )
            END DO
         END DO
      END DO

      ! compute the pressure gradients in the diagonal directions
      DO jk = 2, jpkm1
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1   ! vector opt.
               zmskd1  = tmask(ji+1,jj+1,jk  ) * tmask(ji  ,jj,jk  )      ! level jk   mask in the 1st diagnonal
               zmskd1m = tmask(ji+1,jj+1,jk-1) * tmask(ji  ,jj,jk-1)      ! level jk-1    "               "     
               zmskd2  = tmask(ji  ,jj+1,jk  ) * tmask(ji+1,jj,jk  )      ! level jk   mask in the 2nd diagnonal
               zmskd2m = tmask(ji  ,jj+1,jk-1) * tmask(ji+1,jj,jk-1)      ! level jk-1    "               "     
               ! hydrostatic pressure gradient along s-surfaces
               zhpitra(ji,jj,jk) = zhpitra(ji,jj,jk-1)                                                       &
                  &              + zdistr(ji,jj) * zmskd1  * ( fse3t(ji+1,jj+1,jk  ) * rhd(ji+1,jj+1,jk)     &
                  &                                           -fse3t(ji  ,jj  ,jk  ) * rhd(ji  ,jj  ,jk) )   &
                  &              + zdistr(ji,jj) * zmskd1m * ( fse3t(ji+1,jj+1,jk-1) * rhd(ji+1,jj+1,jk-1)   &
                  &                                           -fse3t(ji  ,jj  ,jk-1) * rhd(ji  ,jj  ,jk-1) )
               zhpjtra(ji,jj,jk) = zhpjtra(ji,jj,jk-1)                                                       &
                  &              + zdistr(ji,jj) * zmskd2  * ( fse3t(ji  ,jj+1,jk  ) * rhd(ji  ,jj+1,jk)     &
                  &                                           -fse3t(ji+1,jj  ,jk  ) * rhd(ji+1,jj,  jk) )   &
                  &              + zdistr(ji,jj) * zmskd2m * ( fse3t(ji  ,jj+1,jk-1) * rhd(ji  ,jj+1,jk-1)   &
                  &                                           -fse3t(ji+1,jj  ,jk-1) * rhd(ji+1,jj,  jk-1) )
               ! s-coordinate pressure gradient correction
               zuap = - zdistr(ji,jj) * zmskd1 * ( rhd   (ji+1,jj+1,jk) + rhd   (ji  ,jj,jk) )   &
                  &                            * ( fsdept(ji+1,jj+1,jk) - fsdept(ji  ,jj,jk) )
               zvap = - zdistr(ji,jj) * zmskd2 * ( rhd   (ji  ,jj+1,jk) + rhd   (ji+1,jj,jk) )   &
                  &                            * ( fsdept(ji  ,jj+1,jk) - fsdept(ji+1,jj,jk) )
               ! back rotation
               zhpine(ji,jj,jk) = zcosa(ji,jj) * ( zhpitra(ji,jj,jk) + zuap )   &
                  &             - zsina(ji,jj) * ( zhpjtra(ji,jj,jk) + zvap )
               zhpjne(ji,jj,jk) = zsina(ji,jj) * ( zhpitra(ji,jj,jk) + zuap )   &
                  &             + zcosa(ji,jj) * ( zhpjtra(ji,jj,jk) + zvap )
            END DO
         END DO
      END DO

      ! interpolate and add to the general trend
      DO jk = 2, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               ! averaging
               zhpirot(ji,jj,jk) = 0.5 * ( zhpine(ji,jj,jk) + zhpine(ji  ,jj-1,jk) )
               zhpjrot(ji,jj,jk) = 0.5 * ( zhpjne(ji,jj,jk) + zhpjne(ji-1,jj  ,jk) )
               ! add to the general momentum trend
               ua(ji,jj,jk) = ua(ji,jj,jk) + zfrot * zhpirot(ji,jj,jk) 
               va(ji,jj,jk) = va(ji,jj,jk) + zfrot * zhpjrot(ji,jj,jk) 
            END DO
         END DO
      END DO
      !
      IF( wrk_not_released(2, 1,2,3)           .OR.   &
          wrk_not_released(3, 1,2,3,4,5,6,7,8) )   CALL ctl_stop('dyn:hpg_rot: failed to release workspace arrays')
      !
   END SUBROUTINE hpg_rot

   !!======================================================================
END MODULE dynhpg
