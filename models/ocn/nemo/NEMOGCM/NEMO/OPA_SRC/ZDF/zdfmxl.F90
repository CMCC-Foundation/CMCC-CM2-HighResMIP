MODULE zdfmxl
   !!======================================================================
   !!                       ***  MODULE  zdfmxl  ***
   !! Ocean physics: mixed layer depth 
   !!======================================================================
   !! History :  1.0  ! 2003-08  (G. Madec)  original code
   !!            3.2  ! 2009-07  (S. Masson, G. Madec)  IOM + merge of DO-loop
   !!----------------------------------------------------------------------
   !!   zdf_mxl      : Compute the turbocline and mixed layer depths.
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain variables
   USE zdf_oce         ! ocean vertical physics
   USE in_out_manager  ! I/O manager
   USE prtctl          ! Print control
   USE iom             ! I/O library
   USE lib_mpp         ! MPP library
   USE trc_oce, ONLY : lk_offline ! offline flag

   IMPLICIT NONE
   PRIVATE

   PUBLIC   zdf_mxl       ! called by step.F90

   INTEGER , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   nmln    !: number of level in the mixed layer (used by TOP)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   hmld    !: mixing layer depth (turbocline)      [m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   hmlp    !: mixed layer depth  (rho=rho0+zdcrit) [m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   hmlpt   !: mixed layer depth at t-points        [m]

   !! * Substitutions
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 4.0 , NEMO Consortium (2011)
   !! $Id: zdfmxl.F90 2758 2011-05-02 14:04:25Z cetlod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION zdf_mxl_alloc()
      !!----------------------------------------------------------------------
      !!               ***  FUNCTION zdf_mxl_alloc  ***
      !!----------------------------------------------------------------------
      IF( .NOT. ALLOCATED( nmln ) ) THEN
         ALLOCATE( nmln(jpi,jpj), hmld(jpi,jpj), hmlp(jpi,jpj), hmlpt(jpi,jpj), STAT= zdf_mxl_alloc )
         !
         IF( lk_mpp             )   CALL mpp_sum ( zdf_mxl_alloc )
         IF( zdf_mxl_alloc /= 0 )   CALL ctl_warn('zdf_mxl_alloc: failed to allocate arrays.')
         !
      ENDIF
   END FUNCTION zdf_mxl_alloc


   SUBROUTINE zdf_mxl( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE zdfmxl  ***
      !!                   
      !! ** Purpose :   Compute the turbocline depth and the mixed layer depth
      !!              with density criteria.
      !!
      !! ** Method  :   The mixed layer depth is the shallowest W depth with 
      !!      the density of the corresponding T point (just bellow) bellow a
      !!      given value defined locally as rho(10m) + zrho_c
      !!               The turbocline depth is the depth at which the vertical
      !!      eddy diffusivity coefficient (resulting from the vertical physics
      !!      alone, not the isopycnal part, see trazdf.F) fall below a given
      !!      value defined locally (avt_c here taken equal to 5 cm/s2)
      !!
      !! ** Action  :   nmln, hmld, hmlp, hmlpt
      !!----------------------------------------------------------------------
      USE wrk_nemo, ONLY:   iwrk_in_use, iwrk_not_released
      USE wrk_nemo, ONLY:   imld => iwrk_2d_1    ! 2D integer workspace
      !!
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      !!
      INTEGER  ::   ji, jj, jk          ! dummy loop indices
      INTEGER  ::   iikn, iiki          ! temporary integer within a do loop
      REAL(wp) ::   zrho_c = 0.01_wp    ! density criterion for mixed layer depth
      REAL(wp) ::   zavt_c = 5.e-4_wp   ! Kz criterion for the turbocline depth
      !!----------------------------------------------------------------------

      IF( iwrk_in_use(2, 1) ) THEN
         CALL ctl_stop('zdf_mxl : requested workspace array unavailable')   ;   RETURN
      ENDIF

      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'zdf_mxl : mixed layer depth'
         IF(lwp) WRITE(numout,*) '~~~~~~~ '
         !                             ! allocate zdfmxl arrays
         IF( zdf_mxl_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'zdf_mxl : unable to allocate arrays' )
      ENDIF

      ! w-level of the mixing and mixed layers
      nmln(:,:) = mbkt(:,:) + 1        ! Initialization to the number of w ocean point
      imld(:,:) = mbkt(:,:) + 1
      DO jk = jpkm1, nlb10, -1         ! from the bottom to nlb10
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( rhop(ji,jj,jk) > rhop(ji,jj,nla10) + zrho_c )   nmln(ji,jj) = jk      ! Mixed layer
               IF( avt (ji,jj,jk) < zavt_c                     )   imld(ji,jj) = jk      ! Turbocline 
            END DO
         END DO
      END DO
      ! depth of the mixing and mixed layers
      DO jj = 1, jpj
         DO ji = 1, jpi
            iiki = imld(ji,jj)
            iikn = nmln(ji,jj)
            hmld (ji,jj) = fsdepw(ji,jj,iiki  ) * tmask(ji,jj,1)    ! Turbocline depth 
            hmlp (ji,jj) = fsdepw(ji,jj,iikn  ) * tmask(ji,jj,1)    ! Mixed layer depth
            hmlpt(ji,jj) = fsdept(ji,jj,iikn-1)                     ! depth of the last T-point inside the mixed layer
         END DO
      END DO
      IF( .NOT.lk_offline ) THEN            ! no need to output in offline mode
         CALL iom_put( "mldr10_1", hmlp )   ! mixed layer depth
         CALL iom_put( "mldkz5"  , hmld )   ! turbocline depth
      ENDIF
      
      IF(ln_ctl)   CALL prt_ctl( tab2d_1=REAL(nmln,wp), clinfo1=' nmln : ', tab2d_2=hmlp, clinfo2=' hmlp : ', ovlap=1 )
      !
      IF( iwrk_not_released(2, 1) )   CALL ctl_stop('zdf_mxl: failed to release workspace array')
      !
   END SUBROUTINE zdf_mxl

   !!======================================================================
END MODULE zdfmxl
