#if defined CCSMCOUPLED
MODULE sbcrnf
   !!======================================================================
   !!                       ***  MODULE  sbcrnf  ***
   !! Ocean forcing:  river runoff
   !!=====================================================================
   !! History :  OPA  ! 2000-11  (R. Hordoir, E. Durand)  NetCDF FORMAT
   !!   NEMO     1.0  ! 2002-09  (G. Madec)  F90: Free form and module
   !!            3.0  ! 2006-07  (G. Madec)  Surface module
   !!            3.2  ! 2009-04  (B. Lemaire)  Introduce iom_put
   !!            3.3  ! 2010-10  (R. Furner, G. Madec) runoff distributed over ocean levels
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   sbc_rnf      : monthly runoffs read in a NetCDF file
   !!   sbc_rnf_init : runoffs initialisation
   !!   rnf_mouth    : set river mouth mask
   !!----------------------------------------------------------------------
   USE dom_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE sbc_oce         ! surface boundary condition variables
   USE closea          ! closed seas
   USE fldread         ! read input field at current time step
   USE restart         ! restart
   USE in_out_manager  ! I/O manager
   USE iom             ! I/O module
   USE lib_mpp         ! MPP library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   sbc_rnf       ! routine call in sbcmod module
   PUBLIC   sbc_rnf_div   ! routine called in sshwzv module
   PUBLIC   sbc_rnf_alloc ! routine call in sbcmod module
   PUBLIC   sbc_rnf_init  ! (PUBLIC for TAM)
   !                                                     !!* namsbc_rnf namelist *
   CHARACTER(len=100), PUBLIC ::   cn_dir       = './'    !: Root directory for location of ssr files
   LOGICAL           , PUBLIC ::   ln_rnf_depth = .false. !: depth       river runoffs attribute specified in a file
   LOGICAL           , PUBLIC ::   ln_rnf_tem   = .false. !: temperature river runoffs attribute specified in a file
   LOGICAL           , PUBLIC ::   ln_rnf_sal   = .false. !: salinity    river runoffs attribute specified in a file
   LOGICAL           , PUBLIC ::   ln_rnf_emp   = .false. !: runoffs into a file to be read or already into precipitation
   TYPE(FLD_N)       , PUBLIC ::   sn_rnf                 !: information about the runoff file to be read
   TYPE(FLD_N)       , PUBLIC ::   sn_cnf                 !: information about the runoff mouth file to be read
   TYPE(FLD_N)                ::   sn_s_rnf               !: information about the salinities of runoff file to be read
   TYPE(FLD_N)                ::   sn_t_rnf               !: information about the temperatures of runoff file to be read
   TYPE(FLD_N)                ::   sn_dep_rnf             !: information about the depth which river inflow affects
   LOGICAL           , PUBLIC ::   ln_rnf_mouth = .false. !: specific treatment in mouths vicinity
   REAL(wp)          , PUBLIC ::   rn_hrnf      = 0._wp   !: runoffs, depth over which enhanced vertical mixing is used
   REAL(wp)          , PUBLIC ::   rn_avt_rnf   = 0._wp   !: runoffs, value of the additional vertical mixing coef. [m2/s]
   REAL(wp)          , PUBLIC ::   rn_rfact     = 1._wp   !: multiplicative factor for runoff
   REAL(wp)          , PUBLIC ::   rn_rnf_bnd   = 0._wp   !: max allowed runoff. [kg/m2/s]

   INTEGER , PUBLIC  ::   nkrnf = 0         !: nb of levels over which Kz is increased at river mouths
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   rnfmsk              !: river mouth mask (hori.)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)     ::   rnfmsk_z            !: river mouth mask (vert.)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   h_rnf               !: depth of runoff in m
   INTEGER,  PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   nk_rnf              !: depth of runoff in model levels
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   rnf_tsc_b, rnf_tsc  !: before and now T & S runoff contents   [K.m/s & PSU.m/s]

   REAL(wp) ::   r1_rau0   ! = 1 / rau0

   TYPE(FLD), PUBLIC, ALLOCATABLE, DIMENSION(:) ::   sf_rnf       ! structure: river runoff (file information, fields read) (PUBLIC for TAM)
   TYPE(FLD), PUBLIC, ALLOCATABLE, DIMENSION(:) ::   sf_s_rnf     ! structure: river runoff salinity (file information, fields read)  (PUBLIC for TAM)
   TYPE(FLD), PUBLIC, ALLOCATABLE, DIMENSION(:) ::   sf_t_rnf     ! structure: river runoff temperature (file information, fields read) (PUBLIC for TAM)

   !! * Substitutions
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION sbc_rnf_alloc()
      !!----------------------------------------------------------------------
      !!                ***  ROUTINE sbc_rnf_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( rnfmsk(jpi,jpj)         , rnfmsk_z(jpk)          ,     &
         &      h_rnf (jpi,jpj)         , nk_rnf  (jpi,jpj)      ,     &
         &      rnf_tsc_b(jpi,jpj,jpts) , rnf_tsc (jpi,jpj,jpts) , STAT=sbc_rnf_alloc )
         !
      IF( lk_mpp            )   CALL mpp_sum ( sbc_rnf_alloc )
      IF( sbc_rnf_alloc > 0 )   CALL ctl_warn('sbc_rnf_alloc: allocation of arrays failed')
   END FUNCTION sbc_rnf_alloc

   SUBROUTINE sbc_rnf( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sbc_rnf CESM coupled case ***
      !!       
      !! ** Purpose :   Handle runoff provided by the CESM coupler in a way
      !!                similar to the forced case
      !!
      !! ** Method  :   Runoff is kept in rnf (i.e. separated from emp) in order
      !!                to take it into account in sbc_rnf_div and optionally:
      !!                - distribute it over more than one level
      !!                - increase diffusivity where rnf>0.0
      !!
      !! ** Action  :   handle runoff received from the coupler
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt          ! ocean time step
      !!
      INTEGER           ::   ji, jj, jk    ! dummy loop indices
      !!----------------------------------------------------------------------
      !
      IF( kt == nit000 )   CALL sbc_rnf_init                           ! Read namelist and allocate structures

      !                                            ! ---------------------------------------- !
      IF( kt /= nit000 ) THEN                      !          Swap of forcing fields          !
         !                                         ! ---------------------------------------- !
         rnf_b    (:,:  ) = rnf    (:,:  )               ! Swap the ocean forcing fields except at nit000
         rnf_tsc_b(:,:,:) = rnf_tsc(:,:,:)               ! where before fields are set at the end of the routine
         !
      ENDIF

      !
      IF( MOD( kt - 1, nn_fsbc ) == 0 ) THEN
         !                                                     ! set temperature & salinity content of runoffs
         rnf_tsc(:,:,jp_tem) = sst_m(:,:) * rnf(:,:) * rau0r   ! assume runoff temperature is the SST
         !                                                     ! use S=0 for runoffs (done one for all in the init)
         IF ( ln_rnf_depth ) THEN
            nk_rnf(:,:) = 1
            DO jj = 1, jpj  
               DO ji = 1, jpi  
                  IF( rnf(ji,jj) /= 0._wp ) THEN  
                     jk = 2  
                     DO WHILE ( jk /= mbkt(ji,jj) .AND. fsdepw(ji,jj,jk+1) < rn_hrnf ) ;  jk = jk + 1 ;  END DO  
                     nk_rnf(ji,jj) = jk  
                  ENDIF  
               END DO  
            END DO  
            h_rnf(:,:) = 0.0_wp
            DO jj = 1, jpj                                ! set the associated depth 
               DO ji = 1, jpi 
                  DO jk = 1, nk_rnf(ji,jj)                        
                     h_rnf(ji,jj) = h_rnf(ji,jj) + fse3t(ji,jj,jk)  
                  END DO
               END DO
            END DO
         END IF

         IF ( ln_rnf_mouth ) THEN
            rnfmsk(:,:) = 0.0_wp
            WHERE (rnf(:,:) /= 0.0_wp)
               rnfmsk(:,:) = 0.5_wp
            END WHERE
         ENDIF
      ENDIF
      !
      IF( kt == nit000 ) THEN                          !   set the forcing field at nit000 - 1    !
         !                                             ! ---------------------------------------- !
         IF( ln_rstart .AND.    &                               !* Restart: read in restart file
            & iom_varid( numror, 'rnf_b', ldstop = .FALSE. ) > 0 ) THEN
            IF(lwp) WRITE(numout,*) '          nit000-1 runoff forcing fields red in the restart file'
            CALL iom_get( numror, jpdom_autoglo, 'rnf_b', rnf_b )     ! before runoff
            CALL iom_get( numror, jpdom_autoglo, 'rnf_hc_b', rnf_tsc_b(:,:,jp_tem) )   ! before heat content of runoff
            IF (ln_rnf_sal) CALL iom_get( numror, jpdom_autoglo, 'rnf_sc_b', rnf_tsc_b(:,:,jp_sal) )   ! before salinity content of runoff
         ELSE                                                   !* no restart: set from nit000 values
            IF(lwp) WRITE(numout,*) '          nit000-1 runoff forcing fields set to nit000'
             rnf_b    (:,:  ) = rnf    (:,:  )
             rnf_tsc_b(:,:,:) = rnf_tsc(:,:,:)
         ENDIF
      ENDIF
      !                                                ! ---------------------------------------- !
      IF( lrst_oce ) THEN                              !      Write in the ocean restart file     !
         !                                             ! ---------------------------------------- !
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'sbcrnf : runoff forcing fields written in ocean restart file ',   &
            &                    'at it= ', kt,' date= ', ndastp
         IF(lwp) WRITE(numout,*) '~~~~'
         CALL iom_rstput( kt, nitrst, numrow, 'rnf_b' , rnf )
         CALL iom_rstput( kt, nitrst, numrow, 'rnf_hc_b', rnf_tsc(:,:,jp_tem) )
         IF (ln_rnf_sal) CALL iom_rstput( kt, nitrst, numrow, 'rnf_sc_b', rnf_tsc(:,:,jp_sal) )
      ENDIF
      !
   END SUBROUTINE sbc_rnf


   SUBROUTINE sbc_rnf_div( phdivn )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sbc_rnf  ***
      !!
      !! ** Purpose :   update the horizontal divergence with the runoff inflow
      !!
      !! ** Method  :
      !!                CAUTION : rnf is positive (inflow) decreasing the
      !!                          divergence and expressed in m/s
      !!
      !! ** Action  :   phdivn   decreased by the runoff inflow
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:,:), INTENT(inout) ::   phdivn   ! horizontal divergence
      !!
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   r1_rau0   ! local scalar
      REAL(wp) ::   zfact     ! local scalar
      !!----------------------------------------------------------------------
      !
      zfact = 0.5_wp
      !
      r1_rau0 = 1._wp / rau0
      IF( ln_rnf_depth ) THEN      !==   runoff distributed over several levels   ==!
         IF( lk_vvl ) THEN             ! variable volume case
            DO jj = 1, jpj                   ! update the depth over which runoffs are distributed
               DO ji = 1, jpi
                  h_rnf(ji,jj) = 0._wp
                  DO jk = 1, nk_rnf(ji,jj)                           ! recalculates h_rnf to be the depth in metres
                     h_rnf(ji,jj) = h_rnf(ji,jj) + fse3t(ji,jj,jk)   ! to the bottom of the relevant grid box
                  END DO
                  !                          ! apply the runoff input flow
                  DO jk = 1, nk_rnf(ji,jj)
                     phdivn(ji,jj,jk) = phdivn(ji,jj,jk) - ( rnf(ji,jj) + rnf_b(ji,jj) ) * zfact * r1_rau0 / h_rnf(ji,jj)
                  END DO
               END DO
            END DO
         ELSE                          ! constant volume case : just apply the runoff input flow
            DO jj = 1, jpj
               DO ji = 1, jpi
                  DO jk = 1, nk_rnf(ji,jj)
                     phdivn(ji,jj,jk) = phdivn(ji,jj,jk) - ( rnf(ji,jj) + rnf_b(ji,jj) ) * zfact * r1_rau0 / h_rnf(ji,jj)
                  END DO
               END DO
            END DO
         ENDIF
      ELSE                       !==   runoff put only at the surface   ==!
         IF( lk_vvl ) THEN              ! variable volume case
            h_rnf(:,:) = fse3t(:,:,1)   ! recalculate h_rnf to be depth of top box
         ENDIF
         phdivn(:,:,1) = phdivn(:,:,1) - ( rnf(:,:) + rnf_b(:,:) ) * zfact * r1_rau0 / fse3t(:,:,1)
      ENDIF
      !
   END SUBROUTINE sbc_rnf_div


   SUBROUTINE sbc_rnf_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sbc_rnf_init  ***
      !!
      !! ** Purpose :   Initialisation of the runoffs if (ln_rnf=T)
      !!                when runoff is provided by the coupler.
      !!                Zero salinity assumed.
      !!                Runoff temperature is the SST if ln_rnf_tem=T
      !!                zero otherwise.
      !!
      !! ** Method  : - read the runoff namsbc_rnf namelist
      !!
      !! ** Action  : - read parameters
      !!----------------------------------------------------------------------
      INTEGER           ::   ji, jj, jk    ! dummy loop indices
      !!
      NAMELIST/namsbc_rnf/ cn_dir, ln_rnf_emp, ln_rnf_depth, ln_rnf_tem, ln_rnf_sal,   &
         &                 sn_rnf, sn_cnf    , sn_s_rnf    , sn_t_rnf  , sn_dep_rnf,   &  
         &                 ln_rnf_mouth      , rn_hrnf     , rn_avt_rnf, rn_rfact,     &
         &                 rn_rnf_bnd  
      !!----------------------------------------------------------------------

      !                                   ! ============
      !                                   !   Namelist
      !                                   ! ============
      !
      REWIND ( numnam )                         ! Read Namelist namsbc_rnf
      READ   ( numnam, namsbc_rnf )

      !                                         ! Control print
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'sbc_rnf : runoff received from coupler (NCAR CESM framework)'
         WRITE(numout,*) '~~~~~~~ '
         WRITE(numout,*) '   Namelist namsbc_rnf'
         WRITE(numout,*) '      runoff in a file to be read                ln_rnf_emp   = ', ln_rnf_emp
         WRITE(numout,*) '      specific river mouths treatment            ln_rnf_mouth = ', ln_rnf_mouth
         WRITE(numout,*) '      river mouth additional Kz                  rn_avt_rnf   = ', rn_avt_rnf
         WRITE(numout,*) '      depth of river mouth additional mixing     rn_hrnf      = ', rn_hrnf
         WRITE(numout,*) '      multiplicative factor for runoff           rn_rfact     = ', rn_rfact    
         WRITE(numout,*) '      runoff temperature assumed to be the SST   ln_rnf_tem   = ', ln_rnf_tem
         WRITE(numout,*) '      runoff salinity                            ln_rnf_sal   = ', ln_rnf_sal
         WRITE(numout,*) '      max allowed runoff [kgm-2s-1]              rn_rnf_bnd   = ', rn_rnf_bnd
      ENDIF

      IF ( ln_rnf_emp ) THEN
         CALL ctl_stop( 'STOP', 'sbc_rnf_init : ln_rnf_emp=.T. is incompatible with CESM coupled framework' )
      ENDIF
      IF ( ln_rnf_sal ) THEN
         CALL ctl_stop( 'STOP', 'sbc_rnf_init : ln_rnf_sal=.T. is incompatible with CESM coupled framework' )
      ENDIF

      !                                         !==  allocate runoff arrays
      IF( sbc_rnf_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'sbc_rnf_alloc : unable to allocate arrays' )
      !
      nk_rnf(:,:) = 1  
      IF ( ln_rnf_depth ) THEN
         IF (lwp) THEN
            WRITE(numout,*) '      runoff vertically distributed over the first rn_hrnf=', rn_hrnf, ' m'
         END IF
         DO jj = 1, jpj  
            DO ji = 1, jpi  
               IF( rnf(ji,jj) /= 0._wp ) THEN  
                  jk = 2  
                  DO WHILE ( jk /= mbkt(ji,jj) .AND. fsdepw(ji,jj,jk+1) < rn_hrnf ) ;  jk = jk + 1 ;  END DO  
                  nk_rnf(ji,jj) = jk  
               ENDIF  
            END DO  
         END DO  
         DO jj = 1, jpj                                ! set the associated depth 
            DO ji = 1, jpi 
               h_rnf(ji,jj) = 0._wp
               DO jk = 1, nk_rnf(ji,jj)                        
                  h_rnf(ji,jj) = h_rnf(ji,jj) + fse3t(ji,jj,jk)  
               END DO
            END DO
         END DO
      ELSE
         h_rnf (:,:) = fse3t(:,:,1)
      ENDIF
      !
      rnf_tsc(:,:,:) = 0._wp                    ! runoffs temperature & salinty contents initilisation
      !
      !                                   ! ========================
      !                                   !   River mouth vicinity
      !                                   ! ========================
      !
      IF( ln_rnf_mouth ) THEN                   ! Specific treatment in vicinity of river mouths :
         !                                      !    - Increase Kz in surface layers ( rn_hrnf > 0 )
         !                                      !    - set to zero SSS damping (ln_ssr=T)
         !                                      !    - mixed upstream-centered (ln_traadv_cen2=T)
         !
         IF ( ln_rnf_depth )   CALL ctl_warn( 'sbc_rnf_init: increased mixing turned on but effects may already',   &
            &                                              'be spread through depth by ln_rnf_depth'               )
         !
         nkrnf = 0                                  ! Number of level over which Kz increase
         IF( rn_hrnf > 0._wp ) THEN
            nkrnf = 2
            DO WHILE( nkrnf /= jpkm1 .AND. gdepw_0(nkrnf+1) < rn_hrnf )   ;   nkrnf = nkrnf + 1   ;   END DO
            IF( ln_sco )   &
               CALL ctl_warn( 'sbc_rnf: number of levels over which Kz is increased is computed for zco...' )
         ENDIF
         !
         rnfmsk(:,:) = 0.0_wp 
         WHERE (rnf(:,:) /= 0.0_wp)
            rnfmsk(:,:) = 0.5_wp 
         END WHERE
         !
         rnfmsk_z(:) = 0._wp
         rnfmsk_z(1:nkrnf) = e3t_0(1:nkrnf)/SUM(e3t_0(1:nkrnf))
         !
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '          Specific treatment used in vicinity of river mouths :'
         IF(lwp) WRITE(numout,*) '             - Increase Kz in surface layers (if rn_hrnf > 0 )'
         IF(lwp) WRITE(numout,*) '               by ', rn_avt_rnf,' m2/s  over ', nkrnf, ' w-levels'
         IF(lwp) WRITE(numout,*) '             - set to zero SSS damping       (if ln_ssr=T)'
         IF(lwp) WRITE(numout,*) '             - mixed upstream-centered       (if ln_traadv_cen2=T)'
         !
      ELSE                                      ! No treatment at river mouths
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '          No specific treatment at river mouths'
         rnfmsk  (:,:) = 0._wp
         rnfmsk_z(:)   = 0._wp
         nkrnf = 0
      ENDIF

   END SUBROUTINE sbc_rnf_init

   !!======================================================================
END MODULE sbcrnf
#endif
