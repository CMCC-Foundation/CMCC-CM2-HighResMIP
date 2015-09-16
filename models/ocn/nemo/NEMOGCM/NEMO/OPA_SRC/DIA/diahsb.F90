MODULE diahsb
   !!======================================================================
   !!                       ***  MODULE  diahsb  ***
   !! Ocean diagnostics: Heat, salt and volume budgets
   !!======================================================================
   !! History :  3.3  ! 2010-09  (M. Leclair)  Original code 
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE sbc_oce         ! surface thermohaline fluxes
   USE in_out_manager  ! I/O manager
   USE domvvl          ! vertical scale factors
   USE traqsr          ! penetrative solar radiation
   USE trabbc          ! bottom boundary condition 
   USE lib_mpp         ! distributed memory computing library
   USE trabbc          ! bottom boundary condition
   USE obc_par         ! (for lk_obc)
   USE bdy_par         ! (for lk_bdy)
   USE timing          ! preformance summary
   USE iom             ! I/O manager
   USE lib_fortran
   USE restart         ! ocean restart
   USE sbcrnf

#if defined CCSMCOUPLED
   USE qflxice         ! frazil ice formation
#endif

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dia_hsb        ! routine called by step.F90
   PUBLIC   dia_hsb_init   ! routine called by nemogcm.F90
   PUBLIC   dia_hsb_rst    ! routine called by step.F90

   LOGICAL, PUBLIC ::   ln_diahsb   !: check the heat and salt budgets

   INTEGER                                 ::   numhsb                           !
   REAL(dp)                                ::   surf_tot   , vol_tot             !
   REAL(dp), SAVE                          ::   frc_t      , frc_s     , frc_v   ! global forcing trends
   REAL(dp), SAVE                          ::   frc_wn_t      , frc_wn_s ! global forcing trends
   REAL(dp)                                ::   fact1                            ! conversion factors
   REAL(dp)                                ::   fact21    , fact22               !     -         -
   REAL(dp)                                ::   fact31    , fact32               !     -         -
   REAL(dp), DIMENSION(:,:)  , ALLOCATABLE ::   surf      , ssh_ini              !
   REAL(dp), DIMENSION(:,:,:), ALLOCATABLE ::   hc_loc_ini, sc_loc_ini, e3t_ini  !
   REAL(dp), DIMENSION(:,:)  , ALLOCATABLE ::   ssh_hc_loc_ini, ssh_sc_loc_ini

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE dia_hsb( kt )
      !!---------------------------------------------------------------------------
      !!                  ***  ROUTINE dia_hsb  ***
      !!     
      !! ** Purpose: Compute the ocean global heat content, salt content and volume conservation
      !!	
      !! ** Method : - Compute the deviation of heat content, salt content and volume
      !!	            at the current time step from their values at nit000
      !!	            - Compute the contribution of forcing and remove it from these deviations
      !!
      !! ** Action : Write the results in the 'heat_salt_volume_budgets.txt' ASCII file
      !!---------------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      !!
      INTEGER    ::   jk                          ! dummy loop indice
      REAL(dp)   ::   zdiff_hc    , zdiff_sc      ! heat and salt content variations
      REAL(dp)   ::   zdiff_hc1   , zdiff_sc1     ! heat and salt content variations of ssh
      REAL(dp)   ::   zdiff_v1    , zdiff_v2      ! volume variation
      REAL(dp)   ::   zerr_hc1    , zerr_sc1      ! Non conservation due to free surface
      REAL(dp)   ::   z1_rau0                     ! local scalars
      REAL(dp)   ::   zdeltat                     !    -     -
      REAL(dp)   ::   z_frc_trd_t , z_frc_trd_s   !    -     -
      REAL(dp)   ::   z_frc_trd_v                 !    -     -
      REAL(dp)   ::   z_wn_trd_t , z_wn_trd_s   !    -     -
      REAL(dp)   ::   z_ssh_hc , z_ssh_sc   !    -     -
      !!---------------------------------------------------------------------------
      IF( nn_timing == 1 )   CALL timing_start('dia_hsb')

      ! ------------------------- !
      ! 1 - Trends due to forcing !
      ! ------------------------- !
      z1_rau0 = 1.e0 / rau0
      z_frc_trd_v = - rau0r * glob_sum( ( emp(:,:) - rnf(:,:) ) * surf(:,:) )     ! volume fluxes
      z_frc_trd_t =           glob_sum( sbc_tsc(:,:,jp_tem) * surf(:,:) )     ! heat fluxes
      z_frc_trd_s =           glob_sum( sbc_tsc(:,:,jp_sal) * surf(:,:) )     ! salt fluxes
      ! Add runoff heat & salt input
#if defined CCSMCOUPLED
      IF( ln_rnf_cpl)   z_frc_trd_t = z_frc_trd_t + glob_sum( rnf_tsc(:,:,jp_tem) * surf(:,:) )
#else
      IF( ln_rnf    )   z_frc_trd_t = z_frc_trd_t + glob_sum( rnf_tsc(:,:,jp_tem) * surf(:,:) )
#endif
      IF( ln_rnf_sal)   z_frc_trd_s = z_frc_trd_s + glob_sum( rnf_tsc(:,:,jp_sal) * surf(:,:) )
      ! Add penetrative solar radiation
      IF( ln_traqsr )   z_frc_trd_t = z_frc_trd_t + ro0cpr * glob_sum( qsr     (:,:) * surf(:,:) )
      ! Add geothermal heat flux
      IF( ln_trabbc )   z_frc_trd_t = z_frc_trd_t +  glob_sum( qgh_trd0(:,:) * surf(:,:) )
      IF( .NOT. lk_vvl ) THEN
         z_wn_trd_t = - glob_sum( surf(:,:) * wn(:,:,1) * tsb(:,:,1,jp_tem) )
         z_wn_trd_s = - glob_sum( surf(:,:) * wn(:,:,1) * tsb(:,:,1,jp_sal) )
      ENDIF

      frc_v = frc_v + z_frc_trd_v * rdt
      frc_t = frc_t + z_frc_trd_t * rdt
      frc_s = frc_s + z_frc_trd_s * rdt

#if defined CCSMCOUPLED
      ! Add frazil ice formation heat flux
      z_frc_trd_t = -0.5_wp * glob_sum( QICE(:,:) * surf(:,:) )
      frc_t = frc_t + z_frc_trd_t 
#endif

      !                                          ! Advection flux through fixed surface (z=0)
      IF( .NOT. lk_vvl ) THEN
         frc_wn_t = frc_wn_t + z_wn_trd_t * rdt
         frc_wn_s = frc_wn_s + z_wn_trd_s * rdt
      ENDIF

      ! ----------------------- !
      ! 2 -  Content variations !
      ! ----------------------- !
      zdiff_v2 = 0.d0
      zdiff_hc = 0.d0
      zdiff_sc = 0.d0

      ! volume variation (calculated with ssh)
      zdiff_v1 = glob_sum( surf(:,:) * ( sshn(:,:) - ssh_ini(:,:) ) )

      ! heat & salt content variation (associated with ssh)
      IF( .NOT. lk_vvl ) THEN
         z_ssh_hc = glob_sum( surf(:,:) * ( tsn(:,:,1,jp_tem) * sshn(:,:) - ssh_hc_loc_ini(:,:) ) )
         z_ssh_sc = glob_sum( surf(:,:) * ( tsn(:,:,1,jp_sal) * sshn(:,:) - ssh_sc_loc_ini(:,:) ) )
      ENDIF

      DO jk = 1, jpkm1
        ! volume variation (calculated with scale factors)
         zdiff_v2 = zdiff_v2 + glob_sum( surf(:,:) * tmask(:,:,jk)   &
            &                       * ( fse3t_n(:,:,jk)         &
            &                           - e3t_ini(:,:,jk) ) )
         ! heat content variation
         zdiff_hc = zdiff_hc + glob_sum( surf(:,:) * tmask(:,:,jk)          &
            &                       * ( fse3t_n(:,:,jk) * tsn(:,:,jk,jp_tem)   &
            &                           - hc_loc_ini(:,:,jk) ) )
         ! salt content variation
         zdiff_sc = zdiff_sc + glob_sum( surf(:,:) * tmask(:,:,jk)          &
            &                       * ( fse3t_n(:,:,jk) * tsn(:,:,jk,jp_sal)   &
            &                           - sc_loc_ini(:,:,jk) ) )
      ENDDO

      ! Substract forcing from heat content, salt content and volume variations
      zdiff_v1 = zdiff_v1 - frc_v
      IF( lk_vvl )   zdiff_v2 = zdiff_v2 - frc_v
      zdiff_hc = zdiff_hc - frc_t
      zdiff_sc = zdiff_sc - frc_s
      IF( .NOT. lk_vvl ) THEN
         zdiff_hc1 = zdiff_hc + z_ssh_hc 
         zdiff_sc1 = zdiff_sc + z_ssh_sc
         zerr_hc1  = z_ssh_hc - frc_wn_t
         zerr_sc1  = z_ssh_sc - frc_wn_s
      ENDIF
      
      ! ----------------------- !
      ! 3 - Diagnostics writing !
      ! ----------------------- !
      zdeltat  = 1.e0 / ( ( kt - nit000 + 1 ) * rdt )
      IF( lk_vvl ) THEN
         vol_tot   = 0.0_wp                                                   ! total ocean volume
         DO jk = 1, jpkm1
            vol_tot  = vol_tot + sum( surf(:,:) * tmask(:,:,jk) * fse3t_n(:,:,jk) )
         END DO
         IF (lk_mpp) THEN
            CALL mpp_sum( vol_tot )
         ENDIF
      ENDIF

      IF( lk_vvl ) THEN
        CALL iom_put( 'bgtemper' , zdiff_hc / vol_tot )              ! Temperature variation (C) 
        CALL iom_put( 'bgsaline' , zdiff_sc / vol_tot )              ! Salinity    variation (psu)
        CALL iom_put( 'bgheatco' , zdiff_hc * 1.e-20_wp * rau0 * rcp )   ! Heat content variation (1.e20 J) 
        CALL iom_put( 'bgsaltco' , zdiff_sc * 1.e-9_wp    )              ! Salt content variation (psu*km3)
        CALL iom_put( 'bgvolssh' , zdiff_v1 * 1.e-9_wp    )              ! volume ssh variation (km3)  
        CALL iom_put( 'bgvole3t' , zdiff_v2 * 1.e-9_wp    )              ! volume e3t variation (km3)  
        CALL iom_put( 'bgfrcvol' , frc_v    * 1.e-9_wp    )              ! vol - surface forcing (km3) 
        CALL iom_put( 'bgfrctem' , frc_t / vol_tot    )              ! hc  - surface forcing (C) 
        CALL iom_put( 'bgfrcsal' , frc_s / vol_tot    )              ! sc  - surface forcing (psu) 
      ELSE
        CALL iom_put( 'bgtemper' , zdiff_hc1 / vol_tot)              ! Heat content variation (C) 
        CALL iom_put( 'bgsaline' , zdiff_sc1 / vol_tot)              ! Salt content variation (psu)
        CALL iom_put( 'bgheatco' , zdiff_hc1 * 1.e-20_wp * rau0 * rcp )  ! Heat content variation (1.e20 J) 
        CALL iom_put( 'bgsaltco' , zdiff_sc1 * 1.e-9_wp    )             ! Salt content variation (psu*km3)
        CALL iom_put( 'bgvolssh' , zdiff_v1 * 1.e-9_wp    )              ! volume ssh variation (km3)  
        CALL iom_put( 'bgfrcvol' , frc_v    * 1.e-9_wp    )              ! vol - surface forcing (km3) 
        CALL iom_put( 'bgfrctem' , frc_t / vol_tot    )              ! hc  - surface forcing (C) 
        CALL iom_put( 'bgfrcsal' , frc_s / vol_tot    )              ! sc  - surface forcing (psu) 
        CALL iom_put( 'bgmistem' , zerr_hc1 / vol_tot )              ! hc  - error due to free surface (C)
        CALL iom_put( 'bgmissal' , zerr_sc1 / vol_tot )              ! sc  - error due to free surface (psu)
      ENDIF
      !
      IF( lrst_oce )   CALL dia_hsb_rst( kt, 'WRITE' )

      IF( nn_timing == 1 )   CALL timing_stop('dia_hsb')

      !
   END SUBROUTINE dia_hsb


   SUBROUTINE dia_hsb_init
      !!---------------------------------------------------------------------------
      !!                  ***  ROUTINE dia_hsb  ***
      !!     
      !! ** Purpose: Initialization for the heat salt volume budgets
      !!	
      !! ** Method : Compute initial heat content, salt content and volume
      !!
      !! ** Action : - Compute initial heat content, salt content and volume
      !!             - Initialize forcing trends
      !!             - Compute coefficients for conversion
      !!---------------------------------------------------------------------------
      CHARACTER (len=32) ::   cl_name  ! output file name
      INTEGER            ::   jk       ! dummy loop indice
      INTEGER            ::   ierror   ! local integer
      !!
      NAMELIST/namhsb/ ln_diahsb
      !!----------------------------------------------------------------------
      !
      REWIND ( numnam )              ! Read Namelist namhsb 
      READ   ( numnam, namhsb )
      !
      IF(lwp) THEN                   ! Control print
         WRITE(numout,*)
         WRITE(numout,*) 'dia_hsb_init : check the heat and salt budgets'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namhsb : set hsb parameters'
         WRITE(numout,*) '      Switch for hsb diagnostic (T) or not (F)  ln_diahsb  = ', ln_diahsb
         WRITE(numout,*)
      ENDIF

      IF( .NOT. ln_diahsb )   RETURN
!      IF( .NOT. lk_mpp_rep ) &
!        CALL ctl_stop (' Your global mpp_sum if performed in single precision - 64 bits -', &
!             &         ' whereas the global sum to be precise must be done in double precision ',&
!             &         ' please add key_mpp_rep')

      ! ------------------- !
      ! 1 - Allocate memory !
      ! ------------------- !
      ALLOCATE( hc_loc_ini(jpi,jpj,jpk), sc_loc_ini(jpi,jpj,jpk), &
         &      ssh_hc_loc_ini(jpi,jpj), ssh_sc_loc_ini(jpi,jpj), &
         &      e3t_ini(jpi,jpj,jpk)                            , &
         &      surf(jpi,jpj),  ssh_ini(jpi,jpj), STAT=ierror )
      IF( ierror > 0 ) THEN
         CALL ctl_stop( 'dia_hsb: unable to allocate hc_loc_ini' )   ;   RETURN
      ENDIF

      ! ----------------------------------------------- !
      ! 2 - Time independant variables and file opening !
      ! ----------------------------------------------- !
      IF(lwp) WRITE(numout,*) "dia_hsb: heat salt volume budgets activated"
      IF(lwp) WRITE(numout,*) "~~~~~~~"
      IF( lk_obc .or. lk_bdy ) THEN
         CALL ctl_warn( 'dia_hsb does not take open boundary fluxes into account' )         
      ENDIF

      surf(:,:) = e1t(:,:) * e2t(:,:) * tmask(:,:,1) * tmask_i(:,:)      ! masked surface grid cell area
      surf_tot  = glob_sum( surf(:,:) )                                       ! total ocean surface area
      vol_tot   = 0.d0                                                   ! total ocean volume
      DO jk = 1, jpkm1
         vol_tot  = vol_tot + glob_sum( surf(:,:) * tmask(:,:,jk)     &
            &                         * fse3t_n(:,:,jk)         )
      END DO

      ! --------------- !
      ! 3 - Conversions ! (factors will be multiplied by duration afterwards)
      ! --------------- !

      ! heat content variation   =>   equivalent heat flux:
      fact1  = rau0 * rcp / surf_tot                                         ! [C*m3]   ->  [W/m2]
      ! salt content variation   =>   equivalent EMP and equivalent "flow": 
      fact21 = 1.e3  / ( soce * surf_tot )                                   ! [psu*m3] ->  [mm/s]
      fact22 = 1.e-6 / soce                                                  ! [psu*m3] ->  [Sv]
      ! volume variation         =>   equivalent EMP and equivalent "flow":
      fact31 = 1.e3  / surf_tot                                              ! [m3]     ->  [mm/s]
      fact32 = 1.e-6                                                         ! [m3]     ->  [SV]

      ! ---------------------------------- !
      ! 4 - initial conservation variables !
      ! ---------------------------------- !
      ssh_ini(:,:) = sshn(:,:)                                       ! initial ssh
      DO jk = 1, jpk
         e3t_ini   (:,:,jk) = fse3t_n(:,:,jk)                        ! initial vertical scale factors
         hc_loc_ini(:,:,jk) = tsn(:,:,jk,jp_tem) * fse3t_n(:,:,jk)   ! initial heat content
         sc_loc_ini(:,:,jk) = tsn(:,:,jk,jp_sal) * fse3t_n(:,:,jk)   ! initial salt content
      END DO
      frc_v = 0.d0                                           ! volume       trend due to forcing
      frc_t = 0.d0                                           ! heat content   -    -   -    -   
      frc_s = 0.d0                                           ! salt content   -    -   -    -         
      IF( .NOT. lk_vvl ) THEN
         ssh_hc_loc_ini(:,:) = tsn(:,:,1,jp_tem) * ssh_ini(:,:)   ! initial heat content associated with ssh
         ssh_sc_loc_ini(:,:) = tsn(:,:,1,jp_sal) * ssh_ini(:,:)   ! initial salt content associated with ssh
         frc_wn_t = 0.d0
         frc_wn_s = 0.d0
      ENDIF
      !
      CALL dia_hsb_rst( nit000, 'READ' )  !* read or initialize all required files
      !
   END SUBROUTINE dia_hsb_init

   SUBROUTINE dia_hsb_rst( kt, cdrw )
     !!---------------------------------------------------------------------
     !!                   ***  ROUTINE limdia_rst  ***
     !!                     
     !! ** Purpose :   Read or write DIA file in restart file
     !!
     !! ** Method  :   use of IOM library
     !!----------------------------------------------------------------------
     INTEGER         , INTENT(in) ::   kt     ! ocean time-step
     CHARACTER(len=*), INTENT(in) ::   cdrw   ! "READ"/"WRITE" flag
     !
     INTEGER ::   jk   ! 
     INTEGER ::   id1   ! local integers
     !!----------------------------------------------------------------------
     !
     IF( TRIM(cdrw) == 'READ' ) THEN        ! Read/initialise 
        IF( ln_rstart ) THEN                   !* Read the restart file
           !id1 = iom_varid( numror, 'frc_vol'  , ldstop = .FALSE. )
           !
           IF(lwp) WRITE(numout,*) '~~~~~~~'
           IF(lwp) WRITE(numout,*) ' dia_hsb_rst at it= ', kt,' date= ', ndastp
           IF(lwp) WRITE(numout,*) '~~~~~~~'
           CALL iom_get( numror, 'frc_v', frc_v )
           CALL iom_get( numror, 'frc_t', frc_t )
           CALL iom_get( numror, 'frc_s', frc_s )
           IF( .NOT. lk_vvl ) THEN
              CALL iom_get( numror, 'frc_wn_t', frc_wn_t )
              CALL iom_get( numror, 'frc_wn_s', frc_wn_s )
           ENDIF
           CALL iom_get( numror, jpdom_autoglo, 'ssh_ini', ssh_ini )
           CALL iom_get( numror, jpdom_autoglo, 'e3t_ini', e3t_ini )
           CALL iom_get( numror, jpdom_autoglo, 'hc_loc_ini', hc_loc_ini )
           CALL iom_get( numror, jpdom_autoglo, 'sc_loc_ini', sc_loc_ini )
           IF( .NOT. lk_vvl ) THEN
              CALL iom_get( numror, jpdom_autoglo, 'ssh_hc_loc_ini', ssh_hc_loc_ini )
              CALL iom_get( numror, jpdom_autoglo, 'ssh_sc_loc_ini', ssh_sc_loc_ini )
           ENDIF
       ELSE
          IF(lwp) WRITE(numout,*) '~~~~~~~'
          IF(lwp) WRITE(numout,*) ' dia_hsb at initial state '
          IF(lwp) WRITE(numout,*) '~~~~~~~'
          ssh_ini(:,:) = sshn(:,:)                                       ! initial ssh
          DO jk = 1, jpk
             e3t_ini   (:,:,jk) = fse3t_n(:,:,jk)                        ! initial vertical scale factors
             hc_loc_ini(:,:,jk) = tsn(:,:,jk,jp_tem) * fse3t_n(:,:,jk)   ! initial heat content
             sc_loc_ini(:,:,jk) = tsn(:,:,jk,jp_sal) * fse3t_n(:,:,jk)   ! initial salt content
          END DO
          frc_v = 0.0_wp                                           ! volume       trend due to forcing
          frc_t = 0.0_wp                                           ! heat content   -    -   -    -   
          frc_s = 0.0_wp                                           ! salt content   -    -   -    -        
          IF( .NOT. lk_vvl ) THEN
             ssh_hc_loc_ini(:,:) = tsn(:,:,1,jp_tem) * sshn(:,:)   ! initial heat content in ssh
             ssh_sc_loc_ini(:,:) = tsn(:,:,1,jp_sal) * sshn(:,:)   ! initial salt content in ssh
             frc_wn_t = 0.0_wp                                       ! initial heat content misfit due to free surface
             frc_wn_s = 0.0_wp                                       ! initial salt content misfit due to free surface
          ENDIF
       ENDIF

     ELSEIF( TRIM(cdrw) == 'WRITE' ) THEN   ! Create restart file
        !                                   ! -------------------
        IF(lwp) WRITE(numout,*) '~~~~~~~'
        IF(lwp) WRITE(numout,*) ' dia_hsb_rst at it= ', kt,' date= ', ndastp
        IF(lwp) WRITE(numout,*) '~~~~~~~'

        CALL iom_rstput( kt, nitrst, numrow, 'frc_v'   , frc_v     )
        CALL iom_rstput( kt, nitrst, numrow, 'frc_t'   , frc_t     )
        CALL iom_rstput( kt, nitrst, numrow, 'frc_s'   , frc_s     )
        IF( .NOT. lk_vvl ) THEN
           CALL iom_rstput( kt, nitrst, numrow, 'frc_wn_t', frc_wn_t )
           CALL iom_rstput( kt, nitrst, numrow, 'frc_wn_s', frc_wn_s )
        ENDIF
        CALL iom_rstput( kt, nitrst, numrow, 'ssh_ini', ssh_ini )
        CALL iom_rstput( kt, nitrst, numrow, 'e3t_ini', e3t_ini )
        CALL iom_rstput( kt, nitrst, numrow, 'hc_loc_ini', hc_loc_ini )
        CALL iom_rstput( kt, nitrst, numrow, 'sc_loc_ini', sc_loc_ini )
        IF( .NOT. lk_vvl ) THEN
           CALL iom_rstput( kt, nitrst, numrow, 'ssh_hc_loc_ini', ssh_hc_loc_ini )
           CALL iom_rstput( kt, nitrst, numrow, 'ssh_sc_loc_ini', ssh_sc_loc_ini )
        ENDIF
        !
     ENDIF
     !
   END SUBROUTINE dia_hsb_rst

   !!======================================================================
END MODULE diahsb
