MODULE diahsb
   !!======================================================================
   !!                       ***  MODULE  diahsb  ***
   !! Ocean diagnostics: Heat, salt and volume budgets
   !!======================================================================
   !! History :  3.3  ! 2010-09  (M. Leclair)  Original code 
   !!                 ! 2012-10  (C. Rousset)  add iom_put
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
   USE bdy_par         ! (for lk_bdy)
   USE iom             ! I/O manager
   USE lib_fortran     ! glob_sum
   USE restart         ! ocean restart
   USE wrk_nemo        ! work arrays
   USE sbcrnf          ! river runoffd
#if defined CCSMCOUPLED
   USE qflxice         ! frazil ice formation
#endif

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dia_hsb        ! routine called by step.F90
   PUBLIC   dia_hsb_init   ! routine called by nemogcm.F90
   PUBLIC   dia_hsb_rst    ! routine called by step.F90

   LOGICAL, PUBLIC ::   ln_diahsb   !: check the heat and salt budgets

   REAL(dp)                                ::   surf_tot         !
   REAL(dp)                                ::   vol_tot          ! volume
   REAL(dp), SAVE                          ::   frc_t      , frc_s     , frc_v   ! global forcing trends
   REAL(dp), SAVE                          ::   frc_wn_t      , frc_wn_s ! global forcing trends
   REAL(dp), DIMENSION(:,:)  , ALLOCATABLE ::   surf      , ssh_ini              !
   REAL(dp), DIMENSION(:,:,:), ALLOCATABLE ::   hc_loc_ini, sc_loc_ini, e3t_ini  !
   REAL(dp), DIMENSION(:,:)  , ALLOCATABLE ::   ssh_hc_loc_ini, ssh_sc_loc_ini

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: diahsb.F90 4624 2014-04-28 12:09:03Z acc $
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
      !!---------------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      !!
      INTEGER    ::   jk                          ! dummy loop indice
      REAL(dp)   ::   zdiff_hc    , zdiff_sc      ! heat and salt content variations
      REAL(dp)   ::   zdiff_hc1   , zdiff_sc1     ! -   -   -   -   -   -   -   - 
      REAL(dp)   ::   zdiff_v1    , zdiff_v2      ! volume variation
      REAL(dp)   ::   zerr_hc1    , zerr_sc1       ! heat and salt content misfit
      REAL(dp)   ::   z_frc_trd_t , z_frc_trd_s   !    -     -
      REAL(dp)   ::   z_frc_trd_v                 !    -     -
      REAL(dp)   ::   z_wn_trd_t , z_wn_trd_s   !    -     -
      REAL(dp)   ::   z_ssh_hc , z_ssh_sc   !    -     -
      !!---------------------------------------------------------------------------

      ! ------------------------- !
      ! 1 - Trends due to forcing !
      ! ------------------------- !
      z_frc_trd_v = - rau0r * glob_sum( ( emp(:,:) - rnf(:,:) ) * surf(:,:) ) ! volume fluxes
      z_frc_trd_t =           glob_sum( sbc_tsc(:,:,jp_tem) * surf(:,:) )       ! heat fluxes
      z_frc_trd_s =           glob_sum( sbc_tsc(:,:,jp_sal) * surf(:,:) )       ! salt fluxes
      ! Add runoff heat & salt input
#if defined CCSMCOUPLED
      IF( ln_rnf_cpl)   z_frc_trd_t = z_frc_trd_t + glob_sum( rnf_tsc(:,:,jp_tem) * surf(:,:) )
#else
      IF( ln_rnf    )   z_frc_trd_t = z_frc_trd_t + glob_sum( rnf_tsc(:,:,jp_tem) * surf(:,:) )
#endif
      IF( ln_rnf_sal)   z_frc_trd_s = z_frc_trd_s + glob_sum( rnf_tsc(:,:,jp_sal) * surf(:,:) )

      ! Add penetrative solar radiation
      IF( ln_traqsr )   z_frc_trd_t = z_frc_trd_t + ro0cpr * glob_sum( qsr(:,:) * surf(:,:) )
      ! Add geothermal heat flux
      IF( ln_trabbc )   z_frc_trd_t = z_frc_trd_t +     glob_sum( qgh_trd0(:,:) * surf(:,:) )
      !
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

      ! ------------------------ !
      ! 2 -  Content variations !
      ! ------------------------ !
      zdiff_v2 = 0.0_wp
      zdiff_hc = 0.0_wp
      zdiff_sc = 0.0_wp

      ! volume variation (calculated with ssh)
      zdiff_v1 = glob_sum( surf(:,:) * ( sshn(:,:) - ssh_ini(:,:) ) )

      ! heat & salt content variation (associated with ssh)
      IF( .NOT. lk_vvl ) THEN
         z_ssh_hc = glob_sum( surf(:,:) * ( tsn(:,:,1,jp_tem) * sshn(:,:) - ssh_hc_loc_ini(:,:) ) )
         z_ssh_sc = glob_sum( surf(:,:) * ( tsn(:,:,1,jp_sal) * sshn(:,:) - ssh_sc_loc_ini(:,:) ) )
      ENDIF

      IF( lk_vvl ) THEN
         DO jk = 1, jpkm1
            ! volume variation (calculated with scale factors)
            zdiff_v2 = zdiff_v2 + sum( surf(:,:) * tmask(:,:,jk) &
               &                       * ( fse3t_n(:,:,jk) - e3t_ini(:,:,jk) ) )
         ENDDO
         IF (lk_mpp) THEN
            CALL mpp_sum( zdiff_v2 )
         ENDIF
      ENDIF

      DO jk = 1, jpkm1
         ! heat content variation
         zdiff_hc = zdiff_hc + sum(  surf(:,:) * tmask(:,:,jk) & 
            &                      * ( fse3t_n(:,:,jk) * tsn(:,:,jk,jp_tem) - hc_loc_ini(:,:,jk) ) )
         ! salt content variation
         zdiff_sc = zdiff_sc + sum(  surf(:,:) * tmask(:,:,jk)   &
            &                      * ( fse3t_n(:,:,jk) * tsn(:,:,jk,jp_sal) - sc_loc_ini(:,:,jk) ) )
      ENDDO

      IF (lk_mpp) THEN
         CALL mpp_sum( zdiff_hc )
         CALL mpp_sum( zdiff_sc )
      ENDIF

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
      INTEGER            ::   jk       ! dummy loop indice
      INTEGER            ::   ierror   ! local integer
      !!
      NAMELIST/namhsb/ ln_diahsb
      !
      !!----------------------------------------------------------------------

      REWIND( numnam)              ! Namelist namhsb in reference namelist
      READ  ( numnam, namhsb, IOSTAT = ierror, ERR = 901)
901   IF( ierror /= 0 ) CALL ctl_stop ( 'dia_hsb: namhsb in namelist' )

      IF(lwp) WRITE ( numout, namhsb )

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
         &      e3t_ini(jpi,jpj,jpk), surf(jpi,jpj),  ssh_ini(jpi,jpj), STAT=ierror )
      IF( ierror > 0 ) THEN
         CALL ctl_stop( 'dia_hsb: unable to allocate hc_loc_ini' )   ;   RETURN
      ENDIF

      IF(.NOT. lk_vvl ) ALLOCATE( ssh_hc_loc_ini(jpi,jpj), ssh_sc_loc_ini(jpi,jpj),STAT=ierror )
      IF( ierror > 0 ) THEN
         CALL ctl_stop( 'dia_hsb: unable to allocate hc_loc_ini' )   ;   RETURN
      ENDIF

      ! ----------------------------------------------- !
      ! 2 - Time independant variables and file opening !
      ! ----------------------------------------------- !
      IF(lwp) WRITE(numout,*) "dia_hsb: heat salt volume budgets activated"
      IF(lwp) WRITE(numout,*) '~~~~~~~'
      surf(:,:) = e1e2t(:,:) * tmask(:,:,1) * tmask_i(:,:)      ! masked surface grid cell area
      surf_tot  = glob_sum( surf(:,:) )                                       ! total ocean surface area

      vol_tot   = 0.0_wp                                                   ! total ocean volume
      DO jk = 1, jpkm1
         vol_tot  = vol_tot + sum( surf(:,:) * tmask(:,:,jk) * fse3t_n(:,:,jk) )
      END DO
      IF (lk_mpp) THEN
         CALL mpp_sum( vol_tot )
      ENDIF

      IF( lk_bdy ) CALL ctl_warn( 'dia_hsb does not take open boundary fluxes into account' )         
      !
      ! ---------------------------------- !
      ! 4 - initial conservation variables !
      ! ---------------------------------- !
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
