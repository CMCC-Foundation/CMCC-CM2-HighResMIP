MODULE istate
   !!======================================================================
   !!                     ***  MODULE  istate  ***
   !! Ocean state   :  initial state setting
   !!=====================================================================
   !! History :  OPA  !  1989-12  (P. Andrich)  Original code
   !!            5.0  !  1991-11  (G. Madec)  rewritting
   !!            6.0  !  1996-01  (G. Madec)  terrain following coordinates
   !!            8.0  !  2001-09  (M. Levy, M. Ben Jelloul)  istate_eel
   !!            8.0  !  2001-09  (M. Levy, M. Ben Jelloul)  istate_uvg
   !!   NEMO     1.0  !  2003-08  (G. Madec, C. Talandier)  F90: Free form, modules + EEL R5
   !!             -   !  2004-05  (A. Koch-Larrouy)  istate_gyre 
   !!            2.0  !  2006-07  (S. Masson)  distributed restart using iom
   !!            3.3  !  2010-10  (C. Ethe) merge TRC-TRA
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   istate_init   : initial state setting
   !!   istate_tem    : analytical profile for initial Temperature
   !!   istate_sal    : analytical profile for initial Salinity
   !!   istate_eel    : initial state setting of EEL R5 configuration
   !!   istate_gyre   : initial state setting of GYRE configuration
   !!   istate_uvg    : initial velocity in geostropic balance
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and active tracers 
   USE dom_oce         ! ocean space and time domain 
   USE daymod          ! calendar
   USE eosbn2          ! eq. of state, Brunt Vaisala frequency (eos     routine)
   USE ldftra_oce      ! ocean active tracers: lateral physics
   USE zdf_oce         ! ocean vertical physics
   USE phycst          ! physical constants
   USE dtatem          ! temperature data                 (dta_tem routine)
   USE dtasal          ! salinity data                    (dta_sal routine)
   USE restart         ! ocean restart                   (rst_read routine)
   USE in_out_manager  ! I/O manager
   USE iom             ! I/O library
   USE zpshde          ! partial step: hor. derivative (zps_hde routine)
   USE eosbn2          ! equation of state            (eos bn2 routine)
   USE domvvl          ! varying vertical mesh
   USE dynspg_oce      ! pressure gradient schemes
   USE dynspg_flt      ! pressure gradient schemes
   USE dynspg_exp      ! pressure gradient schemes
   USE dynspg_ts       ! pressure gradient schemes
   USE traswp          ! Swap arrays                      (tra_swp routine)
   USE lib_mpp         ! MPP library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   istate_init   ! routine called by step.F90

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: istate.F90 2777 2011-06-07 09:55:02Z smasson $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE istate_init
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE istate_init  ***
      !! 
      !! ** Purpose :   Initialization of the dynamics and tracer fields.
      !!----------------------------------------------------------------------
      ! - ML - needed for initialization of e3t_b
      INTEGER  ::  jk     ! dummy loop indice

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'istate_ini : Initialization of the dynamics and tracers'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~'

      rhd  (:,:,:) = 0.e0
      rhop (:,:,:) = 0.e0
      rn2  (:,:,:) = 0.e0 
      ta   (:,:,:) = 0.e0    
      sa   (:,:,:) = 0.e0

      IF( ln_rstart ) THEN                    ! Restart from a file
         !                                    ! -------------------
         neuler = 1                              ! Set time-step indicator at nit000 (leap-frog)
         CALL rst_read                           ! Read the restart file
         CALL tra_swap                           ! swap 3D arrays (t,s)  in a 4D array (ts)
         CALL day_init                           ! model calendar (using both namelist and restart infos)
      ELSE
         !                                    ! Start from rest
         !                                    ! ---------------
         numror = 0                              ! define numror = 0 -> no restart file to read
         neuler = 0                              ! Set time-step indicator at nit000 (euler forward)
         CALL day_init                           ! model calendar (using both namelist and restart infos)
         !                                       ! Initialization of ocean to zero
         !   before fields     !       now fields     
         sshb (:,:)   = 0.e0   ;   sshn (:,:)   = 0.e0
         ub   (:,:,:) = 0.e0   ;   un   (:,:,:) = 0.e0
         vb   (:,:,:) = 0.e0   ;   vn   (:,:,:) = 0.e0  
         rotb (:,:,:) = 0.e0   ;   rotn (:,:,:) = 0.e0
         hdivb(:,:,:) = 0.e0   ;   hdivn(:,:,:) = 0.e0
         !
         IF( cp_cfg == 'eel' ) THEN
            CALL istate_eel                      ! EEL   configuration : start from pre-defined U,V T-S fields
         ELSEIF( cp_cfg == 'gyre' ) THEN         
            CALL istate_gyre                     ! GYRE  configuration : start from pre-defined T-S fields
         ELSE
            !                                    ! Other configurations: Initial T-S fields
#if defined key_dtatem
            CALL dta_tem( nit000 )                  ! read 3D temperature data
            tb(:,:,:) = t_dta(:,:,:)   ;   tn(:,:,:) = t_dta(:,:,:)
            
#else
            IF(lwp) WRITE(numout,*)                 ! analytical temperature profile
            IF(lwp) WRITE(numout,*)'             Temperature initialization using an analytic profile'
            CALL istate_tem
#endif
#if defined key_dtasal
            CALL dta_sal( nit000 )                  ! read 3D salinity data
            sb(:,:,:) = s_dta(:,:,:)   ;   sn(:,:,:) = s_dta(:,:,:)
#else
            ! No salinity data
            IF(lwp)WRITE(numout,*)                  ! analytical salinity profile
            IF(lwp)WRITE(numout,*)'             Salinity initialisation using a constant value'
            CALL istate_sal
#endif
         ENDIF
         !
         CALL tra_swap                     ! swap 3D arrays (tb,sb,tn,sn)  in a 4D array
         CALL eos( tsb, rhd, rhop )        ! before potential and in situ densities
#if ! defined key_c1d
         IF( ln_zps )   CALL zps_hde( nit000, jpts, tsb, gtsu, gtsv,  & ! zps: before hor. gradient
            &                                       rhd, gru , grv  )   ! of t,s,rd at ocean bottom
#endif
         !   
         ! - ML - sshn could be modified by istate_eel, so that initialization of fse3t_b is done here
         IF( lk_vvl ) THEN
            DO jk = 1, jpk
               fse3t_b(:,:,jk) = fse3t_n(:,:,jk)
            ENDDO
         ENDIF
         ! 
      ENDIF
      !
      IF( lk_agrif ) THEN                  ! read free surface arrays in restart file
         IF( ln_rstart ) THEN
            IF( lk_dynspg_flt )   CALL flt_rst( nit000, 'READ' )      ! read or initialize the following fields
            !                                                         ! gcx, gcxb for agrif_opa_init
         ENDIF                                                        ! explicit case not coded yet with AGRIF
      ENDIF
      !
   END SUBROUTINE istate_init


   SUBROUTINE istate_tem
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE istate_tem  ***
      !!   
      !! ** Purpose :   Intialization of the temperature field with an 
      !!      analytical profile or a file (i.e. in EEL configuration)
      !!
      !! ** Method  :   Use Philander analytic profile of temperature
      !!
      !! References :  Philander ???
      !!----------------------------------------------------------------------
      INTEGER :: ji, jj, jk
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'istate_tem : initial temperature profile'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~'
      !
      DO jk = 1, jpk
         DO jj = 1, jpj
            DO ji = 1, jpi
               tn(ji,jj,jk) = (  ( ( 7.5 - 0.*ABS(gphit(ji,jj))/30. )   &
                  &               *( 1.-TANH((fsdept(ji,jj,jk)-80.)/30.) )   &
                  &            + 10.*(5000.-fsdept(ji,jj,jk))/5000.)  ) * tmask(ji,jj,jk)
               tb(ji,jj,jk) = tn(ji,jj,jk)
          END DO
        END DO
      END DO
      !
      IF(lwp) CALL prizre( tn    , jpi   , jpj   , jpk   , jpj/2 ,   &
         &                 1     , jpi   , 5     , 1     , jpk   ,   &
         &                 1     , 1.    , numout                  )
      !
   END SUBROUTINE istate_tem


   SUBROUTINE istate_sal
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE istate_sal  ***
      !!
      !! ** Purpose :   Intialize the salinity field with an analytic profile
      !!
      !! ** Method  :   Use to a constant value 35.5
      !!              
      !! ** Action  :   Initialize sn and sb
      !!----------------------------------------------------------------------
      REAL(wp) ::   zsal = 35.50_wp
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'istate_sal : initial salinity : ', zsal
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~'
      !
      sn(:,:,:) = zsal * tmask(:,:,:)
      sb(:,:,:) = sn(:,:,:)
      !
   END SUBROUTINE istate_sal


   SUBROUTINE istate_eel
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE istate_eel  ***
      !! 
      !! ** Purpose :   Initialization of the dynamics and tracers for EEL R5
      !!      configuration (channel with or without a topographic bump)
      !!
      !! ** Method  : - set temprature field
      !!              - set salinity field
      !!              - set velocity field including horizontal divergence
      !!                and relative vorticity fields
      !!----------------------------------------------------------------------
      USE divcur     ! hor. divergence & rel. vorticity      (div_cur routine)
      USE iom
 
      INTEGER  ::   inum              ! temporary logical unit
      INTEGER  ::   ji, jj, jk        ! dummy loop indices
      INTEGER  ::   ijloc
      REAL(wp) ::   zh1, zh2, zslope, zcst, zfcor   ! temporary scalars
      REAL(wp) ::   zt1  = 15._wp                   ! surface temperature value (EEL R5)
      REAL(wp) ::   zt2  =  5._wp                   ! bottom  temperature value (EEL R5)
      REAL(wp) ::   zsal = 35.0_wp                  ! constant salinity (EEL R2, R5 and R6)
      REAL(wp) ::   zueel = 0.1_wp                  ! constant uniform zonal velocity (EEL R5)
      REAL(wp), DIMENSION(jpiglo,jpjglo) ::   zssh  ! initial ssh over the global domain
      !!----------------------------------------------------------------------

      SELECT CASE ( jp_cfg ) 
         !                                              ! ====================
         CASE ( 5 )                                     ! EEL R5 configuration
            !                                           ! ====================
            !
            ! set temperature field with a linear profile
            ! -------------------------------------------
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'istate_eel : EEL R5: linear temperature profile'
            IF(lwp) WRITE(numout,*) '~~~~~~~~~~'
            !
            zh1 = gdept_0(  1  )
            zh2 = gdept_0(jpkm1)
            !
            zslope = ( zt1 - zt2 ) / ( zh1 - zh2 )
            zcst   = ( zt1 * ( zh1 - zh2) - ( zt1 - zt2 ) * zh1 ) / ( zh1 - zh2 )
            !
            DO jk = 1, jpk
               tn(:,:,jk) = ( zt2 + zt1 * exp( - fsdept(:,:,jk) / 1000 ) ) * tmask(:,:,jk)
               tb(:,:,jk) = tn(:,:,jk)
            END DO
            !
            IF(lwp) CALL prizre( tn    , jpi   , jpj   , jpk   , jpj/2 ,   &
               &                 1     , jpi   , 5     , 1     , jpk   ,   &
               &                 1     , 1.    , numout                  )
            !
            ! set salinity field to a constant value
            ! --------------------------------------
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'istate_eel : EEL R5: constant salinity field, S = ', zsal
            IF(lwp) WRITE(numout,*) '~~~~~~~~~~'
            !
            sn(:,:,:) = zsal * tmask(:,:,:)
            sb(:,:,:) = sn(:,:,:)
            !
            ! set the dynamics: U,V, hdiv, rot (and ssh if necessary)
            ! ----------------
            ! Start EEL5 configuration with barotropic geostrophic velocities 
            ! according the sshb and sshn SSH imposed.
            ! we assume a uniform grid (hence the use of e1t(1,1) for delta_y)
            ! we use the Coriolis frequency at mid-channel.   
            ub(:,:,:) = zueel * umask(:,:,:)
            un(:,:,:) = ub(:,:,:)
            ijloc = mj0(INT(jpjglo-1)/2)
            zfcor = ff(1,ijloc)
            !
            DO jj = 1, jpjglo
               zssh(:,jj) = - (FLOAT(jj)- FLOAT(jpjglo-1)/2.)*zueel*e1t(1,1)*zfcor/grav 
            END DO
            !
            IF(lwp) THEN
               WRITE(numout,*) ' Uniform zonal velocity for EEL R5:',zueel
               WRITE(numout,*) ' Geostrophic SSH profile as a function of y:'
               WRITE(numout,'(12(1x,f6.2))') zssh(1,:)
            ENDIF
            !
            DO jj = 1, nlcj
               DO ji = 1, nlci
                  sshb(ji,jj) = zssh( mig(ji) , mjg(jj) ) * tmask(ji,jj,1)
               END DO
            END DO
            sshb(nlci+1:jpi,      :   ) = 0.e0      ! set to zero extra mpp columns
            sshb(      :   ,nlcj+1:jpj) = 0.e0      ! set to zero extra mpp rows
            !
            sshn(:,:) = sshb(:,:)                   ! set now ssh to the before value
            !
            IF( nn_rstssh /= 0 ) THEN  
               nn_rstssh = 0                        ! hand-made initilization of ssh 
               CALL ctl_warn( 'istate_eel: force nn_rstssh = 0' )
            ENDIF
            !
            CALL div_cur( nit000 )                  ! horizontal divergence and relative vorticity (curl)
            ! N.B. the vertical velocity will be computed from the horizontal divergence field
            ! in istate by a call to wzv routine


            !                                     ! ==========================
         CASE ( 2 , 6 )                           ! EEL R2 or R6 configuration
            !                                     ! ==========================
            !
            ! set temperature field with a NetCDF file
            ! ----------------------------------------
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'istate_eel : EEL R2 or R6: read initial temperature in a NetCDF file'
            IF(lwp) WRITE(numout,*) '~~~~~~~~~~'
            !
            CALL iom_open ( 'eel.initemp', inum )
            CALL iom_get ( inum, jpdom_data, 'initemp', tb ) ! read before temprature (tb)
            CALL iom_close( inum )
            !
            tn(:,:,:) = tb(:,:,:)                            ! set nox temperature to tb
            !
            IF(lwp) CALL prizre( tn    , jpi   , jpj   , jpk   , jpj/2 ,   &
               &                 1     , jpi   , 5     , 1     , jpk   ,   &
               &                 1     , 1.    , numout                  )
            !
            ! set salinity field to a constant value
            ! --------------------------------------
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'istate_eel : EEL R5: constant salinity field, S = ', zsal
            IF(lwp) WRITE(numout,*) '~~~~~~~~~~'
            !
            sn(:,:,:) = zsal * tmask(:,:,:)
            sb(:,:,:) = sn(:,:,:)
            !
            !                                    ! ===========================
         CASE DEFAULT                            ! NONE existing configuration
            !                                    ! ===========================
            WRITE(ctmp1,*) 'EEL with a ', jp_cfg,' km resolution is not coded'
            CALL ctl_stop( ctmp1 )
            !
      END SELECT
      !
   END SUBROUTINE istate_eel


   SUBROUTINE istate_gyre
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE istate_gyre  ***
      !! 
      !! ** Purpose :   Initialization of the dynamics and tracers for GYRE
      !!      configuration (double gyre with rotated domain)
      !!
      !! ** Method  : - set temprature field
      !!              - set salinity field
      !!----------------------------------------------------------------------
      INTEGER :: ji, jj, jk  ! dummy loop indices
      INTEGER            ::   inum          ! temporary logical unit
      INTEGER, PARAMETER ::   ntsinit = 0   ! (0/1) (analytical/input data files) T&S initialization
      !!----------------------------------------------------------------------

      SELECT CASE ( ntsinit)

      CASE ( 0 )                  ! analytical T/S profil deduced from LEVITUS
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'istate_gyre : initial analytical T and S profil deduced from LEVITUS '
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'

         DO jk = 1, jpk
            DO jj = 1, jpj
               DO ji = 1, jpi
                  tn(ji,jj,jk) = (  16. - 12. * TANH( (fsdept(ji,jj,jk) - 400) / 700 )         )   &
                       &           * (-TANH( (500-fsdept(ji,jj,jk)) / 150 ) + 1) / 2               &
                       &       + (      15. * ( 1. - TANH( (fsdept(ji,jj,jk)-50.) / 1500.) )       &
                       &                - 1.4 * TANH((fsdept(ji,jj,jk)-100.) / 100.)               &    
                       &                + 7.  * (1500. - fsdept(ji,jj,jk)) / 1500.             )   & 
                       &           * (-TANH( (fsdept(ji,jj,jk) - 500) / 150) + 1) / 2
                  tn(ji,jj,jk) = tn(ji,jj,jk) * tmask(ji,jj,jk)
                  tb(ji,jj,jk) = tn(ji,jj,jk)

                  sn(ji,jj,jk) =  (  36.25 - 1.13 * TANH( (fsdept(ji,jj,jk) - 305) / 460 )  )  &
                     &              * (-TANH((500 - fsdept(ji,jj,jk)) / 150) + 1) / 2          &
                     &          + (  35.55 + 1.25 * (5000. - fsdept(ji,jj,jk)) / 5000.         &
                     &                - 1.62 * TANH( (fsdept(ji,jj,jk) - 60.  ) / 650. )       &
                     &                + 0.2  * TANH( (fsdept(ji,jj,jk) - 35.  ) / 100. )       &
                     &                + 0.2  * TANH( (fsdept(ji,jj,jk) - 1000.) / 5000.)    )  &
                     &              * (-TANH((fsdept(ji,jj,jk) - 500) / 150) + 1) / 2 
                  sn(ji,jj,jk) = sn(ji,jj,jk) * tmask(ji,jj,jk)
                  sb(ji,jj,jk) = sn(ji,jj,jk)
               END DO
            END DO
         END DO

      CASE ( 1 )                  ! T/S data fields read in dta_tem.nc/data_sal.nc files
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'istate_gyre : initial T and S read from dta_tem.nc/data_sal.nc files'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
         IF(lwp) WRITE(numout,*) '              NetCDF FORMAT'

         ! Read temperature field
         ! ----------------------
         CALL iom_open ( 'data_tem', inum )
         CALL iom_get ( inum, jpdom_data, 'votemper', tn ) 
         CALL iom_close( inum )

         tn(:,:,:) = tn(:,:,:) * tmask(:,:,:) 
         tb(:,:,:) = tn(:,:,:)

         ! Read salinity field
         ! -------------------
         CALL iom_open ( 'data_sal', inum )
         CALL iom_get ( inum, jpdom_data, 'vosaline', sn ) 
         CALL iom_close( inum )

         sn(:,:,:)  = sn(:,:,:) * tmask(:,:,:) 
         sb(:,:,:)  = sn(:,:,:)

      END SELECT

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) '              Initial temperature and salinity profiles:'
         WRITE(numout, "(9x,' level   gdept_0   temperature   salinity   ')" )
         WRITE(numout, "(10x, i4, 3f10.2)" ) ( jk, gdept_0(jk), tn(2,2,jk), sn(2,2,jk), jk = 1, jpk )
      ENDIF

   END SUBROUTINE istate_gyre


   SUBROUTINE istate_uvg
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE istate_uvg  ***
      !!
      !! ** Purpose :   Compute the geostrophic velocities from (tn,sn) fields
      !!
      !! ** Method  :   Using the hydrostatic hypothesis the now hydrostatic 
      !!      pressure is computed by integrating the in-situ density from the
      !!      surface to the bottom.
      !!                 p=integral [ rau*g dz ]
      !!----------------------------------------------------------------------
      USE wrk_nemo, ONLY:   wrk_in_use, wrk_not_released
      USE wrk_nemo, ONLY:   zprn => wrk_3d_1    ! 3D workspace

      USE dynspg          ! surface pressure gradient             (dyn_spg routine)
      USE divcur          ! hor. divergence & rel. vorticity      (div_cur routine)
      USE lbclnk          ! ocean lateral boundary condition (or mpp link)

      INTEGER ::   ji, jj, jk        ! dummy loop indices
      INTEGER ::   indic             ! ???
      REAL(wp) ::   zmsv, zphv, zmsu, zphu, zalfg     ! temporary scalars
      !!----------------------------------------------------------------------

      IF(wrk_in_use(3, 1) ) THEN
         CALL ctl_stop('istate_uvg: requested workspace array unavailable')   ;   RETURN
      ENDIF

      IF(lwp) WRITE(numout,*) 
      IF(lwp) WRITE(numout,*) 'istate_uvg : Start from Geostrophy'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~'

      ! Compute the now hydrostatic pressure
      ! ------------------------------------

      zalfg = 0.5 * grav * rau0
      
      zprn(:,:,1) = zalfg * fse3w(:,:,1) * ( 1 + rhd(:,:,1) )       ! Surface value

      DO jk = 2, jpkm1                                              ! Vertical integration from the surface
         zprn(:,:,jk) = zprn(:,:,jk-1)   &
            &         + zalfg * fse3w(:,:,jk) * ( 2. + rhd(:,:,jk) + rhd(:,:,jk-1) )
      END DO  

      ! Compute geostrophic balance
      ! ---------------------------
      DO jk = 1, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vertor opt.
               zmsv = 1. / MAX(  umask(ji-1,jj+1,jk) + umask(ji  ,jj+1,jk)   &
                               + umask(ji-1,jj  ,jk) + umask(ji  ,jj  ,jk) , 1.  )
               zphv = ( zprn(ji  ,jj+1,jk) - zprn(ji-1,jj+1,jk) ) * umask(ji-1,jj+1,jk) / e1u(ji-1,jj+1)   &
                    + ( zprn(ji+1,jj+1,jk) - zprn(ji  ,jj+1,jk) ) * umask(ji  ,jj+1,jk) / e1u(ji  ,jj+1)   &
                    + ( zprn(ji  ,jj  ,jk) - zprn(ji-1,jj  ,jk) ) * umask(ji-1,jj  ,jk) / e1u(ji-1,jj  )   &
                    + ( zprn(ji+1,jj  ,jk) - zprn(ji  ,jj  ,jk) ) * umask(ji  ,jj  ,jk) / e1u(ji  ,jj  )
               zphv = 1. / rau0 * zphv * zmsv * vmask(ji,jj,jk)

               zmsu = 1. / MAX(  vmask(ji+1,jj  ,jk) + vmask(ji  ,jj  ,jk)   &
                               + vmask(ji+1,jj-1,jk) + vmask(ji  ,jj-1,jk) , 1.  )
               zphu = ( zprn(ji+1,jj+1,jk) - zprn(ji+1,jj  ,jk) ) * vmask(ji+1,jj  ,jk) / e2v(ji+1,jj  )   &
                    + ( zprn(ji  ,jj+1,jk) - zprn(ji  ,jj  ,jk) ) * vmask(ji  ,jj  ,jk) / e2v(ji  ,jj  )   &
                    + ( zprn(ji+1,jj  ,jk) - zprn(ji+1,jj-1,jk) ) * vmask(ji+1,jj-1,jk) / e2v(ji+1,jj-1)   &
                    + ( zprn(ji  ,jj  ,jk) - zprn(ji  ,jj-1,jk) ) * vmask(ji  ,jj-1,jk) / e2v(ji  ,jj-1)
               zphu = 1. / rau0 * zphu * zmsu * umask(ji,jj,jk)

               ! Compute the geostrophic velocities
               un(ji,jj,jk) = -2. * zphu / ( ff(ji,jj) + ff(ji  ,jj-1) )
               vn(ji,jj,jk) =  2. * zphv / ( ff(ji,jj) + ff(ji-1,jj  ) )
            END DO
         END DO
      END DO

      IF(lwp) WRITE(numout,*) '         we force to zero bottom velocity'

      ! Susbtract the bottom velocity (level jpk-1 for flat bottom case)
      ! to have a zero bottom velocity

      DO jk = 1, jpkm1
         un(:,:,jk) = ( un(:,:,jk) - un(:,:,jpkm1) ) * umask(:,:,jk)
         vn(:,:,jk) = ( vn(:,:,jk) - vn(:,:,jpkm1) ) * vmask(:,:,jk)
      END DO

      CALL lbc_lnk( un, 'U', -1. )
      CALL lbc_lnk( vn, 'V', -1. )
      
      ub(:,:,:) = un(:,:,:)
      vb(:,:,:) = vn(:,:,:)
      
      ! WARNING !!!!!
      ! after initializing u and v, we need to calculate the initial streamfunction bsf.
      ! Otherwise, only the trend will be computed and the model will blow up (inconsistency).
      ! to do that, we call dyn_spg with a special trick:
      ! we fill ua and va with the velocities divided by dt, and the streamfunction will be brought to the
      ! right value assuming the velocities have been set up in one time step.
      ! we then set bsfd to zero (first guess for next step is d(psi)/dt = 0.)
      !  sets up s false trend to calculate the barotropic streamfunction.

      ua(:,:,:) = ub(:,:,:) / rdt
      va(:,:,:) = vb(:,:,:) / rdt

      ! calls dyn_spg. we assume euler time step, starting from rest.
      indic = 0
      CALL dyn_spg( nit000, indic )       ! surface pressure gradient

      ! the new velocity is ua*rdt

      CALL lbc_lnk( ua, 'U', -1. )
      CALL lbc_lnk( va, 'V', -1. )

      ub(:,:,:) = ua(:,:,:) * rdt
      vb(:,:,:) = va(:,:,:) * rdt
      ua(:,:,:) = 0.e0
      va(:,:,:) = 0.e0
      un(:,:,:) = ub(:,:,:)
      vn(:,:,:) = vb(:,:,:)
       
      ! Compute the divergence and curl

      CALL div_cur( nit000 )            ! now horizontal divergence and curl

      hdivb(:,:,:) = hdivn(:,:,:)       ! set the before to the now value
      rotb (:,:,:) = rotn (:,:,:)       ! set the before to the now value
      !
      IF( wrk_not_released(3, 1) ) THEN
         CALL ctl_stop('istate_uvg: failed to release workspace array')
      ENDIF
      !
   END SUBROUTINE istate_uvg

   !!=====================================================================
END MODULE istate
