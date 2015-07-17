MODULE p4zsed
   !!======================================================================
   !!                         ***  MODULE p4sed  ***
   !! TOP :   PISCES Compute loss of organic matter in the sediments
   !!======================================================================
   !! History :   1.0  !  2004-03 (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   'key_pisces'                                       PISCES bio-model
   !!----------------------------------------------------------------------
   !!   p4z_sed        :  Compute loss of organic matter in the sediments
   !!   p4z_sbc        :  Read and interpolate time-varying nutrients fluxes
   !!   p4z_sed_init   :  Initialization of p4z_sed
   !!----------------------------------------------------------------------
   USE trc
   USE oce_trc         !
   USE sms_pisces
   USE prtctl_trc
   USE p4zbio
   USE p4zint
   USE p4zopt
   USE p4zsink
   USE p4zrem
   USE p4zlim
   USE iom


   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_sed   
   PUBLIC   p4z_sed_init   
   PUBLIC   p4z_sed_alloc

   !! * Shared module variables
   LOGICAL, PUBLIC :: ln_dustfer  = .FALSE.    !: boolean for dust input from the atmosphere
   LOGICAL, PUBLIC :: ln_river    = .FALSE.    !: boolean for river input of nutrients
   LOGICAL, PUBLIC :: ln_ndepo    = .FALSE.    !: boolean for atmospheric deposition of N
   LOGICAL, PUBLIC :: ln_sedinput = .FALSE.    !: boolean for Fe input from sediments

   REAL(wp), PUBLIC :: sedfeinput = 1.E-9_wp   !: Coastal release of Iron
   REAL(wp), PUBLIC :: dustsolub  = 0.014_wp   !: Solubility of the dust

   !! * Module variables
   REAL(wp) :: ryyss                  !: number of seconds per year 
   REAL(wp) :: ryyss1                 !: inverse of ryyss
   REAL(wp) :: rmtss                  !: number of seconds per month
   REAL(wp) :: rday1                  !: inverse of rday

   INTEGER , PARAMETER :: jpmth = 12  !: number of months per year
   INTEGER , PARAMETER :: jpyr  = 1   !: one year

   INTEGER ::  numdust                !: logical unit for surface fluxes data
   INTEGER ::  nflx1 , nflx2          !: first and second record used
   INTEGER ::  nflx11, nflx12

   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: dustmo    !: set of dust fields
   REAL(wp), ALLOCATABLE, SAVE,   DIMENSION(:,:) :: dust      !: dust fields
   REAL(wp), ALLOCATABLE, SAVE,   DIMENSION(:,:) :: rivinp, cotdep    !: river input fields
   REAL(wp), ALLOCATABLE, SAVE,   DIMENSION(:,:) :: nitdep    !: atmospheric N deposition 
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: ironsed   !: Coastal supply of iron

   REAL(wp) :: sumdepsi, rivalkinput, rivpo4input, nitdepinput

   !!* Substitution
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Header:$ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS


   SUBROUTINE p4z_sed( kt, jnt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_sed  ***
      !!
      !! ** Purpose :   Compute loss of organic matter in the sediments. This
      !!              is by no way a sediment model. The loss is simply 
      !!              computed to balance the inout from rivers and dust
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      USE wrk_nemo, ONLY: wrk_in_use, wrk_not_released
      USE wrk_nemo, ONLY: zsidep => wrk_2d_1, zwork => wrk_2d_2, zwork1 => wrk_2d_3
      USE wrk_nemo, ONLY: znitrpot => wrk_3d_2, zirondep => wrk_3d_3
      !
      INTEGER, INTENT(in) ::   kt, jnt ! ocean time step
      INTEGER  ::   ji, jj, jk, ikt
#if ! defined key_sed
      REAL(wp) ::   zsumsedsi, zsumsedpo4, zsumsedcal
      REAL(wp) ::   zrivalk, zrivsil, zrivpo4
#endif
      REAL(wp) ::   zdenitot, znitrpottot, zlim, zfact
      REAL(wp) ::   zwsbio3, zwsbio4, zwscal
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------

      IF( ( wrk_in_use(2, 1,2,3) ) .OR. ( wrk_in_use(3, 2,3) ) ) THEN
         CALL ctl_stop('p4z_sed: requested workspace arrays unavailable')  ;  RETURN
      END IF

      IF( jnt == 1  .AND.  ln_dustfer  )  CALL p4z_sbc( kt )

      ! Iron and Si deposition at the surface
      ! -------------------------------------

      DO jj = 1, jpj
         DO ji = 1, jpi
            zirondep(ji,jj,1) = ( dustsolub * dust(ji,jj) / ( 55.85 * rmtss ) + 3.e-10 * ryyss1 )   &
               &             * rfact2 / fse3t(ji,jj,1)
            zsidep  (ji,jj)   = 8.8 * 0.075 * dust(ji,jj) * rfact2 / ( fse3t(ji,jj,1) * 28.1 * rmtss )
         END DO
      END DO

      ! Iron solubilization of particles in the water column
      ! ----------------------------------------------------

      DO jk = 2, jpkm1
         zirondep(:,:,jk) = dust(:,:) / ( 10. * 55.85 * rmtss ) * rfact2 * 1.e-4
      END DO

      ! Add the external input of nutrients, carbon and alkalinity
      ! ----------------------------------------------------------

      trn(:,:,1,jppo4) = trn(:,:,1,jppo4) + rivinp(:,:) * rfact2 
      trn(:,:,1,jpno3) = trn(:,:,1,jpno3) + (rivinp(:,:) + nitdep(:,:)) * rfact2
      trn(:,:,1,jpfer) = trn(:,:,1,jpfer) + rivinp(:,:) * 3.e-5 * rfact2
      trn(:,:,1,jpsil) = trn(:,:,1,jpsil) + zsidep (:,:) + cotdep(:,:)   * rfact2 / 6.
      trn(:,:,1,jpdic) = trn(:,:,1,jpdic) + rivinp(:,:) * 2.631 * rfact2
      trn(:,:,1,jptal) = trn(:,:,1,jptal) + (cotdep(:,:) - rno3*(rivinp(:,:) +  nitdep(:,:) ) ) * rfact2


      ! Add the external input of iron which is 3D distributed
      ! (dust, river and sediment mobilization)
      ! ------------------------------------------------------

      DO jk = 1, jpkm1
         trn(:,:,jk,jpfer) = trn(:,:,jk,jpfer) + zirondep(:,:,jk) + ironsed(:,:,jk) * rfact2
      END DO


#if ! defined key_sed
      ! Loss of biogenic silicon, Caco3 organic carbon in the sediments. 
      ! First, the total loss is computed.
      ! The factor for calcite comes from the alkalinity effect
      ! -------------------------------------------------------------
      DO jj = 1, jpj
         DO ji = 1, jpi
            ikt = mbkt(ji,jj) 
# if defined key_kriest
            zwork (ji,jj) = trn(ji,jj,ikt,jpdsi) * wscal (ji,jj,ikt)
            zwork1(ji,jj) = trn(ji,jj,ikt,jppoc) * wsbio3(ji,jj,ikt)
# else
            zwork (ji,jj) = trn(ji,jj,ikt,jpdsi) * wsbio4(ji,jj,ikt)
            zwork1(ji,jj) = trn(ji,jj,ikt,jpgoc) * wsbio4(ji,jj,ikt) + trn(ji,jj,ikt,jppoc) * wsbio3(ji,jj,ikt) 
# endif
         END DO
      END DO
      zsumsedsi  = glob_sum( zwork (:,:) * e1e2t(:,:) ) * rday1
      zsumsedpo4 = glob_sum( zwork1(:,:) * e1e2t(:,:) ) * rday1
      DO jj = 1, jpj
         DO ji = 1, jpi
            ikt = mbkt(ji,jj) 
            zwork (ji,jj) = trn(ji,jj,ikt,jpcal) * wscal (ji,jj,ikt)
         END DO
      END DO
      zsumsedcal = glob_sum( zwork (:,:) * e1e2t(:,:) ) * 2.0 * rday1
#endif

      ! Then this loss is scaled at each bottom grid cell for
      ! equilibrating the total budget of silica in the ocean.
      ! Thus, the amount of silica lost in the sediments equal
      ! the supply at the surface (dust+rivers)
      ! ------------------------------------------------------

      DO jj = 1, jpj
         DO ji = 1, jpi
            ikt = mbkt(ji,jj)
            zfact = xstep / fse3t(ji,jj,ikt)
            zwsbio3 = 1._wp - zfact * wsbio3(ji,jj,ikt)
            zwsbio4 = 1._wp - zfact * wsbio4(ji,jj,ikt)
            zwscal  = 1._wp - zfact * wscal (ji,jj,ikt)
            !
# if defined key_kriest
            trn(ji,jj,ikt,jpdsi) = trn(ji,jj,ikt,jpdsi) * zwsbio4
            trn(ji,jj,ikt,jpnum) = trn(ji,jj,ikt,jpnum) * zwsbio4
            trn(ji,jj,ikt,jppoc) = trn(ji,jj,ikt,jppoc) * zwsbio3
            trn(ji,jj,ikt,jpsfe) = trn(ji,jj,ikt,jpsfe) * zwsbio3
# else
            trn(ji,jj,ikt,jpdsi) = trn(ji,jj,ikt,jpdsi) * zwscal 
            trn(ji,jj,ikt,jpgoc) = trn(ji,jj,ikt,jpgoc) * zwsbio4
            trn(ji,jj,ikt,jppoc) = trn(ji,jj,ikt,jppoc) * zwsbio3
            trn(ji,jj,ikt,jpbfe) = trn(ji,jj,ikt,jpbfe) * zwsbio4
            trn(ji,jj,ikt,jpsfe) = trn(ji,jj,ikt,jpsfe) * zwsbio3
# endif
            trn(ji,jj,ikt,jpcal) = trn(ji,jj,ikt,jpcal) * zwscal
         END DO
      END DO

#if ! defined key_sed
      zrivsil =  1._wp - ( sumdepsi + rivalkinput * ryyss1 / 6. ) / zsumsedsi 
      zrivalk =  1._wp - ( rivalkinput * ryyss1 ) / zsumsedcal 
      zrivpo4 =  1._wp - ( rivpo4input * ryyss1 ) / zsumsedpo4 
      DO jj = 1, jpj
         DO ji = 1, jpi
            ikt = mbkt(ji,jj)
            zfact = xstep / fse3t(ji,jj,ikt)
            zwsbio3 = zfact * wsbio3(ji,jj,ikt)
            zwsbio4 = zfact * wsbio4(ji,jj,ikt)
            zwscal  = zfact * wscal (ji,jj,ikt)
            trn(ji,jj,ikt,jptal) =  trn(ji,jj,ikt,jptal) + trn(ji,jj,ikt,jpcal) * zwscal  * zrivalk * 2.0
            trn(ji,jj,ikt,jpdic) =  trn(ji,jj,ikt,jpdic) + trn(ji,jj,ikt,jpcal) * zwscal  * zrivalk
# if defined key_kriest
            trn(ji,jj,ikt,jpsil) =  trn(ji,jj,ikt,jpsil) + trn(ji,jj,ikt,jpdsi) * zwsbio4 * zrivsil 
            trn(ji,jj,ikt,jpdoc) =  trn(ji,jj,ikt,jpdoc) + trn(ji,jj,ikt,jppoc) * zwsbio3 * zrivpo4 
# else
            trn(ji,jj,ikt,jpsil) =  trn(ji,jj,ikt,jpsil) + trn(ji,jj,ikt,jpdsi) * zwscal  * zrivsil 
            trn(ji,jj,ikt,jpdoc) =  trn(ji,jj,ikt,jpdoc)   &
            &                     + ( trn(ji,jj,ikt,jppoc) * zwsbio3 + trn(ji,jj,ikt,jpgoc) * zwsbio4 ) * zrivpo4
# endif
         END DO
      END DO
# endif

      ! Nitrogen fixation (simple parameterization). The total gain
      ! from nitrogen fixation is scaled to balance the loss by 
      ! denitrification
      ! -------------------------------------------------------------

      zdenitot = glob_sum( denitr(:,:,:)  * cvol(:,:,:) * xnegtr(:,:,:) ) * rdenit

      ! Potential nitrogen fixation dependant on temperature and iron
      ! -------------------------------------------------------------

!CDIR NOVERRCHK
      DO jk = 1, jpk
!CDIR NOVERRCHK
         DO jj = 1, jpj
!CDIR NOVERRCHK
            DO ji = 1, jpi
               zlim = ( 1.- xnanono3(ji,jj,jk) - xnanonh4(ji,jj,jk) )
               IF( zlim <= 0.2 )   zlim = 0.01
               znitrpot(ji,jj,jk) = MAX( 0.e0, ( 0.6 * tgfunc(ji,jj,jk) - 2.15 ) * rday1 )   &
# if defined key_degrad
               &                  * facvol(ji,jj,jk)   &
# endif
               &                  * zlim * rfact2 * trn(ji,jj,jk,jpfer)   &
               &                  / ( conc3 + trn(ji,jj,jk,jpfer) ) * ( 1.- EXP( -etot(ji,jj,jk) / 50.) )
            END DO
         END DO 
      END DO

      znitrpottot = glob_sum( znitrpot(:,:,:) * cvol(:,:,:) )

      ! Nitrogen change due to nitrogen fixation
      ! ----------------------------------------

      DO jk = 1, jpk
         DO jj = 1, jpj
            DO ji = 1, jpi
               zfact = znitrpot(ji,jj,jk) * 1.e-7
               trn(ji,jj,jk,jpnh4) = trn(ji,jj,jk,jpnh4) + zfact
               trn(ji,jj,jk,jpoxy) = trn(ji,jj,jk,jpoxy) + zfact   * o2nit
               trn(ji,jj,jk,jppo4) = trn(ji,jj,jk,jppo4) + 30./ 46.* zfact
            END DO
         END DO
      END DO

#if defined key_diatrc
      zfact = 1.e+3 * rfact2r
#  if  ! defined key_iomput
      trc2d(:,:,jp_pcs0_2d + 11) = zirondep(:,:,1)         * zfact * fse3t(:,:,1) * tmask(:,:,1)
      trc2d(:,:,jp_pcs0_2d + 12) = znitrpot(:,:,1) * 1.e-7 * zfact * fse3t(:,:,1) * tmask(:,:,1)
#  else
      zwork (:,:)  =  ( zirondep(:,:,1) + ironsed(:,:,1) * rfact2 ) * zfact * fse3t(:,:,1) * tmask(:,:,1) 
      zwork1(:,:)  =  znitrpot(:,:,1) * 1.e-7                       * zfact * fse3t(:,:,1) * tmask(:,:,1)
      IF( jnt == nrdttrc ) THEN
         CALL iom_put( "Irondep", zwork  )  ! surface downward net flux of iron
         CALL iom_put( "Nfix"   , zwork1 )  ! nitrogen fixation at surface
      ENDIF
#  endif
#endif
      !
       IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('sed ')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=trn, mask=tmask, clinfo=ctrcnm)
       ENDIF

      IF( ( wrk_not_released(2, 1,2,3) ) .OR. ( wrk_not_released(3, 2,3) ) )   &
        &         CALL ctl_stop('p4z_sed: failed to release workspace arrays')

   END SUBROUTINE p4z_sed

   SUBROUTINE p4z_sbc( kt )

      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_sbc  ***
      !!
      !! ** Purpose :   Read and interpolate the external sources of 
      !!                nutrients
      !!
      !! ** Method  :   Read the files and interpolate the appropriate variables
      !!
      !! ** input   :   external netcdf files
      !!
      !!----------------------------------------------------------------------
      !! * arguments
      INTEGER, INTENT( in  ) ::   kt   ! ocean time step

      !! * Local declarations
      INTEGER :: imois, i15, iman 
      REAL(wp) :: zxy

      !!---------------------------------------------------------------------

      ! Initialization
      ! --------------

      i15 = nday / 16
      iman  = INT( raamo )
      imois = nmonth + i15 - 1
      IF( imois == 0 ) imois = iman

      ! Calendar computation
      IF( kt == nit000 .OR. imois /= nflx1 ) THEN

         IF( kt == nit000 )  nflx1  = 0

         ! nflx1 number of the first file record used in the simulation
         ! nflx2 number of the last  file record

         nflx1 = imois
         nflx2 = nflx1 + 1
         nflx1 = MOD( nflx1, iman )
         nflx2 = MOD( nflx2, iman )
         IF( nflx1 == 0 )   nflx1 = iman
         IF( nflx2 == 0 )   nflx2 = iman
         IF(lwp) WRITE(numout,*) 
         IF(lwp) WRITE(numout,*) ' p4z_sbc : first record file used nflx1 ',nflx1
         IF(lwp) WRITE(numout,*) ' p4z_sbc : last  record file used nflx2 ',nflx2

      ENDIF

      ! 3. at every time step interpolation of fluxes
      ! ---------------------------------------------

      zxy = FLOAT( nday + 15 - 30 * i15 ) / 30
      dust(:,:) = ( (1.-zxy) * dustmo(:,:,nflx1) + zxy * dustmo(:,:,nflx2) )

   END SUBROUTINE p4z_sbc


   SUBROUTINE p4z_sed_init

      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_sed_init  ***
      !!
      !! ** Purpose :   Initialization of the external sources of nutrients
      !!
      !! ** Method  :   Read the files and compute the budget
      !!      called at the first timestep (nit000)
      !!
      !! ** input   :   external netcdf files
      !!
      !!----------------------------------------------------------------------
      USE wrk_nemo, ONLY: wrk_in_use, wrk_not_released
      USE wrk_nemo, ONLY: zriverdoc => wrk_2d_1, zriver => wrk_2d_2, zndepo => wrk_2d_3
      USE wrk_nemo, ONLY: zcmask => wrk_3d_2
      !
      INTEGER :: ji, jj, jk, jm
      INTEGER :: numriv, numbath, numdep
      REAL(wp) ::   zcoef
      REAL(wp) ::   expide, denitide,zmaskt
      !
      NAMELIST/nampissed/ ln_dustfer, ln_river, ln_ndepo, ln_sedinput, sedfeinput, dustsolub
      !!----------------------------------------------------------------------

      IF( ( wrk_in_use(2, 1,2,3) ) .OR. ( wrk_in_use(3, 2) ) ) THEN
         CALL ctl_stop('p4z_sed_init: requested workspace arrays unavailable')  ;  RETURN
      END IF
      !
      REWIND( numnat )                     ! read numnat
      READ  ( numnat, nampissed )

      IF(lwp) THEN
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist : nampissed '
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~ '
         WRITE(numout,*) '    Dust input from the atmosphere           ln_dustfer  = ', ln_dustfer
         WRITE(numout,*) '    River input of nutrients                 ln_river    = ', ln_river
         WRITE(numout,*) '    Atmospheric deposition of N              ln_ndepo    = ', ln_ndepo
         WRITE(numout,*) '    Fe input from sediments                  ln_sedinput = ', ln_sedinput
         WRITE(numout,*) '    Coastal release of Iron                  sedfeinput  =', sedfeinput
         WRITE(numout,*) '    Solubility of the dust                   dustsolub   =', dustsolub
      ENDIF

      ! Dust input from the atmosphere
      ! ------------------------------
      IF( ln_dustfer ) THEN 
         IF(lwp) WRITE(numout,*) '    Initialize dust input from atmosphere '
         IF(lwp) WRITE(numout,*) '    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ '
         CALL iom_open ( 'dust.orca.nc', numdust )
         DO jm = 1, jpmth
            CALL iom_get( numdust, jpdom_data, 'dust', dustmo(:,:,jm), jm )
         END DO
         CALL iom_close( numdust )
      ELSE
         dustmo(:,:,:) = 0.e0
         dust(:,:) = 0.0
      ENDIF

      ! Nutrient input from rivers
      ! --------------------------
      IF( ln_river ) THEN
         IF(lwp) WRITE(numout,*) '    Initialize the nutrient input by rivers from river.orca.nc file'
         IF(lwp) WRITE(numout,*) '    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         CALL iom_open ( 'river.orca.nc', numriv )
         CALL iom_get  ( numriv, jpdom_data, 'riverdic', zriver   (:,:), jpyr )
         CALL iom_get  ( numriv, jpdom_data, 'riverdoc', zriverdoc(:,:), jpyr )
         CALL iom_close( numriv )
      ELSE
         zriver   (:,:) = 0.e0
         zriverdoc(:,:) = 0.e0
      endif

      ! Nutrient input from dust
      ! ------------------------
      IF( ln_ndepo ) THEN
         IF(lwp) WRITE(numout,*) '    Initialize the nutrient input by dust from ndeposition.orca.nc'
         IF(lwp) WRITE(numout,*) '    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         CALL iom_open ( 'ndeposition.orca.nc', numdep )
         CALL iom_get  ( numdep, jpdom_data, 'ndep', zndepo(:,:), jpyr )
         CALL iom_close( numdep )
      ELSE
         zndepo(:,:) = 0.e0
      ENDIF

      ! Coastal and island masks
      ! ------------------------
      IF( ln_sedinput ) THEN     
         IF(lwp) WRITE(numout,*) '    Computation of an island mask to enhance coastal supply of iron'
         IF(lwp) WRITE(numout,*) '    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         IF(lwp) WRITE(numout,*) '       from bathy.orca.nc file '
         CALL iom_open ( 'bathy.orca.nc', numbath )
         CALL iom_get  ( numbath, jpdom_data, 'bathy', zcmask(:,:,:), jpyr )
         CALL iom_close( numbath )
         !
         DO jk = 1, 5
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1
                  IF( tmask(ji,jj,jk) /= 0. ) THEN
                     zmaskt = tmask(ji+1,jj,jk) * tmask(ji-1,jj,jk) * tmask(ji,jj+1,jk)    &
                        &                       * tmask(ji,jj-1,jk) * tmask(ji,jj,jk+1)
                     IF( zmaskt == 0. )   zcmask(ji,jj,jk ) = MAX( 0.1, zcmask(ji,jj,jk) ) 
                  ENDIF
               END DO
            END DO
         END DO
         DO jk = 1, jpk
            DO jj = 1, jpj
               DO ji = 1, jpi
                  expide   = MIN( 8.,( fsdept(ji,jj,jk) / 500. )**(-1.5) )
                  denitide = -0.9543 + 0.7662 * LOG( expide ) - 0.235 * LOG( expide )**2
                  zcmask(ji,jj,jk) = zcmask(ji,jj,jk) * MIN( 1., EXP( denitide ) / 0.5 )
               END DO
            END DO
         END DO
      ELSE
         zcmask(:,:,:) = 0.e0
      ENDIF

      CALL lbc_lnk( zcmask , 'T', 1. )      ! Lateral boundary conditions on zcmask   (sign unchanged)


      !                                    ! Number of seconds per year and per month
      ryyss  = nyear_len(1) * rday
      rmtss  = ryyss / raamo
      rday1  = 1. / rday
      ryyss1 = 1. / ryyss
      !                                    ! ocean surface cell

      ! total atmospheric supply of Si
      ! ------------------------------
      sumdepsi = 0.e0
      DO jm = 1, jpmth
         zcoef = 1. / ( 12. * rmtss ) * 8.8 * 0.075 / 28.1        
         sumdepsi = sumdepsi + glob_sum( dustmo(:,:,jm) * e1e2t(:,:) ) * zcoef
      ENDDO

      ! N/P and Si releases due to coastal rivers
      ! -----------------------------------------
      DO jj = 1, jpj
         DO ji = 1, jpi
            zcoef = ryyss * e1e2t(ji,jj)  * fse3t(ji,jj,1) * tmask(ji,jj,1) 
            cotdep(ji,jj) =  zriver(ji,jj)                  *1E9 / ( 12. * zcoef + rtrn )
            rivinp(ji,jj) = (zriver(ji,jj)+zriverdoc(ji,jj)) *1E9 / ( 31.6* zcoef + rtrn )
            nitdep(ji,jj) = 7.6 * zndepo(ji,jj)                  / ( 14E6*ryyss*fse3t(ji,jj,1) + rtrn )
         END DO
      END DO
      ! Lateral boundary conditions on ( cotdep, rivinp, nitdep )   (sign unchanged)
      CALL lbc_lnk( cotdep , 'T', 1. )  ;  CALL lbc_lnk( rivinp , 'T', 1. )  ;  CALL lbc_lnk( nitdep , 'T', 1. )

      rivpo4input = glob_sum( rivinp(:,:) * cvol(:,:,1) ) * ryyss
      rivalkinput = glob_sum( cotdep(:,:) * cvol(:,:,1) ) * ryyss
      nitdepinput = glob_sum( nitdep(:,:) * cvol(:,:,1) ) * ryyss


      ! Coastal supply of iron
      ! -------------------------
      DO jk = 1, jpkm1
         ironsed(:,:,jk) = sedfeinput * zcmask(:,:,jk) / ( fse3t(:,:,jk) * rday )
      END DO
      CALL lbc_lnk( ironsed , 'T', 1. )      ! Lateral boundary conditions on ( ironsed )   (sign unchanged)

      IF( ( wrk_not_released(2, 1,2,3) ) .OR. ( wrk_not_released(3, 2) ) )   &
        &         CALL ctl_stop('p4z_sed_init: failed to release workspace arrays')

   END SUBROUTINE p4z_sed_init

   INTEGER FUNCTION p4z_sed_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_sed_alloc  ***
      !!----------------------------------------------------------------------

      ALLOCATE( dustmo(jpi,jpj,jpmth), dust(jpi,jpj)       ,     &
        &       rivinp(jpi,jpj)      , cotdep(jpi,jpj)     ,     &
        &       nitdep(jpi,jpj)      , ironsed(jpi,jpj,jpk), STAT=p4z_sed_alloc )  

      IF( p4z_sed_alloc /= 0 ) CALL ctl_warn('p4z_sed_alloc : failed to allocate arrays.')

   END FUNCTION p4z_sed_alloc
#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_sed                         ! Empty routine
   END SUBROUTINE p4z_sed
#endif 

   !!======================================================================
END MODULE  p4zsed
