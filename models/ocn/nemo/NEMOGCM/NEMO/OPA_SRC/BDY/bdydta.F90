MODULE bdydta
   !!======================================================================
   !!                       ***  MODULE bdydta  ***
   !! Open boundary data : read the data for the unstructured open boundaries.
   !!======================================================================
   !! History :  1.0  !  2005-01  (J. Chanut, A. Sellar)  Original code
   !!             -   !  2007-01  (D. Storkey) Update to use IOM module
   !!             -   !  2007-07  (D. Storkey) add bdy_dta_fla
   !!            3.0  !  2008-04  (NEMO team)  add in the reference version
   !!            3.3  !  2010-09  (E.O'Dea) modifications for Shelf configurations 
   !!            3.3  !  2010-09  (D.Storkey) add ice boundary conditions
   !!----------------------------------------------------------------------
#if defined key_bdy
   !!----------------------------------------------------------------------
   !!   'key_bdy'                     Unstructured Open Boundary Conditions
   !!----------------------------------------------------------------------
   !!   bdy_dta_frs    : read u, v, t, s data along open boundaries
   !!   bdy_dta_fla : read depth-mean velocities and elevation along open boundaries        
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE bdy_oce         ! ocean open boundary conditions
   USE bdytides        ! tidal forcing at boundaries
   USE iom
   USE ioipsl
   USE in_out_manager  ! I/O logical units
#if defined key_lim2
   USE ice_2
#endif

   IMPLICIT NONE
   PRIVATE

   PUBLIC   bdy_dta_frs      ! routines called by step.F90
   PUBLIC   bdy_dta_fla 
   PUBLIC   bdy_dta_alloc    ! routine called by bdy_init.F90

   INTEGER ::   numbdyt, numbdyu, numbdyv                      ! logical units for T-, U-, & V-points data file, resp.
   INTEGER ::   ntimes_bdy                                     ! exact number of time dumps in data files
   INTEGER ::   nbdy_b, nbdy_a                                 ! record of bdy data file for before and after time step
   INTEGER ::   numbdyt_bt, numbdyu_bt, numbdyv_bt             ! logical unit for T-, U- & V-points data file, resp.
   INTEGER ::   ntimes_bdy_bt                                  ! exact number of time dumps in data files
   INTEGER ::   nbdy_b_bt, nbdy_a_bt                           ! record of bdy data file for before and after time step

   INTEGER, DIMENSION (jpbtime) ::   istep, istep_bt           ! time array in seconds in each data file

   REAL(wp) ::  zoffset                                        ! time offset between time origin in file & start time of model run

   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   tbdydta, sbdydta   ! time interpolated values of T and S bdy data   
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   ubdydta, vbdydta   ! time interpolated values of U and V bdy data 
   REAL(wp), DIMENSION(jpbdim,2)     ::   ubtbdydta, vbtbdydta ! Arrays used for time interpolation of bdy data   
   REAL(wp), DIMENSION(jpbdim,2)     ::   sshbdydta            ! bdy data of ssh

#if defined key_lim2
   REAL(wp), DIMENSION(jpbdim,2)     ::   frld_bdydta          ! }
   REAL(wp), DIMENSION(jpbdim,2)     ::   hicif_bdydta         ! } Arrays used for time interp. of ice bdy data 
   REAL(wp), DIMENSION(jpbdim,2)     ::   hsnif_bdydta         ! }
#endif

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: bdydta.F90 2715 2011-03-30 15:58:35Z rblod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

  FUNCTION bdy_dta_alloc()
     !!----------------------------------------------------------------------
     USE lib_mpp, ONLY: ctl_warn, mpp_sum
     !
     INTEGER :: bdy_dta_alloc
     !!----------------------------------------------------------------------
     !
     ALLOCATE(tbdydta(jpbdim,jpk,2), sbdydta(jpbdim,jpk,2), &
              ubdydta(jpbdim,jpk,2), vbdydta(jpbdim,jpk,2), Stat=bdy_dta_alloc)

     IF( lk_mpp           ) CALL mpp_sum ( bdy_dta_alloc )
     IF(bdy_dta_alloc /= 0) CALL ctl_warn('bdy_dta_alloc: failed to allocate arrays')

   END FUNCTION bdy_dta_alloc


   SUBROUTINE bdy_dta_frs( kt )
      !!----------------------------------------------------------------------
      !!                   ***  SUBROUTINE bdy_dta_frs  ***
      !!                    
      !! ** Purpose :   Read unstructured boundary data for FRS condition.
      !!
      !! ** Method  :   At the first timestep, read in boundary data for two
      !!                times from the file and time-interpolate. At other 
      !!                timesteps, check to see if we need another time from 
      !!                the file. If so read it in. Time interpolate.
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step index (for timesplitting option, otherwise zero)
      !!
      CHARACTER(LEN=80), DIMENSION(3) ::   clfile               ! names of input files
      CHARACTER(LEN=70 )              ::   clunits              ! units attribute of time coordinate
      LOGICAL ::   lect                                         ! flag for reading
      INTEGER ::   it, ib, ik, igrd                             ! dummy loop indices
      INTEGER ::   igrd_start, igrd_end                         ! start and end of loops on igrd
      INTEGER ::   idvar                                        ! netcdf var ID
      INTEGER ::   iman, i15, imois                             ! Time variables for monthly clim forcing
      INTEGER ::   ntimes_bdyt, ntimes_bdyu, ntimes_bdyv
      INTEGER ::   itimer, totime
      INTEGER ::   ii, ij                                       ! array addresses
      INTEGER ::   ipi, ipj, ipk, inum                          ! local integers (NetCDF read)
      INTEGER ::   iyear0, imonth0, iday0
      INTEGER ::   ihours0, iminutes0, isec0
      INTEGER ::   iyear, imonth, iday, isecs
      INTEGER, DIMENSION(jpbtime) ::   istept, istepu, istepv   ! time arrays from data files
      REAL(wp) ::   dayfrac, zxy, zoffsett
      REAL(wp) ::   zoffsetu, zoffsetv
      REAL(wp) ::   dayjul0, zdayjulini
      REAL(wp), DIMENSION(jpbtime)      ::   zstepr             ! REAL time array from data files
      REAL(wp), DIMENSION(jpbdta,1,jpk) ::   zdta               ! temporary array for data fields
      !!---------------------------------------------------------------------------


      IF( ln_dyn_frs .OR. ln_tra_frs    &
         &               .OR. ln_ice_frs ) THEN  ! If these are both false then this routine does nothing

      ! -------------------- !
      !    Initialization    !
      ! -------------------- !

      lect   = .false.           ! If true, read a time record

      ! Some time variables for monthly climatological forcing:
      ! *******************************************************

!!gm  here  use directely daymod calendar variables
 
      iman = INT( raamo )      ! Number of months in a year

      i15 = INT( 2*REAL( nday, wp ) / ( REAL( nmonth_len(nmonth), wp ) + 0.5 ) )
      ! i15=0 if the current day is in the first half of the month, else i15=1

      imois = nmonth + i15 - 1            ! imois is the first month record
      IF( imois == 0 )   imois = iman

      ! Time variable for non-climatological forcing:
      ! *********************************************
      itimer = (kt-nit000+1)*rdt      ! current time in seconds for interpolation 


      !                                                !-------------------!
      IF( kt == nit000 ) THEN                          !  First call only  !
         !                                             !-------------------!
         istep(:) = 0
         nbdy_b   = 0
         nbdy_a   = 0

         ! Get time information from bdy data file
         ! ***************************************

         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*)    'bdy_dta_frs : Initialize unstructured boundary data'
         IF(lwp) WRITE(numout,*)    '~~~~~~~' 

         IF     ( nn_dtactl == 0 ) THEN
            !
            IF(lwp) WRITE(numout,*) '          Bdy data are taken from initial conditions'
            !
         ELSEIF (nn_dtactl == 1) THEN
            !
            IF(lwp) WRITE(numout,*) '          Bdy data are read in netcdf files'
            !
            dayfrac = adatrj  - REAL( itimer, wp ) / 86400.   ! day fraction at time step kt-1
            dayfrac = dayfrac - INT ( dayfrac )               !
            totime  = ( nitend - nit000 + 1 ) * rdt           ! Total time of the run to verify that all the
            !                                                 ! necessary time dumps in file are included
            !
            clfile(1) = cn_dta_frs_T
            clfile(2) = cn_dta_frs_U
            clfile(3) = cn_dta_frs_V
            !                                                  
            ! how many files are we to read in?
            igrd_start = 1
            igrd_end   = 3
            IF(.NOT. ln_tra_frs .AND. .NOT. ln_ice_frs) THEN       ! No T-grid file.
               igrd_start = 2
            ELSEIF ( .NOT. ln_dyn_frs ) THEN                           ! No U-grid or V-grid file.
               igrd_end   = 1         
            ENDIF

            DO igrd = igrd_start, igrd_end                     !  loop over T, U & V grid  !
               !                                               !---------------------------!
               CALL iom_open( clfile(igrd), inum )
               CALL iom_gettime( inum, zstepr, kntime=ntimes_bdy, cdunits=clunits ) 

               SELECT CASE( igrd )
                  CASE (1)   ;   numbdyt = inum
                  CASE (2)   ;   numbdyu = inum
                  CASE (3)   ;   numbdyv = inum
               END SELECT

               ! Calculate time offset 
               READ(clunits,7000) iyear0, imonth0, iday0, ihours0, iminutes0, isec0
               ! Convert time origin in file to julian days 
               isec0 = isec0 + ihours0*60.*60. + iminutes0*60.
               CALL ymds2ju(iyear0, imonth0, iday0, REAL(isec0, wp), dayjul0)
               ! Compute model initialization time 
               iyear  = ndastp / 10000
               imonth = ( ndastp - iyear * 10000 ) / 100
               iday   = ndastp - iyear * 10000 - imonth * 100
               isecs  = dayfrac * 86400
               CALL ymds2ju(iyear, imonth, iday, REAL(isecs, wp) , zdayjulini)
               ! offset from initialization date:
               zoffset = (dayjul0-zdayjulini)*86400
               !
7000           FORMAT('seconds since ', I4.4,'-',I2.2,'-',I2.2,' ',I2.2,':',I2.2,':',I2.2)

               !! TO BE DONE... Check consistency between calendar from file 
               !! (available optionally from iom_gettime) and calendar in model 
               !! when calendar in model available outside of IOIPSL.

               IF(lwp) WRITE(numout,*) 'number of times: ',ntimes_bdy
               IF(lwp) WRITE(numout,*) 'offset: ',zoffset
               IF(lwp) WRITE(numout,*) 'totime: ',totime
               IF(lwp) WRITE(numout,*) 'zstepr: ',zstepr(1:ntimes_bdy)

               ! Check that there are not too many times in the file. 
               IF( ntimes_bdy > jpbtime ) THEN
                  WRITE(ctmp1,*) 'Check file: ', clfile(igrd), 'jpbtime= ', jpbtime, ' ntimes_bdy= ', ntimes_bdy
                  CALL ctl_stop( 'Number of time dumps in files exceed jpbtime parameter', ctmp1 )
               ENDIF

               ! Check that time array increases:
               it = 1
               DO WHILE( zstepr(it+1) > zstepr(it) .AND. it /= ntimes_bdy - 1 ) 
                  it = it + 1
               END DO
               !
               IF( it /= ntimes_bdy-1 .AND. ntimes_bdy > 1 ) THEN
                     WRITE(ctmp1,*) 'Check file: ', clfile(igrd)
                     CALL ctl_stop( 'Time array in unstructured boundary data files',   &
                        &           'does not continuously increase.'               , ctmp1 )
               ENDIF
               !
               ! Check that times in file span model run time:
               IF( zstepr(1) + zoffset > 0 ) THEN
                     WRITE(ctmp1,*) 'Check file: ', clfile(igrd)
                     CALL ctl_stop( 'First time dump in bdy file is after model initial time', ctmp1 )
               END IF
               IF( zstepr(ntimes_bdy) + zoffset < totime ) THEN
                     WRITE(ctmp1,*) 'Check file: ', clfile(igrd)
                     CALL ctl_stop( 'Last time dump in bdy file is before model final time', ctmp1 )
               END IF
               !
               SELECT CASE( igrd )
                  CASE (1)
                    ntimes_bdyt = ntimes_bdy
                    zoffsett = zoffset
                    istept(:) = INT( zstepr(:) + zoffset )
                    numbdyt = inum
                  CASE (2)
                    ntimes_bdyu = ntimes_bdy
                    zoffsetu = zoffset
                    istepu(:) = INT( zstepr(:) + zoffset )
                    numbdyu = inum
                  CASE (3)
                    ntimes_bdyv = ntimes_bdy
                    zoffsetv = zoffset
                    istepv(:) = INT( zstepr(:) + zoffset )
                    numbdyv = inum
               END SELECT
               !
            END DO                                         ! end loop over T, U & V grid 

            IF (igrd_start == 1 .and. igrd_end == 3) THEN
               ! Only test differences if we are reading in 3 files
               ! Verify time consistency between files  
               IF( ntimes_bdyu /= ntimes_bdyt .OR. ntimes_bdyv /= ntimes_bdyt ) THEN
                  CALL ctl_stop( 'Bdy data files must have the same number of time dumps',   &
                  &           'Multiple time frequencies not implemented yet'  )
               ENDIF
               ntimes_bdy = ntimes_bdyt
               !
               IF( zoffsetu /= zoffsett .OR. zoffsetv /= zoffsett ) THEN
                  CALL ctl_stop( 'Bdy data files must have the same time origin',   &
                  &           'Multiple time frequencies not implemented yet' )
               ENDIF
               zoffset = zoffsett
            ENDIF

            IF( igrd_start == 1 ) THEN   ;   istep(:) = istept(:)
            ELSE                         ;   istep(:) = istepu(:)
            ENDIF

            ! Check number of time dumps:              
            IF( ntimes_bdy == 1 .AND. .NOT. ln_clim ) THEN
              CALL ctl_stop( 'There is only one time dump in data files',   &
                 &           'Choose ln_clim=.true. in namelist for constant bdy forcing.' )
            ENDIF

            IF( ln_clim ) THEN
              IF( ntimes_bdy /= 1 .AND. ntimes_bdy /= 12 ) THEN
                 CALL ctl_stop( 'For climatological boundary forcing (ln_clim=.true.),',   &
                    &           'bdy data files must contain 1 or 12 time dumps.' )
              ELSEIF( ntimes_bdy ==  1 ) THEN
                IF(lwp) WRITE(numout,*)
                IF(lwp) WRITE(numout,*) 'We assume constant boundary forcing from bdy data files'
              ELSEIF( ntimes_bdy == 12 ) THEN
                IF(lwp) WRITE(numout,*)
                IF(lwp) WRITE(numout,*) 'We assume monthly (and cyclic) boundary forcing from bdy data files'
              ENDIF
            ENDIF

            ! Find index of first record to read (before first model time). 
            it = 1
            DO WHILE( istep(it+1) <= 0 .AND. it <= ntimes_bdy - 1 )
               it = it + 1
            END DO
            nbdy_b = it
            !
            IF(lwp) WRITE(numout,*) 'Time offset is ',zoffset
            IF(lwp) WRITE(numout,*) 'First record to read is ',nbdy_b

         ENDIF ! endif (nn_dtactl == 1)


         ! 1.2  Read first record in file if necessary (ie if nn_dtactl == 1)
         ! *****************************************************************

         IF( nn_dtactl == 0 ) THEN      ! boundary data arrays are filled with initial conditions
            !
            IF (ln_tra_frs) THEN
               igrd = 1            ! T-points data 
               DO ib = 1, nblen(igrd)
                  ii = nbi(ib,igrd)
                  ij = nbj(ib,igrd)
                  DO ik = 1, jpkm1
                     tbdy(ib,ik) = tn(ii,ij,ik)
                     sbdy(ib,ik) = sn(ii,ij,ik)
                  END DO
               END DO
            ENDIF

            IF(ln_dyn_frs) THEN
               igrd = 2            ! U-points data 
               DO ib = 1, nblen(igrd)
                  ii = nbi(ib,igrd)
                  ij = nbj(ib,igrd)
                  DO ik = 1, jpkm1
                     ubdy(ib,ik) = un(ii, ij, ik)
                  END DO
               END DO
               !
               igrd = 3            ! V-points data 
               DO ib = 1, nblen(igrd)            
                  ii = nbi(ib,igrd)
                  ij = nbj(ib,igrd)
                  DO ik = 1, jpkm1
                     vbdy(ib,ik) = vn(ii, ij, ik)
                  END DO
               END DO
            ENDIF
            !
#if defined key_lim2
            IF( ln_ice_frs ) THEN
               igrd = 1            ! T-points data
               DO ib = 1, nblen(igrd)
                  frld_bdy (ib) =  frld(nbi(ib,igrd), nbj(ib,igrd))
                  hicif_bdy(ib) = hicif(nbi(ib,igrd), nbj(ib,igrd))
                  hsnif_bdy(ib) = hsnif(nbi(ib,igrd), nbj(ib,igrd))
               END DO
            ENDIF
#endif
         ELSEIF( nn_dtactl == 1 ) THEN    ! Set first record in the climatological case:   
            !
            IF( ln_clim .AND. ntimes_bdy == 1 ) THEN
               nbdy_a = 1
            ELSEIF( ln_clim .AND. ntimes_bdy == iman ) THEN
               nbdy_b = 0
               nbdy_a = imois
            ELSE
               nbdy_a = nbdy_b
            ENDIF
   
            ! Read first record:
            ipj  = 1
            ipk  = jpk
            igrd = 1
            ipi  = nblendta(igrd)

            IF(ln_tra_frs) THEN
               !
               igrd = 1                                           ! Temperature
               IF( nblendta(igrd) <=  0 ) THEN 
                  idvar = iom_varid( numbdyt, 'votemper' )
                  nblendta(igrd) = iom_file(numbdyt)%dimsz(1,idvar)
               ENDIF
               IF(lwp) WRITE(numout,*) 'Dim size for votemper is ', nblendta(igrd)
               ipi = nblendta(igrd)
               CALL iom_get ( numbdyt, jpdom_unknown, 'votemper', zdta(1:ipi,1:ipj,1:ipk), nbdy_a )
               !
               DO ib = 1, nblen(igrd)
                  DO ik = 1, jpkm1
                     tbdydta(ib,ik,2) =  zdta(nbmap(ib,igrd),1,ik)
                  END DO
               END DO
               !
               igrd = 1                                           ! salinity
               IF( nblendta(igrd) .le. 0 ) THEN 
                  idvar = iom_varid( numbdyt, 'vosaline' )
                  nblendta(igrd) = iom_file(numbdyt)%dimsz(1,idvar)
               ENDIF
               IF(lwp) WRITE(numout,*) 'Dim size for vosaline is ', nblendta(igrd)
               ipi = nblendta(igrd)
               CALL iom_get ( numbdyt, jpdom_unknown, 'vosaline', zdta(1:ipi,1:ipj,1:ipk), nbdy_a )
               !
               DO ib = 1, nblen(igrd)
                  DO ik = 1, jpkm1
                     sbdydta(ib,ik,2) =  zdta(nbmap(ib,igrd),1,ik)
                  END DO
               END DO
            ENDIF  ! ln_tra_frs
 
            IF( ln_dyn_frs ) THEN
               !
               igrd = 2                                           ! u-velocity
               IF ( nblendta(igrd) .le. 0 ) THEN 
                 idvar = iom_varid( numbdyu,'vozocrtx' )
                 nblendta(igrd) = iom_file(numbdyu)%dimsz(1,idvar)
               ENDIF
               IF(lwp) WRITE(numout,*) 'Dim size for vozocrtx is ', nblendta(igrd)
               ipi = nblendta(igrd)
               CALL iom_get ( numbdyu, jpdom_unknown,'vozocrtx',zdta(1:ipi,1:ipj,1:ipk),nbdy_a )
               DO ib = 1, nblen(igrd)
                  DO ik = 1, jpkm1
                     ubdydta(ib,ik,2) =  zdta(nbmap(ib,igrd),1,ik)
                  END DO
               END DO
               !
               igrd = 3                                           ! v-velocity
               IF ( nblendta(igrd) .le. 0 ) THEN 
                 idvar = iom_varid( numbdyv,'vomecrty' )
                 nblendta(igrd) = iom_file(numbdyv)%dimsz(1,idvar)
               ENDIF
               IF(lwp) WRITE(numout,*) 'Dim size for vomecrty is ', nblendta(igrd)
               ipi = nblendta(igrd)
               CALL iom_get ( numbdyv, jpdom_unknown,'vomecrty',zdta(1:ipi,1:ipj,1:ipk),nbdy_a )
               DO ib = 1, nblen(igrd)
                  DO ik = 1, jpkm1
                     vbdydta(ib,ik,2) =  zdta(nbmap(ib,igrd),1,ik)
                  END DO
               END DO
            ENDIF ! ln_dyn_frs

#if defined key_lim2
            IF( ln_ice_frs ) THEN
              !
              igrd=1                                              ! leads fraction
              IF(lwp) WRITE(numout,*) 'Dim size for ildsconc is ',nblendta(igrd)
              ipi=nblendta(igrd)
              CALL iom_get ( numbdyt, jpdom_unknown,'ildsconc',zdta(1:ipi,:,1),nbdy_a )
              DO ib=1, nblen(igrd)
                frld_bdydta(ib,2) =  zdta(nbmap(ib,igrd),1,1)
              END DO
              !
              igrd=1                                              ! ice thickness
              IF(lwp) WRITE(numout,*) 'Dim size for iicethic is ',nblendta(igrd)
              ipi=nblendta(igrd)
              CALL iom_get ( numbdyt, jpdom_unknown,'iicethic',zdta(1:ipi,:,1),nbdy_a )
              DO ib=1, nblen(igrd)
                hicif_bdydta(ib,2) =  zdta(nbmap(ib,igrd),1,1)
              END DO
              !
              igrd=1                                              ! snow thickness
              IF(lwp) WRITE(numout,*) 'Dim size for isnowthi is ',nblendta(igrd)
              ipi=nblendta(igrd)
              CALL iom_get ( numbdyt, jpdom_unknown,'isnowthi',zdta(1:ipi,:,1),nbdy_a )
              DO ib=1, nblen(igrd)
                hsnif_bdydta(ib,2) =  zdta(nbmap(ib,igrd),1,1)
              END DO
            ENDIF ! just if ln_ice_frs is set
#endif

            IF( .NOT.ln_clim .AND. istep(1) > 0 ) THEN     ! First data time is after start of run
               nbdy_b = nbdy_a                                 ! Put first value in both time levels
               IF( ln_tra_frs ) THEN
                 tbdydta(:,:,1) = tbdydta(:,:,2)
                 sbdydta(:,:,1) = sbdydta(:,:,2)
               ENDIF
               IF( ln_dyn_frs ) THEN
                 ubdydta(:,:,1) = ubdydta(:,:,2)
                 vbdydta(:,:,1) = vbdydta(:,:,2)
               ENDIF
#if defined key_lim2
               IF( ln_ice_frs ) THEN
                  frld_bdydta (:,1) =  frld_bdydta(:,2)
                  hicif_bdydta(:,1) = hicif_bdydta(:,2)
                  hsnif_bdydta(:,1) = hsnif_bdydta(:,2)
               ENDIF
#endif
            END IF
            !
         END IF   ! nn_dtactl == 0/1
 
         ! In the case of constant boundary forcing fill bdy arrays once for all
         IF( ln_clim .AND. ntimes_bdy == 1 ) THEN
            IF( ln_tra_frs ) THEN
               tbdy  (:,:) = tbdydta  (:,:,2)
               sbdy  (:,:) = sbdydta  (:,:,2)
            ENDIF
            IF( ln_dyn_frs) THEN
               ubdy  (:,:) = ubdydta  (:,:,2)
               vbdy  (:,:) = vbdydta  (:,:,2)
            ENDIF
#if defined key_lim2
            IF( ln_ice_frs ) THEN
               frld_bdy (:) = frld_bdydta (:,2)
               hicif_bdy(:) = hicif_bdydta(:,2)
               hsnif_bdy(:) = hsnif_bdydta(:,2)
            ENDIF
#endif

            IF( ln_tra_frs .OR. ln_ice_frs) CALL iom_close( numbdyt )
            IF( ln_dyn_frs                    ) CALL iom_close( numbdyu )
            IF( ln_dyn_frs                    ) CALL iom_close( numbdyv )
         END IF
         !
      ENDIF                                            ! End if nit000


      !                                                !---------------------!
      IF( nn_dtactl == 1 .AND. ntimes_bdy > 1 ) THEN    !  at each time step  !
         !                                             !---------------------!
         ! Read one more record if necessary
         !**********************************

         IF( ln_clim .AND. imois /= nbdy_b ) THEN      ! remember that nbdy_b=0 for kt=nit000
            nbdy_b = imois
            nbdy_a = imois + 1
            nbdy_b = MOD( nbdy_b, iman )   ;   IF( nbdy_b == 0 ) nbdy_b = iman
            nbdy_a = MOD( nbdy_a, iman )   ;   IF( nbdy_a == 0 ) nbdy_a = iman
            lect=.true.
         ELSEIF( .NOT.ln_clim .AND. itimer >= istep(nbdy_a) ) THEN

            IF( nbdy_a < ntimes_bdy ) THEN
               nbdy_b = nbdy_a
               nbdy_a = nbdy_a + 1
               lect  =.true.
            ELSE
               ! We have reached the end of the file
               ! put the last data time into both time levels
               nbdy_b = nbdy_a
               IF(ln_tra_frs) THEN
                  tbdydta(:,:,1) =  tbdydta(:,:,2)
                  sbdydta(:,:,1) =  sbdydta(:,:,2)
               ENDIF
               IF(ln_dyn_frs) THEN
                  ubdydta(:,:,1) =  ubdydta(:,:,2)
                  vbdydta(:,:,1) =  vbdydta(:,:,2)
               ENDIF
#if defined key_lim2
               IF(ln_ice_frs) THEN
                  frld_bdydta (:,1) =  frld_bdydta (:,2)
                  hicif_bdydta(:,1) =  hicif_bdydta(:,2)
                  hsnif_bdydta(:,1) =  hsnif_bdydta(:,2)
               ENDIF
#endif
            END IF ! nbdy_a < ntimes_bdy
            !
        END IF
         
        IF( lect ) THEN           ! Swap arrays
           IF( ln_tra_frs ) THEN
             tbdydta(:,:,1) =  tbdydta(:,:,2)
             sbdydta(:,:,1) =  sbdydta(:,:,2)
           ENDIF
           IF( ln_dyn_frs ) THEN
             ubdydta(:,:,1) =  ubdydta(:,:,2)
             vbdydta(:,:,1) =  vbdydta(:,:,2)
           ENDIF
#if defined key_lim2
           IF( ln_ice_frs ) THEN
             frld_bdydta (:,1) =  frld_bdydta (:,2)
             hicif_bdydta(:,1) =  hicif_bdydta(:,2)
             hsnif_bdydta(:,1) =  hsnif_bdydta(:,2)
           ENDIF
#endif 
           ! read another set
           ipj  = 1
           ipk  = jpk

           IF( ln_tra_frs ) THEN
              ! 
              igrd = 1                                   ! temperature
              ipi  = nblendta(igrd)
              CALL iom_get ( numbdyt, jpdom_unknown, 'votemper', zdta(1:ipi,1:ipj,1:ipk), nbdy_a )
              DO ib = 1, nblen(igrd)
                 DO ik = 1, jpkm1
                    tbdydta(ib,ik,2) = zdta(nbmap(ib,igrd),1,ik)
                 END DO
              END DO
              !
              igrd = 1                                   ! salinity
              ipi  = nblendta(igrd)
              CALL iom_get ( numbdyt, jpdom_unknown, 'vosaline', zdta(1:ipi,1:ipj,1:ipk), nbdy_a )
              DO ib = 1, nblen(igrd)
                 DO ik = 1, jpkm1
                    sbdydta(ib,ik,2) = zdta(nbmap(ib,igrd),1,ik)
                 END DO
              END DO
           ENDIF ! ln_tra_frs

           IF(ln_dyn_frs) THEN
              !
              igrd = 2                                   ! u-velocity
              ipi  = nblendta(igrd)
              CALL iom_get ( numbdyu, jpdom_unknown,'vozocrtx',zdta(1:ipi,1:ipj,1:ipk),nbdy_a )
              DO ib = 1, nblen(igrd)
                DO ik = 1, jpkm1
                  ubdydta(ib,ik,2) =  zdta(nbmap(ib,igrd),1,ik)
                END DO
              END DO
              !
              igrd = 3                                   ! v-velocity
              ipi  = nblendta(igrd)
              CALL iom_get ( numbdyv, jpdom_unknown,'vomecrty',zdta(1:ipi,1:ipj,1:ipk),nbdy_a )
              DO ib = 1, nblen(igrd)
                 DO ik = 1, jpkm1
                    vbdydta(ib,ik,2) =  zdta(nbmap(ib,igrd),1,ik)
                 END DO
              END DO
           ENDIF ! ln_dyn_frs
           !
#if defined key_lim2
           IF(ln_ice_frs) THEN
             !
             igrd = 1                                    ! ice concentration
             ipi=nblendta(igrd)
             CALL iom_get ( numbdyt, jpdom_unknown,'ildsconc',zdta(1:ipi,:,1),nbdy_a )
             DO ib=1, nblen(igrd)
               frld_bdydta(ib,2) =  zdta( nbmap(ib,igrd), 1, 1 )
             END DO
             !
             igrd=1                                      ! ice thickness
             ipi=nblendta(igrd)
             CALL iom_get ( numbdyt, jpdom_unknown,'iicethic',zdta(1:ipi,:,1),nbdy_a )
             DO ib=1, nblen(igrd)
               hicif_bdydta(ib,2) =  zdta( nbmap(ib,igrd), 1, 1 )
             END DO
             !
             igrd=1                                      ! snow thickness
             ipi=nblendta(igrd)
             CALL iom_get ( numbdyt, jpdom_unknown,'isnowthi',zdta(1:ipi,:,1),nbdy_a )
             DO ib=1, nblen(igrd)
               hsnif_bdydta(ib,2) =  zdta( nbmap(ib,igrd), 1, 1 )
             END DO
           ENDIF ! ln_ice_frs
#endif
           !
           IF(lwp) WRITE(numout,*) 'bdy_dta_frs : first record file used nbdy_b ',nbdy_b
           IF(lwp) WRITE(numout,*) '~~~~~~~~  last  record file used nbdy_a ',nbdy_a
           IF (.NOT.ln_clim) THEN
              IF(lwp) WRITE(numout,*) 'first  record time (s): ', istep(nbdy_b)
              IF(lwp) WRITE(numout,*) 'model time (s)        : ', itimer
              IF(lwp) WRITE(numout,*) 'second record time (s): ', istep(nbdy_a)
           ENDIF
           !
       ENDIF ! end lect=.true.


       ! Interpolate linearly
       ! ********************
       ! 
       IF( ln_clim ) THEN   ;   zxy = REAL( nday                   ) / REAL( nmonth_len(nbdy_b) ) + 0.5 - i15
       ELSEIF( istep(nbdy_b) == istep(nbdy_a) ) THEN 
                                    zxy = 0.0_wp
       ELSE                     ;   zxy = REAL( istep(nbdy_b) - itimer ) / REAL( istep(nbdy_b) - istep(nbdy_a) )
       END IF

          IF(ln_tra_frs) THEN
             igrd = 1                                   ! temperature & salinity
             DO ib = 1, nblen(igrd)
               DO ik = 1, jpkm1
                 tbdy(ib,ik) = zxy * tbdydta(ib,ik,2) + (1.-zxy) * tbdydta(ib,ik,1)
                 sbdy(ib,ik) = zxy * sbdydta(ib,ik,2) + (1.-zxy) * sbdydta(ib,ik,1)
               END DO
             END DO
          ENDIF

          IF(ln_dyn_frs) THEN
             igrd = 2                                   ! u-velocity
             DO ib = 1, nblen(igrd)
               DO ik = 1, jpkm1
                 ubdy(ib,ik) = zxy * ubdydta(ib,ik,2) + (1.-zxy) * ubdydta(ib,ik,1)   
               END DO
             END DO
             !
             igrd = 3                                   ! v-velocity
             DO ib = 1, nblen(igrd)
               DO ik = 1, jpkm1
                 vbdy(ib,ik) = zxy * vbdydta(ib,ik,2) + (1.-zxy) * vbdydta(ib,ik,1)   
               END DO
             END DO
          ENDIF

#if defined key_lim2
          IF(ln_ice_frs) THEN
            igrd=1
            DO ib=1, nblen(igrd)
               frld_bdy(ib) = zxy *  frld_bdydta(ib,2) + (1.-zxy) *  frld_bdydta(ib,1)
              hicif_bdy(ib) = zxy * hicif_bdydta(ib,2) + (1.-zxy) * hicif_bdydta(ib,1)
              hsnif_bdy(ib) = zxy * hsnif_bdydta(ib,2) + (1.-zxy) * hsnif_bdydta(ib,1)
            END DO
          ENDIF ! just if ln_ice_frs is true
#endif

      END IF                       !end if ((nn_dtactl==1).AND.(ntimes_bdy>1))
    

      !                                                !---------------------!
      !                                                !     last call       !
      !                                                !---------------------!
      IF( kt == nitend ) THEN
          IF(ln_tra_frs .or. ln_ice_frs) CALL iom_close( numbdyt )              ! Closing of the 3 files
          IF(ln_dyn_frs) CALL iom_close( numbdyu )
          IF(ln_dyn_frs) CALL iom_close( numbdyv )
      ENDIF
      !
      ENDIF ! ln_dyn_frs .OR. ln_tra_frs
      !
   END SUBROUTINE bdy_dta_frs


   SUBROUTINE bdy_dta_fla( kt, jit, icycl )
      !!---------------------------------------------------------------------------
      !!                      ***  SUBROUTINE bdy_dta_fla  ***
      !!                    
      !! ** Purpose :   Read unstructured boundary data for Flather condition
      !!
      !! ** Method  :  At the first timestep, read in boundary data for two
      !!               times from the file and time-interpolate. At other 
      !!               timesteps, check to see if we need another time from 
      !!               the file. If so read it in. Time interpolate.
      !!---------------------------------------------------------------------------
!!gm DOCTOR names :   argument integer :  start with "k"
      INTEGER, INTENT( in ) ::   kt          ! ocean time-step index
      INTEGER, INTENT( in ) ::   jit         ! barotropic time step index
      INTEGER, INTENT( in ) ::   icycl       ! number of cycles need for final file close
      !                                      ! (for timesplitting option, otherwise zero)
      !!
      LOGICAL ::   lect                      ! flag for reading
      INTEGER ::   it, ib, igrd              ! dummy loop indices
      INTEGER ::   idvar                     ! netcdf var ID
      INTEGER ::   iman, i15, imois          ! Time variables for monthly clim forcing
      INTEGER ::   ntimes_bdyt, ntimes_bdyu, ntimes_bdyv
      INTEGER ::   itimer, totime
      INTEGER ::   ipi, ipj, ipk, inum       ! temporary integers (NetCDF read)
      INTEGER ::   iyear0, imonth0, iday0
      INTEGER ::   ihours0, iminutes0, isec0
      INTEGER ::   iyear, imonth, iday, isecs
      INTEGER, DIMENSION(jpbtime) ::   istept, istepu, istepv   ! time arrays from data files
      REAL(wp) ::   dayfrac, zxy, zoffsett
      REAL(wp) ::   zoffsetu, zoffsetv
      REAL(wp) ::   dayjul0, zdayjulini
      REAL(wp) ::   zinterval_s, zinterval_e                    ! First and last interval in time axis
      REAL(wp), DIMENSION(jpbtime)      ::   zstepr             ! REAL time array from data files
      REAL(wp), DIMENSION(jpbdta,1)     ::   zdta               ! temporary array for data fields
      CHARACTER(LEN=80), DIMENSION(6)   ::   clfile
      CHARACTER(LEN=70 )                ::   clunits            ! units attribute of time coordinate
      !!---------------------------------------------------------------------------

!!gm   add here the same style as in bdy_dta_frs
!!gm      clearly bdy_dta_fla and bdy_dta_frs  can be combined...   
!!gm      too many things duplicated in the read of data...   simplification can be done

      ! -------------------- !
      !    Initialization    !
      ! -------------------- !

      lect   = .false.           ! If true, read a time record

      ! Some time variables for monthly climatological forcing:
      ! *******************************************************
 !!gm  here  use directely daymod variables
 
      iman  = INT( raamo ) ! Number of months in a year

      i15 = INT( 2*REAL( nday, wp ) / ( REAL( nmonth_len(nmonth), wp ) + 0.5 ) )
      ! i15=0 if the current day is in the first half of the month, else i15=1

      imois = nmonth + i15 - 1            ! imois is the first month record
      IF( imois == 0 ) imois = iman

      ! Time variable for non-climatological forcing:
      ! *********************************************

      itimer = ((kt-1)-nit000+1)*rdt                      ! current time in seconds for interpolation 
      itimer = itimer + jit*rdt/REAL(nn_baro,wp)      ! in non-climatological case

      IF ( ln_tides ) THEN

         ! -------------------------------------!
         ! Update BDY fields with tidal forcing !
         ! -------------------------------------!  

         CALL tide_update( kt, jit ) 
  
      ENDIF

      IF ( ln_dyn_fla ) THEN

         ! -------------------------------------!
         ! Update BDY fields with model data    !
         ! -------------------------------------!  

      !                                                !-------------------!
      IF( kt == nit000 .and. jit ==2 ) THEN            !  First call only  !
         !                                             !-------------------!
         istep_bt(:) = 0
         nbdy_b_bt    = 0
         nbdy_a_bt    = 0

         ! Get time information from bdy data file
         ! ***************************************

        IF(lwp) WRITE(numout,*)
        IF(lwp) WRITE(numout,*)    'bdy_dta_fla :Initialize unstructured boundary data for barotropic variables.'
        IF(lwp) WRITE(numout,*)    '~~~~~~~' 

        IF( nn_dtactl == 0 ) THEN
          IF(lwp) WRITE(numout,*)  'Bdy data are taken from initial conditions'

        ELSEIF (nn_dtactl == 1) THEN
          IF(lwp) WRITE(numout,*)  'Bdy data are read in netcdf files'

          dayfrac = adatrj  - REAL(itimer,wp)/86400. ! day fraction at time step kt-1
          dayfrac = dayfrac - INT (dayfrac)          !
          totime = (nitend-nit000+1)*rdt             ! Total time of the run to verify that all the
                                                     ! necessary time dumps in file are included

          clfile(4) = cn_dta_fla_T
          clfile(5) = cn_dta_fla_U
          clfile(6) = cn_dta_fla_V

          DO igrd = 4,6

            CALL iom_open( clfile(igrd), inum )
            CALL iom_gettime( inum, zstepr, kntime=ntimes_bdy_bt, cdunits=clunits ) 

            SELECT CASE( igrd )
               CASE (4) 
                  numbdyt_bt = inum
               CASE (5) 
                  numbdyu_bt = inum
               CASE (6) 
                  numbdyv_bt = inum
            END SELECT

            ! Calculate time offset 
            READ(clunits,7000) iyear0, imonth0, iday0, ihours0, iminutes0, isec0
            ! Convert time origin in file to julian days 
            isec0 = isec0 + ihours0*60.*60. + iminutes0*60.
            CALL ymds2ju(iyear0, imonth0, iday0, REAL(isec0, wp), dayjul0)
            ! Compute model initialization time 
            iyear  = ndastp / 10000
            imonth = ( ndastp - iyear * 10000 ) / 100
            iday   = ndastp - iyear * 10000 - imonth * 100
            isecs  = dayfrac * 86400
            CALL ymds2ju(iyear, imonth, iday, REAL(isecs, wp) , zdayjulini)
            ! zoffset from initialization date:
            zoffset = (dayjul0-zdayjulini)*86400
            !

7000 FORMAT('seconds since ', I4.4,'-',I2.2,'-',I2.2,' ',I2.2,':',I2.2,':',I2.2)

            !! TO BE DONE... Check consistency between calendar from file 
            !! (available optionally from iom_gettime) and calendar in model 
            !! when calendar in model available outside of IOIPSL.

            ! Check that there are not too many times in the file. 
            IF (ntimes_bdy_bt > jpbtime) CALL ctl_stop( &
                 'Number of time dumps in bdy file exceed jpbtime parameter', &
                 'Check file:' // TRIM(clfile(igrd))  )

            ! Check that time array increases (or interp will fail):
            DO it = 2, ntimes_bdy_bt
               IF ( zstepr(it-1) >= zstepr(it) ) THEN
                  CALL ctl_stop('Time array in unstructured boundary data file', &
                       'does not continuously increase.',               &
                       'Check file:' // TRIM(clfile(igrd))  )
                  EXIT
               END IF
            END DO

            IF ( .NOT. ln_clim ) THEN
               ! Check that times in file span model run time:

               ! Note: the fields may be time means, so we allow nit000 to be before
               ! first time in the file, provided that it falls inside the meaning
               ! period of the first field.  Until we can get the meaning period
               ! from the file, use the interval between fields as a proxy.
               ! If nit000 is before the first time, use the value at first time
               ! instead of extrapolating.  This is done by putting time 1 into
               ! both time levels.
               ! The same applies to the last time level: see setting of lect below.

               IF ( ntimes_bdy_bt == 1 ) CALL ctl_stop( &
                    'There is only one time dump in data files', &
                    'Set ln_clim=.true. in namelist for constant bdy forcing.' )

               zinterval_s = zstepr(2) - zstepr(1)
               zinterval_e = zstepr(ntimes_bdy_bt) - zstepr(ntimes_bdy_bt-1)

               IF( zstepr(1) + zoffset > 0 ) THEN
                     WRITE(ctmp1,*) 'Check file: ', clfile(igrd)
                     CALL ctl_stop( 'First time dump in bdy file is after model initial time', ctmp1 )
               END IF
               IF( zstepr(ntimes_bdy_bt) + zoffset < totime ) THEN
                     WRITE(ctmp1,*) 'Check file: ', clfile(igrd)
                     CALL ctl_stop( 'Last time dump in bdy file is before model final time', ctmp1 )
               END IF
            END IF ! .NOT. ln_clim

            IF ( igrd .EQ. 4) THEN
              ntimes_bdyt = ntimes_bdy_bt
              zoffsett = zoffset
              istept(:) = INT( zstepr(:) + zoffset )
            ELSE IF (igrd .EQ. 5) THEN
              ntimes_bdyu = ntimes_bdy_bt
              zoffsetu = zoffset
              istepu(:) = INT( zstepr(:) + zoffset )
            ELSE IF (igrd .EQ. 6) THEN
              ntimes_bdyv = ntimes_bdy_bt
              zoffsetv = zoffset
              istepv(:) = INT( zstepr(:) + zoffset )
            ENDIF

          ENDDO

      ! Verify time consistency between files  

          IF ( ntimes_bdyu /= ntimes_bdyt .OR. ntimes_bdyv /= ntimes_bdyt ) THEN
             CALL ctl_stop( &
             'Time axis lengths differ between bdy data files', &
             'Multiple time frequencies not implemented yet' )
          ELSE
            ntimes_bdy_bt = ntimes_bdyt
          ENDIF

          IF (zoffsetu.NE.zoffsett .OR. zoffsetv.NE.zoffsett) THEN
            CALL ctl_stop( & 
            'Bdy data files must have the same time origin', &
            'Multiple time frequencies not implemented yet'  )
          ENDIF
          zoffset = zoffsett

      !! Check that times are the same in the three files... HERE.
          istep_bt(:) = istept(:)

      ! Check number of time dumps:              
          IF (ln_clim) THEN
            SELECT CASE ( ntimes_bdy_bt )
            CASE( 1 )
              IF(lwp) WRITE(numout,*)
              IF(lwp) WRITE(numout,*) 'We assume constant boundary forcing from bdy data files'
              IF(lwp) WRITE(numout,*)             
            CASE( 12 )
              IF(lwp) WRITE(numout,*)
              IF(lwp) WRITE(numout,*) 'We assume monthly (and cyclic) boundary forcing from bdy data files'
              IF(lwp) WRITE(numout,*) 
            CASE DEFAULT
              CALL ctl_stop( &
                'For climatological boundary forcing (ln_clim=.true.),',&
                'bdy data files must contain 1 or 12 time dumps.' )
            END SELECT
          ENDIF

      ! Find index of first record to read (before first model time). 

          it=1
          DO WHILE ( ((istep_bt(it+1)) <= 0 ).AND.(it.LE.(ntimes_bdy_bt-1)))
            it=it+1
          END DO
          nbdy_b_bt = it

          IF(lwp) WRITE(numout,*) 'Time offset is ',zoffset
          IF(lwp) WRITE(numout,*) 'First record to read is ',nbdy_b_bt

        ENDIF ! endif (nn_dtactl == 1)

      ! 1.2  Read first record in file if necessary (ie if nn_dtactl == 1)
      ! *****************************************************************

        IF ( nn_dtactl == 0) THEN
          ! boundary data arrays are filled with initial conditions
          igrd = 5            ! U-points data 
          DO ib = 1, nblen(igrd)              
            ubtbdy(ib) = un(nbi(ib,igrd), nbj(ib,igrd), 1)
          END DO

          igrd = 6            ! V-points data 
          DO ib = 1, nblen(igrd)              
            vbtbdy(ib) = vn(nbi(ib,igrd), nbj(ib,igrd), 1)
          END DO

          igrd = 4            ! T-points data 
          DO ib = 1, nblen(igrd)              
            sshbdy(ib) = sshn(nbi(ib,igrd), nbj(ib,igrd))
          END DO

        ELSEIF (nn_dtactl == 1) THEN
 
        ! Set first record in the climatological case:   
          IF ((ln_clim).AND.(ntimes_bdy_bt==1)) THEN
            nbdy_a_bt = 1
          ELSEIF ((ln_clim).AND.(ntimes_bdy_bt==iman)) THEN
            nbdy_b_bt = 0
            nbdy_a_bt = imois
          ELSE
            nbdy_a_bt = nbdy_b_bt
          END IF
 
         ! Open Netcdf files:

          CALL iom_open ( cn_dta_fla_T, numbdyt_bt )
          CALL iom_open ( cn_dta_fla_U, numbdyu_bt )
          CALL iom_open ( cn_dta_fla_V, numbdyv_bt )

         ! Read first record:
          ipj=1
          igrd=4
          ipi=nblendta(igrd)

          ! ssh
          igrd=4
          IF ( nblendta(igrd) .le. 0 ) THEN 
            idvar = iom_varid( numbdyt_bt,'sossheig' )
            nblendta(igrd) = iom_file(numbdyt_bt)%dimsz(1,idvar)
          ENDIF
          WRITE(numout,*) 'Dim size for sossheig is ',nblendta(igrd)
          ipi=nblendta(igrd)

          CALL iom_get ( numbdyt_bt, jpdom_unknown,'sossheig',zdta(1:ipi,1:ipj),nbdy_a_bt )

          DO ib=1, nblen(igrd)
            sshbdydta(ib,2) =  zdta(nbmap(ib,igrd),1)
          END DO
 
          ! u-velocity
          igrd=5
          IF ( nblendta(igrd) .le. 0 ) THEN 
            idvar = iom_varid( numbdyu_bt,'vobtcrtx' )
            nblendta(igrd) = iom_file(numbdyu_bt)%dimsz(1,idvar)
          ENDIF
          WRITE(numout,*) 'Dim size for vobtcrtx is ',nblendta(igrd)
          ipi=nblendta(igrd)

          CALL iom_get ( numbdyu_bt, jpdom_unknown,'vobtcrtx',zdta(1:ipi,1:ipj),nbdy_a_bt )

          DO ib=1, nblen(igrd)
            ubtbdydta(ib,2) =  zdta(nbmap(ib,igrd),1)
          END DO

          ! v-velocity
          igrd=6
          IF ( nblendta(igrd) .le. 0 ) THEN 
            idvar = iom_varid( numbdyv_bt,'vobtcrty' )
            nblendta(igrd) = iom_file(numbdyv_bt)%dimsz(1,idvar)
          ENDIF
          WRITE(numout,*) 'Dim size for vobtcrty is ',nblendta(igrd)
          ipi=nblendta(igrd)

          CALL iom_get ( numbdyv_bt, jpdom_unknown,'vobtcrty',zdta(1:ipi,1:ipj),nbdy_a_bt )

          DO ib=1, nblen(igrd)
            vbtbdydta(ib,2) =  zdta(nbmap(ib,igrd),1)
          END DO

        END IF
 
        ! In the case of constant boundary forcing fill bdy arrays once for all
        IF ((ln_clim).AND.(ntimes_bdy_bt==1)) THEN

          ubtbdy  (:) = ubtbdydta  (:,2)
          vbtbdy  (:) = vbtbdydta  (:,2)
          sshbdy  (:) = sshbdydta  (:,2)

          CALL iom_close( numbdyt_bt )
          CALL iom_close( numbdyu_bt )
          CALL iom_close( numbdyv_bt )

        END IF

      ENDIF ! End if nit000

      ! -------------------- !
      ! 2. At each time step !
      ! -------------------- !

      IF ((nn_dtactl==1).AND.(ntimes_bdy_bt>1)) THEN 

      ! 2.1 Read one more record if necessary
      !**************************************

        IF ( (ln_clim).AND.(imois/=nbdy_b_bt) ) THEN ! remember that nbdy_b_bt=0 for kt=nit000
         nbdy_b_bt = imois
         nbdy_a_bt = imois+1
         nbdy_b_bt = MOD( nbdy_b_bt, iman )
         IF( nbdy_b_bt == 0 ) nbdy_b_bt = iman
         nbdy_a_bt = MOD( nbdy_a_bt, iman )
         IF( nbdy_a_bt == 0 ) nbdy_a_bt = iman
         lect=.true.

        ELSEIF ((.NOT.ln_clim).AND.(itimer >= istep_bt(nbdy_a_bt))) THEN
          nbdy_b_bt=nbdy_a_bt
          nbdy_a_bt=nbdy_a_bt+1
          lect=.true.
        END IF
         
        IF (lect) THEN

        ! Swap arrays
          sshbdydta(:,1) =  sshbdydta(:,2)
          ubtbdydta(:,1) =  ubtbdydta(:,2)
          vbtbdydta(:,1) =  vbtbdydta(:,2)
 
        ! read another set

          ipj=1
          ipk=jpk
          igrd=4
          ipi=nblendta(igrd)

          
          ! ssh
          igrd=4
          ipi=nblendta(igrd)

          CALL iom_get ( numbdyt_bt, jpdom_unknown,'sossheig',zdta(1:ipi,1:ipj),nbdy_a_bt )

          DO ib=1, nblen(igrd)
            sshbdydta(ib,2) =  zdta(nbmap(ib,igrd),1)
          END DO

          ! u-velocity
          igrd=5
          ipi=nblendta(igrd)

          CALL iom_get ( numbdyu_bt, jpdom_unknown,'vobtcrtx',zdta(1:ipi,1:ipj),nbdy_a_bt )

          DO ib=1, nblen(igrd)
            ubtbdydta(ib,2) =  zdta(nbmap(ib,igrd),1)
          END DO

          ! v-velocity
          igrd=6
          ipi=nblendta(igrd)

          CALL iom_get ( numbdyv_bt, jpdom_unknown,'vobtcrty',zdta(1:ipi,1:ipj),nbdy_a_bt )

          DO ib=1, nblen(igrd)
            vbtbdydta(ib,2) =  zdta(nbmap(ib,igrd),1)
          END DO


         IF(lwp) WRITE(numout,*) 'bdy_dta_fla : first record file used nbdy_b_bt ',nbdy_b_bt
         IF(lwp) WRITE(numout,*) '~~~~~~~~  last  record file used nbdy_a_bt ',nbdy_a_bt
         IF (.NOT.ln_clim) THEN
           IF(lwp) WRITE(numout,*) 'first  record time (s): ', istep_bt(nbdy_b_bt)
           IF(lwp) WRITE(numout,*) 'model time (s)        : ', itimer
           IF(lwp) WRITE(numout,*) 'second record time (s): ', istep_bt(nbdy_a_bt)
         ENDIF
        END IF ! end lect=.true.


      ! 2.2   Interpolate linearly:
      ! ***************************
    
        IF (ln_clim) THEN
          zxy = REAL( nday, wp ) / REAL( nmonth_len(nbdy_b_bt), wp ) + 0.5 - i15
        ELSE          
          zxy = REAL(istep_bt(nbdy_b_bt)-itimer, wp) / REAL(istep_bt(nbdy_b_bt)-istep_bt(nbdy_a_bt), wp)
        END IF

          igrd=4
          DO ib=1, nblen(igrd)
            sshbdy(ib) = zxy      * sshbdydta(ib,2) + &
                       (1.-zxy) * sshbdydta(ib,1)   
          END DO

          igrd=5
          DO ib=1, nblen(igrd)
            ubtbdy(ib) = zxy      * ubtbdydta(ib,2) + &
                         (1.-zxy) * ubtbdydta(ib,1)   
          END DO

          igrd=6
          DO ib=1, nblen(igrd)
            vbtbdy(ib) = zxy      * vbtbdydta(ib,2) + &
                         (1.-zxy) * vbtbdydta(ib,1)   
          END DO


      END IF !end if ((nn_dtactl==1).AND.(ntimes_bdy_bt>1))
    
      ! ------------------- !
      ! Last call kt=nitend !
      ! ------------------- !

      ! Closing of the 3 files
      IF( kt == nitend   .and. jit == icycl ) THEN
          CALL iom_close( numbdyt_bt )
          CALL iom_close( numbdyu_bt )
          CALL iom_close( numbdyv_bt )
      ENDIF

      ENDIF ! ln_dyn_frs

      END SUBROUTINE bdy_dta_fla


#else
   !!----------------------------------------------------------------------
   !!   Dummy module                   NO Unstruct Open Boundary Conditions
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE bdy_dta_frs( kt )              ! Empty routine
      WRITE(*,*) 'bdy_dta_frs: You should not have seen this print! error?', kt
   END SUBROUTINE bdy_dta_frs
   SUBROUTINE bdy_dta_fla( kt, kit, icycle )      ! Empty routine
      WRITE(*,*) 'bdy_dta_frs: You should not have seen this print! error?', kt, kit
   END SUBROUTINE bdy_dta_fla
#endif

   !!==============================================================================
END MODULE bdydta
