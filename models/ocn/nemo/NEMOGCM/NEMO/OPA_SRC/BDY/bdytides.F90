MODULE bdytides
   !!======================================================================
   !!                       ***  MODULE  bdytides  ***
   !! Ocean dynamics:   Tidal forcing at open boundaries
   !!======================================================================
   !! History :  2.0  !  2007-01  (D.Storkey)  Original code
   !!            2.3  !  2008-01  (J.Holt)  Add date correction. Origins POLCOMS v6.3 2007
   !!            3.0  !  2008-04  (NEMO team)  add in the reference version
   !!            3.3  !  2010-09  (D.Storkey and E.O'Dea)  bug fixes
   !!----------------------------------------------------------------------
#if defined key_bdy
   !!----------------------------------------------------------------------
   !!   'key_bdy'     Unstructured Open Boundary Condition
   !!----------------------------------------------------------------------
   !!   PUBLIC
   !!      tide_init     : read of namelist
   !!      tide_data     : read in and initialisation of tidal constituents at boundary
   !!      tide_update   : calculation of tidal forcing at each timestep
   !!   PRIVATE
   !!      uvset         :\
   !!      vday          : | Routines to correct tidal harmonics forcing for
   !!      shpen         : | start time of integration
   !!      ufset         : |
   !!      vset          :/
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers 
   USE dom_oce         ! ocean space and time domain
   USE iom
   USE in_out_manager  ! I/O units
   USE phycst          ! physical constants
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE bdy_par         ! Unstructured boundary parameters
   USE bdy_oce         ! ocean open boundary conditions
   USE daymod          ! calendar

   IMPLICIT NONE
   PRIVATE

   PUBLIC   tide_init     ! routine called in bdyini
   PUBLIC   tide_data     ! routine called in bdyini
   PUBLIC   tide_update   ! routine called in bdydyn

   LOGICAL, PUBLIC            ::   ln_tide_date          !: =T correct tide phases and amplitude for model start date
   INTEGER, PUBLIC, PARAMETER ::   jptides_max = 15      !: Max number of tidal contituents
   INTEGER, PUBLIC            ::   ntide                 !: Actual number of tidal constituents

   CHARACTER(len=80), PUBLIC                         ::   filtide    !: Filename root for tidal input files
   CHARACTER(len= 4), PUBLIC, DIMENSION(jptides_max) ::   tide_cpt   !: Names of tidal components used.

   INTEGER , PUBLIC, DIMENSION(jptides_max) ::   nindx        !: ???
   REAL(wp), PUBLIC, DIMENSION(jptides_max) ::   tide_speed   !: Phase speed of tidal constituent (deg/hr)
   
   REAL(wp), DIMENSION(jpbdim,jptides_max)  ::   ssh1, ssh2   ! Tidal constituents : SSH
   REAL(wp), DIMENSION(jpbdim,jptides_max)  ::   u1  , u2     ! Tidal constituents : U
   REAL(wp), DIMENSION(jpbdim,jptides_max)  ::   v1  , v2     ! Tidal constituents : V
   
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: bdytides.F90 2528 2010-12-27 17:33:53Z rblod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE tide_init
      !!----------------------------------------------------------------------
      !!                    ***  SUBROUTINE tide_init  ***
      !!                     
      !! ** Purpose : - Read in namelist for tides
      !!
      !!----------------------------------------------------------------------
      INTEGER ::   itide                  ! dummy loop index 
      !!
      NAMELIST/nambdy_tide/ln_tide_date, filtide, tide_cpt, tide_speed
      !!----------------------------------------------------------------------

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'tide_init : initialization of tidal harmonic forcing at open boundaries'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~'

      ! Namelist nambdy_tide : tidal harmonic forcing at open boundaries
      ln_tide_date = .false.
      filtide(:) = ''
      tide_speed(:) = 0.0
      tide_cpt(:) = ''
      REWIND( numnam )                                ! Read namelist parameters
      READ  ( numnam, nambdy_tide )
      !                                               ! Count number of components specified
      ntide = jptides_max
      DO itide = 1, jptides_max
        IF( tide_cpt(itide) == '' ) THEN
           ntide = itide-1
           exit
        ENDIF
      END DO

      !                                               ! find constituents in standard list
      DO itide = 1, ntide
         nindx(itide) = 0
         IF( TRIM( tide_cpt(itide) ) == 'Q1'  )   nindx(itide) =  1
         IF( TRIM( tide_cpt(itide) ) == 'O1'  )   nindx(itide) =  2
         IF( TRIM( tide_cpt(itide) ) == 'P1'  )   nindx(itide) =  3
         IF( TRIM( tide_cpt(itide) ) == 'S1'  )   nindx(itide) =  4
         IF( TRIM( tide_cpt(itide) ) == 'K1'  )   nindx(itide) =  5
         IF( TRIM( tide_cpt(itide) ) == '2N2' )   nindx(itide) =  6
         IF( TRIM( tide_cpt(itide) ) == 'MU2' )   nindx(itide) =  7
         IF( TRIM( tide_cpt(itide) ) == 'N2'  )   nindx(itide) =  8
         IF( TRIM( tide_cpt(itide) ) == 'NU2' )   nindx(itide) =  9
         IF( TRIM( tide_cpt(itide) ) == 'M2'  )   nindx(itide) = 10
         IF( TRIM( tide_cpt(itide) ) == 'L2'  )   nindx(itide) = 11
         IF( TRIM( tide_cpt(itide) ) == 'T2'  )   nindx(itide) = 12
         IF( TRIM( tide_cpt(itide) ) == 'S2'  )   nindx(itide) = 13
         IF( TRIM( tide_cpt(itide) ) == 'K2'  )   nindx(itide) = 14
         IF( TRIM( tide_cpt(itide) ) == 'M4'  )   nindx(itide) = 15
         IF( nindx(itide) == 0  .AND. lwp ) THEN
            WRITE(ctmp1,*) 'constitunent', itide,':', tide_cpt(itide), 'not in standard list'
            CALL ctl_warn( ctmp1 )
         ENDIF
      END DO
      !                                               ! Parameter control and print
      IF( ntide < 1 ) THEN 
         CALL ctl_stop( '          Did not find any tidal components in namelist nambdy_tide' )
      ELSE
         IF(lwp) WRITE(numout,*) '          Namelist nambdy_tide : tidal harmonic forcing at open boundaries'
         IF(lwp) WRITE(numout,*) '             tidal components specified ', ntide
         IF(lwp) WRITE(numout,*) '                ', tide_cpt(1:ntide)
         IF(lwp) WRITE(numout,*) '             associated phase speeds (deg/hr) : '
         IF(lwp) WRITE(numout,*) '                ', tide_speed(1:ntide)
      ENDIF

      ! Initialisation of tidal harmonics arrays
      sshtide(:) = 0.e0
      utide  (:) = 0.e0
      vtide  (:) = 0.e0
      !
   END SUBROUTINE tide_init


   SUBROUTINE tide_data
      !!----------------------------------------------------------------------
      !!                    ***  SUBROUTINE tide_data  ***
      !!                    
      !! ** Purpose : - Read in tidal harmonics data and adjust for the start 
      !!                time of the model run. 
      !!
      !!----------------------------------------------------------------------
      INTEGER ::   itide, igrd, ib        ! dummy loop indices
      CHARACTER(len=80) :: clfile         ! full file name for tidal input file 
      INTEGER ::   ipi, ipj, inum, idvar  ! temporary integers (netcdf read)
      INTEGER, DIMENSION(6) :: lendta=0   ! length of data in the file (note may be different from nblendta!)
      REAL(wp) ::  z_arg, z_atde, z_btde, z1t, z2t           
      REAL(wp), DIMENSION(jpbdta,1) ::   zdta   ! temporary array for data fields
      REAL(wp), DIMENSION(jptides_max) :: z_vplu, z_ftc
      !!------------------------------------------------------------------------------

      ! Open files and read in tidal forcing data
      ! -----------------------------------------

      ipj = 1

      DO itide = 1, ntide
         !                                                              ! SSH fields
         clfile = TRIM(filtide)//TRIM(tide_cpt(itide))//'_grid_T.nc'
         IF(lwp) WRITE(numout,*) 'Reading data from file ', clfile
         CALL iom_open( clfile, inum )
         igrd = 4
         IF( nblendta(igrd) <= 0 ) THEN 
            idvar = iom_varid( inum,'z1' )
            IF(lwp) WRITE(numout,*) 'iom_file(1)%ndims(idvar) : ',iom_file%ndims(idvar)
            nblendta(igrd) = iom_file(inum)%dimsz(1,idvar)
            WRITE(numout,*) 'Dim size for z1 is ', nblendta(igrd)
         ENDIF
         ipi = nblendta(igrd)
         CALL iom_get( inum, jpdom_unknown, 'z1', zdta(1:ipi,1:ipj) )
         DO ib = 1, nblenrim(igrd)
            ssh1(ib,itide) = zdta(nbmap(ib,igrd),1)
         END DO
         CALL iom_get( inum, jpdom_unknown, 'z2', zdta(1:ipi,1:ipj) )
         DO ib = 1, nblenrim(igrd)
            ssh2(ib,itide) = zdta(nbmap(ib,igrd),1)
         END DO
         CALL iom_close( inum )
         !
         !                                                              ! U fields
         clfile = TRIM(filtide)//TRIM(tide_cpt(itide))//'_grid_U.nc'
         IF(lwp) WRITE(numout,*) 'Reading data from file ', clfile
         CALL iom_open( clfile, inum )
         igrd = 5
         IF( lendta(igrd) <= 0 ) THEN 
            idvar = iom_varid( inum,'u1' )
            lendta(igrd) = iom_file(inum)%dimsz(1,idvar)
            WRITE(numout,*) 'Dim size for u1 is ',lendta(igrd)
         ENDIF
         ipi = lendta(igrd)
         CALL iom_get( inum, jpdom_unknown, 'u1', zdta(1:ipi,1:ipj) )
         DO ib = 1, nblenrim(igrd)
            u1(ib,itide) = zdta(nbmap(ib,igrd),1)
         END DO
         CALL iom_get( inum, jpdom_unknown, 'u2', zdta(1:ipi,1:ipj) )
         DO ib = 1, nblenrim(igrd)
            u2(ib,itide) = zdta(nbmap(ib,igrd),1)
         END DO
         CALL iom_close( inum )
         !
         !                                                              ! V fields
         clfile = TRIM(filtide)//TRIM(tide_cpt(itide))//'_grid_V.nc'
         if(lwp) write(numout,*) 'Reading data from file ', clfile
         CALL iom_open( clfile, inum )
         igrd = 6
         IF( lendta(igrd) <= 0 ) THEN 
            idvar = iom_varid( inum,'v1' )
            lendta(igrd) = iom_file(inum)%dimsz(1,idvar)
            WRITE(numout,*) 'Dim size for v1 is ', lendta(igrd)
         ENDIF
         ipi = lendta(igrd)
         CALL iom_get( inum, jpdom_unknown, 'v1', zdta(1:ipi,1:ipj) )
         DO ib = 1, nblenrim(igrd)
            v1(ib,itide) = zdta(nbmap(ib,igrd),1)
         END DO
         CALL iom_get( inum, jpdom_unknown, 'v2', zdta(1:ipi,1:ipj) )
         DO ib=1, nblenrim(igrd)
            v2(ib,itide) = zdta(nbmap(ib,igrd),1)
         END DO
         CALL iom_close( inum )
         !
      END DO ! end loop on tidal components

      IF( ln_tide_date ) THEN      ! correct for date factors

!! used nmonth, nyear and nday from daymod....
         ! Calculate date corrects for 15 standard consituents
         ! This is the initialisation step, so nday, nmonth, nyear are the 
         ! initial date/time of the integration.
           print *, nday,nmonth,nyear
           nyear  = int(ndate0 / 10000  )                           ! initial year
           nmonth = int((ndate0 - nyear * 10000 ) / 100 )          ! initial month
           nday   = int(ndate0 - nyear * 10000 - nmonth * 100)

         CALL uvset( 0, nday, nmonth, nyear, z_ftc, z_vplu )

         IF(lwp) WRITE(numout,*) 'Correcting tide for date:', nday, nmonth, nyear

         DO itide = 1, ntide       ! loop on tidal components
            !
            IF( nindx(itide) /= 0 ) THEN
!!gm use rpi  and rad global variable  
               z_arg = 3.14159265d0 * z_vplu(nindx(itide)) / 180.0d0
               z_atde=z_ftc(nindx(itide))*cos(z_arg)
               z_btde=z_ftc(nindx(itide))*sin(z_arg)
               IF(lwp) WRITE(numout,'(2i5,8f10.6)') itide, nindx(itide), tide_speed(itide),   &
                  &                                 z_ftc(nindx(itide)), z_vplu(nindx(itide))
            ELSE
               z_atde = 1.0_wp
               z_btde = 0.0_wp
            ENDIF
            !                                         !  elevation         
            igrd = 4
            DO ib = 1, nblenrim(igrd)                
               z1t = z_atde * ssh1(ib,itide) + z_btde * ssh2(ib,itide)
               z2t = z_atde * ssh2(ib,itide) - z_btde * ssh1(ib,itide)
               ssh1(ib,itide) = z1t
               ssh2(ib,itide) = z2t
            END DO
            !                                         !  u       
            igrd = 5
            DO ib = 1, nblenrim(igrd)                
               z1t = z_atde * u1(ib,itide) + z_btde * u2(ib,itide)
               z2t = z_atde * u2(ib,itide) - z_btde * u1(ib,itide)
               u1(ib,itide) = z1t
               u2(ib,itide) = z2t
            END DO
            !                                         !  v       
            igrd = 6
            DO ib = 1, nblenrim(igrd)                
               z1t = z_atde * v1(ib,itide) + z_btde * v2(ib,itide)
               z2t = z_atde * v2(ib,itide) - z_btde * v1(ib,itide)
               v1(ib,itide) = z1t
               v2(ib,itide) = z2t
            END DO
            !
         END DO   ! end loop on tidal components
         !
      ENDIF ! date correction
      !
   END SUBROUTINE tide_data


   SUBROUTINE tide_update ( kt, jit )
      !!----------------------------------------------------------------------
      !!                 ***  SUBROUTINE tide_update  ***
      !!                
      !! ** Purpose : - Add tidal forcing to sshbdy, ubtbdy and vbtbdy arrays. 
      !!                
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt      ! Main timestep counter
!!gm doctor jit ==> kit
      INTEGER, INTENT( in ) ::   jit     ! Barotropic timestep counter (for timesplitting option)
      !!
      INTEGER  ::   itide, igrd, ib      ! dummy loop indices
      REAL(wp) ::   z_arg, z_sarg            !            
      REAL(wp), DIMENSION(jptides_max) :: z_sist, z_cost
      !!----------------------------------------------------------------------

      ! Note tide phase speeds are in deg/hour, so we need to convert the
      ! elapsed time in seconds to hours by dividing by 3600.0
      IF( jit == 0 ) THEN  
         z_arg = kt * rdt * rad / 3600.0
      ELSE                              ! we are in a barotropic subcycle (for timesplitting option)
!         z_arg = ( (kt-1) * rdt + jit * rdt / REAL(nn_baro,lwp) ) * rad / 3600.0
         z_arg = ( (kt-1) * rdt + jit * rdt / REAL(nn_baro,wp) ) * rad / 3600.0
      ENDIF

      DO itide = 1, ntide
         z_sarg = z_arg * tide_speed(itide)
         z_cost(itide) = COS( z_sarg )
         z_sist(itide) = SIN( z_sarg )
      END DO

      ! summing of tidal constituents into BDY arrays
      sshtide(:) = 0.0
      utide (:) = 0.0
      vtide (:) = 0.0
      !
      DO itide = 1, ntide
         igrd=4                              ! SSH on tracer grid.
         DO ib = 1, nblenrim(igrd)
            sshtide(ib) =sshtide(ib)+ ssh1(ib,itide)*z_cost(itide) + ssh2(ib,itide)*z_sist(itide)
            !    if(lwp) write(numout,*) 'z',ib,itide,sshtide(ib), ssh1(ib,itide),ssh2(ib,itide)
         END DO
         igrd=5                              ! U grid
         DO ib=1, nblenrim(igrd)
            utide(ib) = utide(ib)+ u1(ib,itide)*z_cost(itide) + u2(ib,itide)*z_sist(itide)
            !    if(lwp) write(numout,*) 'u',ib,itide,utide(ib), u1(ib,itide),u2(ib,itide)
         END DO
         igrd=6                              ! V grid
         DO ib=1, nblenrim(igrd)
            vtide(ib) = vtide(ib)+ v1(ib,itide)*z_cost(itide) + v2(ib,itide)*z_sist(itide)
            !    if(lwp) write(numout,*) 'v',ib,itide,vtide(ib), v1(ib,itide),v2(ib,itide)
         END DO
      END DO
      !
   END SUBROUTINE tide_update

!!gm  doctor naming of dummy argument variables!!!   and all local variables

   SUBROUTINE uvset( ihs, iday, imnth, iyr, f, z_vplu )
      !!----------------------------------------------------------------------
      !!                   ***  SUBROUTINE uvset  ***
      !!                     
      !! ** Purpose : - adjust tidal forcing for date factors
      !!                
      !!----------------------------------------------------------------------
      implicit none
      INTEGER, INTENT( in ) ::   ihs     ! Start time hours
      INTEGER, INTENT( in ) ::   iday    ! start time days
      INTEGER, INTENT( in ) ::   imnth   ! start time month
      INTEGER, INTENT( in ) ::   iyr     ! start time year
      !!
!!gm  nc is jptides_max ????
      INTEGER , PARAMETER     ::  nc =15    ! maximum number of constituents
      CHARACTER(len=8), DIMENSION(nc) :: cname
      INTEGER                 ::   year, vd, ivdy, ndc, i, k
      REAL(wp)                ::   ss, h, p, en, p1, rtd
      REAL(wp), DIMENSION(nc) ::   f                          ! nodal correction
      REAL(wp), DIMENSION(nc) ::   z_vplu                     ! phase correction
      REAL(wp), DIMENSION(nc) ::   u, v, zig
      !!
      DATA cname/  'q1'    ,    'o1'    ,     'p1'   ,    's1'    ,     'k1'    ,   &
         &         '2n2'   ,    'mu2'   ,     'n2'   ,    'nu2'   ,     'm2'    ,   &
         &         'l2'    ,    't2'    ,     's2'   ,    'k2'    ,     'm4'      /
      DATA zig/ .2338507481, .2433518789, .2610826055, .2617993878,  .2625161701,   &
         &      .4868657873, .4881373225, .4963669182, .4976384533,  .5058680490,   &
         &      .5153691799, .5228820265, .5235987756, .5250323419, 1.011736098   /
      !!----------------------------------------------------------------------
!
!     ihs             -  start time gmt on ...
!     iday/imnth/iyr  -  date   e.g.   12/10/87
!
      CALL vday(iday,imnth,iyr,ivdy)
      vd = ivdy
!
!rp   note change of year number for d. blackman shpen
!rp   if(iyr.ge.1000) year=iyr-1900
!rp   if(iyr.lt.1000) year=iyr
      year = iyr
!
!.....year  =  year of required data
!.....vd    =  day of required data..set up for 0000gmt day year
      ndc = nc
!.....ndc   =  number of constituents allowed
!!gm use rpi ?
      rtd = 360.0 / 6.2831852
      DO i = 1, ndc
         zig(i) = zig(i)*rtd
         ! sigo(i)= zig(i)
      END DO

!!gm try to avoid RETURN  in F90
      IF( year == 0 )   RETURN
      
         CALL shpen( year, vd, ss, h , p , en, p1 )
         CALL ufset( p   , en, u , f )
         CALL vset ( ss  , h , p , en, p1, v )
         !
         DO k = 1, nc
            z_vplu(k) = v(k) + u(k)
            z_vplu(k) = z_vplu(k) + dble(ihs) * zig(k)
            DO WHILE( z_vplu(k) < 0    )
               z_vplu(k) = z_vplu(k) + 360.0
            END DO
            DO WHILE( z_vplu(k) > 360. )
               z_vplu(k) = z_vplu(k) - 360.0
            END DO
         END DO
         !
      END SUBROUTINE uvset


      SUBROUTINE vday( iday, imnth, iy, ivdy )
      !!----------------------------------------------------------------------
      !!                   *** SUBROUTINE vday  ***
      !!                  
      !! ** Purpose : - adjust tidal forcing for date factors
      !!                
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in   ) ::   iday, imnth, iy   ! ????
      INTEGER, INTENT(  out) ::   ivdy              ! ???
      !! 
      INTEGER ::   iyr
      !!------------------------------------------------------------------------------

!!gm   nday_year in day mode is the variable compiuted here, no?
!!gm      nday_year ,   &  !: curent day counted from jan 1st of the current year

      !calculate day number in year from day/month/year
       if(imnth.eq.1) ivdy=iday
       if(imnth.eq.2) ivdy=iday+31
       if(imnth.eq.3) ivdy=iday+59
       if(imnth.eq.4) ivdy=iday+90
       if(imnth.eq.5) ivdy=iday+120
       if(imnth.eq.6) ivdy=iday+151
       if(imnth.eq.7) ivdy=iday+181
       if(imnth.eq.8) ivdy=iday+212
       if(imnth.eq.9) ivdy=iday+243
       if(imnth.eq.10) ivdy=iday+273
       if(imnth.eq.11) ivdy=iday+304
       if(imnth.eq.12) ivdy=iday+334
       iyr=iy
       if(mod(iyr,4).eq.0.and.imnth.gt.2) ivdy=ivdy+1
       if(mod(iyr,100).eq.0.and.imnth.gt.2) ivdy=ivdy-1
       if(mod(iyr,400).eq.0.and.imnth.gt.2) ivdy=ivdy+1
      !
   END SUBROUTINE vday

!!doctor norme for dummy arguments

   SUBROUTINE shpen( year, vd, s, h, p, en, p1 )
      !!----------------------------------------------------------------------
      !!                    ***  SUBROUTINE shpen  ***
      !!                   
      !! ** Purpose : - calculate astronomical arguments for tides
      !!                this version from d. blackman 30 nove 1990
      !!
      !!----------------------------------------------------------------------
!!gm add INTENT in, out or inout....    DOCTOR name....
!!gm please do not use variable name with a single letter:  impossible to search in a code
      INTEGER  ::   year,vd
      REAL(wp) ::   s,h,p,en,p1
      !!
      INTEGER  ::   yr,ilc,icent,it,iday,ild,ipos,nn,iyd
      REAL(wp) ::   cycle,t,td,delt(84),delta,deltat
      !!
      DATA delt /-5.04, -3.90, -2.87, -0.58,  0.71,  1.80,           &
         &        3.08,  4.63,  5.86,  7.21,  8.58, 10.50, 12.10,    &
         &       12.49, 14.41, 15.59, 15.81, 17.52, 19.01, 18.39,    &
         &       19.55, 20.36, 21.01, 21.81, 21.76, 22.35, 22.68,    &
         &       22.94, 22.93, 22.69, 22.94, 23.20, 23.31, 23.63,    &
         &       23.47, 23.68, 23.62, 23.53, 23.59, 23.99, 23.80,    &
         &       24.20, 24.99, 24.97, 25.72, 26.21, 26.37, 26.89,    &
         &       27.68, 28.13, 28.94, 29.42, 29.66, 30.29, 30.96,    &
         &       31.09, 31.59, 31.52, 31.92, 32.45, 32.91, 33.39,    &
         &       33.80, 34.23, 34.73, 35.40, 36.14, 36.99, 37.87,    &
         &       38.75, 39.70, 40.70, 41.68, 42.82, 43.96, 45.00,    &
         &       45.98, 47.00, 48.03, 49.10, 50.10, 50.97, 51.81,    &
         &       52.57 /
      !!----------------------------------------------------------------------

      cycle = 360.0
      ilc = 0
      icent = year / 100
      yr = year - icent * 100
      t = icent - 20
!
!     for the following equations
!     time origin is fixed at 00 hr of jan 1st,2000.
!     see notes by cartwright
!
!!gm  old coding style, use CASE instead  and avoid GOTO (obsolescence in fortran 90)
!!gm obsol(   1): Arithmetic IF statement is used   ===>  remove this in Fortran 90
      it = icent - 20
      if (it) 1,2,2
    1 iday = it/4 -it
      go to 3
    2 iday = (it+3)/4 - it
!
!     t is in julian century
!     correction in gegorian calander where only century year divisible
!     by 4 is leap year.
!
    3 continue
!
      td = 0.0
!
!!gm obsol(   1): Arithmetic IF statement is used   ===>  remove this in Fortran 90
      if (yr) 4,5,4
!
    4 iyd = 365*yr
      ild = (yr-1)/4
      if((icent - (icent/4)*4) .eq. 0) ilc = 1
      td = iyd + ild + ilc
!
    5 td = td + iday + vd -1.0 - 0.5
      t = t + (td/36525.0)
!
      ipos=year-1899
      if (ipos .lt. 0) go to 7
      if (ipos .gt. 83) go to 6
!
      delta = (delt(ipos+1)+delt(ipos))/2.0
      go to 7
!
    6 delta= (65.0-50.5)/20.0*(year-1980)+50.5
!
    7 deltat = delta * 1.0e-6
!
!!gm   precision of the computation   :  example for s it should be replace by:
!!gm  s   = 218.3165 + (481267.8813 - 0.0016*t)*t + 152.0*deltat   ==>  more precise  modify the last digits results
      s   = 218.3165 + 481267.8813*t - 0.0016*t*t + 152.0*deltat
      h   = 280.4661 + 36000.7698 *t + 0.0003*t*t +  11.0*deltat
      p   =  83.3535 + 4069.0139  *t - 0.0103*t*t +       deltat
      en  = 234.9555 + 1934.1363  *t - 0.0021*t*t +       deltat
      p1  = 282.9384 + 1.7195     *t + 0.0005*t*t
      !
      nn = s / cycle
      s = s - nn * cycle
      IF( s < 0.e0 )   s = s + cycle
      !
      nn = h / cycle
      h = h - cycle * nn
      IF( h < 0.e0 )   h = h + cycle
      !
      nn = p / cycle
      p = p - cycle * nn
      IF( p < 0.e0)   p = p + cycle
      !
      nn = en / cycle
      en = en - cycle * nn
      IF( en < 0.e0 )   en = en + cycle
      en = cycle - en
      !
      nn = p1 / cycle
      p1 = p1 - nn * cycle
      !
   END SUBROUTINE shpen


   SUBROUTINE ufset( p, cn, b, a )
      !!----------------------------------------------------------------------
      !!                    ***  SUBROUTINE ufset  ***
      !!                    
      !! ** Purpose : - calculate nodal parameters for the tides
      !!                
      !!----------------------------------------------------------------------
!!gm  doctor naming of dummy argument variables!!!   and all local variables
!!gm  nc is jptides_max ????
      integer nc
      parameter (nc=15)
      REAL(wp) p,cn
      !!
      
!!gm  rad is already a public variable defined in phycst.F90 ....   ==> doctor norme  local real start with "z"
      REAL(wp) ::   w1, w2, w3, w4, w5, w6, w7, w8, nw, pw, rad
      REAL(wp) ::   a(nc), b(nc)
      INTEGER  ::   ndc, k
      !!----------------------------------------------------------------------

      ndc = nc

!    a=f       ,  b =u
!    t is  zero as compared to tifa.
!! use rad defined in phycst   (i.e.  add a USE phycst at the begining of the module
      rad = 6.2831852d0/360.0
      pw = p  * rad
      nw = cn * rad
      w1 = cos(   nw )
      w2 = cos( 2*nw )
      w3 = cos( 3*nw )
      w4 = sin(   nw )
      w5 = sin( 2*nw )
      w6 = sin( 3*nw )
      w7 = 1. - 0.2505 * COS( 2*pw      ) - 0.1102 * COS(2*pw-nw )   &
         &    - 0.156  * COS( 2*pw-2*nw ) - 0.037  * COS(     nw )
      w8 =    - 0.2505 * SIN( 2*pw      ) - 0.1102 * SIN(2*pw-nw )   &
         &    - 0.0156 * SIN( 2*pw-2*nw ) - 0.037  * SIN(     nw )
      !
      a(1) = 1.0089 + 0.1871 * w1 - 0.0147 * w2 + 0.0014 * w3
      b(1) =          0.1885 * w4 - 0.0234 * w5 + 0.0033 * w6
      !   q1
      a(2) = a(1)
      b(2) = b(1)
      !   o1
      a(3) = 1.0
      b(3) = 0.0
      !   p1
      a(4) = 1.0
      b(4) = 0.0
      !   s1
      a(5) = 1.0060+0.1150*w1- 0.0088*w2 +0.0006*w3
      b(5) = -0.1546*w4 + 0.0119*w5 -0.0012*w6
      !   k1
      a(6) =1.0004 -0.0373*w1+ 0.0002*w2
      b(6) = -0.0374*w4
      !  2n2
      a(7) = a(6)
      b(7) = b(6)
      !  mu2
      a(8) = a(6)
      b(8) = b(6)
      !   n2
      a(9) = a(6)
      b(9) = b(6)
      !  nu2
      a(10) = a(6)
      b(10) = b(6)
      !   m2
      a(11) = SQRT( w7 * w7 + w8 * w8 )
      b(11) = ATAN( w8 / w7 )
!!gmuse rpi  instead of 3.141992 ???   true pi is rpi=3.141592653589793_wp  .....   ????
      IF( w7 < 0.e0 )   b(11) = b(11) + 3.141992
      !   l2
      a(12) = 1.0
      b(12) = 0.0
      !   t2
      a(13)= a(12)
      b(13)= b(12)
      !   s2
      a(14) = 1.0241+0.2863*w1+0.0083*w2 -0.0015*w3
      b(14) = -0.3096*w4 + 0.0119*w5 - 0.0007*w6
      !   k2
      a(15) = a(6)*a(6)
      b(15) = 2*b(6)
      !   m4
!!gm  old coding,  remove GOTO and label of lines
!!gm obsol(   1): Arithmetic IF statement is used   ===>  remove this in Fortran 90
      DO 40 k = 1,ndc
         b(k) = b(k)/rad
32       if (b(k)) 34,35,35
34       b(k) = b(k) + 360.0
         go to 32
35       if (b(k)-360.0) 40,37,37
37       b(k) = b(k)-360.0
         go to 35
40    continue
      !
   END SUBROUTINE ufset
   

   SUBROUTINE vset( s,h,p,en,p1,v)
      !!----------------------------------------------------------------------
      !!                    ***  SUBROUTINE vset  ***
      !!                   
      !! ** Purpose : - calculate tidal phases for 0000gmt on start day of run
      !!                
      !!----------------------------------------------------------------------
!!gm  doctor naming of dummy argument variables!!!   and all local variables
!!gm  nc is jptides_max ????
!!gm  en argument is not used: suppress it ?
      integer nc
      parameter (nc=15)
      real(wp) s,h,p,en,p1
      real(wp)   v(nc)
      !!
      integer ndc, k
      !!----------------------------------------------------------------------

      ndc = nc
      !   v s  are computed here.
      v(1) =-3*s +h +p +270      ! Q1
      v(2) =-2*s +h +270.0       ! O1
      v(3) =-h +270              ! P1
      v(4) =180                  ! S1
      v(5) =h +90.0              ! K1
      v(6) =-4*s +2*h +2*p       ! 2N2
      v(7) =-4*(s-h)             ! MU2
      v(8) =-3*s +2*h +p         ! N2
      v(9) =-3*s +4*h -p         ! MU2
      v(10) =-2*s +2*h           ! M2
      v(11) =-s +2*h -p +180     ! L2 
      v(12) =-h +p1              ! T2
      v(13) =0                   ! S2
      v(14) =h+h                 ! K2
      v(15) =2*v(10)             ! M4
!
!!gm  old coding,  remove GOTO and label of lines
!!gm obsol(   1): Arithmetic IF statement is used   ===>  remove this in Fortran 90
      do 72 k = 1, ndc
69       if( v(k) )   70,71,71
70       v(k) = v(k)+360.0
         go to 69
71       if( v(k) - 360.0 )   72,73,73
73       v(k) = v(k)-360.0
         go to 71
72    continue
      !
   END SUBROUTINE vset

#else
   !!----------------------------------------------------------------------
   !!   Dummy module         NO Unstruct Open Boundary Conditions for tides
   !!----------------------------------------------------------------------
!!gm  are you sure we need to define filtide and tide_cpt ?
   CHARACTER(len=80), PUBLIC               ::   filtide                !: Filename root for tidal input files
   CHARACTER(len=4 ), PUBLIC, DIMENSION(1) ::   tide_cpt               !: Names of tidal components used.

CONTAINS
   SUBROUTINE tide_init                ! Empty routine
   END SUBROUTINE tide_init
   SUBROUTINE tide_data                ! Empty routine
   END SUBROUTINE tide_data
   SUBROUTINE tide_update( kt, kit )   ! Empty routine
      WRITE(*,*) 'tide_update: You should not have seen this print! error?', kt, kit
   END SUBROUTINE tide_update
#endif

   !!======================================================================
END MODULE bdytides
