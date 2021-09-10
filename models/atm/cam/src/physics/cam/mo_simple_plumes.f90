!>
!!
!! @brief Module MO_SIMPLE_PLUMES: provides anthropogenic aerorol optial properties as a function of lat, lon
!!   height, time, and wavelength
!!
!! @remarks
!!
!! @author Bjorn Stevens & Karsten Peters MPI-M, Hamburg (v1-beta release 2015-12-05)
!!
!! $ID: n/a$
!!
!! @par Origin
!!   Based on code originally developed at the MPI by Karsten Peters, Bjorn Stevens, Stephanie Fiedler
!!   and Stefan Kinne with input from Thorsten Mauritsen and Robert Pincus
!!
!! @par Copyright
!! 
!
MODULE MO_SIMPLE_PLUMES
!EB+
use shr_kind_mod     , only: r8 => shr_kind_r8
use shr_const_mod, only: pi => shr_const_pi
use abortutils,      only: endrun
use cam_logfile,     only: iulog
!EB-
  USE netcdf

  IMPLICIT NONE

  INTEGER, PARAMETER ::                        &
       nplumes   = 9                          ,& !< Number of plumes
       nfeatures = 2                          ,& !< Number of features per plume
       ntimes    = 52                         ,& !< Number of times resolved per year (52 => weekly resolution)
       nyears    = 251                           !< Number of years of available forcing

  REAL(r8), PARAMETER    ::                        &
!EB       pi      = 2.*ASIN(1.),                  & !< half the unit circle (radians)
!EB       deg2rad = pi/180.                         !< conversion factor from degrees to radians
       deg2rad = pi/180._r8                         !< conversion factor from degrees to radians

  LOGICAL, SAVE ::                             &
       sp_initialized = .FALSE.                  !< parameter determining whether input needs to be read

  REAL(r8) ::                                      &
       plume_lat   (nplumes)                  ,& !< latitude where plume maximizes
       plume_lon   (nplumes)                  ,& !< longitude where plume maximizes
       beta_a      (nplumes)                  ,& !< parameter a for beta function vertical profile
       beta_b      (nplumes)                  ,& !< parameter b for beta function vertical profile
       aod_spmx    (nplumes)                  ,& !< aod at 550 for simple plume (maximum)
       aod_fmbg    (nplumes)                  ,& !< aod at 550 for fine mode background (for twomey effect)
       asy550      (nplumes)                  ,& !< asymmetry parameter for plume at 550nm
       ssa550      (nplumes)                  ,& !< single scattering albedo for plume at 550nm
       angstrom    (nplumes)                  ,& !< angstrom parameter for plume 
       sig_lon_E   (nfeatures,nplumes)        ,& !< Eastward extent of plume feature
       sig_lon_W   (nfeatures,nplumes)        ,& !< Westward extent of plume feature
       sig_lat_E   (nfeatures,nplumes)        ,& !< Southward extent of plume feature
       sig_lat_W   (nfeatures,nplumes)        ,& !< Northward extent of plume feature
       theta       (nfeatures,nplumes)        ,& !< Rotation angle of feature
       ftr_weight  (nfeatures,nplumes)        ,& !< Feature weights = (nfeatures + 1) to account for BB background
       time_weight (nfeatures,nplumes)        ,& !< Time-weights = (nfeatures +1) to account for BB background
       year_weight (nyears,nplumes)           ,& !< Yearly weight for plume
       ann_cycle   (nfeatures,ntimes,nplumes)    !< annual cycle for feature

  PUBLIC sp_aop_profile

CONTAINS
  !
  ! ------------------------------------------------------------------------------------------------------------------------
  ! SP_SETUP:  This subroutine should be called at initialization to read the netcdf data that describes the simple plume
  ! climatology.  The information needs to be either read by each processor or distributed to processors.
  !
  SUBROUTINE sp_setup
    !
    ! ---------- 
    !
    INTEGER :: iret, ncid, DimID, VarID, xdmy
    !
    ! ---------- 
    !    
!ED+ SSP5 RCP5 scenario instead of Historical
    iret = nf90_open("./MACv2.0-SP_v1-beta.nc", NF90_NOWRITE, ncid)
!ED- 
!EB+    IF (iret /= NF90_NOERR) STOP 'NetCDF File not opened'
    IF (iret /= NF90_NOERR) call endrun('NetCDF File not opened')
!EB-
    !
    ! read dimensions and make sure file conforms to expected size
    !
    iret = nf90_inq_dimid(ncid, "plume_number"  , DimId)
    iret = nf90_inquire_dimension(ncid, DimId, len = xdmy)
!EB+    IF (xdmy /= nplumes) STOP 'NetCDF improperly dimensioned -- plume_number'
    IF (xdmy /= nplumes) call endrun('NetCDF improperly dimensioned -- plume_number')
!EB-

    iret = nf90_inq_dimid(ncid, "plume_feature", DimId)
    iret = nf90_inquire_dimension(ncid, DimId, len = xdmy)
!EB+    IF (xdmy /= nfeatures) STOP 'NetCDF improperly dimensioned -- plume_feature'
    IF (xdmy /= nfeatures) call endrun('NetCDF improperly dimensioned -- plume_feature')
!EB-

    iret = nf90_inq_dimid(ncid, "year_fr"   , DimId)
    iret = nf90_inquire_dimension(ncid, DimID, len = xdmy)
!EB+    IF (xdmy /= ntimes) STOP 'NetCDF improperly dimensioned -- year_fr'
    IF (xdmy /= ntimes) call endrun('NetCDF improperly dimensioned -- year_fr')
!EB-
    iret = nf90_inq_dimid(ncid, "years"   , DimId)
    iret = nf90_inquire_dimension(ncid, DimID, len = xdmy)
!EB+    IF (xdmy /= nyears) STOP 'NetCDF improperly dimensioned -- years'
    IF (xdmy /= nyears) call endrun('NetCDF improperly dimensioned -- years')
!EB-
    !
    ! read variables that define the simple plume climatology
    !
    iret = nf90_inq_varid(ncid, "plume_lat", VarId)
    iret = nf90_get_var(ncid, VarID, plume_lat(:), start=(/1/),count=(/nplumes/))
!EB+    IF (iret /= NF90_NOERR) STOP 'NetCDF Error reading plume_lat'
    IF (iret /= NF90_NOERR) call endrun( 'NetCDF Error reading plume_lat')
    iret = nf90_inq_varid(ncid, "plume_lon", VarId)
    iret = nf90_get_var(ncid, VarID, plume_lon(:), start=(/1/),count=(/nplumes/))
!EB+    IF (iret /= NF90_NOERR) STOP 'NetCDF Error reading plume_lon'
    IF (iret /= NF90_NOERR) call endrun( 'NetCDF Error reading plume_lon')
    iret = nf90_inq_varid(ncid, "beta_a"   , VarId)
    iret = nf90_get_var(ncid, VarID, beta_a(:)   , start=(/1/),count=(/nplumes/))
!EB+    IF (iret /= NF90_NOERR) STOP 'NetCDF Error reading beta_a'
    IF (iret /= NF90_NOERR) call endrun( 'NetCDF Error reading beta_a')
    iret = nf90_inq_varid(ncid, "beta_b"   , VarId)
    iret = nf90_get_var(ncid, VarID, beta_b(:)   , start=(/1/),count=(/nplumes/))
!EB+    IF (iret /= NF90_NOERR) STOP 'NetCDF Error reading beta_b'
    IF (iret /= NF90_NOERR) call endrun( 'NetCDF Error reading beta_b')
    iret = nf90_inq_varid(ncid, "aod_spmx" , VarId)
    iret = nf90_get_var(ncid, VarID, aod_spmx(:)  , start=(/1/),count=(/nplumes/))
!EB+    IF (iret /= NF90_NOERR) STOP 'NetCDF Error reading aod_spmx'
    IF (iret /= NF90_NOERR) call endrun( 'NetCDF Error reading aod_spmx')
    iret = nf90_inq_varid(ncid, "aod_fmbg" , VarId)
    iret = nf90_get_var(ncid, VarID, aod_fmbg(:)  , start=(/1/),count=(/nplumes/))
!EB+    IF (iret /= NF90_NOERR) STOP 'NetCDF Error reading aod_fmbg'
    IF (iret /= NF90_NOERR) call endrun( 'NetCDF Error reading aod_fmbg')
    iret = nf90_inq_varid(ncid, "ssa550"   , VarId)
    iret = nf90_get_var(ncid, VarID, ssa550(:)  , start=(/1/),count=(/nplumes/))
!EB+    IF (iret /= NF90_NOERR) STOP 'NetCDF Error reading ssa550'
    IF (iret /= NF90_NOERR) call endrun( 'NetCDF Error reading ssa550')
    iret = nf90_inq_varid(ncid, "asy550"   , VarId)
    iret = nf90_get_var(ncid, VarID, asy550(:)  , start=(/1/),count=(/nplumes/))
!EB+    IF (iret /= NF90_NOERR) STOP 'NetCDF Error reading asy550'
    IF (iret /= NF90_NOERR) call endrun( 'NetCDF Error reading asy550')
    iret = nf90_inq_varid(ncid, "angstrom" , VarId)
    iret = nf90_get_var(ncid, VarID, angstrom(:), start=(/1/),count=(/nplumes/))
!EB+    IF (iret /= NF90_NOERR) STOP 'NetCDF Error reading anstrom'
    IF (iret /= NF90_NOERR) call endrun( 'NetCDF Error reading anstrom')

    iret = nf90_inq_varid(ncid, "sig_lat_W"     , VarId)
    iret = nf90_get_var(ncid, VarID, sig_lat_W(:,:)    , start=(/1,1/),count=(/nfeatures,nplumes/))
!EB+    IF (iret /= NF90_NOERR) STOP 'NetCDF Error reading sig_lat_W'
    IF (iret /= NF90_NOERR) call endrun( 'NetCDF Error reading sig_lat_W')
    iret = nf90_inq_varid(ncid, "sig_lat_E"     , VarId)
    iret = nf90_get_var(ncid, VarID, sig_lat_E(:,:)    , start=(/1,1/),count=(/nfeatures,nplumes/))
!EB+    IF (iret /= NF90_NOERR) STOP 'NetCDF Error reading sig_lat_E'
    IF (iret /= NF90_NOERR) call endrun( 'NetCDF Error reading sig_lat_E')
    iret = nf90_inq_varid(ncid, "sig_lon_E"     , VarId)
    iret = nf90_get_var(ncid, VarID, sig_lon_E(:,:)    , start=(/1,1/),count=(/nfeatures,nplumes/))
!EB+    IF (iret /= NF90_NOERR) STOP 'NetCDF Error reading sig_lon_E'
    IF (iret /= NF90_NOERR) call endrun( 'NetCDF Error reading sig_lon_E')
    iret = nf90_inq_varid(ncid, "sig_lon_W"     , VarId)
    iret = nf90_get_var(ncid, VarID, sig_lon_W(:,:)    , start=(/1,1/),count=(/nfeatures,nplumes/))
!EB+    IF (iret /= NF90_NOERR) STOP 'NetCDF Error reading sig_lon_W'
    IF (iret /= NF90_NOERR) call endrun( 'NetCDF Error reading sig_lon_W')
    iret = nf90_inq_varid(ncid, "theta"         , VarId)
    iret = nf90_get_var(ncid, VarID, theta(:,:)        , start=(/1,1/),count=(/nfeatures,nplumes/))
!EB+    IF (iret /= NF90_NOERR) STOP 'NetCDF Error reading theta'
    IF (iret /= NF90_NOERR) call endrun( 'NetCDF Error reading theta')
    iret = nf90_inq_varid(ncid, "ftr_weight"    , VarId)
    iret = nf90_get_var(ncid, VarID, ftr_weight(:,:)   , start=(/1,1/),count=(/nfeatures,nplumes/))
!EB+    IF (iret /= NF90_NOERR) STOP 'NetCDF Error reading plume_lat'
    IF (iret /= NF90_NOERR) call endrun( 'NetCDF Error reading plume_lat')
    iret = nf90_inq_varid(ncid, "year_weight"   , VarId)
    iret = nf90_get_var(ncid, VarID, year_weight(:,:)  , start=(/1,1/),count=(/nyears,nplumes   /))
!EB+    IF (iret /= NF90_NOERR) STOP 'NetCDF Error reading year_weight'
    IF (iret /= NF90_NOERR) call endrun( 'NetCDF Error reading year_weight')
    iret = nf90_inq_varid(ncid, "ann_cycle"     , VarId)
    iret = nf90_get_var(ncid, VarID, ann_cycle(:,:,:)  , start=(/1,1,1/),count=(/nfeatures,ntimes,nplumes/))
!EB+    IF (iret /= NF90_NOERR) STOP 'NetCDF Error reading ann_cycle'
    IF (iret /= NF90_NOERR) call endrun( 'NetCDF Error reading ann_cycle')

    iret = nf90_close(ncid)

    sp_initialized = .TRUE.

    RETURN
  END SUBROUTINE sp_setup
  !
  ! ------------------------------------------------------------------------------------------------------------------------
  ! SET_TIME_WEIGHT:  The simple plume model assumes that meteorology constrains plume shape and that only source strength
  ! influences the amplitude of a plume associated with a given source region.   This routine retrieves the temporal weights
  ! for the plumes.  Each plume feature has its own temporal weights which varies yearly.  The annual cycle is indexed by
  ! week in the year and superimposed on the yearly mean value of the weight. 
  !
  SUBROUTINE set_time_weight(year_fr)
    !
    ! ---------- 
    !
    REAL(r8), INTENT(IN) ::  &
         year_fr           !< Fractional Year (1850.0 - 2100.99)

    INTEGER          ::  &
         idate(8)       ,& !< integer array of system clock date yyyy, mm, dd
         iyear          ,& !< Integer year values between 1 and 156 (1850-2100) 
         iweek          ,& !< Integer index (between 1 and ntimes); for ntimes=52 this corresponds to weeks (roughly)
         iplume            ! plume number
    !
    ! ---------- 
    !
    iyear = FLOOR(year_fr) - 1850 + 1
    iweek = FLOOR((year_fr - FLOOR(year_fr)) * ntimes) + 1

!EB+    IF ((iweek > ntimes) .OR. (iweek < 1) .OR. (iyear > nyears) .OR. (iyear < 1)) STOP 'Time out of bounds in set_time_weight'
!    write(iulog,*)' EB+    iweek = ',iweek,', iyear = ',iyear
    IF ((iweek > ntimes) .OR. (iweek < 1) .OR. (iyear > nyears) .OR. (iyear < 1)) call endrun( 'Time out of bounds in set_time_weight')
    DO iplume=1,nplumes
      time_weight(1,iplume) = year_weight(iyear,iplume) * ann_cycle(1,iweek,iplume)
      time_weight(2,iplume) = year_weight(iyear,iplume) * ann_cycle(2,iweek,iplume)
    END DO
    !
    ! check if code is being used in beta version and stop if it is after June 1, 2016
    !
!EB+    CALL DATE_AND_TIME(VALUES=idate)
!EB+    IF (idate(1) >= 2016 .AND. idate(2) >= 6) THEN
!EB+      PRINT '(A102)', 'System clock says thde date is after June 2016 at which time this beta release should be depreciated.'
!EB+      PRINT '(A85)' , 'Terminating; please check with developers for a finalized CMIP6 version of the code.'
!EB+      STOP
!EB+      call endrun
!EB+    END IF
    
    RETURN
  END SUBROUTINE set_time_weight
  !
  ! ------------------------------------------------------------------------------------------------------------------------
  ! SP_AOP_PROFILE:  This subroutine calculates the simple plume aerosol and cloud active optical properites based on the
  ! the simple plume fit to the MPI Aerosol Climatology (Version 2).  It sums over nplumes to provide a profile of aerosol
  ! optical properties on a host models vertical grid. 
  !
  SUBROUTINE sp_aop_profile                                                                           ( &
       nlevels        ,ncol           ,lambda         ,oro            ,lon            ,lat            , &
       year_fr        ,z              ,dz             ,xre            ,aod_prof       ,ssa_prof       , &
       asy_prof       )
    !
    ! ---------- 
    !
    INTEGER, INTENT(IN)        :: &
         nlevels,                 & !< number of levels
         ncol                       !< number of columns

    REAL(r8), INTENT(IN)           :: &
         lambda,                  & !< wavelength
         year_fr,                 & !< Fractional Year (1903.0 is the 0Z on the first of January 1903, Gregorian)
         oro(ncol),               & !< orographic height (m)
         lon(ncol),               & !< longitude 
         lat(ncol),               & !< latitude
         z (ncol,nlevels),        & !< height above sea-level (m)
         dz(ncol,nlevels)           !< level thickness (difference between half levels)

    REAL(r8), INTENT(OUT)          :: &
         xre(ncol)              , & !< multiplicative change factor for effective radius  (1.0 implies no change)
         aod_prof(ncol,nlevels) , & !< profile of aerosol optical depth
         ssa_prof(ncol,nlevels) , & !< profile of single scattering albedo
         asy_prof(ncol,nlevels)     !< profile of asymmetry parameter

    INTEGER                    :: iplume, icol, k

    REAL(r8)                       ::  &
         eta(ncol,nlevels),        & !< normalized height (by 15 km)
         z_beta(ncol,nlevels),     & !< profile for scaling column optical depth
         prof(ncol,nlevels),       & !< scaled profile (by beta function)
         beta_sum(ncol),           & !< vertical sum of beta function
         ssa(ncol),                & !< aerosol optical depth 
         asy(ncol),                & !< aerosol optical depth 
         cw_an(ncol),              & !< column weight for simple plume (anthropogenic) aod at 550 nm
         cw_bg(ncol),              & !< column weight for fine-mode background aod at 550 nm
         caod_sp(ncol),            & !< column simple plume (anthropogenic) aod at 550 nm
         caod_bg(ncol),            & !< column fine-mode background aod at 550 nm
         a_plume1,                 & !< gaussian longitude factor for feature 1
         a_plume2,                 & !< gaussian longitude factor for feature 2
         b_plume1,                 & !< gaussian latitude factor for feature 1
         b_plume2,                 & !< gaussian latitude factor for feature 2
         delta_lat,                & !< latitude offset
         delta_lon,                & !< longitude offset
         lon1,                     & !< rotated longitude for feature 1
         lat1,                     & !< rotated latitude for feature 2
         lon2,                     & !< rotated longitude for feature 1
         lat2,                     & !< rotated latitude for feature 2
         f1,                       & !< contribution from feature 1
         f2,                       & !< contribution from feature 2
         aod_550,                  & !< aerosol optical depth at 550nm
         aod_lmd,                  & !< aerosol optical depth at input wavelength
         lfactor                     !< factor to compute wavelength dependence of optical properties
    !
    ! ---------- 
    !
    ! initialize input data (by calling setup at first instance) 
    !
    IF (.NOT.sp_initialized) CALL sp_setup
    !
    ! get time weights
    !
    CALL set_time_weight(year_fr)
    !
    ! initialize variables, including output
    !
    DO k=1,nlevels
      DO icol=1,ncol
        aod_prof(icol,k) = 0.0
        ssa_prof(icol,k) = 0.0
        asy_prof(icol,k) = 0.0
        z_beta(icol,k)   = MERGE(1.0, 0.0, z(icol,k) >= oro(icol))
        eta(icol,k)      = MAX(0.0,MIN(1.0,z(icol,k)/15000.))
      END DO
    END DO
    DO icol=1,ncol
      xre(icol)      = 1.0_r8
      caod_sp(icol)  = 0.00_r8
      caod_bg(icol)  = 0.02_r8
    END DO
    !
    ! sum contribution from plumes to construct composite profiles of aerosol otpical properties
    !
    DO iplume=1,nplumes
      !
      ! calculate vertical distribution function from parameters of beta distribution
      !
      DO icol=1,ncol
        beta_sum(icol) = 0._r8
      END DO
      DO k=1,nlevels
        DO icol=1,ncol
          prof(icol,k)   = (eta(icol,k)**(beta_a(iplume)-1.) * (1.-eta(icol,k))**(beta_b(iplume)-1.))*dz(icol,k)
          beta_sum(icol) = beta_sum(icol) + prof(icol,k)
        END DO
      END DO
      DO k=1,nlevels
        DO icol=1,ncol
          prof(icol,k)   = prof(icol,k) / beta_sum(icol) * z_beta(icol,k)
        END DO
      END DO
      !
      ! calculate plume weights
      !
      DO icol=1,ncol
        !
        ! get plume-center relative spatial parameters for specifying amplitude of plume at given lat and lon
        !
        delta_lat = lat(icol) - plume_lat(iplume)
        delta_lon = lon(icol) - plume_lon(iplume)
!EB        delta_lon = MERGE ( delta_lon-SIGN(360.,delta_lon) , delta_lon , ABS(delta_lon) > 180. )
        delta_lon = MERGE ( delta_lon-SIGN(360.,delta_lon) , delta_lon , ABS(delta_lon) > 180._r8 )
        a_plume1  = 0.5_r8/ (MERGE(sig_lon_E(1,iplume), sig_lon_W(1,iplume), delta_lon > 0)**2)
        b_plume1  = 0.5_r8 / (MERGE(sig_lat_E(1,iplume), sig_lat_W(1,iplume), delta_lon > 0)**2)
        a_plume2  = 0.5_r8 / (MERGE(sig_lon_E(2,iplume), sig_lon_W(2,iplume), delta_lon > 0)**2)
        b_plume2  = 0.5_r8 / (MERGE(sig_lat_E(2,iplume), sig_lat_W(2,iplume), delta_lon > 0)**2)
        !
        ! adjust for a plume specific rotation which helps match plume state to climatology.
        !
        lon1 =   COS(theta(1,iplume))*(delta_lon) + SIN(theta(1,iplume))*(delta_lat)
        lat1 = - SIN(theta(1,iplume))*(delta_lon) + COS(theta(1,iplume))*(delta_lat)
        lon2 =   COS(theta(2,iplume))*(delta_lon) + SIN(theta(2,iplume))*(delta_lat)
        lat2 = - SIN(theta(2,iplume))*(delta_lon) + COS(theta(2,iplume))*(delta_lat)
        !
        ! calculate contribution to plume from its different features, to get a column weight for the anthropogenic
        ! (cw_an) and the fine-mode background aerosol (cw_bg)
        !
        f1 = time_weight(1,iplume) * ftr_weight(1,iplume) * EXP(-1.* (a_plume1 * ((lon1)**2) + (b_plume1 * ((lat1)**2)))) 
        f2 = time_weight(2,iplume) * ftr_weight(2,iplume) * EXP(-1.* (a_plume2 * ((lon2)**2) + (b_plume2 * ((lat2)**2)))) 

        cw_an(icol) = f1 * aod_spmx(iplume) + f2 * aod_spmx(iplume)  
        cw_bg(icol) = f1 * aod_fmbg(iplume) + f2 * aod_fmbg(iplume) 
        !
        ! calculate wavelength-dependent scattering properties
        !
        lfactor   = MIN(1.0_r8,700.0_r8/lambda)
        ssa(icol) = (ssa550(iplume) * lfactor**4) / ((ssa550(iplume) * lfactor**4) + ((1-ssa550(iplume)) * lfactor))
        asy(icol) =  asy550(iplume) * SQRT(lfactor)
      END DO
      !
      ! distribute plume optical properties across its vertical profile weighting by optical depth and scaling for
      ! wavelength using the anstrom parameter. 
      !      
      lfactor = EXP(-angstrom(iplume) * LOG(lambda/550.0_r8))
      DO k=1,nlevels
        DO icol = 1,ncol
          aod_550          = prof(icol,k)     * cw_an(icol)
          aod_lmd          = aod_550          * lfactor
          caod_sp(icol)    = caod_sp(icol)    + prof(icol,k) * cw_an(icol)
          caod_bg(icol)    = caod_bg(icol)    + prof(icol,k) * cw_bg(icol)
          aod_prof(icol,k) = aod_prof(icol,k) + aod_lmd
          ssa_prof(icol,k) = ssa_prof(icol,k) + aod_lmd * ssa(iplume)
          asy_prof(icol,k) = asy_prof(icol,k) + aod_lmd * ssa(iplume) * asy(iplume)
        END DO
      END DO
    END DO
    !
    ! complete optical depth weighting
    !
    DO k=1,nlevels
      DO icol = 1,ncol
        asy_prof(icol,k) = MERGE(asy_prof(icol,k)/ssa_prof(icol,k), 0.0_r8, ssa_prof(icol,k) > 0._r8)
        ssa_prof(icol,k) = MERGE(ssa_prof(icol,k)/aod_prof(icol,k), 0.0_r8, aod_prof(icol,k) > 0._r8)
      END DO
    END DO
    !
    ! calcuate effective radius normalization (divisor) factor
    !
    DO icol=1,ncol
      xre(icol) = (20 * LOG((1000.0_r8 * caod_bg(icol)) + 3.0_r8))/(20 * LOG((1000.0_r8 * (caod_sp(icol) + caod_bg(icol))) + 3.0_r8))
    END DO

    RETURN
  END SUBROUTINE sp_aop_profile

END MODULE MO_SIMPLE_PLUMES
