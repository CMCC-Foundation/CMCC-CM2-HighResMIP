MODULE bdy_oce
   !!======================================================================
   !!                       ***  MODULE bdy_oce   ***
   !! Unstructured Open Boundary Cond. :   define related variables
   !!======================================================================
   !! History :  1.0  !  2001-05  (J. Chanut, A. Sellar)  Original code
   !!            3.0  !  2008-04  (NEMO team)  add in the reference version     
   !!            3.3  !  2010-09  (D. Storkey) add ice boundary conditions
   !!----------------------------------------------------------------------
#if defined key_bdy 
   !!----------------------------------------------------------------------
   !!   'key_bdy'                      Unstructured Open Boundary Condition
   !!----------------------------------------------------------------------
   USE par_oce         ! ocean parameters
   USE bdy_par         ! Unstructured boundary parameters
   USE lib_mpp         ! distributed memory computing

   IMPLICIT NONE
   PUBLIC

   !!----------------------------------------------------------------------
   !! Namelist variables
   !!----------------------------------------------------------------------
   CHARACTER(len=80) ::   cn_mask        !: Name of unstruct. bdy mask file
   CHARACTER(len=80) ::   cn_dta_frs_T   !: Name of unstruct. bdy data file at T points for FRS conditions
   CHARACTER(len=80) ::   cn_dta_frs_U   !: Name of unstruct. bdy data file at U points for FRS conditions
   CHARACTER(len=80) ::   cn_dta_frs_V   !: Name of unstruct. bdy data file at V points for FRS conditions
   CHARACTER(len=80) ::   cn_dta_fla_T   !: Name of unstruct. bdy data file at T points for Flather scheme
   CHARACTER(len=80) ::   cn_dta_fla_U   !: Name of unstruct. bdy data file at U points for Flather scheme
   CHARACTER(len=80) ::   cn_dta_fla_V   !: Name of unstruct. bdy data file at V points for Flather scheme
   !
   LOGICAL ::   ln_tides = .false.    !: =T apply tidal harmonic forcing along open boundaries
   LOGICAL ::   ln_vol  = .false.     !: =T volume correction             
   LOGICAL ::   ln_mask = .false.     !: =T read bdymask from file
   LOGICAL ::   ln_clim = .false.     !: =T bdy data files contain  1 time dump  (-->bdy forcing will be constant) 
   !                                  !                         or 12 months     (-->bdy forcing will be cyclic) 
   LOGICAL ::   ln_dyn_fla  = .false. !: =T Flather boundary conditions on barotropic velocities
   LOGICAL ::   ln_dyn_frs  = .false. !: =T FRS boundary conditions on velocities
   LOGICAL ::   ln_tra_frs  = .false. !: =T FRS boundary conditions on tracers (T and S)
   LOGICAL ::   ln_ice_frs  = .false. !: =T FRS boundary conditions on seaice (leads fraction, ice depth, snow depth)
   !
   INTEGER ::   nn_rimwidth = 7       !: boundary rim width
   INTEGER ::   nn_dtactl   = 1	     !: = 0 use the initial state as bdy dta ; = 1 read it in a NetCDF file
   INTEGER ::   nn_volctl   = 1       !: = 0 the total volume will have the variability of the surface Flux E-P 
   !                                  !  = 1 the volume will be constant during all the integration.

   !!----------------------------------------------------------------------
   !! Global variables
   !!----------------------------------------------------------------------
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   bdytmask   !: Mask defining computational domain at T-points
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   bdyumask   !: Mask defining computational domain at U-points
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   bdyvmask   !: Mask defining computational domain at V-points

   !!----------------------------------------------------------------------
   !! Unstructured open boundary data variables
   !!----------------------------------------------------------------------
   INTEGER, DIMENSION(jpbgrd) ::   nblen    = 0           !: Size of bdy data on a proc for each grid type
   INTEGER, DIMENSION(jpbgrd) ::   nblenrim = 0           !: Size of bdy data on a proc for first rim ind
   INTEGER, DIMENSION(jpbgrd) ::   nblendta = 0           !: Size of bdy data in file

   INTEGER, DIMENSION(jpbdim,jpbgrd) ::   nbi, nbj        !: i and j indices of bdy dta
   INTEGER, DIMENSION(jpbdim,jpbgrd) ::   nbr             !: Discrete distance from rim points
   INTEGER, DIMENSION(jpbdim,jpbgrd) ::   nbmap           !: Indices of data in file for data in memory 
    
   REAL(wp) ::   bdysurftot	                            !: Lateral surface of unstructured open boundary

   REAL(wp), DIMENSION(jpbdim)        ::   flagu, flagv   !: Flag for normal velocity compnt for velocity components
   REAL(wp), DIMENSION(jpbdim,jpbgrd) ::   nbw            !: Rim weights of bdy data

   REAL(wp), DIMENSION(jpbdim)     ::   sshbdy            !: Now clim of bdy sea surface height (Flather)
   REAL(wp), DIMENSION(jpbdim)     ::   ubtbdy, vbtbdy    !: Now clim of bdy barotropic velocity components
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   tbdy  , sbdy      !: Now clim of bdy temperature and salinity  
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   ubdy  , vbdy	  !: Now clim of bdy velocity components
   REAL(wp), DIMENSION(jpbdim) ::   sshtide               !: Tidal boundary array : SSH
   REAL(wp), DIMENSION(jpbdim) ::   utide, vtide          !: Tidal boundary array : U and V
#if defined key_lim2
   REAL(wp), DIMENSION(jpbdim) ::   frld_bdy    !: now ice leads fraction climatology   
   REAL(wp), DIMENSION(jpbdim) ::   hicif_bdy   !: Now ice  thickness climatology
   REAL(wp), DIMENSION(jpbdim) ::   hsnif_bdy   !: now snow thickness
#endif

   !!----------------------------------------------------------------------
   !! NEMO/OPA 4.0 , NEMO Consortium (2011)
   !! $Id: bdy_oce.F90 2715 2011-03-30 15:58:35Z rblod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   FUNCTION bdy_oce_alloc()
      !!----------------------------------------------------------------------
      USE lib_mpp, ONLY: ctl_warn, mpp_sum
      !
      INTEGER :: bdy_oce_alloc
      !!----------------------------------------------------------------------
      !
      ALLOCATE( bdytmask(jpi,jpj) , tbdy(jpbdim,jpk) , sbdy(jpbdim,jpk) ,     &
         &      bdyumask(jpi,jpj) , ubdy(jpbdim,jpk) ,                        &
         &      bdyvmask(jpi,jpj) , vbdy(jpbdim,jpk) ,                    STAT=bdy_oce_alloc )
         !
      IF( lk_mpp             )   CALL mpp_sum ( bdy_oce_alloc )
      IF( bdy_oce_alloc /= 0 )   CALL ctl_warn('bdy_oce_alloc: failed to allocate arrays.')
      !
   END FUNCTION bdy_oce_alloc

#else
   !!----------------------------------------------------------------------
   !!   Dummy module                NO Unstructured Open Boundary Condition
   !!----------------------------------------------------------------------
   LOGICAL ::   ln_tides = .false.  !: =T apply tidal harmonic forcing along open boundaries
#endif

   !!======================================================================
END MODULE bdy_oce
