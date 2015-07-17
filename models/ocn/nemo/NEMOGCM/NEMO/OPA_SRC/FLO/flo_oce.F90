MODULE flo_oce
   !!======================================================================
   !!                     ***  MODULE flo_oce  ***
   !! lagrangian floats :   define in memory all floats parameters and variables
   !!======================================================================
   !! History :   OPA  ! 1999-10  (CLIPPER projet)
   !!   NEMO      1.0  ! 2002-11  (G. Madec, A. Bozec)  F90: Free form and module
   !!----------------------------------------------------------------------
#if   defined key_floats   ||   defined key_esopa
   !!----------------------------------------------------------------------
   !!   'key_floats'                                        drifting floats
   !!----------------------------------------------------------------------
   USE par_oce         ! ocean parameters
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! MPP library

   IMPLICIT NONE
   PUBLIC

   PUBLIC   flo_oce_alloc   ! Routine called in floats.F90

   LOGICAL, PUBLIC, PARAMETER ::   lk_floats = .TRUE.    !: float flag

   !! float parameters
   !! ----------------
   INTEGER, PUBLIC, PARAMETER ::   jpnfl     = 23                  !: total number of floats during the run
   INTEGER, PUBLIC, PARAMETER ::   jpnnewflo =  0                  !: number of floats added in a new run
   INTEGER, PUBLIC, PARAMETER ::   jpnrstflo = jpnfl - jpnnewflo   !: number of floats for the restart

   !! float variables
   !! ---------------
   INTEGER , PUBLIC, DIMENSION(jpnfl) ::   nisobfl   !: =0 for a isobar float , =1 for a float following the w velocity
   INTEGER , PUBLIC, DIMENSION(jpnfl) ::   ngrpfl    !: number to identify searcher group

   REAL(wp), PUBLIC, DIMENSION(jpnfl) ::   flxx , flyy , flzz    !: long, lat, depth of float (decimal degree, m >0)
   REAL(wp), PUBLIC, DIMENSION(jpnfl) ::   tpifl, tpjfl, tpkfl   !: (i,j,k) indices of float position

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   wb   !: vertical velocity at previous time step (m s-1).
   
   !                                            !!! * namelist namflo : langrangian floats *
   LOGICAL, PUBLIC  ::   ln_rstflo  = .FALSE.    !: T/F float restart 
   LOGICAL, PUBLIC  ::   ln_argo    = .FALSE.    !: T/F argo type floats
   LOGICAL, PUBLIC  ::   ln_flork4  = .FALSE.    !: T/F 4th order Runge-Kutta
   INTEGER, PUBLIC  ::   nn_writefl = 150        !: frequency of float output file 
   INTEGER, PUBLIC  ::   nn_stockfl = 450        !: frequency of float restart file

   !!----------------------------------------------------------------------
   !! NEMO/OPA 4.0 , NEMO Consortium (2011)
   !! $Id: flo_oce.F90 2715 2011-03-30 15:58:35Z rblod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION flo_oce_alloc()
      !!----------------------------------------------------------------------
      !!                 ***  FUNCTION flo_oce_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( wb(jpi,jpj,jpk)   , STAT=flo_oce_alloc )
      !
      IF( lk_mpp             )   CALL mpp_sum ( flo_oce_alloc )
      IF( flo_oce_alloc /= 0 )   CALL ctl_warn('flo_oce_alloc: failed to allocate arrays')
   END FUNCTION flo_oce_alloc

#else
   !!----------------------------------------------------------------------
   !!   Default option :                                 NO drifting floats
   !!----------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_floats = .FALSE.   !: float flag
#endif

   !!======================================================================
END MODULE flo_oce
