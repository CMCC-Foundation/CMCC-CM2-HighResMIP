MODULE trc
   !!======================================================================
   !!                      ***  MODULE  trc  ***
   !! Passive tracers   :  module for tracers defined
   !!======================================================================
   !! History :   OPA  !  1996-01  (M. Levy)  Original code
   !!              -   !  1999-07  (M. Levy)  for LOBSTER1 or NPZD model
   !!              -   !  2000-04  (O. Aumont, M.A. Foujols)  HAMOCC3 and P3ZD
   !!   NEMO      1.0  !  2004-03  (C. Ethe)  Free form and module
   !!----------------------------------------------------------------------
#if defined key_top
   !!----------------------------------------------------------------------
   !!   'key_top'                                                TOP models
   !!----------------------------------------------------------------------
   USE par_oce
   USE par_trc
   
   IMPLICIT NONE
   PUBLIC

   PUBLIC   trc_alloc   ! called by nemogcm.F90

   !! passive tracers names and units (read in namelist)
   !! --------------------------------------------------
   CHARACTER(len=12), PUBLIC, DIMENSION(jptra) ::   ctrcnm     !: tracer name 
   CHARACTER(len=12), PUBLIC, DIMENSION(jptra) ::   ctrcun     !: tracer unit
   CHARACTER(len=80), PUBLIC, DIMENSION(jptra) ::   ctrcnl     !: tracer long name 
   
   
   !! parameters for the control of passive tracers
   !! --------------------------------------------------
   INTEGER, PUBLIC                   ::   numnat   !: the number of the passive tracer NAMELIST
   LOGICAL, PUBLIC, DIMENSION(jptra) ::   lutini   !:  initialisation from FILE or not (NAMELIST)
   LOGICAL, PUBLIC, DIMENSION(jptra) ::   lutsav   !:  save the tracer or not

   !! passive tracers fields (before,now,after)
   !! --------------------------------------------------
   REAL(wp), PUBLIC ::   trai                          !: initial total tracer
   REAL(wp), PUBLIC ::   areatot                       !: total volume 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION (:,:,:)   ::   cvol   !: volume correction -degrad option- 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION (:,:,:,:) ::   trn    !: traceur concentration for now time step
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION (:,:,:,:) ::   tra    !: traceur concentration for next time step
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION (:,:,:,:) ::   trb    !: traceur concentration for before time step

   !! interpolated gradient
   !!--------------------------------------------------  
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION (:,:,:) ::   gtru   !: hor. gradient at u-points at bottom ocean level
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION (:,:,:) ::   gtrv   !: hor. gradient at v-points at bottom ocean level
   
   !! passive tracers restart (input and output)
   !! ------------------------------------------  
   LOGICAL          , PUBLIC ::  ln_rsttr        !: boolean term for restart i/o for passive tracers (namelist)
   LOGICAL          , PUBLIC ::  lrst_trc        !: logical to control the trc restart write
   INTEGER          , PUBLIC ::  nn_dttrc        !: frequency of step on passive tracers
   INTEGER          , PUBLIC ::  nutwrs          !: output FILE for passive tracers restart
   INTEGER          , PUBLIC ::  nutrst          !: logical unit for restart FILE for passive tracers
   INTEGER          , PUBLIC ::  nn_rsttr        !: control of the time step ( 0 or 1 ) for pass. tr.
   CHARACTER(len=50), PUBLIC ::  cn_trcrst_in    !: suffix of pass. tracer restart name (input)
   CHARACTER(len=50), PUBLIC ::  cn_trcrst_out   !: suffix of pass. tracer restart name (output)
   
   !! information for outputs
   !! --------------------------------------------------
   INTEGER , PUBLIC ::   nn_writetrc   !: time step frequency for concentration outputs (namelist)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   rdttrc        !: vertical profile of passive tracer time step
   
# if defined key_diatrc && ! defined key_iomput
   !! additional 2D/3D outputs namelist
   !! --------------------------------------------------
   INTEGER         , PUBLIC                      ::   nn_writedia   !: frequency of additional arrays outputs(namelist)
   CHARACTER(len= 8), PUBLIC, DIMENSION(jpdia2d) ::   ctrc2d      !: 2d output field name
   CHARACTER(len= 8), PUBLIC, DIMENSION(jpdia2d) ::   ctrc2u      !: 2d output field unit   
   CHARACTER(len= 8), PUBLIC, DIMENSION(jpdia3d) ::   ctrc3d      !: 3d output field name
   CHARACTER(len= 8), PUBLIC, DIMENSION(jpdia3d) ::   ctrc3u      !: 3d output field unit
   CHARACTER(len=80), PUBLIC, DIMENSION(jpdia2d) ::   ctrc2l      !: 2d output field long name
   CHARACTER(len=80), PUBLIC, DIMENSION(jpdia3d) ::   ctrc3l      !: 3d output field long name

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,  :) ::   trc2d    !:  additional 2d outputs  
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) ::   trc3d    !:  additional 3d outputs  
# endif

# if defined key_diabio || defined key_trdmld_trc
   !                                                              !!*  namtop_XXX namelist *
   INTEGER , PUBLIC                               ::   nn_writebio   !: time step frequency for biological outputs 
   CHARACTER(len=8 ), PUBLIC, DIMENSION(jpdiabio) ::   ctrbio      !: biological trends name      
   CHARACTER(len=20), PUBLIC, DIMENSION(jpdiabio) ::   ctrbiu      !: biological trends unit   
   CHARACTER(len=80), PUBLIC, DIMENSION(jpdiabio) ::   ctrbil      !: biological trends long name
# endif
# if defined key_diabio
   !! Biological trends
   !! -----------------
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) ::   trbio   !: biological trends
# endif

   
   !! passive tracers data read and at given time_step
   !! --------------------------------------------------
# if defined key_dtatrc
   INTEGER , PUBLIC, DIMENSION(jptra) ::   numtr   !: logical unit for passive tracers data
# endif

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3.1 , NEMO Consortium (2010)
   !! $Id: trc.F90 2715 2011-03-30 15:58:35Z rblod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION trc_alloc()
      !!-------------------------------------------------------------------
      !!                    *** ROUTINE trc_alloc ***
      !!-------------------------------------------------------------------
      USE lib_mpp, ONLY: ctl_warn
      !!-------------------------------------------------------------------
      !
      ALLOCATE( cvol(jpi,jpj,jpk      ) ,                           &
         &      trn (jpi,jpj,jpk,jptra) ,                           &
         &      tra (jpi,jpj,jpk,jptra) ,                           &
         &      trb (jpi,jpj,jpk,jptra) ,                           &
         &      gtru(jpi,jpj    ,jptra) , gtrv(jpi,jpj,jptra) ,     &
# if defined key_diatrc && ! defined key_iomput
         &      trc2d(jpi,jpj,jpdia2d), trc3d(jpi,jpj,jpk,jpdia3d), &
# endif
# if defined key_diabio
         &      trbio(jpi,jpj,jpk,jpdiabio),                        &
#endif
               rdttrc(jpk) ,  STAT=trc_alloc )      

      IF( trc_alloc /= 0 )   CALL ctl_warn('trc_alloc: failed to allocate arrays')
      !
   END FUNCTION trc_alloc

#else
   !!----------------------------------------------------------------------
   !!  Empty module :                                     No passive tracer
   !!----------------------------------------------------------------------
#endif

   !!======================================================================
END MODULE trc
