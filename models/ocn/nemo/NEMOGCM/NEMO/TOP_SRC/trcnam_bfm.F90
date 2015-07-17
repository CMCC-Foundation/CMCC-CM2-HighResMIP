MODULE trcnam_bfm
   !!======================================================================
   !!                      ***  MODULE trcnam_bfm  ***
   !! TOP :   initialisation of some run parameters for BFM bio-model
   !!======================================================================
   !! History :   2.0  !  2007-12  (C. Ethe, G. Madec) Original code
   !!----------------------------------------------------------------------
#if defined key_bfm
   !!----------------------------------------------------------------------
   !!   'key_bfm'   :                                       BFM model
   !!----------------------------------------------------------------------
   !! trc_nam_bfm      : BFM model initialisation
   !!----------------------------------------------------------------------
   USE oce_trc         ! Ocean variables
   USE par_trc         ! TOP parameters
   USE trc             ! TOP variables

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_nam_bfm   ! called by trcnam.F90 module

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcnam_my_trc.F90 2528 2010-12-27 17:33:53Z rblod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trc_nam_bfm
      !!----------------------------------------------------------------------
      !!                     ***  trc_nam_bfm  ***  
      !!
      !! ** Purpose :   read BFM namelist
      !!
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_nam_bfm : read BFM namelists'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~~'
      !
   END SUBROUTINE trc_nam_bfm
   
#else
   !!----------------------------------------------------------------------
   !!  Dummy module :                                             No MY_TRC
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_nam_bfm                      ! Empty routine

   END  SUBROUTINE  trc_nam_bfm
#endif  

   !!======================================================================
END MODULE trcnam_bfm
