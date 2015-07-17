MODULE par_bfm
   !!======================================================================
   !!                        ***  par_bfm  ***
   !! TOP :   set the BFM parameters
   !!======================================================================
   !! History :   2.0  !  2007-12  (C. Ethe, G. Madec)  revised architecture
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: par_my_trc.F90 2528 2010-12-27 17:33:53Z rblod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

   IMPLICIT NONE

#if defined key_bfm
   !!---------------------------------------------------------------------
   !!   'key_bfm'                     user defined tracers (BFM)
   !!---------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_bfm     = .TRUE.   !: BFM flag 
   INTEGER, PUBLIC, PARAMETER ::   jp_bfm     =  1       !: number of BFM tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_bfm_2d  =  0       !: additional 2d output arrays ('key_trc_diaadd')
   INTEGER, PUBLIC, PARAMETER ::   jp_bfm_3d  =  0       !: additional 3d output arrays ('key_trc_diaadd')
   INTEGER, PUBLIC, PARAMETER ::   jp_bfm_trd =  0       !: number of sms trends for BFM

#else
   !!---------------------------------------------------------------------
   !!   Default                           No user defined tracers (MY_TRC)
   !!---------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_bfm     = .FALSE.   !: BFM flag 
   INTEGER, PUBLIC, PARAMETER ::   jp_bfm     =  0       !: number of BFM tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_bfm_2d  =  0       !: additional 2d output arrays ('key_trc_diaadd')
   INTEGER, PUBLIC, PARAMETER ::   jp_bfm_3d  =  0       !: additional 3d output arrays ('key_trc_diaadd')
   INTEGER, PUBLIC, PARAMETER ::   jp_bfm_trd =  0       !: number of sms trends for BFM
#endif
   !!======================================================================
END MODULE par_bfm
