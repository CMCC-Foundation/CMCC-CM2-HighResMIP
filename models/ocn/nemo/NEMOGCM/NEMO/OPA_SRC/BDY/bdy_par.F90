MODULE bdy_par
   !!======================================================================
   !!                      ***  MODULE bdy_par   ***
   !! Unstructured Open Boundary Cond. :   define related parameters
   !!======================================================================
   !! History :  1.0  !  2005-01  (J. Chanut, A. Sellar)  Original code
   !!            3.0  !  2008-04  (NEMO team)  add in the reference version
   !!            3.3  !  2010-09  (D. Storkey and E. O'Dea) update for Shelf configurations
   !!----------------------------------------------------------------------
#if defined   key_bdy
   !!----------------------------------------------------------------------
   !!   'key_bdy' :                    Unstructured Open Boundary Condition
   !!----------------------------------------------------------------------

   IMPLICIT NONE
   PUBLIC

   LOGICAL, PUBLIC, PARAMETER ::   lk_bdy  = .TRUE.   !: Unstructured Ocean Boundary Condition flag
   INTEGER, PUBLIC, PARAMETER ::   jpbdta  = 20000    !: Max length of bdy field in file
   INTEGER, PUBLIC, PARAMETER ::   jpbdim  = 20000    !: Max length of bdy field on a processor
   INTEGER, PUBLIC, PARAMETER ::   jpbtime = 1000     !: Max number of time dumps per file
   INTEGER, PUBLIC, PARAMETER ::   jpbgrd  = 6	      !: Number of horizontal grid types used  (T, u, v, f)
#else
   !!----------------------------------------------------------------------
   !!   Default option :            NO Unstructured open boundary condition
   !!----------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_bdy = .FALSE.   !: Unstructured Ocean Boundary Condition flag
#endif

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: bdy_par.F90 2528 2010-12-27 17:33:53Z rblod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!======================================================================
END MODULE bdy_par
