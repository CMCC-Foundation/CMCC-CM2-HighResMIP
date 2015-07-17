MODULE traswp
   !!==============================================================================
   !!                       ***  MODULE  traswp  ***
   !! Ocean active tracers: swapping array 
   !!==============================================================================
   USE par_oce         ! ocean parameters
   USE oce             ! ocean dynamics and active tracers

   IMPLICIT NONE
   PRIVATE

   PUBLIC   tra_swap     ! routine called by step.F90
   PUBLIC   tra_unswap   ! routine called by step.F90

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: traswp.F90 2715 2011-03-30 15:58:35Z rblod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE tra_swap 
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_swp  ***
      !!                   
      !! ** Purpose : Store temperature and salinity aaray into a 4D array 
      !!
      !!----------------------------------------------------------------------
      !
      tsn(:,:,:,jp_tem) = tn(:,:,:)      ;      tsn(:,:,:,jp_sal) = sn(:,:,:)
      tsb(:,:,:,jp_tem) = tb(:,:,:)      ;      tsb(:,:,:,jp_sal) = sb(:,:,:)
      tsa(:,:,:,jp_tem) = ta(:,:,:)      ;      tsa(:,:,:,jp_sal) = sa(:,:,:)
      !
   END SUBROUTINE tra_swap

   SUBROUTINE tra_unswap 
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_unswap  ***
      !!                   
      !! ** Purpose : Store temperature and salinity aaray into a 4D array 
      !!
      !!----------------------------------------------------------------------
      !
      tn(:,:,:) = tsn(:,:,:,jp_tem)      ;      sn(:,:,:) = tsn(:,:,:,jp_sal)
      tb(:,:,:) = tsb(:,:,:,jp_tem)      ;      sb(:,:,:) = tsb(:,:,:,jp_sal)
      ta(:,:,:) = tsa(:,:,:,jp_tem)      ;      sa(:,:,:) = tsa(:,:,:,jp_sal)
      !
   END SUBROUTINE tra_unswap

   !!======================================================================
END MODULE traswp
