MODULE bdyice
   !!======================================================================
   !!                       ***  MODULE  bdyice  ***
   !! Unstructured Open Boundary Cond. :  Flow Relaxation Scheme applied sea-ice
   !!======================================================================
   !!  History :  3.3  !  2010-09 (D. Storkey)  Original code
   !!----------------------------------------------------------------------
#if defined   key_bdy   &&   defined key_lim2
   !!----------------------------------------------------------------------
   !!   'key_bdy'            and                 Unstructured Open Boundary Conditions
   !!   'key_lim2'                                                 LIM-2 sea ice model
   !!----------------------------------------------------------------------
   !!   bdy_ice_frs        : Relaxation of tracers on unstructured open boundaries
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers variables
   USE ice_2           ! LIM_2 ice variables
   USE dom_oce         ! ocean space and time domain variables 
   USE bdy_oce         ! ocean open boundary conditions
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE in_out_manager  ! write to numout file
   USE lib_mpp         ! distributed memory computing
   
   IMPLICIT NONE
   PRIVATE

   PUBLIC   bdy_ice_frs    ! routine called in sbcmod

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: bdyice.F90 2715 2011-03-30 15:58:35Z rblod $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE bdy_ice_frs( kt )
      !!------------------------------------------------------------------------------
      !!                 ***  SUBROUTINE bdy_ice_frs  ***
      !!                    
      !! ** Purpose : Apply the Flow Relaxation Scheme for sea-ice fields in the case 
      !!              of unstructured open boundaries. Currently only tested for LIM2.
      !! 
      !! Reference : Engedahl H., 1995: Use of the flow relaxation scheme in a three-
      !!             dimensional baroclinic ocean model with realistic topography. Tellus, 365-382.
      !!------------------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt   ! model time step index
      !!
      INTEGER  ::   jb, jk, jgrd   ! dummy loop indices
      INTEGER  ::   ii, ij         ! local scalar
      REAL(wp) ::   zwgt, zwgt1    ! local scalar
      !!------------------------------------------------------------------------------
      !
      jgrd = 1      ! Everything is at T-points here
      !
      IF( ln_ice_frs ) THEN     ! update ice fields by relaxation at the boundary
         DO jb = 1, nblen(jgrd)
            DO jk = 1, jpkm1
               ii    = nbi(jb,jgrd)
               ij    = nbj(jb,jgrd)
               zwgt  = nbw(jb,jgrd)
               zwgt1 = 1.e0 - nbw(jb,jgrd)
               frld (ii,ij) = ( frld (ii,ij) * zwgt1 + frld_bdy (jb) * zwgt ) * tmask(ii,ij,1)     ! Leads fraction 
               hicif(ii,ij) = ( hicif(ii,ij) * zwgt1 + hicif_bdy(jb) * zwgt ) * tmask(ii,ij,1)     ! Ice depth 
               hsnif(ii,ij) = ( hsnif(ii,ij) * zwgt1 + hsnif_bdy(jb) * zwgt ) * tmask(ii,ij,1)     ! Snow depth
            END DO
         END DO 
         CALL lbc_lnk( frld, 'T', 1. )                                         ! lateral boundary conditions
         CALL lbc_lnk( hicif, 'T', 1. )   ;   CALL lbc_lnk( hsnif, 'T', 1. )
         !
      ELSE                          ! we have called this routine without ln_ice_frs not set
         IF( kt == nit000 )   CALL ctl_warn( 'E R R O R (possible) called bdy_ice_frs when ln_ice_frs is false?' )
      ENDIF
      !      
   END SUBROUTINE bdy_ice_frs
#else
   !!---------------------------------------------------------------------------------
   !!   Default option                                                    Empty module
   !!---------------------------------------------------------------------------------
CONTAINS
   SUBROUTINE bdy_ice_frs( kt )      ! Empty routine
      WRITE(*,*) 'bdy_ice_frs: You should not have seen this print! error?', kt
   END SUBROUTINE bdy_ice_frs
#endif

   !!=================================================================================
END MODULE bdyice
