MODULE bdytra
   !!======================================================================
   !!                       ***  MODULE  bdytra  ***
   !! Ocean tracers:   Flow Relaxation Scheme of tracers on each open boundary
   !!======================================================================
   !! History :  1.0  !  2005-01  (J. Chanut, A. Sellar)  Original code
   !!            3.0  !  2008-04  (NEMO team)  add in the reference version
   !!----------------------------------------------------------------------
#if defined key_bdy
   !!----------------------------------------------------------------------
   !!   'key_bdy'                     Unstructured Open Boundary Conditions
   !!----------------------------------------------------------------------
   !!   bdy_tra_frs        : Relaxation of tracers on unstructured open boundaries
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain variables 
   USE bdy_oce         ! ocean open boundary conditions
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE in_out_manager  ! I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC bdy_tra_frs     ! routine called in tranxt.F90 

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: bdytra.F90 2528 2010-12-27 17:33:53Z rblod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE bdy_tra_frs( kt )
      !!----------------------------------------------------------------------
      !!                 ***  SUBROUTINE bdy_tra_frs  ***
      !!                    
      !! ** Purpose : Apply the Flow Relaxation Scheme for tracers in the  
      !!              case of unstructured open boundaries.
      !! 
      !! Reference : Engedahl H., 1995, Tellus, 365-382.
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt
      !! 
      REAL(wp) ::   zwgt           ! boundary weight
      INTEGER  ::   ib, ik, igrd   ! dummy loop indices
      INTEGER  ::   ii, ij         ! 2D addresses
      !!----------------------------------------------------------------------
      !
      IF(ln_tra_frs) THEN       ! If this is false, then this routine does nothing. 
         !
         IF( kt == nit000 ) THEN
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'bdy_tra_frs : Flow Relaxation Scheme for tracers'
            IF(lwp) WRITE(numout,*) '~~~~~~~'
         ENDIF
         !
         igrd = 1                       ! Everything is at T-points here
         DO ib = 1, nblen(igrd)
            DO ik = 1, jpkm1
               ii = nbi(ib,igrd)
               ij = nbj(ib,igrd)
               zwgt = nbw(ib,igrd)
               ta(ii,ij,ik) = ( ta(ii,ij,ik) * (1.-zwgt) + tbdy(ib,ik) * zwgt ) * tmask(ii,ij,ik)         
               sa(ii,ij,ik) = ( sa(ii,ij,ik) * (1.-zwgt) + sbdy(ib,ik) * zwgt ) * tmask(ii,ij,ik)
            END DO
         END DO 
         !
         CALL lbc_lnk( ta, 'T', 1. )   ; CALL lbc_lnk( sa, 'T', 1. )    ! Boundary points should be updated
         !
      ENDIF ! ln_tra_frs
      !
   END SUBROUTINE bdy_tra_frs
   
#else
   !!----------------------------------------------------------------------
   !!   Dummy module                   NO Unstruct Open Boundary Conditions
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE bdy_tra_frs(kt)      ! Empty routine
      WRITE(*,*) 'bdy_tra_frs: You should not have seen this print! error?', kt
   END SUBROUTINE bdy_tra_frs
#endif

   !!======================================================================
END MODULE bdytra
