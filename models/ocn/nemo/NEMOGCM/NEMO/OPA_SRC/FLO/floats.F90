MODULE floats
   !!======================================================================
   !!                       ***  MODULE  floats  ***
   !! Ocean floats : floats
   !!======================================================================
   !! History :  OPA  !          (CLIPPER)   original Code
   !!   NEMO     1.0  ! 2002-06  (A. Bozec)  F90, Free form and module
   !!----------------------------------------------------------------------
#if   defined key_floats   ||   defined key_esopa
   !!----------------------------------------------------------------------
   !!   'key_floats'                                     float trajectories
   !!----------------------------------------------------------------------
   !!   flo_stp   : float trajectories computation
   !!   flo_init  : initialization of float trajectories computation
   !!----------------------------------------------------------------------
   USE oce             ! ocean variables
   USE flo_oce         ! floats variables
   USE lib_mpp         ! distributed memory computing
   USE flodom          ! initialisation Module 
   USE flowri          ! float output                     (flo_wri routine)
   USE flo4rk          ! Trajectories, Runge Kutta scheme (flo_4rk routine)
   USE floblk          ! Trajectories, Blanke scheme      (flo_blk routine)
   USE in_out_manager  ! I/O manager

   IMPLICIT NONE
   PRIVATE  

   PUBLIC   flo_stp    ! routine called by step.F90
   PUBLIC   flo_init   ! routine called by opa.F90

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: floats.F90 2715 2011-03-30 15:58:35Z rblod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE flo_stp( kt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE flo_stp  ***
      !!                    
      !! ** Purpose :   Compute the geographical position (lat., long., depth)
      !!      of each float at each time step with one of the algorithm.
      !! 
      !! ** Method  :   The position of a float is computed with Bruno Blanke 
      !!        algorithm by default and with a 4th order Runge-Kutta scheme
      !!        if ln_flork4 =T
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in  ) ::   kt   ! ocean time step
      !!----------------------------------------------------------------------
      !
      IF( ln_flork4 ) THEN   ;   CALL flo_4rk( kt )        ! Trajectories using a 4th order Runge Kutta scheme
      ELSE                   ;   CALL flo_blk( kt )        ! Trajectories using Blanke' algorithme
      ENDIF
      !
      IF( lk_mpp )   CALL mppsync   ! synchronization of all the processor
      !
      IF( kt == nit000 .OR. MOD( kt, nn_writefl ) == 0 )   CALL flo_wri( kt )      ! trajectories file 
      IF( kt == nitend .OR. MOD( kt, nn_stockfl ) == 0 )   CALL flo_wri( kt )      ! restart file 
      !
      wb(:,:,:) = wn(:,:,:)         ! Save the old vertical velocity field
      !
   END SUBROUTINE flo_stp


   SUBROUTINE flo_init
      !!----------------------------------------------------------------
      !!                 ***  ROUTINE flo_init  ***
      !!                   
      !! ** Purpose :   Read the namelist of floats
      !!----------------------------------------------------------------------
      NAMELIST/namflo/ ln_rstflo, nn_writefl, nn_stockfl, ln_argo, ln_flork4 
      !!---------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'flo_stp : call floats routine '
      IF(lwp) WRITE(numout,*) '~~~~~~~'

      REWIND( numnam )              ! Namelist namflo : floats
      READ  ( numnam, namflo )
      !
      IF(lwp) THEN                  ! control print
         WRITE(numout,*)
         WRITE(numout,*) '         Namelist floats :'
         WRITE(numout,*) '            restart                          ln_rstflo  = ', ln_rstflo
         WRITE(numout,*) '            frequency of float output file   nn_writefl = ', nn_writefl
         WRITE(numout,*) '            frequency of float restart file  nn_stockfl = ', nn_stockfl
         WRITE(numout,*) '            Argo type floats                 ln_argo    = ', ln_argo
         WRITE(numout,*) '            Computation of T trajectories    ln_flork4  = ', ln_flork4
      ENDIF
      !
      !                             ! allocate floats arrays
      IF( flo_oce_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'flo_init : unable to allocate arrays' )
      !
      !                             ! allocate flowri arrays
      IF( flo_wri_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'flo_wri : unable to allocate arrays' )
      !
      CALL flo_dom                  ! compute/read initial position of floats

      wb(:,:,:) = wn(:,:,:)         ! set wb for computation of floats trajectories at the first time step
      !
   END SUBROUTINE flo_init

#  else
   !!----------------------------------------------------------------------
   !!   Default option :                                       Empty module
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE flo_stp( kt )          ! Empty routine
      WRITE(*,*) 'flo_stp: You should not have seen this print! error?', kt
   END SUBROUTINE flo_stp
   SUBROUTINE flo_init          ! Empty routine
   END SUBROUTINE flo_init
#endif

   !!======================================================================
 END MODULE floats
