MODULE trcini_pisces
   !!======================================================================
   !!                         ***  MODULE trcini_pisces  ***
   !! TOP :   initialisation of the PISCES biochemical model
   !!======================================================================
   !! History :    -   !  1988-07  (E. Maier-Reiner) Original code
   !!              -   !  1999-10  (O. Aumont, C. Le Quere)
   !!              -   !  2002     (O. Aumont)  PISCES
   !!             1.0  !  2005-03  (O. Aumont, A. El Moussaoui) F90
   !!             2.0  !  2007-12  (C. Ethe, G. Madec) from trcini.pisces.h90
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   'key_pisces'                                       PISCES bio-model
   !!----------------------------------------------------------------------
   !! trc_ini_pisces   : PISCES biochemical model initialisation
   !!----------------------------------------------------------------------
   USE par_trc         ! TOP parameters
   USE sms_pisces      ! Source Minus Sink variables
   USE trc
   USE oce_trc         ! ocean variables
   USE p4zche 
   USE p4zche          ! 
   USE p4zsink         ! 
   USE p4zopt          ! 
   USE p4zprod         !
   USE p4zrem          ! 
   USE p4zsed          ! 
   USE p4zflx          ! 

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_ini_pisces   ! called by trcini.F90 module

   REAL(wp) :: sco2   =  2.312e-3_wp
   REAL(wp) :: alka0  =  2.423e-3_wp
   REAL(wp) :: oxyg0  =  177.6e-6_wp 
   REAL(wp) :: po4    =  2.174e-6_wp 
   REAL(wp) :: bioma0 =  1.000e-8_wp  
   REAL(wp) :: silic1 =  91.65e-6_wp  
   REAL(wp) :: no3    =  31.04e-6_wp * 7.6_wp

#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcini_pisces.F90 2715 2011-03-30 15:58:35Z rblod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_ini_pisces
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE trc_ini_pisces ***
      !!
      !! ** Purpose :   Initialisation of the PISCES biochemical model
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_ini_pisces :   PISCES biochemical model initialisation'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'

      CALL pisces_alloc()                          ! Allocate PISCES arrays

      !                                            ! Time-step
      rfact   = rdttrc(1)                          ! ---------
      rfactr  = 1. / rfact
      rfact2  = rfact / FLOAT( nrdttrc )
      rfact2r = 1. / rfact2

      IF(lwp) WRITE(numout,*) '    Passive Tracer  time step    rfact  = ', rfact, ' rdt = ', rdttra(1)
      IF(lwp) write(numout,*) '    PISCES  Biology time step    rfact2 = ', rfact2



      ! Set biological ratios
      ! ---------------------
      rno3   = (16.+2.) / 122.
      po4r   =   1.e0   / 122.
      o2nit  =  32.     / 122.
      rdenit =  97.6    /  16.
      o2ut   = 140.     / 122.

      CALL p4z_che        ! initialize the chemical constants

      ! Initialization of tracer concentration in case of  no restart 
      !--------------------------------------------------------------
      IF( .NOT. ln_rsttr ) THEN  
         
         trn(:,:,:,jpdic) = sco2
         trn(:,:,:,jpdoc) = bioma0
         trn(:,:,:,jptal) = alka0
         trn(:,:,:,jpoxy) = oxyg0
         trn(:,:,:,jpcal) = bioma0
         trn(:,:,:,jppo4) = po4 / po4r
         trn(:,:,:,jppoc) = bioma0
#  if ! defined key_kriest
         trn(:,:,:,jpgoc) = bioma0
         trn(:,:,:,jpbfe) = bioma0 * 5.e-6
#  else
         trn(:,:,:,jpnum) = bioma0 / ( 6. * xkr_massp )
#  endif
         trn(:,:,:,jpsil) = silic1
         trn(:,:,:,jpbsi) = bioma0 * 0.15
         trn(:,:,:,jpdsi) = bioma0 * 5.e-6
         trn(:,:,:,jpphy) = bioma0
         trn(:,:,:,jpdia) = bioma0
         trn(:,:,:,jpzoo) = bioma0
         trn(:,:,:,jpmes) = bioma0
         trn(:,:,:,jpfer) = 0.6E-9
         trn(:,:,:,jpsfe) = bioma0 * 5.e-6
         trn(:,:,:,jpdfe) = bioma0 * 5.e-6
         trn(:,:,:,jpnfe) = bioma0 * 5.e-6
         trn(:,:,:,jpnch) = bioma0 * 12. / 55.
         trn(:,:,:,jpdch) = bioma0 * 12. / 55.
         trn(:,:,:,jpno3) = no3
         trn(:,:,:,jpnh4) = bioma0

         ! initialize the half saturation constant for silicate
         ! ----------------------------------------------------
         xksi(:,:)    = 2.e-6
         xksimax(:,:) = xksi(:,:)

      ENDIF

      IF(lwp) WRITE(numout,*) 'Initialization of PISCES tracers done'
      IF(lwp) WRITE(numout,*) ' '
      !
   END SUBROUTINE trc_ini_pisces


   SUBROUTINE pisces_alloc
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE pisces_alloc  ***
      !!
      !! ** Purpose :   Allocate all the dynamic arrays of PISCES 
      !!----------------------------------------------------------------------
      USE p4zint , ONLY : p4z_int_alloc      
      USE p4zsink, ONLY : p4z_sink_alloc      
      USE p4zopt , ONLY : p4z_opt_alloc           
      USE p4zprod, ONLY : p4z_prod_alloc         
      USE p4zrem , ONLY : p4z_rem_alloc           
      USE p4zsed , ONLY : p4z_sed_alloc          
      USE p4zflx , ONLY : p4z_flx_alloc
      !
      INTEGER :: ierr
      !!----------------------------------------------------------------------
      !
      ierr =         sms_pisces_alloc()          ! Start of PISCES-related alloc routines...
      ierr = ierr +     p4z_che_alloc()
      ierr = ierr +     p4z_int_alloc()
      ierr = ierr +    p4z_sink_alloc()
      ierr = ierr +     p4z_opt_alloc()
      ierr = ierr +    p4z_prod_alloc()
      ierr = ierr +     p4z_rem_alloc()
      ierr = ierr +     p4z_sed_alloc()
      ierr = ierr +     p4z_flx_alloc()
      !
      IF( lk_mpp    )   CALL mpp_sum( ierr )
      IF( ierr /= 0 )   CALL ctl_stop( 'STOP', 'pisces_alloc: unable to allocate PISCES arrays' )
      !
   END SUBROUTINE pisces_alloc

#else
   !!----------------------------------------------------------------------
   !!   Dummy module                            No PISCES biochemical model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_ini_pisces             ! Empty routine
   END SUBROUTINE trc_ini_pisces
#endif

   !!======================================================================
END MODULE trcini_pisces
