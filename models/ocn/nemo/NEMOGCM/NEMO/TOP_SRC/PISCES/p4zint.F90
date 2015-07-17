MODULE p4zint
   !!======================================================================
   !!                         ***  MODULE p4zint  ***
   !! TOP :   PISCES interpolation and computation of various accessory fields
   !!======================================================================
   !! History :   1.0  !  2004-03 (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   'key_pisces'                                       PISCES bio-model
   !!----------------------------------------------------------------------
   !!   p4z_int        :  interpolation and computation of various accessory fields
   !!----------------------------------------------------------------------
   USE oce_trc         !
   USE trc
   USE sms_pisces

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_int  
   PUBLIC   p4z_int_alloc

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   tgfunc    !: Temp. dependancy of various biological rates
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   tgfunc2   !: Temp. dependancy of mesozooplankton rates

   REAL(wp) ::   xksilim = 16.5e-6_wp   ! Half-saturation constant for the Si half-saturation constant computation

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: p4zint.F90 2715 2011-03-30 15:58:35Z rblod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_int
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_int  ***
      !!
      !! ** Purpose :   interpolation and computation of various accessory fields
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER  ::   ji, jj
      REAL(wp) ::   zdum
      !!---------------------------------------------------------------------

      ! Computation of phyto and zoo metabolic rate
      ! -------------------------------------------
      tgfunc (:,:,:) = EXP( 0.063913 * tsn(:,:,:,jp_tem) )
      tgfunc2(:,:,:) = EXP( 0.07608  * tsn(:,:,:,jp_tem) )

      ! Computation of the silicon dependant half saturation
      ! constant for silica uptake
      ! ---------------------------------------------------
      DO ji = 1, jpi
         DO jj = 1, jpj
            zdum = trn(ji,jj,1,jpsil) * trn(ji,jj,1,jpsil)
            xksimax(ji,jj) = MAX( xksimax(ji,jj), ( 1.+ 7.* zdum / ( xksilim * xksilim + zdum ) ) * 1e-6 )
         END DO
      END DO
      !
      IF( nday_year == nyear_len(1) ) THEN
         xksi    = xksimax
         xksimax = 0._wp
      ENDIF
      !
   END SUBROUTINE p4z_int


   INTEGER FUNCTION p4z_int_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_int_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( tgfunc(jpi,jpj,jpk), tgfunc2(jpi,jpj,jpk), STAT=p4z_int_alloc )
      !
      IF( p4z_int_alloc /= 0 )   CALL ctl_warn('p4z_int_alloc : failed to allocate arrays.')
      !
   END FUNCTION p4z_int_alloc

#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_int                   ! Empty routine
      WRITE(*,*) 'p4z_int: You should not have seen this print! error?'
   END SUBROUTINE p4z_int
#endif 

   !!======================================================================
END MODULE  p4zint
