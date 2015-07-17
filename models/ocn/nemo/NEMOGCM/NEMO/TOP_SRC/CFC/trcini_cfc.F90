MODULE trcini_cfc
   !!======================================================================
   !!                         ***  MODULE trcini_cfc  ***
   !! TOP :   initialisation of the CFC tracers
   !!======================================================================
   !! History :   2.0  !  2007-12  (C. Ethe, G. Madec) from trcini.cfc.h90
   !!----------------------------------------------------------------------
#if defined key_cfc
   !!----------------------------------------------------------------------
   !!   'key_cfc'                                               CFC tracers
   !!----------------------------------------------------------------------
   !! trc_ini_cfc      : CFC model initialisation
   !!----------------------------------------------------------------------
   USE oce_trc         ! Ocean variables
   USE par_trc         ! TOP parameters
   USE trc             ! TOP variables
   USE trcsms_cfc      ! CFC sms trends

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_ini_cfc   ! called by trcini.F90 module

   CHARACTER (len=34) ::   clname = 'cfc1112.atm'   ! ???

   INTEGER  ::   inum                   ! unit number
   REAL(wp) ::   ylats = -10.           ! 10 degrees south
   REAL(wp) ::   ylatn =  10.           ! 10 degrees north

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcini_cfc.F90 2715 2011-03-30 15:58:35Z rblod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_ini_cfc
      !!----------------------------------------------------------------------
      !!                     ***  trc_ini_cfc  ***  
      !!
      !! ** Purpose :   initialization for cfc model
      !!
      !! ** Method  : - Read the namcfc namelist and check the parameter values
      !!----------------------------------------------------------------------
      INTEGER  ::  ji, jj, jn, jl, jm, js
      REAL(wp) ::  zyy, zyd
      !!----------------------------------------------------------------------

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_ini_cfc: initialisation of CFC chemical model'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~'

      !                                ! Allocate CFC arrays
      IF( trc_sms_cfc_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'trc_ini_cfc: unable to allocate CFC arrays' )


      ! Initialization of boundaries conditions
      ! --------------------------------------- 
      xphem (:,:)    = 0._wp
      p_cfc(:,:,:)   = 0._wp
      
      ! Initialization of qint in case of  no restart 
      !----------------------------------------------
      qtr_cfc(:,:,:) = 0._wp
      IF( .NOT. ln_rsttr ) THEN    
         IF(lwp) THEN
            WRITE(numout,*)
            WRITE(numout,*) 'Initialization de qint ; No restart : qint equal zero '
         ENDIF
         qint_cfc(:,:,:) = 0._wp
         DO jl = 1, jp_cfc
            jn = jp_cfc0 + jl - 1
            trn(:,:,:,jn) = 0._wp
         END DO
      ENDIF


      !   READ CFC partial pressure atmospheric value :
      !     p11(year,nt) = PCFC11  in northern (1) and southern (2) hemisphere 
      !     p12(year,nt) = PCFC12  in northern (1) and southern (2) hemisphere 
      !--------------------------------------------------------------------

      IF(lwp) WRITE(numout,*) 'read of formatted file cfc1112atm'
      
      CALL ctl_opn( inum, clname, 'OLD', 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )
      REWIND(inum)
      
      DO jm = 1, 6        ! Skip over 1st six descriptor lines
         READ(inum,'(1x)')
      END DO
   
      ! file starts in 1931 do jn represent the year in the century.jhh
      ! Read file till the end
      jn = 31
      DO WHILE ( 1 /= 2 )
         READ(inum,*,END=100) zyy, p_cfc(jn,1,1), p_cfc(jn,1,2), p_cfc(jn,2,1), p_cfc(jn,2,2)
         IF ( lwp) THEN
           WRITE(numout,'(f7.2, 4f8.2)' ) &
            &         zyy, p_cfc(jn,1,1), p_cfc(jn,1,2), p_cfc(jn,2,1), p_cfc(jn,2,2)
         ENDIF
         jn = jn + 1
      END DO
 100  npyear = jn - 1
      IF ( lwp) WRITE(numout,*) '    ', npyear ,' years read'

!      p_cfc(32,1:2,1) = 5.e-4      ! modify the values of the first years
!      p_cfc(33,1:2,1) = 8.e-4
!      p_cfc(34,1:2,1) = 1.e-6
!      p_cfc(35,1:2,1) = 2.e-3
!      p_cfc(36,1:2,1) = 4.e-3
!      p_cfc(37,1:2,1) = 6.e-3
!      p_cfc(38,1:2,1) = 8.e-3
!      p_cfc(39,1:2,1) = 1.e-2
      
      IF(lwp) THEN        ! Control print
         WRITE(numout,*)
         WRITE(numout,*) ' Year   p11NH    p11SH    p12NH    p12SH '
         DO jn = 30, 100
            WRITE(numout, '( 1I4, 4F9.2)') jn, p_cfc(jn,1,1), p_cfc(jn,2,1), p_cfc(jn,1,2), p_cfc(jn,2,2)
         END DO
      ENDIF


      ! Interpolation factor of atmospheric partial pressure
      ! Linear interpolation between 2 hemispheric function of latitud between ylats and ylatn
      !---------------------------------------------------------------------------------------
      zyd = ylatn - ylats      
      DO jj = 1 , jpj
         DO ji = 1 , jpi
            IF(     gphit(ji,jj) >= ylatn ) THEN   ;   xphem(ji,jj) = 1.e0
            ELSEIF( gphit(ji,jj) <= ylats ) THEN   ;   xphem(ji,jj) = 0.e0
            ELSE                                   ;   xphem(ji,jj) = ( gphit(ji,jj) - ylats) / zyd
            ENDIF
         END DO
      END DO
      !
      IF(lwp) WRITE(numout,*) 'Initialization of CFC tracers done'
      IF(lwp) WRITE(numout,*) ' '
      !
   END SUBROUTINE trc_ini_cfc
   
#else
   !!----------------------------------------------------------------------
   !!   Dummy module                                         No CFC tracers
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_ini_cfc             ! Empty routine
   END SUBROUTINE trc_ini_cfc
#endif

   !!======================================================================
END MODULE trcini_cfc
