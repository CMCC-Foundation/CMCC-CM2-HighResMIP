MODULE trcrst_pisces
   !!======================================================================
   !!                       ***  MODULE trcrst_pisces  ***
   !! TOP :   create, write, read the restart files of PISCES tracer
   !!======================================================================
   !! History :   1.0  !  2010-01 (C. Ethe) Original
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   'key_pisces'                                               pisces tracers
   !!----------------------------------------------------------------------
   !!   trc_rst_read_pisces   : read  restart file
   !!   trc_rst_wri_pisces    : write restart file
   !!----------------------------------------------------------------------
   USE oce_trc         ! Ocean variables
   USE par_trc         ! TOP parameters
   USE trc             ! TOP variables
   USE trcsms_pisces          ! pisces sms trends
   USE sms_pisces          ! pisces sms variables
   USE iom
   USE trcdta

   IMPLICIT NONE
   PRIVATE

   PUBLIC  trc_rst_read_pisces   ! called by trcini.F90 module
   PUBLIC  trc_rst_wri_pisces   ! called by trcini.F90 module

CONTAINS
   
   SUBROUTINE trc_rst_read_pisces( knum ) 
      !!----------------------------------------------------------------------
      !!                     ***  trc_rst_read_pisces  ***  
      !!
      !! ** Purpose : Read in restart file specific variables from pisces model
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in)  :: knum  ! unit of the restart file
      INTEGER  ::  ji, jj, jk
      REAL(wp) ::  zcaralk, zbicarb, zco3
      REAL(wp) ::  ztmas, ztmas1
      !!----------------------------------------------------------------------

      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_rst_read_pisces : Read specific variables from pisces model '
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'
      ! 
      IF( iom_varid( knum, 'PH', ldstop = .FALSE. ) > 0 ) THEN
         CALL iom_get( knum, jpdom_autoglo, 'PH' , hi(:,:,:)  )
      ELSE
!         hi(:,:,:) = 1.e-9 
         ! Set PH from  total alkalinity, borat (???), akb3 (???) and ak23 (???)
         ! --------------------------------------------------------
         DO jk = 1, jpk
            DO jj = 1, jpj
               DO ji = 1, jpi
                  ztmas   = tmask(ji,jj,jk)
                  ztmas1  = 1. - tmask(ji,jj,jk)
                  zcaralk = trn(ji,jj,jk,jptal) - borat(ji,jj,jk) / (  1. + 1.E-8 / ( rtrn + akb3(ji,jj,jk) )  )
                  zco3    = ( zcaralk - trn(ji,jj,jk,jpdic) ) * ztmas + 0.5e-3 * ztmas1
                  zbicarb = ( 2. * trn(ji,jj,jk,jpdic) - zcaralk )
                 hi(ji,jj,jk) = ( ak23(ji,jj,jk) * zbicarb / zco3 ) * ztmas + 1.e-9 * ztmas1
               END DO
            END DO
         END DO
      ENDIF
      CALL iom_get( knum, jpdom_autoglo, 'Silicalim', xksi(:,:) ) 
      IF( iom_varid( knum, 'Silicamax', ldstop = .FALSE. ) > 0 ) THEN
         CALL iom_get( knum, jpdom_autoglo, 'Silicamax' , xksimax(:,:)  )
      ELSE
         xksimax(:,:) = xksi(:,:)
      ENDIF

   END SUBROUTINE trc_rst_read_pisces

   SUBROUTINE trc_rst_wri_pisces( kt, kitrst, knum )
      !!----------------------------------------------------------------------
      !!                     ***  trc_rst_read_pisces  ***
      !!
      !! ** Purpose : Read in restart file specific variables from pisces model
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in)  :: kt      ! time step
      INTEGER, INTENT(in)  :: kitrst  ! time step of restart write
      INTEGER, INTENT(in)  :: knum    ! unit of the restart file
      !!----------------------------------------------------------------------

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_rst_wri_pisces : Write specific variables from pisces model '
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'

      CALL iom_rstput( kt, kitrst, knum, 'PH', hi(:,:,:) )
      CALL iom_rstput( kt, kitrst, knum, 'Silicalim', xksi(:,:) ) 
      CALL iom_rstput( kt, kitrst, knum, 'Silicamax', xksimax(:,:) )

   END SUBROUTINE trc_rst_wri_pisces


#else
   !!----------------------------------------------------------------------
   !!  Dummy module :                                     No passive tracer
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_rst_read_pisces( knum )
      INTEGER, INTENT(in)  :: knum
      WRITE(*,*) 'trc_rst_read_pisces: You should not have seen this print! error?', knum
   END SUBROUTINE trc_rst_read_pisces

   SUBROUTINE trc_rst_wri_pisces( kt, kitrst, knum )
     INTEGER, INTENT(in)  :: kt, kitrst, knum
     WRITE(*,*) 'trc_rst_wri_pisces: You should not have seen this print! error?', kt, kitrst, knum
   END SUBROUTINE trc_rst_wri_pisces
#endif

   !!======================================================================
END MODULE trcrst_pisces
