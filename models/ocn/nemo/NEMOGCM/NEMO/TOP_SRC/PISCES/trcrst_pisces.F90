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
      IF( lk_dtatrc .AND. ln_pisclo ) CALL pis_dmp_clo  ! restoring of nutrients on close seas
      IF( ln_pisdmp )                 CALL pis_dmp_ini  ! relaxation of some tracers
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_rst_read_pisces : Read specific variables from pisces model '
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'
      ! 
      IF( iom_varid( knum, 'PH', ldstop = .FALSE. ) > 0 ) THEN
         CALL iom_get( knum, jpdom_autoglo, 'PH' , hi(:,:,:)  )
      ELSE
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

   SUBROUTINE pis_dmp_ini
      !!----------------------------------------------------------------------
      !!                    ***  pis_dmp_ini  ***
      !!
      !! ** purpose  : Relaxation of some tracers
      !!----------------------------------------------------------------------
      REAL(wp) ::  alkmean = 2426.     ! mean value of alkalinity ( Glodap ; for Goyet 2391. )
      REAL(wp) ::  po4mean = 2.165     ! mean value of phosphates
      REAL(wp) ::  no3mean = 30.90     ! mean value of nitrate
      REAL(wp) ::  silmean = 91.51     ! mean value of silicate

      REAL(wp) :: zarea, zalksum, zpo4sum, zno3sum, zsilsum


      IF(lwp)  WRITE(numout,*)

      IF( cp_cfg == "orca" .AND. .NOT. lk_c1d ) THEN      ! ORCA condiguration (not 1D) !
         !                                                    ! --------------------------- !
         ! set total alkalinity, phosphate, nitrate & silicate

         zarea   = 1. / areatot * 1.e6
# if defined key_degrad
         zalksum = glob_sum( trn(:,:,:,jptal) * cvol(:,:,:) * facvol(:,:,:) ) * zarea
         zpo4sum = glob_sum( trn(:,:,:,jppo4) * cvol(:,:,:) * facvol(:,:,:) ) * zarea / 122.
         zno3sum = glob_sum( trn(:,:,:,jpno3) * cvol(:,:,:) * facvol(:,:,:) ) * zarea / 7.6
         zsilsum = glob_sum( trn(:,:,:,jpsil) * cvol(:,:,:) * facvol(:,:,:) ) * zarea
# else
         zalksum = glob_sum( trn(:,:,:,jptal) * cvol(:,:,:)  ) * zarea
         zpo4sum = glob_sum( trn(:,:,:,jppo4) * cvol(:,:,:)  ) * zarea / 122.
         zno3sum = glob_sum( trn(:,:,:,jpno3) * cvol(:,:,:)  ) * zarea / 7.6
         zsilsum = glob_sum( trn(:,:,:,jpsil) * cvol(:,:,:)  ) * zarea
# endif

         IF(lwp) WRITE(numout,*) '       TALK mean : ', zalksum
         trn(:,:,:,jptal) = trn(:,:,:,jptal) * alkmean / zalksum
            
         IF(lwp) WRITE(numout,*) '       PO4  mean : ', zpo4sum
         trn(:,:,:,jppo4) = trn(:,:,:,jppo4) * po4mean / zpo4sum

         IF(lwp) WRITE(numout,*) '       NO3  mean : ', zno3sum
         trn(:,:,:,jpno3) = trn(:,:,:,jpno3) * no3mean / zno3sum

         IF(lwp) WRITE(numout,*) '       SiO3 mean : ', zsilsum
         trn(:,:,:,jpsil) = MIN( 400.e-6,trn(:,:,:,jpsil) * silmean / zsilsum )
         !
      ENDIF

!#if defined key_kriest
!     !! Initialize number of particles from a standart restart file
!     !! The name of big organic particles jpgoc has been only change
!     !! and replace by jpnum but the values here are concentration
!     trn(:,:,:,jppoc) = trn(:,:,:,jppoc) + trn(:,:,:,jpnum) 
!     trn(:,:,:,jpnum) = trn(:,:,:,jppoc) / ( 6. * xkr_massp )
!#endif

   END SUBROUTINE pis_dmp_ini

   SUBROUTINE pis_dmp_clo   
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE pis_dmp_clo  ***
      !!
      !! ** Purpose :   Closed sea domain initialization
      !!
      !! ** Method  :   if a closed sea is located only in a model grid point
      !!                we restore to initial data
      !!
      !! ** Action  :   ictsi1(), ictsj1() : south-west closed sea limits (i,j)
      !!                ictsi2(), ictsj2() : north-east Closed sea limits (i,j)
      !!----------------------------------------------------------------------
      INTEGER, PARAMETER           ::   npicts   = 4       !: number of closed sea
      INTEGER, DIMENSION(npicts)   ::   ictsi1, ictsj1     !: south-west closed sea limits (i,j)
      INTEGER, DIMENSION(npicts)   ::   ictsi2, ictsj2     !: north-east closed sea limits (i,j)
      INTEGER :: ji, jj, jk, jn, jc            ! dummy loop indices
      !!----------------------------------------------------------------------

      IF(lwp) WRITE(numout,*) 
      IF(lwp) WRITE(numout,*)' pis_dmp_clo : closed seas '
      IF(lwp) WRITE(numout,*)'~~~~~~~'

      ! initial values
      ictsi1(:) = 1  ;  ictsi2(:) = 1
      ictsj1(:) = 1  ;  ictsj2(:) = 1

      ! set the closed seas (in data domain indices)
      ! -------------------

      IF( cp_cfg == "orca" ) THEN
         !
         SELECT CASE ( jp_cfg )
         !                                           ! =======================
         CASE ( 2 )                                  !  ORCA_R2 configuration
            !                                        ! =======================
            !                                            ! Caspian Sea
            ictsi1(1)   =  11  ;  ictsj1(1)   = 103
            ictsi2(1)   =  17  ;  ictsj2(1)   = 112
            !                                            ! Great North American Lakes
            ictsi1(2)   =  97  ;  ictsj1(2)   = 107
            ictsi2(2)   = 103  ;  ictsj2(2)   = 111
            !                                            ! Black Sea 1 : west part of the Black Sea
            ictsi1(3)   = 174  ; ictsj1(3)   = 107
            ictsi2(3)   = 181  ; ictsj2(3)   = 112
            !                                            ! Black Sea 2 : est part of the Black Sea
            ictsi1(4)   =   2  ;  ictsj1(4)   = 107
            ictsi2(4)   =   6  ;  ictsj2(4)   = 112
            !                                        ! =======================
         CASE ( 4 )                                  !  ORCA_R4 configuration
            !                                        ! =======================
            !                                            ! Caspian Sea
            ictsi1(1)   =  4  ;  ictsj1(1)   = 53
            ictsi2(1)   =  4  ;  ictsj2(1)   = 56
            !                                            ! Great North American Lakes
            ictsi1(2)   = 49  ;  ictsj1(2)   = 55
            ictsi2(2)   = 51  ;  ictsj2(2)   = 56
            !                                            ! Black Sea
            ictsi1(3)   = 88  ;  ictsj1(3)   = 55
            ictsi2(3)   = 91  ;  ictsj2(3)   = 56
            !                                            ! Baltic Sea
            ictsi1(4)   = 75  ;  ictsj1(4)   = 59
            ictsi2(4)   = 76  ;  ictsj2(4)   = 61
            !                                        ! =======================
            !                                        ! =======================
         CASE ( 025 )                                ! ORCA_R025 configuration
            !                                        ! =======================
                                                     ! Caspian + Aral sea
            ictsi1(1)   = 1330 ; ictsj1(1)   = 645
            ictsi2(1)   = 1400 ; ictsj2(1)   = 795
            !                                        ! Azov Sea
            ictsi1(2)   = 1284 ; ictsj1(2)   = 722
            ictsi2(2)   = 1304 ; ictsj2(2)   = 747
            !
         END SELECT
         !
      ENDIF

      ! convert the position in local domain indices
      ! --------------------------------------------
      DO jc = 1, npicts
         ictsi1(jc)   = mi0( ictsi1(jc) )
         ictsj1(jc)   = mj0( ictsj1(jc) )

         ictsi2(jc)   = mi1( ictsi2(jc) )
         ictsj2(jc)   = mj1( ictsj2(jc) )
      END DO

#if defined key_dtatrc
      ! Restore close seas values to initial data
      CALL trc_dta( nit000 ) 
      DO jn = 1, jptra
         IF( lutini(jn) ) THEN
            DO jc = 1, npicts
               DO jk = 1, jpkm1
                  DO jj = ictsj1(jc), ictsj2(jc)
                     DO ji = ictsi1(jc), ictsi2(jc)
                        trn(ji,jj,jk,jn) = trdta(ji,jj,jk,jn) * tmask(ji,jj,jk) 
                        trb(ji,jj,jk,jn) = trn(ji,jj,jk,jn)
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
      ENDDO
#endif
   !
   END SUBROUTINE pis_dmp_clo

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
