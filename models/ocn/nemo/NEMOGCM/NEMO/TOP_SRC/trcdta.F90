MODULE trcdta
   !!======================================================================
   !!                     ***  MODULE  trcdta  ***
   !! TOP :  reads passive tracer data 
   !!=====================================================================
   !! History :   1.0  !  2002-04  (O. Aumont)  original code
   !!              -   !  2004-03  (C. Ethe)  module
   !!              -   !  2005-03  (O. Aumont, A. El Moussaoui) F90
   !!----------------------------------------------------------------------
#if  defined key_top  &&  defined key_dtatrc
   !!----------------------------------------------------------------------
   !!   'key_top'  and  'key_dtatrc'        TOP model + passive tracer data
   !!----------------------------------------------------------------------
   !!   trc_dta      : read ocean passive tracer data
   !!----------------------------------------------------------------------
   USE oce_trc
   USE par_trc
   USE trc
   USE lib_print
   USE iom

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_dta         ! called in trcini.F90 and trcdmp.F90
   PUBLIC   trc_dta_alloc   ! called in nemogcm.F90

   LOGICAL , PUBLIC, PARAMETER ::   lk_dtatrc = .TRUE.   !: temperature data flag
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) ::   trdta   !: tracer data at given time-step

   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:,:,:) ::   tracdta       ! tracer data at two consecutive times
   INTEGER , ALLOCATABLE, SAVE, DIMENSION(:) ::   nlectr      !: switch for reading once
   INTEGER , ALLOCATABLE, SAVE, DIMENSION(:) ::   ntrc1       !: number of 1st month when reading 12 monthly value
   INTEGER , ALLOCATABLE, SAVE, DIMENSION(:) ::   ntrc2       !: number of 2nd month when reading 12 monthly value

   !! * Substitutions
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcdta.F90 2715 2011-03-30 15:58:35Z rblod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_dta( kt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE trc_dta  ***
      !!
      !! ** Purpose :   Reads passive tracer data (Levitus monthly data)
      !!
      !! ** Method  :   Read on unit numtr the interpolated tracer concentra-
      !!      tion onto the global grid. Data begin at january. 
      !!      The value is centered at the middle of month. 
      !!      In the opa model, kt=1 agree with january 1. 
      !!      At each time step, a linear interpolation is applied between 
      !!      two monthly values.
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt     ! ocean time-step
      !!
      CHARACTER (len=39) ::   clname(jptra)
      INTEGER, PARAMETER ::   jpmonth = 12    ! number of months
      INTEGER ::   ji, jj, jn, jl 
      INTEGER ::   imois, iman, i15, ik  ! temporary integers 
      REAL(wp) ::   zxy, zl
!!gm HERE the daymod should be used instead of computation of month and co !!
!!gm      better in case of real calandar and leap-years !
      !!----------------------------------------------------------------------

      DO jn = 1, jptra

         IF( lutini(jn) ) THEN 

            IF ( kt == nit000 ) THEN
               !! 3D tracer data
               IF(lwp)WRITE(numout,*)
               IF(lwp)WRITE(numout,*) ' dta_trc: reading tracer' 
               IF(lwp)WRITE(numout,*) ' data file ', jn, ctrcnm(jn)
               IF(lwp)WRITE(numout,*)
               nlectr(jn) = 0
            ENDIF
            ! Initialization
            iman = jpmonth
            i15  = nday / 16
            imois = nmonth + i15 -1
            IF( imois == 0 ) imois = iman


            ! First call kt=nit000
            ! --------------------

            IF ( kt == nit000 .AND. nlectr(jn) == 0 ) THEN
               ntrc1(jn) = 0
               IF(lwp) WRITE(numout,*) ' trc_dta : Levitus tracer data monthly fields'
               ! open file 
# if defined key_pisces
               clname(jn) = 'data_1m_'//TRIM(ctrcnm(jn))//'_nomask'
# else
               clname(jn) = TRIM(ctrcnm(jn))
# endif
               CALL iom_open ( clname(jn), numtr(jn) )              

            ENDIF

# if defined key_pisces
            ! Read montly file
            IF( ( kt == nit000 .AND. nlectr(jn) == 0)  .OR. imois /= ntrc1(jn) ) THEN
               nlectr(jn) = 1

               ! Calendar computation

               ! ntrc1 number of the first file record used in the simulation
               ! ntrc2 number of the last  file record

               ntrc1(jn) = imois
               ntrc2(jn) = ntrc1(jn) + 1
               ntrc1(jn) = MOD( ntrc1(jn), iman )
               IF ( ntrc1(jn) == 0 ) ntrc1(jn) = iman
               ntrc2(jn) = MOD( ntrc2(jn), iman )
               IF ( ntrc2(jn) == 0 ) ntrc2(jn) = iman
               IF(lwp) WRITE(numout,*) 'first record file used ntrc1 ', ntrc1(jn) 
               IF(lwp) WRITE(numout,*) 'last  record file used ntrc2 ', ntrc2(jn)

               ! Read montly passive tracer data Levitus 

               CALL iom_get ( numtr(jn), jpdom_data, ctrcnm(jn), tracdta(:,:,:,jn,1), ntrc1(jn) )
               CALL iom_get ( numtr(jn), jpdom_data, ctrcnm(jn), tracdta(:,:,:,jn,2), ntrc2(jn) )

               IF(lwp) THEN
                  WRITE(numout,*)
                  WRITE(numout,*) ' read tracer data ', ctrcnm(jn),' ok'
                  WRITE(numout,*)
               ENDIF

               ! Apply Mask
               DO jl = 1, 2
                  tracdta(:,:,:  ,jn,jl) = tracdta(:,:,:,jn,jl) * tmask(:,:,:) 
                  tracdta(:,:,jpk,jn,jl) = 0.
                  IF( ln_zps ) THEN                ! z-coord. with partial steps
                     DO jj = 1, jpj                ! interpolation of temperature at the last level
                        DO ji = 1, jpi
                           ik = mbkt(ji,jj)
                           IF( ik > 2 ) THEN
                              zl = ( gdept_0(ik) - fsdept_0(ji,jj,ik) ) / ( gdept_0(ik) - gdept_0(ik-1) )
                              tracdta(ji,jj,ik,jn,jl) = (1.-zl) * tracdta(ji,jj,ik  ,jn,jl)    &
                                 &                    +     zl  * tracdta(ji,jj,ik-1,jn,jl)
                           ENDIF
                        END DO
                     END DO
                  ENDIF

               END DO

            ENDIF

            IF(lwp) THEN
               WRITE(numout,*) ctrcnm(jn), 'Levitus month ', ntrc1(jn), ntrc2(jn)
               WRITE(numout,*)
               WRITE(numout,*) ' Levitus month = ', ntrc1(jn), '  level = 1'
               CALL prihre( tracdta(1,1,1,jn,1), jpi, jpj, 1, jpi, 20, 1   &
                  &        ,jpj, 20, 1., numout )
               WRITE(numout,*) ' Levitus month = ', ntrc1(jn), '  level = ',jpk/2
               CALL prihre( tracdta(1,1,jpk/2,jn,1), jpi, jpj, 1, jpi,    &
                  &         20, 1, jpj, 20, 1., numout )
               WRITE(numout,*) ' Levitus month = ',ntrc1(jn),'  level = ',jpkm1
               CALL prihre( tracdta(1,1,jpkm1,jn,1), jpi, jpj, 1, jpi,     &
                  &         20, 1, jpj, 20, 1., numout )
            ENDIF

            ! At every time step compute temperature data
            zxy = FLOAT( nday + 15 - 30 * i15 ) / 30.
            trdta(:,:,:,jn) =  ( 1. - zxy ) * tracdta(:,:,:,jn,1)    &
               &              +       zxy   * tracdta(:,:,:,jn,2) 

            IF( jn == jpno3 )   trdta(:,:,:,jn) = trdta(:,:,:,jn) *   7.6e-6
            IF( jn == jpdic )   trdta(:,:,:,jn) = trdta(:,:,:,jn) *   1.0e-6
            IF( jn == jptal )   trdta(:,:,:,jn) = trdta(:,:,:,jn) *   1.0e-6
            IF( jn == jpoxy )   trdta(:,:,:,jn) = trdta(:,:,:,jn) *  44.6e-6
            IF( jn == jpsil )   trdta(:,:,:,jn) = trdta(:,:,:,jn) *   1.0e-6
            IF( jn == jppo4 )   trdta(:,:,:,jn) = trdta(:,:,:,jn) * 122.0e-6

            ! Close the file
            ! --------------
            
            IF( kt == nitend )   CALL iom_close( numtr(jn) )

# else
            ! Read init file only
            IF( kt == nit000  ) THEN
               ntrc1(jn) = 1
               CALL iom_get ( numtr(jn), jpdom_data, ctrcnm(jn), trdta(:,:,:,jn), ntrc1(jn) )
               trdta(:,:,:,jn) = trdta(:,:,:,jn) * tmask(:,:,:)
               CALL iom_close ( numtr(jn) )
            ENDIF 
# endif
         ENDIF

      END DO
      !
   END SUBROUTINE trc_dta


   INTEGER FUNCTION trc_dta_alloc()
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE trc_dta_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( trdta  (jpi,jpj,jpk,jptra  ) ,                    &
         &      tracdta(jpi,jpj,jpk,jptra,2) ,                    &
         &      nlectr(jptra) , ntrc1(jptra) , ntrc2(jptra) , STAT=trc_dta_alloc)
         !
      IF( trc_dta_alloc /= 0 )   CALL ctl_warn('trc_dta_alloc : failed to allocate arrays')
      !
   END FUNCTION trc_dta_alloc

#else
   !!----------------------------------------------------------------------
   !!   Dummy module                              NO 3D passive tracer data
   !!----------------------------------------------------------------------
   LOGICAL , PUBLIC, PARAMETER ::   lk_dtatrc = .FALSE.   !: temperature data flag
CONTAINS
   SUBROUTINE trc_dta( kt )        ! Empty routine
      WRITE(*,*) 'trc_dta: You should not have seen this print! error?', kt
   END SUBROUTINE trc_dta
#endif

   !!======================================================================
END MODULE trcdta
