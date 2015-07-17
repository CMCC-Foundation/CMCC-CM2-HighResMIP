MODULE dtatem
   !!======================================================================
   !!                     ***  MODULE  dtatem  ***
   !! Ocean data  :  read ocean temperature data from monthly atlas data
   !!=====================================================================
   !! History :  OPA  ! 1991-03  ()  Original code
   !!             -   ! 1992-07  (M. Imbard)
   !!            8.0  ! 1999-10  (M.A. Foujols, M. Imbard)  NetCDF FORMAT 
   !!   NEMO     1.0  ! 2002-06  (G. Madec)  F90: Free form and module 
   !!            3.3  ! 2010-10  (C. Bricaud, S. Masson)  use of fldread
   !!----------------------------------------------------------------------
#if defined key_dtatem   ||   defined key_esopa
   !!----------------------------------------------------------------------
   !!   'key_dtatem'                              3D temperature data field
   !!----------------------------------------------------------------------
   !!   dta_tem      : read ocean temperature data
   !!---l-------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE fldread         ! read input fields
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! MPP library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dta_tem    ! called by step.F90 and inidta.F90

   LOGICAL , PUBLIC, PARAMETER                     ::   lk_dtatem = .TRUE. !: temperature data flag
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) ::   t_dta              !: temperature data at given time-step

   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_tem      ! structure of input SST (file informations, fields read)

   !! * Substitutions
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 4.0 , NEMO Consortium (2011)
   !! $Id: dtatem.F90 2715 2011-03-30 15:58:35Z rblod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dta_tem( kt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE dta_tem  ***
      !!                    
      !! ** Purpose :   Reads monthly temperature data 
      !! 
      !! ** Method  :   Read on unit numtdt the interpolated temperature 
      !!      onto the model grid.
      !!      Data begin at january. 
      !!      The value is centered at the middle of month.
      !!      In the opa model, kt=1 agree with january 1.
      !!      At each time step, a linear interpolation is applied between 
      !!      two monthly values.
      !!      Read on unit numtdt
      !!
      !! ** Action  :   define t_dta array at time-step kt
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step
      !
      INTEGER ::   ji, jj, jk, jl, jkk       ! dummy loop indicies
      INTEGER ::   ik, ierr, ierr0, ierr1, ierr2   ! local integers
#if defined key_tradmp
      INTEGER ::   il0, il1, ii0, ii1, ij0, ij1   ! local integers
#endif
      REAL(wp)::   zl
      REAL(wp), DIMENSION(jpk) ::   ztemdta            ! auxiliary array for interpolation
      !
      CHARACTER(len=100)       ::   cn_dir             ! Root directory for location of ssr files
      TYPE(FLD_N)              ::   sn_tem
      LOGICAL , SAVE           ::   linit_tem = .FALSE.
      !!
      NAMELIST/namdta_tem/   cn_dir, sn_tem
      !!----------------------------------------------------------------------
 
      ! 1. Initialization 
      ! -----------------------
      
      IF( kt == nit000 .AND. (.NOT. linit_tem ) ) THEN

         !                   ! set file information
         cn_dir = './'       ! directory in which the model is executed
         ! ... default values (NB: frequency positive => hours, negative => months)
         !            !   file    ! frequency ! variable  ! time intep !  clim   ! 'yearly' or ! weights  ! rotation !
         !            !   name    !  (hours)  !  name     !   (T/F)    !  (T/F)  !  'monthly'  ! filename ! pairs    !
         sn_tem = FLD_N( 'temperature',  -1.  , 'votemper',  .false.   , .true.  ,  'yearly'   , ''       , ''       )

         REWIND( numnam )          ! read in namlist namdta_tem 
         READ( numnam, namdta_tem ) 

         IF(lwp) THEN              ! control print
            WRITE(numout,*)
            WRITE(numout,*) 'dta_tem : Temperature Climatology '
            WRITE(numout,*) '~~~~~~~ '
         ENDIF

                                   ! Allocate temperature data array 
                                ALLOCATE( t_dta(jpi,jpj,jpk)           , STAT=ierr  )
         IF( ierr > 0              )   CALL ctl_stop( 'STOP', 'dta_tem: unable to allocate t_dta array' )
                                   ! Allocate sf_tem structure
                                ierr2 = 0
                                ALLOCATE( sf_tem(1)                    , STAT=ierr0 )
                                ALLOCATE( sf_tem(1)%fnow(jpi,jpj,jpk)  , STAT=ierr1 )
         IF( sn_tem%ln_tint )   ALLOCATE( sf_tem(1)%fdta(jpi,jpj,jpk,2), STAT=ierr2 )
         IF( ierr0+ierr1+ierr2 > 0 )   CALL ctl_stop( 'STOP', 'dta_tem: unable to allocate sf_tem structure' )
         !                         ! fill sf_tem with sn_tem and control print
         CALL fld_fill( sf_tem, (/ sn_tem /), cn_dir, 'dta_tem', 'Temperature data', 'namdta_tem' )
         linit_tem = .TRUE.
         !
      ENDIF
      
      ! 2. Read monthly file
      ! -------------------
         
      CALL fld_read( kt, 1, sf_tem )
       
      IF( lwp .AND. kt == nit000 )THEN 
         WRITE(numout,*)
         WRITE(numout,*) ' read Levitus temperature ok'
         WRITE(numout,*)
      ENDIF
         
#if defined key_tradmp
      IF( cp_cfg == "orca" .AND. jp_cfg == 2 ) THEN      !  ORCA_R2 configuration
         !
         ij0 = 101   ;   ij1 = 109
         ii0 = 141   ;   ii1 = 155
         DO jj = mj0(ij0), mj1(ij1)                      ! Reduced temperature in the Alboran Sea
            DO ji = mi0(ii0), mi1(ii1)
               sf_tem(1)%fnow(ji,jj, 13:13 ) = sf_tem(1)%fnow(ji,jj, 13:13 ) - 0.20
               sf_tem(1)%fnow(ji,jj, 14:15 ) = sf_tem(1)%fnow(ji,jj, 14:15 ) - 0.35  
               sf_tem(1)%fnow(ji,jj, 16:25 ) = sf_tem(1)%fnow(ji,jj, 16:25 ) - 0.40
            END DO
         END DO
         !
         IF( nn_cla == 1 ) THEN 
            !                                         ! New temperature profile at Gibraltar
            il0 = 138   ;   il1 = 138
            ij0 = 101   ;   ij1 = 102
            ii0 = 139   ;   ii1 = 139
            DO jl = mi0(il0), mi1(il1)
               DO jj = mj0(ij0), mj1(ij1)
                  DO ji = mi0(ii0), mi1(ii1)
                     sf_tem(1)%fnow(ji,jj,:) = sf_tem(1)%fnow(jl,jj,:)
                  END DO
               END DO
            END DO
            !                                         ! New temperature profile at Bab el Mandeb
            il0 = 164   ;   il1 = 164
            ij0 =  87   ;   ij1 =  88
            ii0 = 161   ;   ii1 = 163
            DO jl = mi0(il0), mi1(il1)
               DO jj = mj0(ij0), mj1(ij1)
                  DO ji = mi0(ii0), mi1(ii1)
                     sf_tem(1)%fnow(ji,jj,:) = sf_tem(1)%fnow(jl,jj,:)
                  END DO
               END DO
            END DO
         ELSE
            !                                         ! Reduced temperature at Red Sea
            ij0 =  87   ;   ij1 =  96
            ii0 = 148   ;   ii1 = 160
            sf_tem(1)%fnow( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) ,  4:10 ) = 7.0
            sf_tem(1)%fnow( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , 11:13 ) = 6.5
            sf_tem(1)%fnow( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , 14:20 ) = 6.0
         ENDIF
            !
      ENDIF
#endif
         
      t_dta(:,:,:) = sf_tem(1)%fnow(:,:,:) 
         
      IF( ln_sco ) THEN
         DO jj = 1, jpj                  ! interpolation of temperatures
            DO ji = 1, jpi
               DO jk = 1, jpk
                  zl=fsdept_0(ji,jj,jk)
                  IF(zl < gdept_0(1))   ztemdta(jk) =  t_dta(ji,jj,1)
                  IF(zl > gdept_0(jpk)) ztemdta(jk) =  t_dta(ji,jj,jpkm1) 
                  DO jkk = 1, jpkm1
                     IF((zl-gdept_0(jkk))*(zl-gdept_0(jkk+1)).le.0.0) THEN
                        ztemdta(jk) = t_dta(ji,jj,jkk)                                 &
                                  &    + (zl-gdept_0(jkk))/(gdept_0(jkk+1)-gdept_0(jkk))  &
                                  &    * (t_dta(ji,jj,jkk+1) - t_dta(ji,jj,jkk))
                     ENDIF
                  END DO
               END DO
               DO jk = 1, jpkm1
                  t_dta(ji,jj,jk) = ztemdta(jk)
               END DO
               t_dta(ji,jj,jpk) = 0.0
            END DO
         END DO
            
         IF( lwp .AND. kt == nit000 )THEN
            WRITE(numout,*)
            WRITE(numout,*) ' Levitus temperature data interpolated to s-coordinate'
            WRITE(numout,*)
         ENDIF
            
      ELSE
         !                                  ! Mask
         t_dta(:,:,:  ) = t_dta(:,:,:) * tmask(:,:,:)
         t_dta(:,:,jpk) = 0.
         IF( ln_zps ) THEN                ! z-coord. with partial steps
            DO jj = 1, jpj                ! interpolation of temperature at the last level
               DO ji = 1, jpi
                  ik = mbkt(ji,jj)
                  IF( ik > 1 ) THEN
                     zl = ( gdept_0(ik) - fsdept_0(ji,jj,ik) ) / ( gdept_0(ik) - gdept_0(ik-1) )
                     t_dta(ji,jj,ik) = (1.-zl) * t_dta(ji,jj,ik) + zl * t_dta(ji,jj,ik-1)
                  ENDIF
               END DO
            END DO
         ENDIF
         !
      ENDIF
         
      IF( lwp .AND. kt == nit000 ) THEN
         WRITE(numout,*) ' temperature Levitus '
         WRITE(numout,*)
         WRITE(numout,*)'  level = 1'
         CALL prihre( t_dta(:,:,1    ), jpi, jpj, 1, jpi, 20, 1, jpj, 20, 1., numout )
         WRITE(numout,*)'  level = ', jpk/2
         CALL prihre( t_dta(:,:,jpk/2), jpi, jpj, 1, jpi, 20, 1, jpj, 20, 1., numout )
         WRITE(numout,*)'  level = ', jpkm1
         CALL prihre( t_dta(:,:,jpkm1), jpi, jpj, 1, jpi, 20, 1, jpj, 20, 1., numout )
      ENDIF
      !
   END SUBROUTINE dta_tem

#else
   !!----------------------------------------------------------------------
   !!   Default case                           NO 3D temperature data field
   !!----------------------------------------------------------------------
   LOGICAL , PUBLIC, PARAMETER ::   lk_dtatem = .FALSE.   !: temperature data flag
CONTAINS
   SUBROUTINE dta_tem( kt )        ! Empty routine
      WRITE(*,*) 'dta_tem: You should not have seen this print! error?', kt
   END SUBROUTINE dta_tem
#endif
   !!======================================================================
END MODULE dtatem
