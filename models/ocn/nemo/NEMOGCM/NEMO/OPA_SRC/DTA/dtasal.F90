MODULE dtasal
   !!======================================================================
   !!                     ***  MODULE  dtasal  ***
   !! Ocean data  :  read ocean salinity data from monthly atlas data
   !!=====================================================================
   !! History :  OPA  ! 1991-03  ()  Original code
   !!             -   ! 1992-07  (M. Imbard)
   !!            8.0  ! 1999-10  (M.A. Foujols, M. Imbard)  NetCDF FORMAT 
   !!   NEMO     1.0  ! 2002-06  (G. Madec)  F90: Free form and module 
   !!            3.3  ! 2010-10  (C. Bricaud, S. Masson)  use of fldread
   !!----------------------------------------------------------------------
#if defined key_dtasal   ||   defined key_esopa
   !!----------------------------------------------------------------------
   !!   'key_dtasal'                                          salinity data
   !!----------------------------------------------------------------------
   !!   dta_sal      : read ocean salinity data
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE fldread         ! read input fields
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! MPP library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dta_sal    ! called by step.F90 and inidta.F90

   LOGICAL , PUBLIC, PARAMETER                     ::   lk_dtasal = .TRUE. !: salinity data flag
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) ::   s_dta              !: salinity data at given time-step

   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_sal   ! structure of input SST (file informations, fields read)

   !! * Substitutions
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 4.0 , NEMO Consortium (2011)
   !! $Id: dtasal.F90 2715 2011-03-30 15:58:35Z rblod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dta_sal( kt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE dta_sal  ***
      !!        
      !! ** Purpose :   Reads monthly salinity data
      !!             
      !! ** Method  : - Read on unit numsdt the monthly salinity data interpo-
      !!                lated onto the model grid.
      !!              - At each time step, a linear interpolation is applied
      !!                between two monthly values.
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time step
      !
      INTEGER ::   ji, jj, jk, jl, jkk       ! local loop indicies
      INTEGER ::   ik, ierr, ierr0, ierr1, ierr2   ! local integers
#if defined key_tradmp
      INTEGER ::   il0, il1, ii0, ii1, ij0, ij1   ! local integers
#endif
      REAL(wp)::   zl
      REAL(wp), DIMENSION(jpk) :: zsaldta         ! auxiliary array for interpolation
      CHARACTER(len=100)       :: cn_dir          ! Root directory for location of ssr files
      TYPE(FLD_N)              :: sn_sal
      LOGICAL , SAVE           :: linit_sal = .FALSE.
      !!
      NAMELIST/namdta_sal/   cn_dir, sn_sal
      !!----------------------------------------------------------------------
     
      ! 1. Initialization
      ! -----------------------
     
      IF( kt == nit000 .AND. ( .NOT. linit_sal ) ) THEN
        
         !                         ! set file information
         cn_dir = './'             ! directory in which the model is executed
         ! ... default values (NB: frequency positive => hours, negative => months)
         !            !   file    ! frequency ! variable ! time intep !  clim   ! 'yearly' or ! weights  ! rotation   !
         !            !   name    !  (hours)  !  name    !   (T/F)    !  (T/F)  !  'monthly'  ! filename ! pairs      !
         sn_sal = FLD_N( 'salinity',  -1.     ,'vosaline',  .false.   , .true.  ,  'monthly'  , ''       , ''         )

         REWIND ( numnam )         ! read in namlist namdta_sal 
         READ( numnam, namdta_sal ) 

         IF(lwp) THEN              ! control print
            WRITE(numout,*)
            WRITE(numout,*) 'dta_sal : Salinity Climatology '
            WRITE(numout,*) '~~~~~~~ '
         ENDIF

                                   ! Allocate salinity data array 
                                ALLOCATE( s_dta(jpi,jpj,jpk)           , STAT=ierr  )
         IF( ierr > 0              )   CALL ctl_stop( 'STOP', 'dta_sal: unable to allocate s_dta array' )
                                   ! Allocate sf_tem structure
                                ierr2 = 0
                                ALLOCATE( sf_sal(1)                    , STAT=ierr0 )
                                ALLOCATE( sf_sal(1)%fnow(jpi,jpj,jpk)  , STAT=ierr1 )
         IF( sn_sal%ln_tint )   ALLOCATE( sf_sal(1)%fdta(jpi,jpj,jpk,2), STAT=ierr2 )
         IF( ierr0+ierr1+ierr2 > 0 )   CALL ctl_stop( 'STOP', 'dta_sal: unable to allocate sf_sal structure' )
         !                         ! fill sf_sal with sn_sal and control print
         CALL fld_fill( sf_sal, (/ sn_sal /), cn_dir, 'dta_sal', 'Salinity data', 'namdta_sal' )
         linit_sal = .TRUE.        
      ENDIF
     
      ! 2. Read monthly file
      ! -------------------
     
      CALL fld_read( kt, 1, sf_sal )

      IF( lwp .AND. kt == nit000 ) THEN
         WRITE(numout,*)
         WRITE(numout,*) ' read Levitus salinity ok'
         WRITE(numout,*)
      ENDIF

#if defined key_tradmp
      IF( cp_cfg == "orca"  .AND. jp_cfg == 2 ) THEN     !  ORCA_R2 configuration
         !
         ij0 = 101   ;   ij1 = 109
         ii0 = 141   ;   ii1 = 155   
         DO jj = mj0(ij0), mj1(ij1)                  ! Reduced salinity in the Alboran Sea
            DO ji = mi0(ii0), mi1(ii1)
               sf_sal(1)%fnow(ji,jj,13:13) = sf_sal(1)%fnow(ji,jj,13:13) - 0.15
               sf_sal(1)%fnow(ji,jj,14:15) = sf_sal(1)%fnow(ji,jj,14:15) - 0.25
               sf_sal(1)%fnow(ji,jj,16:17) = sf_sal(1)%fnow(ji,jj,16:17) - 0.30
               sf_sal(1)%fnow(ji,jj,18:25) = sf_sal(1)%fnow(ji,jj,18:25) - 0.35
            END DO
         END DO
         !
         IF( nn_cla == 1 ) THEN 
            !                                         ! New salinity profile at Gibraltar
            il0 = 138   ;   il1 = 138   
            ij0 = 101   ;   ij1 = 102
            ii0 = 139   ;   ii1 = 139   
            DO jl = mi0(il0), mi1(il1)
               DO jj = mj0(ij0), mj1(ij1)
                  DO ji = mi0(ii0), mi1(ii1)
                        sf_sal(1)%fnow(ji,jj,:) = sf_sal(1)%fnow(jl,jj,:)
                  END DO
               END DO
            END DO
            !                                         ! New salinity profile at Bab el Mandeb
            il0 = 164   ;   il1 = 164   
            ij0 =  87   ;   ij1 =  88
            ii0 = 161   ;   ii1 = 163   
            DO jl = mi0(il0), mi1(il1)
               DO jj = mj0(ij0), mj1(ij1)
                  DO ji = mi0(ii0), mi1(ii1)
                     sf_sal(1)%fnow(ji,jj,:) = sf_sal(1)%fnow(jl,jj,:)
                  END DO
               END DO
            END DO
            !
         ENDIF
            !
      ENDIF
#endif   
        
      s_dta(:,:,:)=sf_sal(1)%fnow(:,:,:)
        
      IF( ln_sco ) THEN
         DO jj = 1, jpj                  ! interpolation of salinites
            DO ji = 1, jpi
               DO jk = 1, jpk
                  zl=fsdept_0(ji,jj,jk)
                  IF(zl < gdept_0(1)  ) zsaldta(jk) =  s_dta(ji,jj,1    ) 
                  IF(zl > gdept_0(jpk)) zsaldta(jk) =  s_dta(ji,jj,jpkm1) 
                  DO jkk = 1, jpkm1
                     IF((zl-gdept_0(jkk))*(zl-gdept_0(jkk+1)).le.0.0) THEN
                          zsaldta(jk) = s_dta(ji,jj,jkk)                                 &
                                     &           + (zl-gdept_0(jkk))/(gdept_0(jkk+1)-gdept_0(jkk))      &
                                     &                              *(s_dta(ji,jj,jkk+1) - s_dta(ji,jj,jkk))
                     ENDIF
                  END DO
               END DO
               DO jk = 1, jpkm1
                  s_dta(ji,jj,jk) = zsaldta(jk) 
               END DO
               s_dta(ji,jj,jpk) = 0.0 
            END DO
         END DO
           
         IF( lwp .AND. kt==nn_it000 ) THEN
            WRITE(numout,*)
            WRITE(numout,*) ' Levitus salinity data interpolated to s-coordinate'
            WRITE(numout,*)
         ENDIF

      ELSE
         !                                  ! Mask
         s_dta(:,:,:) = s_dta(:,:,:) * tmask(:,:,:)
         s_dta(:,:,jpk) = 0.e0
         IF( ln_zps ) THEN               ! z-coord. partial steps
            DO jj = 1, jpj               ! interpolation of salinity at the last ocean level (i.e. the partial step)
               DO ji = 1, jpi
                  ik = mbkt(ji,jj)
                  IF( ik > 1 ) THEN
                     zl = ( gdept_0(ik) - fsdept_0(ji,jj,ik) ) / ( gdept_0(ik) - gdept_0(ik-1) )
                     s_dta(ji,jj,ik) = (1.-zl) * s_dta(ji,jj,ik) + zl * s_dta(ji,jj,ik-1)
                  ENDIF
               END DO
            END DO
         ENDIF
      ENDIF
        
      IF( lwp .AND. kt == nit000 ) THEN
         WRITE(numout,*)' salinity Levitus '
         WRITE(numout,*)
         WRITE(numout,*)'  level = 1'
         CALL prihre(s_dta(:,:,1),    jpi,jpj,1,jpi,20,1,jpj,20,1.,numout)
         WRITE(numout,*)'  level = ',jpk/2
         CALL prihre(s_dta(:,:,jpk/2),jpi,jpj,1,jpi,20,1,jpj,20,1.,numout)           
         WRITE(numout,*) '  level = ',jpkm1
         CALL prihre(s_dta(:,:,jpkm1),jpi,jpj,1,jpi,20,1,jpj,20,1.,numout)
      ENDIF
      !
   END SUBROUTINE dta_sal

#else
   !!----------------------------------------------------------------------
   !!   Default option:                                    NO salinity data
   !!----------------------------------------------------------------------
   LOGICAL , PUBLIC, PARAMETER ::   lk_dtasal = .FALSE.   !: salinity data flag
CONTAINS
   SUBROUTINE dta_sal( kt )        ! Empty routine
      WRITE(*,*) 'dta_sal: You should not have seen this print! error?', kt
   END SUBROUTINE dta_sal
#endif
   !!======================================================================
END MODULE dtasal
