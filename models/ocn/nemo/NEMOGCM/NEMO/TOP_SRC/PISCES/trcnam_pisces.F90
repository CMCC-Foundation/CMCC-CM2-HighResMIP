MODULE trcnam_pisces
   !!======================================================================
   !!                      ***  MODULE trcnam_lobster  ***
   !! TOP :   initialisation of some run parameters for PISCES bio-model
   !!======================================================================
   !! History :    -   !  1999-10 (M.A. Foujols, M. Levy) original code
   !!              -   !  2000-01 (L. Bopp) hamocc3, p3zd
   !!             1.0  !  2003-08 (C. Ethe)  module F90
   !!             2.0  !  2007-12  (C. Ethe, G. Madec) from trcnam.pisces.h90
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   'key_pisces'   :                                   PISCES bio-model
   !!----------------------------------------------------------------------
   !! trc_nam_pisces       : PISCES model namelist read
   !!----------------------------------------------------------------------
   USE oce_trc         ! Ocean variables
   USE par_trc         ! TOP parameters
   USE trc             ! TOP variables
   USE sms_pisces      ! sms trends


   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_nam_pisces   ! called by trcnam.F90 module


   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcnam_pisces.F90 2715 2011-03-30 15:58:35Z rblod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trc_nam_pisces
      !!----------------------------------------------------------------------
      !!                     ***  trc_nam_pisces  ***  
      !!
      !! ** Purpose :   read PISCES namelist
      !!
      !! ** input   :   file 'namelist.trc.sms' containing the following
      !!             namelist: natext, natbio, natsms
      !!                       natkriest ("key_kriest")
      !!----------------------------------------------------------------------
      !!
#if defined key_diatrc && ! defined key_iomput
      INTEGER ::  jl, jn
      ! definition of additional diagnostic as a structure
      TYPE DIAG
         CHARACTER(len = 20)  :: snamedia   !: short name
         CHARACTER(len = 80 ) :: lnamedia   !: long name
         CHARACTER(len = 20 ) :: unitdia    !: unit
      END TYPE DIAG

      TYPE(DIAG) , DIMENSION(jp_pisces_2d) :: pisdia2d
      TYPE(DIAG) , DIMENSION(jp_pisces_3d) :: pisdia3d
#endif

      NAMELIST/nampisbio/ part, nrdttrc, wsbio, xkmort, ferat3, wsbio2
#if defined key_kriest
      NAMELIST/nampiskrp/ xkr_eta, xkr_zeta, xkr_mass_min, xkr_mass_max
#endif
#if defined key_diatrc && ! defined key_iomput
      NAMELIST/nampisdia/ nn_writedia, pisdia3d, pisdia2d     ! additional diagnostics
#endif
      NAMELIST/nampisdmp/ ln_pisdmp, ln_pisclo

      !!----------------------------------------------------------------------

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_nam_pisces : read PISCES namelists'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'


      !                               ! Open the namelist file
      !                               ! ----------------------
      CALL ctl_opn( numnat, 'namelist_pisces', 'OLD', 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )

      REWIND( numnat )                    
      READ  ( numnat, nampisbio )

      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' Namelist : nampisbio'
         WRITE(numout,*) '    part of calcite not dissolved in guts     part      =', part
         WRITE(numout,*) '    frequence pour la biologie                nrdttrc   =', nrdttrc
         WRITE(numout,*) '    POC sinking speed                         wsbio     =', wsbio
         WRITE(numout,*) '    half saturation constant for mortality    xkmort    =', xkmort
         WRITE(numout,*) '    Fe/C in zooplankton                       ferat3    =', ferat3
         WRITE(numout,*) '    Big particles sinking speed               wsbio2    =', wsbio2
      ENDIF

#if defined key_kriest

      !                               ! nampiskrp : kriest parameters
      !                               ! -----------------------------
      xkr_eta      = 0.62        
      xkr_zeta     = 1.62        
      xkr_mass_min = 0.0002     
      xkr_mass_max = 1.      

      REWIND( numnat )                     ! read natkriest
      READ  ( numnat, nampiskrp )

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) ' Namelist : nampiskrp'
         WRITE(numout,*) '    Sinking  exponent                        xkr_eta      = ', xkr_eta
         WRITE(numout,*) '    N content exponent                       xkr_zeta     = ', xkr_zeta
         WRITE(numout,*) '    Minimum mass for Aggregates              xkr_mass_min = ', xkr_mass_min
         WRITE(numout,*) '    Maximum mass for Aggregates              xkr_mass_max = ', xkr_mass_max
         WRITE(numout,*)
     ENDIF


     ! Computation of some variables
     xkr_massp = 5.7E-6 * 7.6 * xkr_mass_min**xkr_zeta

#endif
      !
#if defined key_diatrc && ! defined key_iomput

      ! Namelist namlobdia
      ! -------------------
      nn_writedia = 10                   ! default values

      DO jl = 1, jp_pisces_2d
         jn = jp_pcs0_2d + jl - 1
         WRITE(ctrc2d(jn),'("2D_",I1)') jn                      ! short name
         WRITE(ctrc2l(jn),'("2D DIAGNOSTIC NUMBER ",I2)') jn    ! long name
         ctrc2u(jn) = ' '                                       ! units
      END DO
      !                                 ! 3D output arrays
      DO jl = 1, jp_pisces_3d
         jn = jp_pcs0_3d + jl - 1
         WRITE(ctrc3d(jn),'("3D_",I1)') jn                      ! short name
         WRITE(ctrc3l(jn),'("3D DIAGNOSTIC NUMBER ",I2)') jn    ! long name
         ctrc3u(jn) = ' '                                       ! units
      END DO

      REWIND( numnat )               ! read natrtd
      READ  ( numnat, nampisdia )

      DO jl = 1, jp_pisces_2d
         jn = jp_pcs0_2d + jl - 1
         ctrc2d(jn) = pisdia2d(jl)%snamedia
         ctrc2l(jn) = pisdia2d(jl)%lnamedia
         ctrc2u(jn) = pisdia2d(jl)%unitdia
      END DO

      DO jl = 1, jp_pisces_3d
         jn = jp_pcs0_3d + jl - 1
         ctrc3d(jn) = pisdia3d(jl)%snamedia
         ctrc3l(jn) = pisdia3d(jl)%lnamedia
         ctrc3u(jn) = pisdia3d(jl)%unitdia
      END DO

      IF(lwp) THEN                   ! control print
         WRITE(numout,*)
         WRITE(numout,*) ' Namelist : natadd'
         WRITE(numout,*) '    frequency of outputs for additional arrays nn_writedia = ', nn_writedia
         DO jl = 1, jp_pisces_3d
            jn = jp_pcs0_3d + jl - 1
            WRITE(numout,*) '   3d output field No : ',jn
            WRITE(numout,*) '   short name         : ', TRIM(ctrc3d(jn))
            WRITE(numout,*) '   long name          : ', TRIM(ctrc3l(jn))
            WRITE(numout,*) '   unit               : ', TRIM(ctrc3u(jn))
            WRITE(numout,*) ' '
         END DO

         DO jl = 1, jp_pisces_2d
            jn = jp_pcs0_2d + jl - 1
            WRITE(numout,*) '   2d output field No : ',jn
            WRITE(numout,*) '   short name         : ', TRIM(ctrc2d(jn))
            WRITE(numout,*) '   long name          : ', TRIM(ctrc2l(jn))
            WRITE(numout,*) '   unit               : ', TRIM(ctrc2u(jn))
            WRITE(numout,*) ' '
         END DO
      ENDIF
#endif

      REWIND( numnat )
      READ  ( numnat, nampisdmp )

      IF(lwp) THEN                         ! control print
         WRITE(numout,*)
         WRITE(numout,*) ' Namelist : nampisdmp'
         WRITE(numout,*) '    Relaxation of tracer to glodap mean value            ln_pisdmp      =', ln_pisdmp
         WRITE(numout,*) '    Restoring of tracer to initial value  on closed seas  ln_pisclo      =', ln_pisclo
         WRITE(numout,*) ' '
      ENDIF

   END SUBROUTINE trc_nam_pisces

#else
   !!----------------------------------------------------------------------
   !!  Dummy module :                                   No PISCES bio-model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_nam_pisces                      ! Empty routine
   END  SUBROUTINE  trc_nam_pisces
#endif  

   !!======================================================================
END MODULE trcnam_pisces
