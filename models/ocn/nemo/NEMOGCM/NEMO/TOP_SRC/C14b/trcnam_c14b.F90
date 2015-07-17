MODULE trcnam_c14b
   !!======================================================================
   !!                         ***  MODULE trcnam_c14b  ***
   !! TOP :   initialisation of some run parameters for C14 chemical model
   !!======================================================================
   !! History :   2.0  !  2007-12  (C. Ethe, G. Madec) from trcnam.cfc.h90
   !!----------------------------------------------------------------------
#if defined key_c14b
   !!----------------------------------------------------------------------
   !!   'key_c14b'                                         C14 bomb tracer
   !!----------------------------------------------------------------------
   !! trc_nam_c14b      : C14 model initialisation
   !!----------------------------------------------------------------------
   USE oce_trc         ! Ocean variables
   USE par_trc         ! TOP parameters
   USE trc             ! TOP variables
   USE trcsms_c14b     ! C14b specific variable

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_nam_c14b   ! called by trcnam.F90 module

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcnam_c14b.F90 2715 2011-03-30 15:58:35Z rblod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trc_nam_c14b
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE trc_nam_c14b  ***
      !!                 
      !! ** Purpose :   Definition some run parameter for C14 model
      !!
      !! ** Method  :   Read the namc14 namelist and check the parameter 
      !!       values called at the first timestep (nit000)
      !!
      !! ** input   :   Namelist namelist_c14b
      !!----------------------------------------------------------------------
      INTEGER ::   numnatb

#if defined key_diatrc && ! defined key_iomput
      ! definition of additional diagnostic as a structure
      INTEGER ::   jl, jn
      TYPE DIAG
         CHARACTER(len = 20)  :: snamedia   !: short name
         CHARACTER(len = 80 ) :: lnamedia   !: long name
         CHARACTER(len = 20 ) :: unitdia    !: unit
      END TYPE DIAG

      TYPE(DIAG) , DIMENSION(jp_c14b_2d) :: c14dia2d
      TYPE(DIAG) , DIMENSION(jp_c14b_3d) :: c14dia3d
#endif
      !!
      NAMELIST/namc14date/ ndate_beg_b, nyear_res_b
#if defined key_diatrc && ! defined key_iomput
      NAMELIST/namc14dia/nn_writedia, c14dia2d, c14dia3d     ! additional diagnostics
#endif
      !!-------------------------------------------------------------------

      ndate_beg_b = 650101            ! default namelist value
      nyear_res_b = 1955

      !                             ! Open namelist file
      CALL ctl_opn( numnatb, 'namelist_c14b', 'OLD', 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )
         
      READ( numnatb , namc14date )     ! read namelist

      IF(lwp) THEN                  ! control print
         WRITE(numout,*)
         WRITE(numout,*) ' trc_nam: Read namdates, namelist for C14 chemical model'
         WRITE(numout,*) ' ~~~~~~~'
         WRITE(numout,*) '    initial calendar date (aammjj) for C14  ndate_beg_b = ', ndate_beg_b
         WRITE(numout,*) '    restoring time constant (year)          nyear_res_b = ', nyear_res_b
      ENDIF
      nyear_beg_b = ndate_beg_b / 10000
      IF(lwp) WRITE(numout,*) '    initial year (aa)                  nyear_beg_b = ', nyear_beg_b
      !
#if defined key_diatrc && ! defined key_iomput

      ! Namelist namc14dia
      ! -------------------
      nn_writedia = 10                   ! default values

      DO jl = 1, jp_c14b_2d
         jn = jp_c14b0_2d + jl - 1
         WRITE(ctrc2d(jn),'("2D_",I1)') jn                      ! short name
         WRITE(ctrc2l(jn),'("2D DIAGNOSTIC NUMBER ",I2)') jn    ! long name
         ctrc2u(jn) = ' '                                       ! units
      END DO
      !                                 ! 3D output arrays
      DO jl = 1, jp_c14b_3d
         jn = jp_c14b0_3d + jl - 1
         WRITE(ctrc3d(jn),'("3D_",I1)') jn                      ! short name
         WRITE(ctrc3l(jn),'("3D DIAGNOSTIC NUMBER ",I2)') jn    ! long name
         ctrc3u(jn) = ' '                                       ! units
      END DO

      REWIND( numnatb )               ! read natrtd
      READ  ( numnatb, namc14dia )

      DO jl = 1, jp_c14b_2d
         jn = jp_c14b0_2d + jl - 1
         ctrc2d(jn) = c14dia2d(jl)%snamedia
         ctrc2l(jn) = c14dia2d(jl)%lnamedia
         ctrc2u(jn) = c14dia2d(jl)%unitdia
      END DO

      DO jl = 1, jp_c14b_3d
         jn = jp_c14b0_3d + jl - 1
         ctrc3d(jn) = c14dia3d(jl)%snamedia
         ctrc3l(jn) = c14dia3d(jl)%lnamedia
         ctrc3u(jn) = c14dia3d(jl)%unitdia
      END DO

      IF(lwp) THEN                   ! control print
         WRITE(numout,*)
         WRITE(numout,*) ' Namelist : natadd'
         WRITE(numout,*) '    frequency of outputs for additional arrays nn_writedia = ', nn_writedia
         DO jl = 1, jp_c14b_3d
            jn = jp_c14b0_3d + jl - 1
            WRITE(numout,*) '   3d output field No : ',jn
            WRITE(numout,*) '   short name         : ', TRIM(ctrc3d(jn))
            WRITE(numout,*) '   long name          : ', TRIM(ctrc3l(jn))
            WRITE(numout,*) '   unit               : ', TRIM(ctrc3u(jn))
            WRITE(numout,*) ' '
         END DO

         DO jl = 1, jp_c14b_2d
            jn = jp_c14b0_2d + jl - 1
            WRITE(numout,*) '   2d output field No : ',jn
            WRITE(numout,*) '   short name         : ', TRIM(ctrc2d(jn))
            WRITE(numout,*) '   long name          : ', TRIM(ctrc2l(jn))
            WRITE(numout,*) '   unit               : ', TRIM(ctrc2u(jn))
            WRITE(numout,*) ' '
         END DO
      ENDIF

#endif

   END SUBROUTINE trc_nam_c14b
   
#else
   !!----------------------------------------------------------------------
   !!  Dummy module :                                                No 14C
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_nam_c14b                      ! Empty routine
   END  SUBROUTINE  trc_nam_c14b
#endif  

   !!======================================================================
END MODULE trcnam_c14b
