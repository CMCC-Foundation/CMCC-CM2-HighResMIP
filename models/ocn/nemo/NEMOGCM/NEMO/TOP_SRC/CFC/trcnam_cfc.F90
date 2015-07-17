MODULE trcnam_cfc
   !!======================================================================
   !!                         ***  MODULE trcnam_cfc  ***
   !! TOP :   initialisation of some run parameters for CFC chemical model
   !!======================================================================
   !! History :   2.0  !  2007-12  (C. Ethe, G. Madec) from trcnam.cfc.h90
   !!----------------------------------------------------------------------
#if defined key_cfc
   !!----------------------------------------------------------------------
   !!   'key_cfc'                                               CFC tracers
   !!----------------------------------------------------------------------
   !! trc_nam_cfc      : CFC model initialisation
   !!----------------------------------------------------------------------
   USE oce_trc         ! Ocean variables
   USE par_trc         ! TOP parameters
   USE trc             ! TOP variables
   USE trcsms_cfc      ! CFC specific variable

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_nam_cfc   ! called by trcnam.F90 module

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcnam_cfc.F90 2715 2011-03-30 15:58:35Z rblod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trc_nam_cfc
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE trc_nam_cfc  ***
      !!                 
      !! ** Purpose :   Definition some run parameter for CFC model
      !!
      !! ** Method  :   Read the namcfc namelist and check the parameter 
      !!       values called at the first timestep (nit000)
      !!
      !! ** input   :   Namelist namcfc
      !!----------------------------------------------------------------------
      INTEGER ::   numnatc
#if defined key_diatrc && ! defined key_iomput
      ! definition of additional diagnostic as a structure
      INTEGER :: jl, jn
      TYPE DIAG
         CHARACTER(len = 20)  :: snamedia   !: short name
         CHARACTER(len = 80 ) :: lnamedia   !: long name
         CHARACTER(len = 20 ) :: unitdia    !: unit
      END TYPE DIAG

      TYPE(DIAG) , DIMENSION(jp_cfc_2d) :: cfcdia2d
#endif
      !!
      NAMELIST/namcfcdate/ ndate_beg, nyear_res
#if defined key_diatrc && ! defined key_iomput
      NAMELIST/namcfcdia/nn_writedia, cfcdia2d     ! additional diagnostics
#endif
      !!-------------------------------------------------------------------

      ndate_beg = 300101            ! default namelist value
      nyear_res = 1950

      !                             ! Open namelist file
      CALL ctl_opn( numnatc, 'namelist_cfc', 'OLD', 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )
         
      READ( numnatc , namcfcdate )     ! read namelist

      IF(lwp) THEN                  ! control print
         WRITE(numout,*)
         WRITE(numout,*) ' trc_nam: Read namdates, namelist for CFC chemical model'
         WRITE(numout,*) ' ~~~~~~~'
         WRITE(numout,*) '    initial calendar date (aammjj) for CFC  ndate_beg = ', ndate_beg
         WRITE(numout,*) '    restoring time constant (year)          nyear_res = ', nyear_res
      ENDIF
      nyear_beg = ndate_beg / 10000
      IF(lwp) WRITE(numout,*) '    initial year (aa)                       nyear_beg = ', nyear_beg
      !
#if defined key_diatrc && ! defined key_iomput

      ! Namelist namcfcdia
      ! -------------------
      nn_writedia = 10                   ! default values

      DO jl = 1, jp_cfc_2d
         jn = jp_cfc0_2d + jl - 1 
         WRITE(ctrc2d(jn),'("2D_",I1)') jn                      ! short name
         WRITE(ctrc2l(jn),'("2D DIAGNOSTIC NUMBER ",I2)') jn    ! long name
         ctrc2u(jn) = ' '                                       ! units
      END DO

      REWIND( numnatc )               ! read natrtd
      READ  ( numnatc, namcfcdia )

      DO jl = 1, jp_cfc_2d
         jn = jp_cfc0_2d + jl - 1
         ctrc2d(jn) = cfcdia2d(jl)%snamedia
         ctrc2l(jn) = cfcdia2d(jl)%lnamedia
         ctrc2u(jn) = cfcdia2d(jl)%unitdia
      END DO


      IF(lwp) THEN                   ! control print
         WRITE(numout,*)
         WRITE(numout,*) ' Namelist : natadd'
         WRITE(numout,*) '    frequency of outputs for additional arrays nn_writedia = ', nn_writedia
         DO jl = 1, jp_cfc_2d
            jn = jp_cfc0_2d + jl - 1
            WRITE(numout,*) '   2d output field No : ',jn
            WRITE(numout,*) '   short name         : ', TRIM(ctrc2d(jn))
            WRITE(numout,*) '   long name          : ', TRIM(ctrc2l(jn))
            WRITE(numout,*) '   unit               : ', TRIM(ctrc2u(jn))
            WRITE(numout,*) ' '
         END DO
      ENDIF
#endif

   END SUBROUTINE trc_nam_cfc
   
#else
   !!----------------------------------------------------------------------
   !!  Dummy module :                                                No CFC
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_nam_cfc                      ! Empty routine
   END  SUBROUTINE  trc_nam_cfc
#endif  

   !!======================================================================
END MODULE trcnam_cfc
