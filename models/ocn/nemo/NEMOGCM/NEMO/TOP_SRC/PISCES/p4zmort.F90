MODULE p4zmort
   !!======================================================================
   !!                         ***  MODULE p4zmort  ***
   !! TOP :   PISCES Compute the mortality terms for phytoplankton
   !!======================================================================
   !! History :   1.0  !  2002     (O. Aumont)  Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   'key_pisces'                                       PISCES bio-model
   !!----------------------------------------------------------------------
   !!   p4z_mort       :   Compute the mortality terms for phytoplankton
   !!   p4z_mort_init  :   Initialize the mortality params for phytoplankton
   !!----------------------------------------------------------------------
   USE trc
   USE oce_trc         !
   USE trc         ! 
   USE sms_pisces      ! 
   USE p4zsink
   USE prtctl_trc

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_mort    
   PUBLIC   p4z_mort_init    


   !! * Shared module variables
   REAL(wp), PUBLIC ::   &
     wchl   = 0.001_wp    ,  &  !:
     wchld  = 0.02_wp     ,  &  !:
     mprat  = 0.01_wp     ,  &  !:
     mprat2 = 0.01_wp     ,  &  !:
     mpratm = 0.01_wp           !:


   !!* Substitution
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: p4zmort.F90 2528 2010-12-27 17:33:53Z rblod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE p4z_mort( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_mort  ***
      !!
      !! ** Purpose :   Calls the different subroutine to initialize and compute
      !!                the different phytoplankton mortality terms
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt ! ocean time step
      !!---------------------------------------------------------------------

      CALL p4z_nano            ! nanophytoplankton

      CALL p4z_diat            ! diatoms

   END SUBROUTINE p4z_mort


   SUBROUTINE p4z_nano
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_nano  ***
      !!
      !! ** Purpose :   Compute the mortality terms for nanophytoplankton
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER  :: ji, jj, jk
      REAL(wp) :: zcompaph
      REAL(wp) :: zfactfe, zfactch, zprcaca, zfracal
      REAL(wp) :: ztortp , zrespp , zmortp , zstep
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------


#if defined key_diatrc
     prodcal(:,:,:) = 0.  !: Initialisation of calcite production variable
#endif

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi

               zcompaph = MAX( ( trn(ji,jj,jk,jpphy) - 1e-8 ), 0.e0 )

# if defined key_degrad
               zstep =  xstep * facvol(ji,jj,jk)  
# else
               zstep =  xstep  
# endif
               !     Squared mortality of Phyto similar to a sedimentation term during
               !     blooms (Doney et al. 1996)
               zrespp = wchl * 1.e6 * zstep * xdiss(ji,jj,jk) * zcompaph * trn(ji,jj,jk,jpphy) 

               !     Phytoplankton mortality. This mortality loss is slightly
               !     increased when nutrients are limiting phytoplankton growth
               !     as observed for instance in case of iron limitation.
               ztortp = mprat * xstep * trn(ji,jj,jk,jpphy) / ( xkmort + trn(ji,jj,jk,jpphy) ) * zcompaph

               zmortp = zrespp + ztortp

               !   Update the arrays TRA which contains the biological sources and sinks

               zfactfe = trn(ji,jj,jk,jpnfe)/(trn(ji,jj,jk,jpphy)+rtrn)
               zfactch = trn(ji,jj,jk,jpnch)/(trn(ji,jj,jk,jpphy)+rtrn)

               tra(ji,jj,jk,jpphy) = tra(ji,jj,jk,jpphy) - zmortp
               tra(ji,jj,jk,jpnch) = tra(ji,jj,jk,jpnch) - zmortp * zfactch
               tra(ji,jj,jk,jpnfe) = tra(ji,jj,jk,jpnfe) - zmortp * zfactfe
               zprcaca = xfracal(ji,jj,jk) * zmortp
#if defined key_diatrc
               prodcal(ji,jj,jk) = prodcal(ji,jj,jk) + zprcaca  ! prodcal=prodcal(nanophy)+prodcal(microzoo)+prodcal(mesozoo)
#endif
               zfracal = 0.5 * xfracal(ji,jj,jk)
               tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) - zprcaca
               tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) - 2. * zprcaca
               tra(ji,jj,jk,jpcal) = tra(ji,jj,jk,jpcal) + zprcaca
#if defined key_kriest
               tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) + zmortp
               tra(ji,jj,jk,jpnum) = tra(ji,jj,jk,jpnum) + ztortp * xkr_dnano + zrespp * xkr_ddiat
               tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) + zmortp * zfactfe
#else
               tra(ji,jj,jk,jpgoc) = tra(ji,jj,jk,jpgoc) + zfracal * zmortp
               tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) + ( 1. - zfracal ) * zmortp
               tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) + ( 1. - zfracal ) * zmortp * zfactfe
               tra(ji,jj,jk,jpbfe) = tra(ji,jj,jk,jpbfe) + zfracal * zmortp * zfactfe
#endif
            END DO
         END DO
      END DO
      !
       IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('nano')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
       ENDIF

   END SUBROUTINE p4z_nano

   SUBROUTINE p4z_diat
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_diat  ***
      !!
      !! ** Purpose :   Compute the mortality terms for diatoms
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER  ::  ji, jj, jk
      REAL(wp) ::  zfactfe,zfactsi,zfactch, zcompadi
      REAL(wp) ::  zrespp2, ztortp2, zmortp2, zstep
      CHARACTER (len=25) :: charout
 
      !!---------------------------------------------------------------------


      !    Aggregation term for diatoms is increased in case of nutrient
      !    stress as observed in reality. The stressed cells become more
      !    sticky and coagulate to sink quickly out of the euphotic zone
      !     ------------------------------------------------------------

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi

               zcompadi = MAX( ( trn(ji,jj,jk,jpdia) - 1e-8), 0. )

               !    Aggregation term for diatoms is increased in case of nutrient
               !    stress as observed in reality. The stressed cells become more
               !    sticky and coagulate to sink quickly out of the euphotic zone
               !     ------------------------------------------------------------

# if defined key_degrad
               zstep =  xstep * facvol(ji,jj,jk)  
# else
               zstep =  xstep  
# endif
               !  Phytoplankton respiration 
               !     ------------------------
               zrespp2  = 1.e6 * zstep * (  wchl + wchld * ( 1.- xlimdia(ji,jj,jk) )  )    &
                  &       * xdiss(ji,jj,jk) * zcompadi * trn(ji,jj,jk,jpdia)

               !     Phytoplankton mortality. 
               !     ------------------------
               ztortp2  = mprat2 * zstep * trn(ji,jj,jk,jpdia)  / ( xkmort + trn(ji,jj,jk,jpdia) ) * zcompadi 

               zmortp2 = zrespp2 + ztortp2

               !   Update the arrays tra which contains the biological sources and sinks
               !   ---------------------------------------------------------------------
               zfactch = trn(ji,jj,jk,jpdch) / ( trn(ji,jj,jk,jpdia) + rtrn )
               zfactfe = trn(ji,jj,jk,jpdfe) / ( trn(ji,jj,jk,jpdia) + rtrn )
               zfactsi = trn(ji,jj,jk,jpbsi) / ( trn(ji,jj,jk,jpdia) + rtrn )

               tra(ji,jj,jk,jpdia) = tra(ji,jj,jk,jpdia) - zmortp2 
               tra(ji,jj,jk,jpdch) = tra(ji,jj,jk,jpdch) - zmortp2 * zfactch
               tra(ji,jj,jk,jpdfe) = tra(ji,jj,jk,jpdfe) - zmortp2 * zfactfe
               tra(ji,jj,jk,jpbsi) = tra(ji,jj,jk,jpbsi) - zmortp2 * zfactsi
               tra(ji,jj,jk,jpdsi) = tra(ji,jj,jk,jpdsi) + zmortp2 * zfactsi
#if defined key_kriest
               tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) + zmortp2  
               tra(ji,jj,jk,jpnum) = tra(ji,jj,jk,jpnum) + ztortp2 * xkr_ddiat + zrespp2 * xkr_daggr
               tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) + zmortp2 * zfactfe
#else
               tra(ji,jj,jk,jpgoc) = tra(ji,jj,jk,jpgoc) + zrespp2 + 0.5 * ztortp2
               tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) + 0.5 * ztortp2
               tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) + 0.5 * ztortp2 * zfactfe
               tra(ji,jj,jk,jpbfe) = tra(ji,jj,jk,jpbfe) + ( zrespp2 + 0.5 * ztortp2 ) * zfactfe
#endif
            END DO
         END DO
      END DO
      !
        IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('diat')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
       ENDIF
             
   END SUBROUTINE p4z_diat

   SUBROUTINE p4z_mort_init

      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_mort_init  ***
      !!
      !! ** Purpose :   Initialization of phytoplankton parameters
      !!
      !! ** Method  :   Read the nampismort namelist and check the parameters
      !!      called at the first timestep
      !!
      !! ** input   :   Namelist nampismort
      !!
      !!----------------------------------------------------------------------

      NAMELIST/nampismort/ wchl, wchld, mprat, mprat2, mpratm

      REWIND( numnat )                     ! read numnat
      READ  ( numnat, nampismort )

      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for phytoplankton mortality, nampismort'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    quadratic mortality of phytoplankton      wchl      =', wchl
         WRITE(numout,*) '    maximum quadratic mortality of diatoms    wchld     =', wchld
         WRITE(numout,*) '    phytoplankton mortality rate              mprat     =', mprat
         WRITE(numout,*) '    Diatoms mortality rate                    mprat2    =', mprat2
         WRITE(numout,*) '    Phytoplankton minimum mortality rate      mpratm    =', mpratm
      ENDIF

   END SUBROUTINE p4z_mort_init

#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_mort                    ! Empty routine
   END SUBROUTINE p4z_mort
#endif 

   !!======================================================================
END MODULE  p4zmort
