MODULE p4zmeso
   !!======================================================================
   !!                         ***  MODULE p4zmeso  ***
   !! TOP :   PISCES Compute the sources/sinks for mesozooplankton
   !!======================================================================
   !! History :   1.0  !  2002     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   'key_pisces'                                       PISCES bio-model
   !!----------------------------------------------------------------------
   !!   p4z_meso       :   Compute the sources/sinks for mesozooplankton
   !!   p4z_meso_init  :   Initialization of the parameters for mesozooplankton
   !!----------------------------------------------------------------------
   USE trc
   USE oce_trc         !
   USE trc         ! 
   USE sms_pisces      ! 
   USE prtctl_trc
   USE p4zint
   USE p4zsink
   USE iom

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_meso              ! called in p4zbio.F90
   PUBLIC   p4z_meso_init         ! called in trcsms_pisces.F90

   !! * Shared module variables
   REAL(wp), PUBLIC ::   &
      xprefc   = 1.0_wp     ,  &  !: 
      xprefp   = 0.2_wp     ,  &  !:
      xprefz   = 1.0_wp     ,  &  !:
      xprefpoc = 0.0_wp     ,  &  !:
      resrat2  = 0.005_wp   ,  &  !:
      mzrat2   = 0.03_wp    ,  &  !:
      grazrat2 = 0.7_wp     ,  &  !:
      xkgraz2  = 20E-6_wp   ,  &  !:
      unass2   = 0.3_wp     ,  &  !:
      sigma2   = 0.6_wp     ,  &  !:
      epsher2  = 0.33_wp    ,  &  !:   
      grazflux = 5.E3_wp 


   !!* Substitution
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: p4zmeso.F90 2528 2010-12-27 17:33:53Z rblod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE p4z_meso( kt, jnt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_meso  ***
      !!
      !! ** Purpose :   Compute the sources/sinks for mesozooplankton
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt, jnt ! ocean time step
      INTEGER  :: ji, jj, jk
      REAL(wp) :: zcompadi, zcompaph, zcompapoc, zcompaz
      REAL(wp) :: zfact, zcompam, zdenom, zgraze2, zstep
      REAL(wp) :: zgrarem2, zgrafer2, zgrapoc2, zprcaca, zmortz2
#if defined key_kriest
      REAL znumpoc
#endif
      REAL(wp) :: zrespz2,ztortz2,zgrazd,zgrazz,zgrazpof
      REAL(wp) :: zgrazn,zgrazpoc,zgraznf,zgrazf
      REAL(wp) :: zgrazfff,zgrazffe
      CHARACTER (len=25) :: charout
#if defined key_diatrc && defined key_iomput
      REAL(wp) :: zrfact2
#endif

      !!---------------------------------------------------------------------

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi

               zcompam = MAX( ( trn(ji,jj,jk,jpmes) - 1.e-9 ), 0.e0 )
# if defined key_degrad
               zstep   = xstep * facvol(ji,jj,jk)
# else
               zstep   = xstep
# endif
               zfact   = zstep * tgfunc(ji,jj,jk) * zcompam

               !  Respiration rates of both zooplankton
               !  -------------------------------------
               zrespz2  = resrat2 * zfact * ( 1. + 3. * nitrfac(ji,jj,jk) )        &
                  &     * trn(ji,jj,jk,jpmes) / ( xkmort + trn(ji,jj,jk,jpmes) )

               !  Zooplankton mortality. A square function has been selected with
               !  no real reason except that it seems to be more stable and may mimic predation
               !  ---------------------------------------------------------------
               ztortz2 = mzrat2 * 1.e6 * zfact * trn(ji,jj,jk,jpmes)
               !

               zcompadi  = MAX( ( trn(ji,jj,jk,jpdia) - 1.e-8 ), 0.e0 )
               zcompaz   = MAX( ( trn(ji,jj,jk,jpzoo) - 1.e-8 ), 0.e0 )
               zcompaph  = MAX( ( trn(ji,jj,jk,jpphy) - 2.e-7 ), 0.e0 )
               zcompapoc = MAX( ( trn(ji,jj,jk,jppoc) - 1.e-8 ), 0.e0 )

               !  Microzooplankton grazing
               !     ------------------------
               zdenom = 1. / (  xkgraz2 + xprefc   * trn(ji,jj,jk,jpdia)   &
                  &                     + xprefz   * trn(ji,jj,jk,jpzoo)   &
                  &                     + xprefp   * trn(ji,jj,jk,jpphy)   &
                  &                     + xprefpoc * trn(ji,jj,jk,jppoc)  )

               zgraze2 = grazrat2 * zstep * Tgfunc2(ji,jj,jk) * zdenom * trn(ji,jj,jk,jpmes) 

               zgrazd   = zgraze2  * xprefc   * zcompadi
               zgrazz   = zgraze2  * xprefz   * zcompaz
               zgrazn   = zgraze2  * xprefp   * zcompaph
               zgrazpoc = zgraze2  * xprefpoc * zcompapoc

               zgraznf  = zgrazn   * trn(ji,jj,jk,jpnfe) / (trn(ji,jj,jk,jpphy) + rtrn)
               zgrazf   = zgrazd   * trn(ji,jj,jk,jpdfe) / (trn(ji,jj,jk,jpdia) + rtrn)
               zgrazpof = zgrazpoc * trn(ji,jj,jk,jpsfe) / (trn(ji,jj,jk,jppoc) + rtrn)
               
               !  Mesozooplankton flux feeding on GOC
               !  ----------------------------------
# if ! defined key_kriest
               zgrazffe = grazflux * zstep * wsbio4(ji,jj,jk)          &
                  &                 * tgfunc2(ji,jj,jk) * trn(ji,jj,jk,jpgoc) * trn(ji,jj,jk,jpmes)
               zgrazfff = zgrazffe * trn(ji,jj,jk,jpbfe) / (trn(ji,jj,jk,jpgoc) + rtrn)
# else
               !!--------------------------- KRIEST3 -------------------------------------------
               !!               zgrazffe = 0.5 * 1.3e-2 / 5.5e-7 * 0.3 * zstep * wsbio3(ji,jj,jk)     &
               !!                  &     * tgfunc(ji,jj,jk) * trn(ji,jj,jk,jppoc) * trn(ji,jj,jk,jpmes)    &
               !! #  if defined key_degrad
               !!                  &     * facvol(ji,jj,jk)          &
               !! #  endif
               !!                  &     /  (trn(ji,jj,jk,jppoc) * 1.e7 + 0.1)
               !!--------------------------- KRIEST3 -------------------------------------------

              zgrazffe = grazflux * zstep * wsbio3(ji,jj,jk)     &
                  &                * tgfunc2(ji,jj,jk) * trn(ji,jj,jk,jppoc) * trn(ji,jj,jk,jpmes)
              zgrazfff = zgrazffe * trn(ji,jj,jk,jpsfe) / (trn(ji,jj,jk,jppoc) + rtrn)
# endif
      
#if defined key_diatrc
              ! Total grazing ( grazing by microzoo is already computed in p4zmicro ) 
              grazing(ji,jj,jk) = grazing(ji,jj,jk) + (  zgrazd + zgrazz + zgrazn + zgrazpoc + zgrazffe )
#endif

              !    Mesozooplankton efficiency
              !    --------------------------
              zgrarem2 = ( zgrazd + zgrazz + zgrazn + zgrazpoc + zgrazffe ) * ( 1. - epsher2 - unass2 )
#if ! defined key_kriest
              zgrafer2 = ( zgrazf + zgraznf + zgrazz * ferat3 + zgrazpof + zgrazfff ) * ( 1.- epsher2 - unass2 ) & 
                  &     + epsher2 * ( zgrazd   * MAX((trn(ji,jj,jk,jpdfe) / (trn(ji,jj,jk,jpdia) + rtrn)-ferat3),0.) &
                  &                 + zgrazn   * MAX((trn(ji,jj,jk,jpnfe) / (trn(ji,jj,jk,jpphy) + rtrn)-ferat3),0.) &
                  &                 + zgrazpoc * MAX((trn(ji,jj,jk,jpsfe) / (trn(ji,jj,jk,jppoc) + rtrn)-ferat3),0.) &
                  &                 + zgrazffe * MAX((trn(ji,jj,jk,jpbfe) / (trn(ji,jj,jk,jpgoc) + rtrn)-ferat3),0.)  )
#else
              zgrafer2 = ( zgrazf + zgraznf + zgrazz * ferat3 + zgrazpof + zgrazfff ) * ( 1. - epsher2 - unass2 ) &
                  &    + epsher2 * ( zgrazd   * MAX((trn(ji,jj,jk,jpdfe) / (trn(ji,jj,jk,jpdia) + rtrn)-ferat3),0.) &
                  &                + zgrazn   * MAX((trn(ji,jj,jk,jpnfe) / (trn(ji,jj,jk,jpphy) + rtrn)-ferat3),0.) &
                  &                + zgrazpoc * MAX((trn(ji,jj,jk,jpsfe) / (trn(ji,jj,jk,jppoc) + rtrn)-ferat3),0.) &
                  &                + zgrazffe * MAX((trn(ji,jj,jk,jpsfe) / (trn(ji,jj,jk,jppoc) + rtrn)-ferat3),0.)  )

#endif
               !   Update the arrays TRA which contain the biological sources and sinks

               zgrapoc2 =  zgrazd + zgrazz  + zgrazn + zgrazpoc + zgrazffe

               tra(ji,jj,jk,jppo4) = tra(ji,jj,jk,jppo4) + zgrarem2 * sigma2
               tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) + zgrarem2 * sigma2
               tra(ji,jj,jk,jpdoc) = tra(ji,jj,jk,jpdoc) + zgrarem2 * ( 1. - sigma2 )
               tra(ji,jj,jk,jpoxy) = tra(ji,jj,jk,jpoxy) - o2ut * zgrarem2 * sigma2
               tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) + zgrafer2
               tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) + zgrarem2 * sigma2
               
#if defined key_kriest
               tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) + zgrapoc2 * unass2
               tra(ji,jj,jk,jpnum) = tra(ji,jj,jk,jpnum) + zgrapoc2 * unass2 * xkr_dmeso
#else
               tra(ji,jj,jk,jpgoc) = tra(ji,jj,jk,jpgoc) + zgrapoc2 * unass2
#endif
               zmortz2 = ztortz2 + zrespz2
               tra(ji,jj,jk,jpmes) = tra(ji,jj,jk,jpmes) - zmortz2 + epsher2 * zgrapoc2
               tra(ji,jj,jk,jpdia) = tra(ji,jj,jk,jpdia) - zgrazd
               tra(ji,jj,jk,jpzoo) = tra(ji,jj,jk,jpzoo) - zgrazz
               tra(ji,jj,jk,jpphy) = tra(ji,jj,jk,jpphy) - zgrazn
               tra(ji,jj,jk,jpnch) = tra(ji,jj,jk,jpnch) - zgrazn * trn(ji,jj,jk,jpnch) / ( trn(ji,jj,jk,jpphy) + rtrn )
               tra(ji,jj,jk,jpdch) = tra(ji,jj,jk,jpdch) - zgrazd * trn(ji,jj,jk,jpdch) / ( trn(ji,jj,jk,jpdia) + rtrn )
               tra(ji,jj,jk,jpbsi) = tra(ji,jj,jk,jpbsi) - zgrazd * trn(ji,jj,jk,jpbsi) / ( trn(ji,jj,jk,jpdia) + rtrn )
               tra(ji,jj,jk,jpdsi) = tra(ji,jj,jk,jpdsi) + zgrazd * trn(ji,jj,jk,jpbsi) / ( trn(ji,jj,jk,jpdia) + rtrn )
               tra(ji,jj,jk,jpnfe) = tra(ji,jj,jk,jpnfe) - zgraznf
               tra(ji,jj,jk,jpdfe) = tra(ji,jj,jk,jpdfe) - zgrazf

               zprcaca = xfracal(ji,jj,jk) * unass2 * zgrazn
#if defined key_diatrc
               prodcal(ji,jj,jk) = prodcal(ji,jj,jk) + zprcaca  ! prodcal=prodcal(nanophy)+prodcal(microzoo)+prodcal(mesozoo)
#endif
               zprcaca = part * zprcaca
               tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) - zprcaca
               tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) - 2. * zprcaca
               tra(ji,jj,jk,jpcal) = tra(ji,jj,jk,jpcal) + zprcaca
#if defined key_kriest
               znumpoc = trn(ji,jj,jk,jpnum) / ( trn(ji,jj,jk,jppoc) + rtrn )
               tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) + zmortz2 - zgrazpoc - zgrazffe
               tra(ji,jj,jk,jpnum) = tra(ji,jj,jk,jpnum) - zgrazpoc * znumpoc &
                  &    + zmortz2  * xkr_dmeso - zgrazffe * znumpoc * wsbio4(ji,jj,jk) / ( wsbio3(ji,jj,jk) + rtrn )
               tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) + ferat3 * zmortz2 &
               &       + unass2 * ( ferat3 * zgrazz + zgraznf + zgrazf + zgrazpof + zgrazfff ) - zgrazfff - zgrazpof
#else
               tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) - zgrazpoc
               tra(ji,jj,jk,jpgoc) = tra(ji,jj,jk,jpgoc) + zmortz2 - zgrazffe
               tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) - zgrazpof
               tra(ji,jj,jk,jpbfe) = tra(ji,jj,jk,jpbfe) + ferat3 * zmortz2 &
               &       + unass2 * ( ferat3 * zgrazz + zgraznf + zgrazf + zgrazpof + zgrazfff ) - zgrazfff
#endif

            END DO
         END DO
      END DO
      !
#if defined key_diatrc && defined key_iomput
      zrfact2 = 1.e3 * rfact2r
      ! Total grazing of phyto by zoo
      grazing(:,:,:) = grazing(:,:,:) * zrfact2 * tmask(:,:,:)
      ! Calcite production
      prodcal(:,:,:) = prodcal(:,:,:) * zrfact2 * tmask(:,:,:)
      IF( jnt == nrdttrc ) then 
         CALL iom_put( "GRAZ" , grazing  )  ! Total grazing of phyto by zooplankton
         CALL iom_put( "PCAL" , prodcal  )  ! Calcite production
      ENDIF
#endif

       IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('meso')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
       ENDIF

   END SUBROUTINE p4z_meso

   SUBROUTINE p4z_meso_init

      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_meso_init  ***
      !!
      !! ** Purpose :   Initialization of mesozooplankton parameters
      !!
      !! ** Method  :   Read the nampismes namelist and check the parameters
      !!      called at the first timestep (nit000)
      !!
      !! ** input   :   Namelist nampismes
      !!
      !!----------------------------------------------------------------------

      NAMELIST/nampismes/ grazrat2,resrat2,mzrat2,xprefc, xprefp, &
         &             xprefz, xprefpoc, xkgraz2, epsher2, sigma2, unass2, grazflux

      REWIND( numnat )                     ! read numnat
      READ  ( numnat, nampismes )


      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' ' 
         WRITE(numout,*) ' Namelist parameters for mesozooplankton, nampismes'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    zoo preference for phyto                  xprefc    =', xprefc
         WRITE(numout,*) '    zoo preference for POC                    xprefp    =', xprefp
         WRITE(numout,*) '    zoo preference for zoo                    xprefz    =', xprefz
         WRITE(numout,*) '    zoo preference for poc                    xprefpoc  =', xprefpoc
         WRITE(numout,*) '    exsudation rate of mesozooplankton        resrat2   =', resrat2
         WRITE(numout,*) '    mesozooplankton mortality rate            mzrat2    =', mzrat2
         WRITE(numout,*) '    maximal mesozoo grazing rate              grazrat2  =', grazrat2
         WRITE(numout,*) '    mesozoo flux feeding rate                 grazflux  =', grazflux
         WRITE(numout,*) '    non assimilated fraction of P by mesozoo  unass2    =', unass2
         WRITE(numout,*) '    Efficicency of Mesozoo growth             epsher2   =', epsher2
         WRITE(numout,*) '    Fraction of mesozoo excretion as DOM      sigma2    =', sigma2
         WRITE(numout,*) '    half sturation constant for grazing 2     xkgraz2   =', xkgraz2
      ENDIF


   END SUBROUTINE p4z_meso_init


#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_meso                    ! Empty routine
   END SUBROUTINE p4z_meso
#endif 

   !!======================================================================
END MODULE  p4zmeso
