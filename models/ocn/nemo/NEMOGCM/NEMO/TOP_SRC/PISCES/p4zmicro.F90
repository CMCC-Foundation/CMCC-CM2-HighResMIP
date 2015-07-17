MODULE p4zmicro
   !!======================================================================
   !!                         ***  MODULE p4zmicro  ***
   !! TOP :   PISCES Compute the sources/sinks for microzooplankton
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   'key_pisces'                                       PISCES bio-model
   !!----------------------------------------------------------------------
   !!   p4z_micro       :   Compute the sources/sinks for microzooplankton
   !!   p4z_micro_init  :   Initialize and read the appropriate namelist
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

   PUBLIC   p4z_micro         ! called in p4zbio.F90
   PUBLIC   p4z_micro_init    ! called in trcsms_pisces.F90

   !! * Shared module variables
   REAL(wp), PUBLIC ::   &
      xpref2c = 0.0_wp       ,  &  !:
      xpref2p = 0.5_wp       ,  &  !:
      xpref2d = 0.5_wp       ,  &  !:
      resrat  = 0.03_wp      ,  &  !:
      mzrat   = 0.0_wp       ,  &  !:
      grazrat = 4.0_wp       ,  &  !:
      xkgraz  = 20E-6_wp     ,  &  !:
      unass   = 0.3_wp       ,  &  !:
      sigma1  = 0.6_wp       ,  &  !:
      epsher  = 0.33_wp


   !!* Substitution
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: p4zmicro.F90 2528 2010-12-27 17:33:53Z rblod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE p4z_micro( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_micro  ***
      !!
      !! ** Purpose :   Compute the sources/sinks for microzooplankton
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt ! ocean time step
      INTEGER  :: ji, jj, jk
      REAL(wp) :: zcompadi, zcompadi2, zcompaz , zcompaph, zcompapoc
      REAL(wp) :: zgraze  , zdenom  , zdenom2, zstep
      REAL(wp) :: zfact   , zinano , zidiat, zipoc
      REAL(wp) :: zgrarem, zgrafer, zgrapoc, zprcaca, zmortz
      REAL(wp) :: zrespz, ztortz
      REAL(wp) :: zgrazp, zgrazm, zgrazsd
      REAL(wp) :: zgrazmf, zgrazsf, zgrazpf
      CHARACTER (len=25) :: charout

      !!---------------------------------------------------------------------


#if defined key_diatrc
      grazing(:,:,:) = 0.  !: Initialisation of  grazing
#endif

      zstep = rfact2 / rday      ! Time step duration for biology

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               zcompaz = MAX( ( trn(ji,jj,jk,jpzoo) - 1.e-9 ), 0.e0 )
# if defined key_degrad
               zstep   = xstep * facvol(ji,jj,jk)
# else
               zstep   = xstep
# endif
               zfact   = zstep * tgfunc(ji,jj,jk) * zcompaz

               !  Respiration rates of both zooplankton
               !  -------------------------------------
               zrespz = resrat * zfact  * ( 1.+ 3.* nitrfac(ji,jj,jk) )     &
                  &            * trn(ji,jj,jk,jpzoo) / ( xkmort + trn(ji,jj,jk,jpzoo) )

               !  Zooplankton mortality. A square function has been selected with
               !  no real reason except that it seems to be more stable and may mimic predation.
               !  ---------------------------------------------------------------
               ztortz = mzrat * 1.e6 * zfact * trn(ji,jj,jk,jpzoo)

               zcompadi  = MAX( ( trn(ji,jj,jk,jpdia) - 1.e-8 ), 0.e0 )
               zcompadi2 = MIN( zcompadi, 5.e-7 )
               zcompaph  = MAX( ( trn(ji,jj,jk,jpphy) - 2.e-7 ), 0.e0 )
               zcompapoc = MAX( ( trn(ji,jj,jk,jppoc) - 1.e-8 ), 0.e0 )
               
               !     Microzooplankton grazing
               !     ------------------------
               zdenom2 = 1./ ( xpref2p * zcompaph + xpref2c * zcompapoc + xpref2d * zcompadi2 + rtrn )

               zgraze = grazrat * zstep * tgfunc(ji,jj,jk) * trn(ji,jj,jk,jpzoo)

               zinano = xpref2p * zcompaph  * zdenom2
               zipoc  = xpref2c * zcompapoc * zdenom2
               zidiat = xpref2d * zcompadi2 * zdenom2

               zdenom = 1./ ( xkgraz + zinano * zcompaph + zipoc * zcompapoc + zidiat * zcompadi2 )

               zgrazp  = zgraze * zinano * zcompaph * zdenom
               zgrazm  = zgraze * zipoc  * zcompapoc * zdenom
               zgrazsd = zgraze * zidiat * zcompadi2 * zdenom

               zgrazpf = zgrazp  * trn(ji,jj,jk,jpnfe) / (trn(ji,jj,jk,jpphy) + rtrn)
               zgrazmf = zgrazm  * trn(ji,jj,jk,jpsfe) / (trn(ji,jj,jk,jppoc) + rtrn)
               zgrazsf = zgrazsd * trn(ji,jj,jk,jpdfe) / (trn(ji,jj,jk,jpdia) + rtrn)
#if defined key_diatrc
               ! Grazing by microzooplankton
               grazing(ji,jj,jk) = grazing(ji,jj,jk) + zgrazp + zgrazm + zgrazsd 
#endif

               !    Various remineralization and excretion terms
               !    --------------------------------------------
               zgrarem = ( zgrazp + zgrazm + zgrazsd ) * ( 1.- epsher - unass )
               zgrafer = ( zgrazpf + zgrazsf + zgrazmf ) * ( 1.- epsher - unass ) &
                  &      + epsher * ( zgrazm  * MAX((trn(ji,jj,jk,jpsfe) / (trn(ji,jj,jk,jppoc)+ rtrn)-ferat3),0.e0) & 
                  &                 + zgrazp  * MAX((trn(ji,jj,jk,jpnfe) / (trn(ji,jj,jk,jpphy)+ rtrn)-ferat3),0.e0) &
                  &                 + zgrazsd * MAX((trn(ji,jj,jk,jpdfe) / (trn(ji,jj,jk,jpdia)+ rtrn)-ferat3),0.e0 )  )

               zgrapoc = (  zgrazp + zgrazm + zgrazsd ) 

               !  Update of the TRA arrays
               !  ------------------------

               tra(ji,jj,jk,jppo4) = tra(ji,jj,jk,jppo4) + zgrarem * sigma1
               tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) + zgrarem * sigma1
               tra(ji,jj,jk,jpdoc) = tra(ji,jj,jk,jpdoc) + zgrarem * (1.-sigma1)
               tra(ji,jj,jk,jpoxy) = tra(ji,jj,jk,jpoxy) - o2ut * zgrarem * sigma1
               tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) + zgrafer
               tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) + zgrapoc * unass
               tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) + zgrarem * sigma1
#if defined key_kriest
               tra(ji,jj,jk,jpnum) = tra(ji,jj,jk,jpnum) + zgrapoc * unass * xkr_ddiat
#endif

               !
               !   Update the arrays TRA which contain the biological sources and sinks
               !   --------------------------------------------------------------------

               zmortz = ztortz + zrespz
               tra(ji,jj,jk,jpzoo) = tra(ji,jj,jk,jpzoo) - zmortz + epsher * zgrapoc 
               tra(ji,jj,jk,jpphy) = tra(ji,jj,jk,jpphy) - zgrazp
               tra(ji,jj,jk,jpdia) = tra(ji,jj,jk,jpdia) - zgrazsd
               tra(ji,jj,jk,jpnch) = tra(ji,jj,jk,jpnch) - zgrazp  * trn(ji,jj,jk,jpnch)/(trn(ji,jj,jk,jpphy)+rtrn)
               tra(ji,jj,jk,jpdch) = tra(ji,jj,jk,jpdch) - zgrazsd * trn(ji,jj,jk,jpdch)/(trn(ji,jj,jk,jpdia)+rtrn)
               tra(ji,jj,jk,jpbsi) = tra(ji,jj,jk,jpbsi) - zgrazsd * trn(ji,jj,jk,jpbsi)/(trn(ji,jj,jk,jpdia)+rtrn)
               tra(ji,jj,jk,jpdsi) = tra(ji,jj,jk,jpdsi) + zgrazsd * trn(ji,jj,jk,jpbsi)/(trn(ji,jj,jk,jpdia)+rtrn)
               tra(ji,jj,jk,jpnfe) = tra(ji,jj,jk,jpnfe) - zgrazpf
               tra(ji,jj,jk,jpdfe) = tra(ji,jj,jk,jpdfe) - zgrazsf
               tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) + zmortz - zgrazm
               tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) + ferat3 * zmortz + unass * ( zgrazpf + zgrazsf ) - (1.-unass) * zgrazmf
               zprcaca = xfracal(ji,jj,jk) * unass * zgrazp
#if defined key_diatrc
               prodcal(ji,jj,jk) = prodcal(ji,jj,jk) + zprcaca  ! prodcal=prodcal(nanophy)+prodcal(microzoo)+prodcal(mesozoo)
#endif
               zprcaca = part * zprcaca
               tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) - zprcaca
               tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) - 2. * zprcaca
               tra(ji,jj,jk,jpcal) = tra(ji,jj,jk,jpcal) + zprcaca
#if defined key_kriest
               tra(ji,jj,jk,jpnum) = tra(ji,jj,jk,jpnum) + ( zmortz - zgrazm ) * xkr_ddiat
#endif
            END DO
         END DO
      END DO
      !
      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('micro')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF

   END SUBROUTINE p4z_micro


   SUBROUTINE p4z_micro_init

      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_micro_init  ***
      !!
      !! ** Purpose :   Initialization of microzooplankton parameters
      !!
      !! ** Method  :   Read the nampiszoo namelist and check the parameters
      !!      called at the first timestep (nit000)
      !!
      !! ** input   :   Namelist nampiszoo
      !!
      !!----------------------------------------------------------------------

      NAMELIST/nampiszoo/ grazrat,resrat,mzrat,xpref2c, xpref2p, &
         &             xpref2d, xkgraz, epsher, sigma1, unass

      REWIND( numnat )                     ! read numnat
      READ  ( numnat, nampiszoo )

      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for microzooplankton, nampiszoo'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    zoo preference for POC                    xpref2c    =', xpref2c
         WRITE(numout,*) '    zoo preference for nano                   xpref2p    =', xpref2p
         WRITE(numout,*) '    zoo preference for diatoms                xpref2d    =', xpref2d
         WRITE(numout,*) '    exsudation rate of microzooplankton       resrat    =', resrat
         WRITE(numout,*) '    microzooplankton mortality rate           mzrat     =', mzrat
         WRITE(numout,*) '    maximal microzoo grazing rate             grazrat   =', grazrat
         WRITE(numout,*) '    non assimilated fraction of P by microzoo unass     =', unass
         WRITE(numout,*) '    Efficicency of microzoo growth            epsher    =', epsher
         WRITE(numout,*) '    Fraction of microzoo excretion as DOM     sigma1    =', sigma1
         WRITE(numout,*) '    half sturation constant for grazing 1     xkgraz    =', xkgraz
      ENDIF

   END SUBROUTINE p4z_micro_init

#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_micro                    ! Empty routine
   END SUBROUTINE p4z_micro
#endif 

   !!======================================================================
END MODULE  p4zmicro
