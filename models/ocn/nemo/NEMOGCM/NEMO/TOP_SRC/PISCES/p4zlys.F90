MODULE p4zlys
   !!======================================================================
   !!                         ***  MODULE p4zlys  ***
   !! TOP :   PISCES 
   !!======================================================================
   !! History :    -   !  1988-07  (E. MAIER-REIMER) Original code
   !!              -   !  1998     (O. Aumont) additions
   !!              -   !  1999     (C. Le Quere) modifications
   !!             1.0  !  2004     (O. Aumont) modifications
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   'key_pisces'                                       PISCES bio-model
   !!----------------------------------------------------------------------
   !!   p4z_lys        :   Compute the CaCO3 dissolution 
   !!   p4z_lys_init   :   Read the namelist parameters
   !!----------------------------------------------------------------------
   USE trc
   USE oce_trc         !
   USE trc
   USE sms_pisces
   USE prtctl_trc
   USE iom

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_lys         ! called in trcsms_pisces.F90
   PUBLIC   p4z_lys_init    ! called in trcsms_pisces.F90

   !! * Shared module variables
   REAL(wp), PUBLIC :: kdca = 0.327e3_wp  !: diss. rate constant calcite
   REAL(wp), PUBLIC :: nca  = 1.0_wp      !: order of reaction for calcite dissolution

   !! * Module variables
   REAL(wp) :: calcon = 1.03E-2           !: mean calcite concentration [Ca2+] in sea water [mole/kg solution]
 
   INTEGER  :: rmtss                      !: number of seconds per month 

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: p4zlys.F90 2715 2011-03-30 15:58:35Z rblod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE p4z_lys( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_lys  ***
      !!
      !! ** Purpose :   CALCULATES DEGREE OF CACO3 SATURATION IN THE WATER
      !!                COLUMN, DISSOLUTION/PRECIPITATION OF CACO3 AND LOSS
      !!                OF CACO3 TO THE CACO3 SEDIMENT POOL.
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      USE wrk_nemo, ONLY: wrk_in_use, wrk_not_released
      USE wrk_nemo, ONLY: zco3 => wrk_3d_2, zcaldiss => wrk_3d_3 
      !
      INTEGER, INTENT(in) ::   kt ! ocean time step
      INTEGER  ::   ji, jj, jk, jn
      REAL(wp) ::   zbot, zalk, zdic, zph, zremco3, zah2
      REAL(wp) ::   zdispot, zfact, zalka
      REAL(wp) ::   zomegaca, zexcess, zexcess0
#if defined key_diatrc && defined key_iomput
      REAL(wp) ::   zrfact2
#endif
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------

      IF(  wrk_in_use(3, 2,3) ) THEN
         CALL ctl_stop('p4z_lys: requested workspace arrays unavailable')  ;  RETURN
      END IF

      zco3(:,:,:) = 0.
# if defined key_diatrc && defined key_iomput
      zcaldiss(:,:,:) = 0.
# endif
      !     -------------------------------------------
      !     COMPUTE [CO3--] and [H+] CONCENTRATIONS
      !     -------------------------------------------
      
      DO jn = 1, 5                               !  BEGIN OF ITERATION
         !
!CDIR NOVERRCHK
         DO jk = 1, jpkm1
!CDIR NOVERRCHK
            DO jj = 1, jpj
!CDIR NOVERRCHK
               DO ji = 1, jpi

                  ! SET DUMMY VARIABLE FOR TOTAL BORATE
                  zbot  = borat(ji,jj,jk)

                  ! SET DUMMY VARIABLE FOR TOTAL BORATE
                  zbot  = borat(ji,jj,jk)
                  zfact = rhop (ji,jj,jk) / 1000. + rtrn

                  ! SET DUMMY VARIABLE FOR [H+]
                  zph   = hi(ji,jj,jk) * tmask(ji,jj,jk) / zfact + ( 1.-tmask(ji,jj,jk) ) * 1.e-9

                  ! SET DUMMY VARIABLE FOR [SUM(CO2)]GIVEN 
                  zdic  = trn(ji,jj,jk,jpdic) / zfact
                  zalka = trn(ji,jj,jk,jptal) / zfact

                  ! CALCULATE [ALK]([CO3--], [HCO3-])
                  zalk  = zalka - (  akw3(ji,jj,jk) / zph - zph   &
                     &             + zbot / (1.+ zph / akb3(ji,jj,jk) )  )

                  ! CALCULATE [H+] and [CO3--]
                  zah2 = SQRT( (zdic-zalk)*(zdic-zalk)+   &
                     &     4.*(zalk*ak23(ji,jj,jk)/ak13(ji,jj,jk))   &
                     &     *(2*zdic-zalk))

                  zah2=0.5*ak13(ji,jj,jk)/zalk*((zdic-zalk)+zah2)
                  zco3(ji,jj,jk) = zalk/(2.+zah2/ak23(ji,jj,jk))*zfact

                  hi(ji,jj,jk)  = zah2*zfact

               END DO
            END DO
         END DO
         !
      END DO 

      !     ---------------------------------------------------------
      !        CALCULATE DEGREE OF CACO3 SATURATION AND CORRESPONDING
      !        DISSOLOUTION AND PRECIPITATION OF CACO3 (BE AWARE OF
      !        MGCO3)
      !     ---------------------------------------------------------

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi

               ! DEVIATION OF [CO3--] FROM SATURATION VALUE
               zomegaca = ( calcon * zco3(ji,jj,jk) ) / aksp(ji,jj,jk)

               ! SET DEGREE OF UNDER-/SUPERSATURATION
               zexcess0 = MAX( 0., ( 1.- zomegaca ) )
               zexcess  = zexcess0**nca

               ! AMOUNT CACO3 (12C) THAT RE-ENTERS SOLUTION
               !       (ACCORDING TO THIS FORMULATION ALSO SOME PARTICULATE
               !       CACO3 GETS DISSOLVED EVEN IN THE CASE OF OVERSATURATION)
# if defined key_degrad
              zdispot = kdca * zexcess * trn(ji,jj,jk,jpcal) * facvol(ji,jj,jk)
# else
              zdispot = kdca * zexcess * trn(ji,jj,jk,jpcal)
# endif

              !  CHANGE OF [CO3--] , [ALK], PARTICULATE [CACO3],
              !       AND [SUM(CO2)] DUE TO CACO3 DISSOLUTION/PRECIPITATION
              zremco3 = zdispot / rmtss
              zco3(ji,jj,jk) = zco3(ji,jj,jk) + zremco3 * rfact
              tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) + 2. * zremco3
              tra(ji,jj,jk,jpcal) = tra(ji,jj,jk,jpcal) -      zremco3
              tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) +      zremco3

# if defined key_diatrc && defined key_iomput
              zcaldiss(ji,jj,jk) = zremco3  ! calcite dissolution
# endif
            END DO
         END DO
      END DO

# if defined key_diatrc
#  if ! defined key_iomput
      trc3d(:,:,:,jp_pcs0_3d    ) = hi  (:,:,:)          * tmask(:,:,:)
      trc3d(:,:,:,jp_pcs0_3d + 1) = zco3(:,:,:)          * tmask(:,:,:)
      trc3d(:,:,:,jp_pcs0_3d + 2) = aksp(:,:,:) / calcon * tmask(:,:,:)
#  else
      zrfact2 = 1.e3 * rfact2r
      CALL iom_put( "PH"    , hi      (:,:,:)           * tmask(:,:,:) )
      CALL iom_put( "CO3"   , zco3    (:,:,:)           * tmask(:,:,:) )
      CALL iom_put( "CO3sat", aksp    (:,:,:) / calcon  * tmask(:,:,:) )
      CALL iom_put( "DCAL"  , zcaldiss(:,:,:) * zrfact2 * tmask(:,:,:) )
#  endif
# endif
      !
       IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('lys ')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
       ENDIF

      IF( wrk_not_released(3, 2,3) ) CALL ctl_stop('p4z_lys: failed to release workspace arrays')
      !
   END SUBROUTINE p4z_lys

   SUBROUTINE p4z_lys_init

      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_lys_init  ***
      !!
      !! ** Purpose :   Initialization of CaCO3 dissolution parameters
      !!
      !! ** Method  :   Read the nampiscal namelist and check the parameters
      !!      called at the first timestep (nit000)
      !!
      !! ** input   :   Namelist nampiscal
      !!
      !!----------------------------------------------------------------------

      NAMELIST/nampiscal/ kdca, nca

      REWIND( numnat )                     ! read numnat
      READ  ( numnat, nampiscal )

      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for CaCO3 dissolution, nampiscal'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    diss. rate constant calcite (per month)   kdca      =', kdca
         WRITE(numout,*) '    order of reaction for calcite dissolution nca       =', nca
      ENDIF

      ! Number of seconds per month 
      rmtss =  nyear_len(1) * rday / raamo

   END SUBROUTINE p4z_lys_init

#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_lys( kt )                   ! Empty routine
      INTEGER, INTENT( in ) ::   kt
      WRITE(*,*) 'p4z_lys: You should not have seen this print! error?', kt
   END SUBROUTINE p4z_lys
#endif 
   !!======================================================================
END MODULE  p4zlys
