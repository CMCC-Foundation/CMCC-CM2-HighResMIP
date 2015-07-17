MODULE p4zflx
   !!======================================================================
   !!                         ***  MODULE p4zflx  ***
   !! TOP :   PISCES CALCULATES GAS EXCHANGE AND CHEMISTRY AT SEA SURFACE
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
   !!   p4z_flx       :   CALCULATES GAS EXCHANGE AND CHEMISTRY AT SEA SURFACE
   !!   p4z_flx_init  :   Read the namelist
   !!----------------------------------------------------------------------
   USE trc
   USE oce_trc         !
   USE trc
   USE sms_pisces
   USE prtctl_trc
   USE p4zche
   USE iom
#if defined key_cpl_carbon_cycle
   USE sbc_oce , ONLY :  atm_co2
#endif

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_flx  
   PUBLIC   p4z_flx_init  
   PUBLIC   p4z_flx_alloc  

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: oce_co2   !: ocean carbon flux 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: satmco2   !: atmospheric pco2 

   REAL(wp) ::  t_oce_co2_flx               !: Total ocean carbon flux 
   REAL(wp) ::  t_atm_co2_flx               !: global mean of atmospheric pco2
   REAL(wp) ::  area                        !: ocean surface
   REAL(wp) ::  atcco2 = 278._wp            !: pre-industrial atmospheric [co2] (ppm) 	
   REAL(wp) ::  atcox  = 0.20946_wp         !:
   REAL(wp) ::  xconv  = 0.01_wp / 3600._wp !: coefficients for conversion 

   !!* Substitution
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: p4zflx.F90 2715 2011-03-30 15:58:35Z rblod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_flx ( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_flx  ***
      !!
      !! ** Purpose :   CALCULATES GAS EXCHANGE AND CHEMISTRY AT SEA SURFACE
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      USE wrk_nemo, ONLY:   wrk_in_use, wrk_not_released
      USE wrk_nemo, ONLY:   zkgco2 => wrk_2d_1 , zkgo2 => wrk_2d_2 , zh2co3 => wrk_2d_3 
      USE wrk_nemo, ONLY:   zoflx  => wrk_2d_4 , zkg   => wrk_2d_5
      USE wrk_nemo, ONLY:   zdpco2 => wrk_2d_6 , zdpo2 => wrk_2d_7
      !
      INTEGER, INTENT(in) ::   kt   !
      !
      INTEGER  ::   ji, jj, jrorr
      REAL(wp) ::   ztc, ztc2, ztc3, zws, zkgwan
      REAL(wp) ::   zfld, zflu, zfld16, zflu16, zfact
      REAL(wp) ::   zph, zah2, zbot, zdic, zalk, zsch_o2, zalka, zsch_co2
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------

      IF( wrk_in_use(2, 1,2,3,4,5,6,7) ) THEN
         CALL ctl_stop('p4z_flx: requested workspace arrays unavailable')   ;   RETURN
      ENDIF

      ! SURFACE CHEMISTRY (PCO2 AND [H+] IN
      !     SURFACE LAYER); THE RESULT OF THIS CALCULATION
      !     IS USED TO COMPUTE AIR-SEA FLUX OF CO2

#if defined key_cpl_carbon_cycle
      satmco2(:,:) = atm_co2(:,:)
#endif

      DO jrorr = 1, 10

!CDIR NOVERRCHK
         DO jj = 1, jpj
!CDIR NOVERRCHK
            DO ji = 1, jpi

               ! DUMMY VARIABLES FOR DIC, H+, AND BORATE
               zbot  = borat(ji,jj,1)
               zfact = rhop(ji,jj,1) / 1000. + rtrn
               zdic  = trn(ji,jj,1,jpdic) / zfact
               zph   = MAX( hi(ji,jj,1), 1.e-10 ) / zfact
               zalka = trn(ji,jj,1,jptal) / zfact

               ! CALCULATE [ALK]([CO3--], [HCO3-])
               zalk  = zalka - (  akw3(ji,jj,1) / zph - zph + zbot / ( 1.+ zph / akb3(ji,jj,1) )  )

               ! CALCULATE [H+] AND [H2CO3]
               zah2   = SQRT(  (zdic-zalk)**2 + 4.* ( zalk * ak23(ji,jj,1)   &
                  &                                        / ak13(ji,jj,1) ) * ( 2.* zdic - zalk )  )
               zah2   = 0.5 * ak13(ji,jj,1) / zalk * ( ( zdic - zalk ) + zah2 )
               zh2co3(ji,jj) = ( 2.* zdic - zalk ) / ( 2.+ ak13(ji,jj,1) / zah2 ) * zfact
               hi(ji,jj,1)   = zah2 * zfact
            END DO
         END DO
      END DO


      ! --------------
      ! COMPUTE FLUXES
      ! --------------

      ! FIRST COMPUTE GAS EXCHANGE COEFFICIENTS
      ! -------------------------------------------

!CDIR NOVERRCHK
      DO jj = 1, jpj
!CDIR NOVERRCHK
         DO ji = 1, jpi
            ztc  = MIN( 35., tsn(ji,jj,1,jp_tem) )
            ztc2 = ztc * ztc
            ztc3 = ztc * ztc2 
            ! Compute the schmidt Number both O2 and CO2
            zsch_co2 = 2073.1 - 125.62 * ztc + 3.6276 * ztc2 - 0.043126 * ztc3
            zsch_o2  = 1953.4 - 128.0  * ztc + 3.9918 * ztc2 - 0.050091 * ztc3
            !  wind speed 
            zws  = wndm(ji,jj) * wndm(ji,jj)
            ! Compute the piston velocity for O2 and CO2
            zkgwan = 0.3 * zws  + 2.5 * ( 0.5246 + 0.016256 * ztc + 0.00049946  * ztc2 )
# if defined key_degrad
            zkgwan = zkgwan * xconv * ( 1.- fr_i(ji,jj) ) * tmask(ji,jj,1) * facvol(ji,jj,1)
#else
            zkgwan = zkgwan * xconv * ( 1.- fr_i(ji,jj) ) * tmask(ji,jj,1)
#endif 
            ! compute gas exchange for CO2 and O2
            zkgco2(ji,jj) = zkgwan * SQRT( 660./ zsch_co2 )
            zkgo2 (ji,jj) = zkgwan * SQRT( 660./ zsch_o2 )
         END DO
      END DO

      DO jj = 1, jpj
         DO ji = 1, jpi
            ! Compute CO2 flux for the sea and air
            zfld = satmco2(ji,jj) * tmask(ji,jj,1) * chemc(ji,jj,1) * zkgco2(ji,jj)
            zflu = zh2co3(ji,jj) * tmask(ji,jj,1) * zkgco2(ji,jj)
            oce_co2(ji,jj) = ( zfld - zflu ) * rfact * e1e2t(ji,jj) * tmask(ji,jj,1) * 1000.
            ! compute the trend
            tra(ji,jj,1,jpdic) = tra(ji,jj,1,jpdic) + ( zfld - zflu ) / fse3t(ji,jj,1)

            ! Compute O2 flux 
            zfld16 = atcox * chemc(ji,jj,2) * tmask(ji,jj,1) * zkgo2(ji,jj)
            zflu16 = trn(ji,jj,1,jpoxy) * tmask(ji,jj,1) * zkgo2(ji,jj)
            tra(ji,jj,1,jpoxy) = tra(ji,jj,1,jpoxy) + ( zfld16 - zflu16 ) / fse3t(ji,jj,1)

#if defined key_diatrc 
            ! Save diagnostics
#  if ! defined key_iomput
            zfact = 1. / e1e2t(ji,jj) / rfact
            trc2d(ji,jj,jp_pcs0_2d    ) = oce_co2(ji,jj) * zfact
            trc2d(ji,jj,jp_pcs0_2d + 1) = ( zfld16 - zflu16 ) * 1000. * tmask(ji,jj,1)
            trc2d(ji,jj,jp_pcs0_2d + 2) = zkgco2(ji,jj) * tmask(ji,jj,1)
            trc2d(ji,jj,jp_pcs0_2d + 3) = ( satmco2(ji,jj) - zh2co3(ji,jj) / ( chemc(ji,jj,1) + rtrn ) ) &
               &                            * tmask(ji,jj,1)
#  else
            zoflx(ji,jj) = ( zfld16 - zflu16 ) * 1000. * tmask(ji,jj,1)
            zkg  (ji,jj) = zkgco2(ji,jj) * tmask(ji,jj,1)
            zdpco2(ji,jj) = ( satmco2(ji,jj) - zh2co3(ji,jj) / ( chemc(ji,jj,1) + rtrn ) ) * tmask(ji,jj,1)
            zdpo2 (ji,jj) = ( atcox  - trn(ji,jj,1,jpoxy) / ( chemc(ji,jj,2) + rtrn ) ) * tmask(ji,jj,1)
#  endif
#endif
         END DO
      END DO

      t_oce_co2_flx = t_oce_co2_flx + glob_sum( oce_co2(:,:) )                     ! Cumulative Total Flux of Carbon
      IF( kt == nitend ) THEN
         t_atm_co2_flx = glob_sum( satmco2(:,:) * e1e2t(:,:) )            ! Total atmospheric pCO2
         !
         t_oce_co2_flx = (-1.) * t_oce_co2_flx  * 12. / 1.e15                      ! Conversion in PgC ; negative for out of the ocean
         t_atm_co2_flx = t_atm_co2_flx  / area                                     ! global mean of atmospheric pCO2
         !
         IF( lwp) THEN
            WRITE(numout,*)
            WRITE(numout,*) ' Global mean of atmospheric pCO2 (ppm) at it= ', kt, ' date= ', ndastp
            WRITE(numout,*) '------------------------------------------------------- :  ',t_atm_co2_flx
            WRITE(numout,*)
            WRITE(numout,*) ' Cumulative total Flux of Carbon out of the ocean (PgC) :'
            WRITE(numout,*) '-------------------------------------------------------  ',t_oce_co2_flx
         ENDIF
         !
      ENDIF

      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('flx ')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF

# if defined key_diatrc && defined key_iomput
      CALL iom_put( "Cflx" , oce_co2(:,:) / e1e2t(:,:) / rfact ) 
      CALL iom_put( "Oflx" , zoflx  )
      CALL iom_put( "Kg"   , zkg    )
      CALL iom_put( "Dpco2", zdpco2 )
      CALL iom_put( "Dpo2" , zdpo2  )
#endif
      !
      IF( wrk_not_released(2, 1,2,3,4,5,6,7) )   CALL ctl_stop('p4z_flx: failed to release workspace arrays')
      !
   END SUBROUTINE p4z_flx


   SUBROUTINE p4z_flx_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_flx_init  ***
      !!
      !! ** Purpose :   Initialization of atmospheric conditions
      !!
      !! ** Method  :   Read the nampisext namelist and check the parameters
      !!      called at the first timestep (nit000)
      !! ** input   :   Namelist nampisext
      !!----------------------------------------------------------------------
      NAMELIST/nampisext/ atcco2
      !!----------------------------------------------------------------------
      !
      REWIND( numnat )                     ! read numnat
      READ  ( numnat, nampisext )
      !
      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for air-sea exchange, nampisext'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    Atmospheric pCO2      atcco2      =', atcco2
      ENDIF
      !
      area = glob_sum( e1e2t(:,:) )        ! interior global domain surface
      !
      oce_co2(:,:)  = 0._wp                ! Initialization of Flux of Carbon
      t_atm_co2_flx = 0._wp
      !
      satmco2(:,:)  = atcco2      ! Initialisation of atmospheric pco2
      t_oce_co2_flx = 0._wp
      !
   END SUBROUTINE p4z_flx_init


   INTEGER FUNCTION p4z_flx_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_flx_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( oce_co2(jpi,jpj), satmco2(jpi,jpj), STAT=p4z_flx_alloc )
      !
      IF( p4z_flx_alloc /= 0 )   CALL ctl_warn('p4z_flx_alloc : failed to allocate arrays')
      !
   END FUNCTION p4z_flx_alloc

#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_flx( kt )                   ! Empty routine
      INTEGER, INTENT( in ) ::   kt
      WRITE(*,*) 'p4z_flx: You should not have seen this print! error?', kt
   END SUBROUTINE p4z_flx
#endif 

   !!======================================================================
END MODULE  p4zflx
