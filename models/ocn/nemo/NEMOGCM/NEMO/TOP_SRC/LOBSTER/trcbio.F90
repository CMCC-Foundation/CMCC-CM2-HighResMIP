MODULE trcbio
   !!======================================================================
   !!                         ***  MODULE trcbio  ***
   !! TOP :   LOBSTER
   !!======================================================================
   !! History :    -   !  1999-07  (M. Levy) Original code
   !!              -   !  2000-12  (E. Kestenare) assign a parameter to name individual tracers
   !!              -   !  2001-03  (M. Levy)  LNO3 + dia2d 
   !!             2.0  !  2007-12  (C. Deltel, G. Madec)  F90
   !!----------------------------------------------------------------------
#if defined key_lobster
   !!----------------------------------------------------------------------
   !!   'key_lobster'                                     LOBSTER bio-model
   !!----------------------------------------------------------------------
   !!   trc_bio        :  
   !!----------------------------------------------------------------------
   USE oce_trc         !
   USE trc             ! 
   USE sms_lobster     ! 
   USE lbclnk          ! 
   USE prtctl_trc      ! Print control for debbuging
   USE trdmod_oce
   USE trdmod_trc
   USE iom
   
   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_bio    ! called in ???

   !!* Substitution
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcbio.F90 2715 2011-03-30 15:58:35Z rblod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trc_bio( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_bio  ***
      !!
      !! ** Purpose :   compute the now trend due to biogeochemical processes
      !!              and add it to the general trend of passive tracers equations
      !!
      !! ** Method  :   each now biological flux is calculated in function of now
      !!              concentrations of tracers.
      !!              depending on the tracer, these fluxes are sources or sinks.
      !!              the total of the sources and sinks for each tracer
      !!              is added to the general trend.
      !!        
      !!                      tra = tra + zf...tra - zftra...
      !!                                     |         |
      !!                                     |         |
      !!                                  source      sink
      !!        
      !!              IF 'key_diabio' defined , the biogeochemical trends
      !!              for passive tracers are saved for futher diagnostics.
      !!---------------------------------------------------------------------
      USE wrk_nemo, ONLY: wrk_in_use,  wrk_not_released
      USE wrk_nemo, ONLY: wrk_3d_2, wrk_4d_1
      !!
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index      
      !!
      INTEGER  ::   ji, jj, jk, jl
      REAL(wp) ::   zdet, zzoo, zphy, zno3, znh4, zdom      ! now concentrations
      REAL(wp) ::   zlno3, zlnh4, zle, zlt                  ! limitation terms for phyto
      REAL(wp) ::   zno3phy, znh4phy, zphynh4, zphydom
      REAL(wp) ::   zphydet, zphyzoo, zdetzoo
      REAL(wp) ::   zzoonh4, zzoodom, zzoodet, zdetnh4, zdetdom
      REAL(wp) ::   znh4no3, zdomnh4, zppz, zpdz, zpppz, zppdz, zfood
      REAL(wp) ::   zfilpz, zfildz, zphya, zzooa, zno3a
      REAL(wp) ::   znh4a, zdeta, zdoma, zzoobod, zboddet, zdomaju
#if defined key_diatrc
      REAL(wp) ::   ze3t
#endif
#if defined key_diatrc && defined key_iomput
      REAL(wp), POINTER,   DIMENSION(:,:,:) :: zw2d
      REAL(wp), POINTER, DIMENSION(:,:,:,:) :: zw3d
#endif
      REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::   ztrbio
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------

#if defined key_diatrc && defined key_iomput
      IF( ( wrk_in_use(3, 2) ) .OR. ( wrk_in_use(4, 1) ) ) THEN
         CALL ctl_stop('trc_bio : requested workspace arrays unavailable.')
         RETURN
      END IF
      ! Set-up pointers into sub-arrays of workspaces
      zw2d => wrk_3d_2(:,:,1:17)
      zw3d => wrk_4d_1(:,:,:,1:3)
#endif

      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) ' trc_bio: LOBSTER bio-model'
         IF(lwp) WRITE(numout,*) ' ~~~~~~~'
      ENDIF

      fbod(:,:) = 0.e0
#if defined key_diatrc && ! defined key_iomput
#  if defined key_iomput
      zw2d  (:,:,:) = 0.e0
      zw3d(:,:,:,:) = 0.e0
#  else
      DO jl = jp_lob0_2d, jp_lob1_2d
         trc2d(:,:,jl) = 0.e0
      END DO 
#  endif
#endif

      IF( l_trdtrc )THEN
         ALLOCATE( ztrbio(jpi,jpj,jpk,jp_lobster_trd) )
         ztrbio(:,:,:,:) = 0.
      ENDIF

      !                                      ! -------------------------- !
      DO jk = 1, jpkbm1                      !  Upper ocean (bio-layers)  !
         !                                   ! -------------------------- !
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1 
               ! trophic variables( det, zoo, phy, no3, nh4, dom)
               ! ------------------------------------------------

               ! negative trophic variables DO not contribute to the fluxes
               zdet = MAX( 0.e0, trn(ji,jj,jk,jp_lob_det) )
               zzoo = MAX( 0.e0, trn(ji,jj,jk,jp_lob_zoo) )
               zphy = MAX( 0.e0, trn(ji,jj,jk,jp_lob_phy) )
               zno3 = MAX( 0.e0, trn(ji,jj,jk,jp_lob_no3) )
               znh4 = MAX( 0.e0, trn(ji,jj,jk,jp_lob_nh4) )
               zdom = MAX( 0.e0, trn(ji,jj,jk,jp_lob_dom) )

               ! Limitations
               zlt   = 1.
               zle   = 1. - EXP( -xpar(ji,jj,jk) / aki / zlt )
               ! psinut,akno3,aknh4 added by asklod AS Kremeur 2005-03
               zlno3 = zno3 * EXP( -psinut * znh4 ) / ( akno3 + zno3 )
               zlnh4 = znh4 / (znh4+aknh4) 

               ! sinks and sources
               !    phytoplankton production and exsudation
               zno3phy = tmumax * zle * zlt * zlno3 * zphy
               znh4phy = tmumax * zle * zlt * zlnh4 * zphy

               !    fphylab added by asklod AS Kremeur 2005-03
               zphydom = rgamma * (1 - fphylab) * (zno3phy + znh4phy)
               zphynh4 = rgamma * fphylab * (zno3phy + znh4phy)

               ! zooplankton production
               !    preferences
               zppz = rppz
               zpdz = 1. - rppz
               zpppz = ( zppz * zphy ) / ( ( zppz * zphy + zpdz * zdet ) + 1.e-13 )
               zppdz = ( zpdz * zdet ) / ( ( zppz * zphy + zpdz * zdet ) + 1.e-13 )
               zfood = zpppz * zphy + zppdz * zdet
               !    filtration
               zfilpz = taus * zpppz / (aks + zfood)
               zfildz = taus * zppdz / (aks + zfood)
               !    grazing
               zphyzoo = zfilpz * zphy * zzoo
               zdetzoo = zfildz * zdet * zzoo

               ! fecal pellets production
               zzoodet = rpnaz * zphyzoo + rdnaz * zdetzoo
 
               ! zooplankton liquide excretion
               zzoonh4 = tauzn * fzoolab * zzoo 
               zzoodom = tauzn * (1 - fzoolab) * zzoo

               ! mortality
               !    phytoplankton mortality 
               zphydet = tmminp * zphy

               !    zooplankton mortality
               !    closure : flux fbod is redistributed below level jpkbio
               zzoobod = tmminz * zzoo * zzoo
               fbod(ji,jj) = fbod(ji,jj) + (1-fdbod) * zzoobod * fse3t(ji,jj,jk)
               zboddet = fdbod * zzoobod

               ! detritus and dom breakdown
               zdetnh4 = taudn * fdetlab * zdet
               zdetdom = taudn * (1 - fdetlab) * zdet 

               zdomnh4 = taudomn * zdom

               ! flux added to express how the excess of nitrogen from
               ! PHY, ZOO and DET to DOM goes directly to NH4 (flux of ajustment)
               zdomaju = (1 - redf/reddom) * (zphydom + zzoodom + zdetdom)

               ! Nitrification
               znh4no3 = taunn * znh4

               ! determination of trends
               !    total trend for each biological tracer
               zphya =   zno3phy + znh4phy - zphynh4 - zphydom - zphyzoo - zphydet
               zzooa =   zphyzoo + zdetzoo - zzoodet - zzoodom - zzoonh4 - zzoobod
               zno3a = - zno3phy + znh4no3
               znh4a = - znh4phy - znh4no3 + zphynh4 + zzoonh4 + zdomnh4 + zdetnh4 + zdomaju
               zdeta =   zphydet + zzoodet - zdetzoo - zdetnh4 - zdetdom + zboddet
               zdoma =   zphydom + zzoodom + zdetdom - zdomnh4 - zdomaju

               ! tracer flux at totox-point added to the general trend
               tra(ji,jj,jk,jp_lob_det) = tra(ji,jj,jk,jp_lob_det) + zdeta
               tra(ji,jj,jk,jp_lob_zoo) = tra(ji,jj,jk,jp_lob_zoo) + zzooa
               tra(ji,jj,jk,jp_lob_phy) = tra(ji,jj,jk,jp_lob_phy) + zphya
               tra(ji,jj,jk,jp_lob_no3) = tra(ji,jj,jk,jp_lob_no3) + zno3a
               tra(ji,jj,jk,jp_lob_nh4) = tra(ji,jj,jk,jp_lob_nh4) + znh4a
               tra(ji,jj,jk,jp_lob_dom) = tra(ji,jj,jk,jp_lob_dom) + zdoma

#if defined key_diabio
               trbio(ji,jj,jk,jp_lob0_trd     ) = zno3phy
               trbio(ji,jj,jk,jp_lob0_trd +  1) = znh4phy
               trbio(ji,jj,jk,jp_lob0_trd +  2) = zphynh4
               trbio(ji,jj,jk,jp_lob0_trd +  3) = zphydom
               trbio(ji,jj,jk,jp_lob0_trd +  4) = zphyzoo
               trbio(ji,jj,jk,jp_lob0_trd +  5) = zphydet
               trbio(ji,jj,jk,jp_lob0_trd +  6) = zdetzoo
               trbio(ji,jj,jk,jp_lob0_trd +  8) = zzoodet
               trbio(ji,jj,jk,jp_lob0_trd +  9) = zzoobod
               trbio(ji,jj,jk,jp_lob0_trd + 10) = zzoonh4
               trbio(ji,jj,jk,jp_lob0_trd + 11) = zzoodom
               trbio(ji,jj,jk,jp_lob0_trd + 12) = znh4no3
               trbio(ji,jj,jk,jp_lob0_trd + 13) = zdomnh4
               trbio(ji,jj,jk,jp_lob0_trd + 14) = zdetnh4
               trbio(ji,jj,jk,jp_lob0_trd + 15) = zdetdom
#endif
               IF( l_trdtrc ) THEN
                  ztrbio(ji,jj,jk,jp_lob0_trd     ) = zno3phy
                  ztrbio(ji,jj,jk,jp_lob0_trd +  1) = znh4phy
                  ztrbio(ji,jj,jk,jp_lob0_trd +  2) = zphynh4
                  ztrbio(ji,jj,jk,jp_lob0_trd +  3) = zphydom
                  ztrbio(ji,jj,jk,jp_lob0_trd +  4) = zphyzoo
                  ztrbio(ji,jj,jk,jp_lob0_trd +  5) = zphydet
                  ztrbio(ji,jj,jk,jp_lob0_trd +  6) = zdetzoo
                  !  trend number 8 in trcsed
                  ztrbio(ji,jj,jk,jp_lob0_trd +  8) = zzoodet
                  ztrbio(ji,jj,jk,jp_lob0_trd +  9) = zzoobod
                  ztrbio(ji,jj,jk,jp_lob0_trd + 10) = zzoonh4
                  ztrbio(ji,jj,jk,jp_lob0_trd + 11) = zzoodom
                  ztrbio(ji,jj,jk,jp_lob0_trd + 12) = znh4no3
                  ztrbio(ji,jj,jk,jp_lob0_trd + 13) = zdomnh4
                  ztrbio(ji,jj,jk,jp_lob0_trd + 14) = zdetnh4
                  ztrbio(ji,jj,jk,jp_lob0_trd + 15) = zdetdom
                  !  trend number 17 in trcexp
                ENDIF

#if defined key_diatrc
               ! convert fluxes in per day
               ze3t = fse3t(ji,jj,jk) * 86400.
#if ! defined key_iomput
               trc2d(ji,jj,jp_lob0_2d    ) = trc2d(ji,jj, jp_lob0_2d    ) + zno3phy * ze3t 
               trc2d(ji,jj,jp_lob0_2d + 1) = trc2d(ji,jj, jp_lob0_2d + 1) + znh4phy * ze3t
               trc2d(ji,jj,jp_lob0_2d + 2) = trc2d(ji,jj, jp_lob0_2d + 2) + zphydom * ze3t
               trc2d(ji,jj,jp_lob0_2d + 3) = trc2d(ji,jj, jp_lob0_2d + 3) + zphynh4 * ze3t
               trc2d(ji,jj,jp_lob0_2d + 4) = trc2d(ji,jj, jp_lob0_2d + 4) + zphyzoo * ze3t
               trc2d(ji,jj,jp_lob0_2d + 5) = trc2d(ji,jj, jp_lob0_2d + 5) + zphydet * ze3t
               trc2d(ji,jj,jp_lob0_2d + 6) = trc2d(ji,jj, jp_lob0_2d + 6) + zdetzoo * ze3t
               ! trend number 8 is in trcsed.F            
               trc2d(ji,jj,jp_lob0_2d +  8) = trc2d(ji,jj,jp_lob0_2d +  8) + zzoodet * ze3t
               trc2d(ji,jj,jp_lob0_2d +  9) = trc2d(ji,jj,jp_lob0_2d +  9) + zzoobod * ze3t
               trc2d(ji,jj,jp_lob0_2d + 10) = trc2d(ji,jj,jp_lob0_2d + 10) + zzoonh4 * ze3t
               trc2d(ji,jj,jp_lob0_2d + 11) = trc2d(ji,jj,jp_lob0_2d + 11) + zzoodom * ze3t
               trc2d(ji,jj,jp_lob0_2d + 12) = trc2d(ji,jj,jp_lob0_2d + 12) + znh4no3 * ze3t
               trc2d(ji,jj,jp_lob0_2d + 13) = trc2d(ji,jj,jp_lob0_2d + 13) + zdomnh4 * ze3t
               trc2d(ji,jj,jp_lob0_2d + 14) = trc2d(ji,jj,jp_lob0_2d + 14) + zdetnh4 * ze3t             
               trc2d(ji,jj,jp_lob0_2d + 15) = trc2d(ji,jj,jp_lob0_2d + 15) + (  zno3phy + znh4phy - zphynh4   &
                  &                                 - zphydom - zphyzoo - zphydet ) * ze3t
               trc2d(ji,jj,jp_lob0_2d + 16) = trc2d(ji,jj,jp_lob0_2d + 16) + (  zphyzoo + zdetzoo - zzoodet   &
                  &                                 - zzoobod - zzoonh4 - zzoodom ) * ze3t
               trc2d(ji,jj,jp_lob0_2d + 17) = trc2d(ji,jj,jp_lob0_2d + 17) + zdetdom * ze3t
               ! trend number 19 is in trcexp.F
#else
               zw2d(ji,jj,1)  = zw2d(ji,jj,1)  + zno3phy * ze3t 
               zw2d(ji,jj,2)  = zw2d(ji,jj,2)  + znh4phy * ze3t
               zw2d(ji,jj,3)  = zw2d(ji,jj,3)  + zphydom * ze3t
               zw2d(ji,jj,4)  = zw2d(ji,jj,4)  + zphynh4 * ze3t
               zw2d(ji,jj,5)  = zw2d(ji,jj,5)  + zphyzoo * ze3t
               zw2d(ji,jj,6)  = zw2d(ji,jj,6)  + zphydet * ze3t
               zw2d(ji,jj,7)  = zw2d(ji,jj,7)  + zdetzoo * ze3t
               zw2d(ji,jj,8)  = zw2d(ji,jj,8)  + zzoodet * ze3t
               zw2d(ji,jj,9)  = zw2d(ji,jj,9)  + zzoobod * ze3t
               zw2d(ji,jj,10) = zw2d(ji,jj,10) + zzoonh4 * ze3t
               zw2d(ji,jj,11) = zw2d(ji,jj,11) + zzoodom * ze3t
               zw2d(ji,jj,12) = zw2d(ji,jj,12) + znh4no3 * ze3t
               zw2d(ji,jj,13) = zw2d(ji,jj,13) + zdomnh4 * ze3t
               zw2d(ji,jj,14) = zw2d(ji,jj,14) + zdetnh4 * ze3t             
               zw2d(ji,jj,15) = zw2d(ji,jj,15) + ( zno3phy + znh4phy - zphynh4 - zphydom - zphyzoo - zphydet ) * ze3t
               zw2d(ji,jj,16) = zw2d(ji,jj,16) + ( zphyzoo + zdetzoo - zzoodet - zzoobod - zzoonh4 - zzoodom ) * ze3t
               zw2d(ji,jj,17) = zw2d(ji,jj,17) + zdetdom * ze3t
#endif
#if defined key_diatrc 
# if ! defined key_iomput
               trc3d(ji,jj,jk,jp_lob0_3d    ) = zno3phy * 86400     
               trc3d(ji,jj,jk,jp_lob0_3d + 1) = znh4phy * 86400     
               trc3d(ji,jj,jk,jp_lob0_3d + 2) = znh4no3 * 86400   
# else
               zw3d(ji,jj,jk,1) = zno3phy * 86400     
               zw3d(ji,jj,jk,2) = znh4phy * 86400     
               zw3d(ji,jj,jk,3) = znh4no3 * 86400   
# endif
#endif  
#endif
            END DO
         END DO
      END DO

      !                                      ! -------------------------- !
      DO jk = jpkb, jpkm1                    !  Upper ocean (bio-layers)  !
         !                                   ! -------------------------- !
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1 
               ! remineralisation of all quantities towards nitrate 

               !    trophic variables( det, zoo, phy, no3, nh4, dom)
               !       negative trophic variables DO not contribute to the fluxes
               zdet = MAX( 0.e0, trn(ji,jj,jk,jp_lob_det) )
               zzoo = MAX( 0.e0, trn(ji,jj,jk,jp_lob_zoo) )
               zphy = MAX( 0.e0, trn(ji,jj,jk,jp_lob_phy) )
               zno3 = MAX( 0.e0, trn(ji,jj,jk,jp_lob_no3) )
               znh4 = MAX( 0.e0, trn(ji,jj,jk,jp_lob_nh4) )
               zdom = MAX( 0.e0, trn(ji,jj,jk,jp_lob_dom) )

               !    Limitations
               zlt   = 0.e0
               zle   = 0.e0
               zlno3 = 0.e0
               zlnh4 = 0.e0

               !    sinks and sources
               !       phytoplankton production and exsudation
               zno3phy = 0.e0
               znh4phy = 0.e0
               zphydom = 0.e0
               zphynh4 = 0.e0

               !    zooplankton production
               zphyzoo = 0.e0      ! grazing
               zdetzoo = 0.e0

               zzoodet = 0.e0      ! fecal pellets production

               zzoonh4 = tauzn * fzoolab * zzoo         ! zooplankton liquide excretion
               zzoodom = tauzn * (1 - fzoolab) * zzoo

               !    mortality
               zphydet = tmminp * zphy      ! phytoplankton mortality 

               zzoobod = 0.e0               ! zooplankton mortality
               zboddet = 0.e0               ! closure : flux fbod is redistributed below level jpkbio

               !    detritus and dom breakdown
               zdetnh4 = taudn * fdetlab * zdet
               zdetdom = taudn * (1 - fdetlab) * zdet 

               zdomnh4 = taudomn * zdom
               zdomaju = (1 - redf/reddom) * (zphydom + zzoodom + zdetdom)

               !    Nitrification
               znh4no3 = taunn * znh4


               ! determination of trends
               !     total trend for each biological tracer
               zphya =   zno3phy + znh4phy - zphynh4 - zphydom - zphyzoo - zphydet
               zzooa =   zphyzoo + zdetzoo - zzoodet - zzoodom - zzoonh4 - zzoobod
               zno3a = - zno3phy + znh4no3
               znh4a = - znh4phy - znh4no3 + zphynh4 + zzoonh4 + zdomnh4 + zdetnh4 + zdomaju
               zdeta = zphydet + zzoodet  - zdetzoo - zdetnh4 - zdetdom + zboddet
               zdoma = zphydom + zzoodom + zdetdom - zdomnh4 - zdomaju

               ! tracer flux at totox-point added to the general trend
               tra(ji,jj,jk,jp_lob_det) = tra(ji,jj,jk,jp_lob_det) + zdeta
               tra(ji,jj,jk,jp_lob_zoo) = tra(ji,jj,jk,jp_lob_zoo) + zzooa
               tra(ji,jj,jk,jp_lob_phy) = tra(ji,jj,jk,jp_lob_phy) + zphya
               tra(ji,jj,jk,jp_lob_no3) = tra(ji,jj,jk,jp_lob_no3) + zno3a
               tra(ji,jj,jk,jp_lob_nh4) = tra(ji,jj,jk,jp_lob_nh4) + znh4a
               tra(ji,jj,jk,jp_lob_dom) = tra(ji,jj,jk,jp_lob_dom) + zdoma
               !
#if defined key_diabio
               trbio(ji,jj,jk,jp_lob0_trd     ) = zno3phy
               trbio(ji,jj,jk,jp_lob0_trd +  1) = znh4phy
               trbio(ji,jj,jk,jp_lob0_trd +  2) = zphynh4
               trbio(ji,jj,jk,jp_lob0_trd +  3) = zphydom
               trbio(ji,jj,jk,jp_lob0_trd +  4) = zphyzoo
               trbio(ji,jj,jk,jp_lob0_trd +  5) = zphydet
               trbio(ji,jj,jk,jp_lob0_trd +  6) = zdetzoo
               trbio(ji,jj,jk,jp_lob0_trd +  8) = zzoodet
               trbio(ji,jj,jk,jp_lob0_trd +  9) = zzoobod
               trbio(ji,jj,jk,jp_lob0_trd + 10) = zzoonh4
               trbio(ji,jj,jk,jp_lob0_trd + 11) = zzoodom
               trbio(ji,jj,jk,jp_lob0_trd + 12) = znh4no3
               trbio(ji,jj,jk,jp_lob0_trd + 13) = zdomnh4
               trbio(ji,jj,jk,jp_lob0_trd + 14) = zdetnh4
               trbio(ji,jj,jk,jp_lob0_trd + 15) = zdetdom
#endif
               IF( l_trdtrc ) THEN
                  ztrbio(ji,jj,jk,jp_lob0_trd     ) = zno3phy
                  ztrbio(ji,jj,jk,jp_lob0_trd +  1) = znh4phy
                  ztrbio(ji,jj,jk,jp_lob0_trd +  2) = zphynh4
                  ztrbio(ji,jj,jk,jp_lob0_trd +  3) = zphydom
                  ztrbio(ji,jj,jk,jp_lob0_trd +  4) = zphyzoo
                  ztrbio(ji,jj,jk,jp_lob0_trd +  5) = zphydet
                  ztrbio(ji,jj,jk,jp_lob0_trd +  6) = zdetzoo
                  !  trend number 8 in trcsed
                  ztrbio(ji,jj,jk,jp_lob0_trd +  8) = zzoodet
                  ztrbio(ji,jj,jk,jp_lob0_trd +  9) = zzoobod
                  ztrbio(ji,jj,jk,jp_lob0_trd + 10) = zzoonh4
                  ztrbio(ji,jj,jk,jp_lob0_trd + 11) = zzoodom
                  ztrbio(ji,jj,jk,jp_lob0_trd + 12) = znh4no3
                  ztrbio(ji,jj,jk,jp_lob0_trd + 13) = zdomnh4
                  ztrbio(ji,jj,jk,jp_lob0_trd + 14) = zdetnh4
                  ztrbio(ji,jj,jk,jp_lob0_trd + 15) = zdetdom
                  !  trend number 17 in trcexp
                ENDIF
#if defined key_diatrc
# if ! defined key_iomput
               trc3d(ji,jj,jk,jp_lob0_3d    ) =  zno3phy * 86400     
               trc3d(ji,jj,jk,jp_lob0_3d + 1) =  znh4phy * 86400     
               trc3d(ji,jj,jk,jp_lob0_3d + 2) =  znh4no3 * 86400     
# else
               zw3d(ji,jj,jk,1) = zno3phy * 86400     
               zw3d(ji,jj,jk,2) = znh4phy * 86400     
               zw3d(ji,jj,jk,3) = znh4no3 * 86400   
# endif
#endif
            END DO
         END DO
      END DO

#if defined key_diatrc
      ! Lateral boundary conditions 
# if ! defined key_iomput
      DO jl = jp_lob0_2d, jp_lob1_2d
          CALL lbc_lnk( trc2d(:,:,jl),'T', 1. )
      END DO 
# else
      DO jl = 1, 17 
          CALL lbc_lnk( zw2d(:,:,jl),'T', 1. )
      END DO
      ! Save diagnostics
      CALL iom_put( "TNO3PHY", zw2d(:,:,1) )
      CALL iom_put( "TNH4PHY", zw2d(:,:,2) )
      CALL iom_put( "TPHYDOM", zw2d(:,:,3) )
      CALL iom_put( "TPHYNH4", zw2d(:,:,4) )
      CALL iom_put( "TPHYZOO", zw2d(:,:,5) )
      CALL iom_put( "TPHYDET", zw2d(:,:,6) )
      CALL iom_put( "TDETZOO", zw2d(:,:,7) )
      CALL iom_put( "TZOODET", zw2d(:,:,8) )
      CALL iom_put( "TZOOBOD", zw2d(:,:,9) )
      CALL iom_put( "TZOONH4", zw2d(:,:,10) )
      CALL iom_put( "TZOODOM", zw2d(:,:,11) )
      CALL iom_put( "TNH4NO3", zw2d(:,:,12) )
      CALL iom_put( "TDOMNH4", zw2d(:,:,13) )
      CALL iom_put( "TDETNH4", zw2d(:,:,14) )
      CALL iom_put( "TPHYTOT", zw2d(:,:,15) )
      CALL iom_put( "TZOOTOT", zw2d(:,:,16) )
      CALL iom_put( "TDETDOM", zw2d(:,:,17) )
# endif
#endif

#if defined key_diatrc
      ! Lateral boundary conditions 
# if ! defined key_iomput
      DO jl = jp_lob0_3d, jp_lob1_3d
          CALL lbc_lnk( trc3d(:,:,1,jl),'T', 1. )
      END DO 
# else
      DO jl = 1, 3
          CALL lbc_lnk( zw3d(:,:,:,jl),'T', 1. )
      END DO
      ! save diagnostics
      CALL iom_put( "FNO3PHY", zw3d(:,:,:,1) )
      CALL iom_put( "FNH4PHY", zw3d(:,:,:,2) )
      CALL iom_put( "FNH4NO3", zw3d(:,:,:,3) )
# endif 
#endif

#if defined key_diabio
      ! Lateral boundary conditions on trcbio
      DO jl = jp_lob0_trd, jp_lob1_trd
          CALL lbc_lnk( trbio(:,:,1,jl),'T', 1. )
      END DO 
#endif
      !
      IF( l_trdtrc ) THEN
         DO jl = jp_lob0_trd, jp_lob1_trd
            CALL trd_mod_trc( ztrbio(:,:,:,jl), jl, kt )   ! handle the trend
         END DO
      ENDIF

      IF( l_trdtrc ) DEALLOCATE( ztrbio )

      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('bio')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
#if defined key_diatrc && defined key_iomput
      IF( ( wrk_not_released(3, 2) ) .OR. ( wrk_not_released(4, 1) ) )  &
        &   CALL ctl_stop('trc_bio : failed to release workspace arrays.')
#endif
      !
   END SUBROUTINE trc_bio

#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE trc_bio( kt )                   ! Empty routine
      INTEGER, INTENT( in ) ::   kt
      WRITE(*,*) 'trc_bio: You should not have seen this print! error?', kt
   END SUBROUTINE trc_bio
#endif 

   !!======================================================================
END MODULE  trcbio
