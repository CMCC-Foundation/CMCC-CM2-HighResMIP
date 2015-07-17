MODULE trcnam_lobster
   !!======================================================================
   !!                      ***  MODULE trcnam_lobster  ***
   !! TOP :   initialisation of some run parameters for LOBSTER bio-model
   !!======================================================================
   !! History :   2.0  !  2007-12  (C. Ethe, G. Madec) from trcnam.lobster1.h90
   !!----------------------------------------------------------------------
#if defined key_lobster
   !!----------------------------------------------------------------------
   !!   'key_lobster'   :                                 LOBSTER bio-model
   !!----------------------------------------------------------------------
   !! trc_nam_lobster   : LOBSTER model namelist read
   !!----------------------------------------------------------------------
   USE oce_trc          ! Ocean variables
   USE par_trc          ! TOP parameters
   USE trc              ! TOP variables
   USE sms_lobster      ! sms trends

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_nam_lobster   ! called by trcnam.F90 module

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcnam_lobster.F90 2715 2011-03-30 15:58:35Z rblod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trc_nam_lobster
      !!----------------------------------------------------------------------
      !!                     ***  trc_nam_lobster  ***  
      !!
      !! ** Purpose :   read LOBSTER namelist
      !!
      !! ** input   :   file 'namelist.trc.sms' containing the following
      !!             namelist: natbio, natopt, and natdbi ("key_diabio")
      !!----------------------------------------------------------------------
      INTEGER ::   numnatl
      !!
#if defined key_diatrc && ! defined key_iomput
      INTEGER :: jl, jn
      ! definition of additional diagnostic as a structure
      TYPE DIAG
         CHARACTER(len = 20)  :: snamedia   !: short name
         CHARACTER(len = 80 ) :: lnamedia   !: long name
         CHARACTER(len = 20 ) :: unitdia    !: unit
      END TYPE DIAG

      TYPE(DIAG) , DIMENSION(jp_lobster_2d) :: lobdia2d
      TYPE(DIAG) , DIMENSION(jp_lobster_3d) :: lobdia3d
#endif
#if defined key_diabio || defined key_trdmld_trc
      INTEGER :: js, jd
      ! definition of additional diagnostic as a structure
      TYPE DIABIO
         CHARACTER(len = 20)  :: snamebio   !: short name
         CHARACTER(len = 80 ) :: lnamebio   !: long name
         CHARACTER(len = 20 ) :: unitbio    !: unit
      END TYPE DIABIO

      TYPE(DIABIO) , DIMENSION(jp_lobster_trd) :: lobdiabio
#endif

      NAMELIST/namlobphy/ apmin, tmumax, rgamma, fphylab, tmmaxp, tmminp, &
         &                rcchl, aki, toptp 
      NAMELIST/namlobnut/ anmin, akno3, aknh4, taunn, psinut
      NAMELIST/namlobzoo/ azmin, eggzoo, rgz, rppz, taus, aks, rpnaz,     &
         &                rdnaz, tauzn, fzoolab, fdbod, tmmaxz, tmminz
      NAMELIST/namlobdet/ admin, taudn, fdetlab, vsed
      NAMELIST/namlobdom/ taudomn
      NAMELIST/namlobsed/ sedlam, sedlostpoc
      NAMELIST/namlobrat/ redf, reddom, slopet, tmaxr, tminr, xhr,        &
         &                filmax, toptgz, tmaxgz, anumin, afdmin

      NAMELIST/namlobopt/ xkg0, xkr0, xkgp, xkrp, xlg, xlr, rpig
#if defined key_diatrc && ! defined key_iomput
      NAMELIST/namlobdia/nn_writedia, lobdia3d, lobdia2d     ! additional diagnostics
#endif
#if defined key_diabio || defined key_trdmld_trc
      NAMELIST/namlobdbi/nwritebio, lobdiabio
#endif
      !!----------------------------------------------------------------------

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_nam_lobster : read LOBSTER namelists'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~~'

      !                               ! Open the namelist file
      !                               ! ----------------------
      CALL ctl_opn( numnatl, 'namelist_lobster', 'OLD', 'FORMATTED', 'SEQUENTIAL', 1, numout, .FALSE. )

      ! namlobphy : parameters for phytoplankton
      apmin   = 0. 
      tmumax  = 0. 
      rgamma  = 0. 
      fphylab = 0.  
      tmmaxp  = 0.
      tmminp  = 0.
      rcchl   = 0. 
      aki     = 0. 
      toptp   = 0. 

      REWIND( numnatl )
      READ  ( numnatl, namlobphy )

      IF(lwp) THEN
          WRITE(numout,*) ' Namelist namlobphy'
          WRITE(numout,*) '    minimum phytoplancton concentration                  apmin     =', apmin
          WRITE(numout,*) '    phyto max growth rate                                tmumax    =', 86400 * tmumax, ' d'
          WRITE(numout,*) '    phytoplankton exudation fraction                     rgamma    =', rgamma
          WRITE(numout,*) '    NH4 fraction of phytoplankton exsudation             fphylab   =', fphylab
          WRITE(numout,*) '    maximal phyto mortality rate                         tmmaxp    =', 86400 * tmmaxp
          WRITE(numout,*) '    minimal phyto mortality rate                         tmminp    =', 86400 * tmminp
          WRITE(numout,*) '    carbone/chlorophyl ratio                             rcchl     =', rcchl
          WRITE(numout,*) '    light hlaf saturation constant                       aki       =', aki
          WRITE(numout,*) '    optimal photosynthesis temperature                   toptp     =', toptp
          WRITE(numout,*) ' '
      ENDIF

      ! namlobnut : parameters for nutrients
      anmin  = 0.
      psinut = 0.
      akno3  = 0.
      aknh4  = 0.
      taunn  = 0.

      REWIND( numnatl )
      READ  ( numnatl, namlobnut )
      IF(lwp) THEN
          WRITE(numout,*) ' Namelist namlobnut'
          WRITE(numout,*) '    minimum nutrients     concentration                  anmin     =', anmin
          WRITE(numout,*) '    half-saturation nutrient for no3 uptake              akno3     =', akno3
          WRITE(numout,*) '    half-saturation nutrient for nh4 uptake              aknh4     =', aknh4
          WRITE(numout,*) '    nitrification rate                                   taunn     =', taunn
          WRITE(numout,*) '    inhibition of no3 uptake by nh4                      psinut    =', psinut
          WRITE(numout,*) ' '
      ENDIF

      ! namlobzoo : parameters for zooplankton
      azmin   = 0.
      rgz     = 0.
      rppz    = 0.
      taus    = 0.
      aks     = 0.
      rpnaz   = 0.
      rdnaz   = 0.
      eggzoo  = 0.
      tauzn   = 0.
      tmmaxz  = 0.
      tmminz  = 0.
      fzoolab = 0.
      fdbod   = 0.

      REWIND( numnatl )
      READ  ( numnatl, namlobzoo )

      IF(lwp) THEN
          WRITE(numout,*) ' Namelist namlobzoo'
          WRITE(numout,*) '    minimum zooplancton   concentration                  azmin     =', azmin
          WRITE(numout,*) '    minimum  for zoo concentration                       eggzoo    =', eggzoo
          WRITE(numout,*) '    widtht of zoo temperature FUNCTION                   rgz       =', rgz
          WRITE(numout,*) '    zoo preference for phyto                             rppz      =', rppz
          WRITE(numout,*) '    maximal zoo grazing rate                             taus      =', 86400 * taus, ' d'
          WRITE(numout,*) '    half saturation constant for zoo food                aks       =', aks
          WRITE(numout,*) '    non-assimilated phyto by zoo                         rpnaz     =', rpnaz
          WRITE(numout,*) '    non-assimilated detritus by zoo                      rdnaz     =', rdnaz
          WRITE(numout,*) '    zoo specific excretion rate                          tauzn     =', 86400 * tauzn
          WRITE(numout,*) '    maximal zoo mortality rate                           tmmaxz    =', 86400 * tmmaxz
          WRITE(numout,*) '    minimal zoo mortality rate                           tmminz    =', 86400 * tmminz
          WRITE(numout,*) '    NH4 fraction of zooplankton excretion                fzoolab   =', fzoolab
          WRITE(numout,*) '    Zooplankton mortality fraction that goes to detritus fdbod     =', fdbod
          WRITE(numout,*) ' '
      ENDIF

      ! namlobdet : parameters for detritus
      admin   = 0.
      taudn   = 0.
      vsed    = 0.
      fdetlab = 0.

      REWIND( numnatl )
      READ  ( numnatl, namlobdet )

      IF(lwp) THEN
          WRITE(numout,*) ' Namelist namlobdet'
          WRITE(numout,*) '    minimum detritus      concentration                  admin     =', admin
          WRITE(numout,*) '    detrital breakdown rate                              taudn     =', 86400 * taudn , ' d'
          WRITE(numout,*) '    detritus sedimentation speed                         vsed      =', 86400 * vsed  , ' d'
          WRITE(numout,*) '    NH4 fraction of detritus dissolution                 fdetlab   =', fdetlab
          WRITE(numout,*) ' '
      ENDIF

      ! namlobdom : parameters for DOM
      taudomn = 0.

      REWIND( numnatl ) 
      READ  ( numnatl, namlobdom )

      IF(lwp) THEN
          WRITE(numout,*) ' Namelist namlobdom'
          WRITE(numout,*) '    dom remineralisation rate                            taudomn   =', taudomn
          WRITE(numout,*) ' '
      ENDIF

      ! namlobsed : parameters from aphotic layers to sediment
      sedlam     = 0.
      sedlostpoc = 0.

      REWIND( numnatl )
      READ  ( numnatl, namlobsed )

      IF(lwp) THEN
          WRITE(numout,*) ' Namelist namlobsed'
          WRITE(numout,*) '    time coeff of POC in sediments                       sedlam    =', sedlam
          WRITE(numout,*) '    Sediment geol loss for POC                           sedlostpoc=', sedlostpoc
          WRITE(numout,*) ' '
      ENDIF

      ! namlobrat : general coefficient
      redf   = 0.
      reddom = 0.
      slopet = 0.
      tmaxr  = 1./(     4.*rday)*0.
      tminr  = 1./(24.*30.*rday)*0.
      xhr    = 0.
      filmax = 0.
      toptgz = 0.
      tmaxgz = 0.
      anumin = 0.
      afdmin = 0.

      REWIND( numnatl )
      READ  ( numnatl, namlobrat )

      IF(lwp) THEN
          WRITE(numout,*) ' Namelist namlobrat'
          WRITE(numout,*) '    redfield ratio  c:n for phyto                        redf      =', redf 
          WRITE(numout,*) '    redfield ratio  c:n for DOM                          reddom    =', reddom 
          WRITE(numout,*) '    van t hoff coefficient                               slopet    =', slopet
          WRITE(numout,*) '    maximum damping for d z or p                         tmaxr     =', tmaxr
          WRITE(numout,*) '    damping-remineralisation rate                        tminr     =', tminr
          WRITE(numout,*) '    coeff for martin''s remineralistion                  xhr       =', xhr
          WRITE(numout,*) '    maximal mass clearance rate for zoo                  filmax    =', filmax
          WRITE(numout,*) '    optimal temperature for zoo growth                   toptgz    =', toptgz
          WRITE(numout,*) '    maximal temperature for zoo growth                   tmaxgz    =', tmaxgz
          WRITE(numout,*) '    nutrient threshold for phyto mort                    anumin    =', anumin
          WRITE(numout,*) '    food threshold for zoo mort                          afdmin    =', afdmin
          WRITE(numout,*) ' '
      ENDIF


      ! namlobopt : optical parameters
      xkg0  = 0. 
      xkr0  = 0.
      xkgp  = 0.
      xkrp  = 0.
      xlg   = 0.
      xlr   = 0.
      rpig  = 0.

      REWIND( numnatl )
      READ  ( numnatl, namlobopt )

      IF(lwp) THEN                         
         WRITE(numout,*)
         WRITE(numout,*) ' Namelist namlobopt'
         WRITE(numout,*) '    green   water absorption coeff                       xkg0  = ', xkg0
         WRITE(numout,*) '    red water absorption coeff                           xkr0  = ', xkr0
         WRITE(numout,*) '    pigment red absorption coeff                         xkrp  = ', xkrp
         WRITE(numout,*) '    pigment green absorption coeff                       xkgp  = ', xkgp
         WRITE(numout,*) '    green chl exposant                                   xlg   = ', xlg
         WRITE(numout,*) '    red   chl exposant                                   xlr   = ', xlr
         WRITE(numout,*) '    chla/chla+phea ratio                                 rpig  = ', rpig
         WRITE(numout,*) ' '
      ENDIF

#if defined key_diatrc && ! defined key_iomput

      ! Namelist namlobdia
      ! -------------------
      nn_writedia = 10                   ! default values

      DO jl = 1, jp_lobster_2d
         jn = jp_lob0_2d + jl - 1
         WRITE(ctrc2d(jn),'("2D_",I1)') jn                      ! short name
         WRITE(ctrc2l(jn),'("2D DIAGNOSTIC NUMBER ",I2)') jn    ! long name
         ctrc2u(jn) = ' '                                       ! units
      END DO
      !                                 ! 3D output arrays
      DO jl = 1, jp_lobster_3d
         jn = jp_lob0_3d + jl - 1
         WRITE(ctrc3d(jn),'("3D_",I1)') jn                      ! short name
         WRITE(ctrc3l(jn),'("3D DIAGNOSTIC NUMBER ",I2)') jn    ! long name
         ctrc3u(jn) = ' '                                       ! units
      END DO

      REWIND( numnatl )               ! read natrtd
      READ  ( numnatl, namlobdia )

      DO jl = 1, jp_lobster_2d
         jn = jp_lob0_2d + jl - 1
         ctrc2d(jn) = lobdia2d(jl)%snamedia
         ctrc2l(jn) = lobdia2d(jl)%lnamedia
         ctrc2u(jn) = lobdia2d(jl)%unitdia
      END DO

      DO jl = 1, jp_lobster_3d
         jn = jp_lob0_3d + jl - 1
         ctrc3d(jn) = lobdia3d(jl)%snamedia
         ctrc3l(jn) = lobdia3d(jl)%lnamedia
         ctrc3u(jn) = lobdia3d(jl)%unitdia
      END DO

      IF(lwp) THEN                   ! control print
         WRITE(numout,*)
         WRITE(numout,*) ' Namelist : natadd'
         WRITE(numout,*) '    frequency of outputs for additional arrays nn_writedia = ', nn_writedia
         DO jl = 1, jp_lobster_3d
            jn = jp_lob0_3d + jl - 1
            WRITE(numout,*) '   3d output field No : ',jn
            WRITE(numout,*) '   short name         : ', TRIM(ctrc3d(jn))
            WRITE(numout,*) '   long name          : ', TRIM(ctrc3l(jn))
            WRITE(numout,*) '   unit               : ', TRIM(ctrc3u(jn))
            WRITE(numout,*) ' '
         END DO

         DO jl = 1, jp_lobster_2d
            jn = jp_lob0_2d + jl - 1
            WRITE(numout,*) '   2d output field No : ',jn
            WRITE(numout,*) '   short name         : ', TRIM(ctrc2d(jn))
            WRITE(numout,*) '   long name          : ', TRIM(ctrc2l(jn))
            WRITE(numout,*) '   unit               : ', TRIM(ctrc2u(jn))
            WRITE(numout,*) ' '
         END DO
      ENDIF
#endif

#if defined key_diabio || defined key_trdmld_trc
      ! namlobdbi : bio diagnostics
      nwritebio = 10                     ! default values

      DO js = 1, jp_lobster_trd
         jd = jp_lob0_trd + js - 1
         IF(     jd <  10 ) THEN   ;   WRITE (ctrbio(jd),'("BIO_",I1)') jd      ! short name
         ELSEIF (jd < 100 ) THEN   ;   WRITE (ctrbio(jd),'("BIO_",I2)') jd   
         ELSE                      ;   WRITE (ctrbio(jd),'("BIO_",I3)') jd
         ENDIF
         WRITE(ctrbil(jd),'("BIOLOGICAL TREND NUMBER ",I2)') jd                 ! long name
         ctrbiu(jd) = 'mmoleN/m3/s '                                            ! units
      END DO

      REWIND( numnatl )
      READ  ( numnatl, namlobdbi ) 
 
      DO js = 1, jp_lobster_trd
         jd = jp_lob0_trd + js - 1
         ctrbio(jd) = lobdiabio(js)%snamebio
         ctrbil(jd) = lobdiabio(js)%lnamebio
         ctrbiu(jd) = lobdiabio(js)%unitbio
      END DO

      IF(lwp) THEN                   ! control print
         WRITE(numout,*)
         WRITE(numout,*) ' Namelist : namlobdbi'
         WRITE(numout,*) '    frequency of outputs for biological trends nwritebio = ', nwritebio
         DO js = 1, jp_lobster_trd
            jd = jp_lob0_trd + js - 1
            WRITE(numout,*) '   biological trend No : ',jd
            WRITE(numout,*) '   short name         : ', TRIM(ctrbio(jd))
            WRITE(numout,*) '   long name          : ', TRIM(ctrbil(jd))
            WRITE(numout,*) '   unit               : ', TRIM(ctrbiu(jd))
            WRITE(numout,*) ' '
         END DO
      END IF
#endif
      !
   END SUBROUTINE trc_nam_lobster
   
#else
   !!----------------------------------------------------------------------
   !!  Dummy module :                                            No LOBSTER
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_nam_lobster                      ! Empty routine
   END  SUBROUTINE  trc_nam_lobster
#endif  

   !!======================================================================
END MODULE trcnam_lobster
