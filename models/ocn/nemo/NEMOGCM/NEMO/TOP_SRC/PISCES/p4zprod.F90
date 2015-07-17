MODULE p4zprod
   !!======================================================================
   !!                         ***  MODULE p4zprod  ***
   !! TOP :   PISCES 
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   'key_pisces'                                       PISCES bio-model
   !!----------------------------------------------------------------------
   !!   p4z_prod       :  
   !!----------------------------------------------------------------------
   USE trc
   USE oce_trc         !
   USE sms_pisces      ! 
   USE prtctl_trc
   USE p4zopt
   USE p4zint
   USE p4zlim
   USE iom

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_prod         ! called in p4zbio.F90
   PUBLIC   p4z_prod_init    ! called in trcsms_pisces.F90
   PUBLIC   p4z_prod_alloc

   REAL(wp), PUBLIC ::   &
     pislope   = 3.0_wp          ,  &  !:
     pislope2  = 3.0_wp          ,  &  !:
     excret    = 10.e-5_wp       , &   !:
     excret2   = 0.05_wp         , &   !:
     chlcnm    = 0.033_wp        , &   !:
     chlcdm    = 0.05_wp         , &   !:
     fecnm     = 10.E-6_wp       , &   !:
     fecdm     = 15.E-6_wp       , &   !:
     grosip    = 0.151_wp

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   prmax   !:
   
   REAL(wp) ::   &
      rday1                      ,  &  !: 0.6 / rday
      texcret                    ,  &  !: 1 - excret 
      texcret2                   ,  &  !: 1 - excret2        
      tpp                              !: Total primary production

   !!* Substitution
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: p4zprod.F90 2730 2011-04-08 11:14:00Z rblod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_prod( kt , jnt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_prod  ***
      !!
      !! ** Purpose :   Compute the phytoplankton production depending on
      !!              light, temperature and nutrient availability
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      USE wrk_nemo, ONLY:   wrk_in_use, wrk_not_released
      USE wrk_nemo, ONLY:   zmixnano   => wrk_2d_1  , zmixdiat    => wrk_2d_2  , zstrn  => wrk_2d_3
      USE wrk_nemo, ONLY:   zpislopead => wrk_3d_2  , zpislopead2 => wrk_3d_3
      USE wrk_nemo, ONLY:   zprdia     => wrk_3d_4  , zprbio      => wrk_3d_5  , zysopt => wrk_3d_6
      USE wrk_nemo, ONLY:   zprorca    => wrk_3d_7  , zprorcad    => wrk_3d_8
      USE wrk_nemo, ONLY:   zprofed    => wrk_3d_9  , zprofen     => wrk_3d_10
      USE wrk_nemo, ONLY:   zprochln   => wrk_3d_11 , zprochld    => wrk_3d_12
      USE wrk_nemo, ONLY:   zpronew    => wrk_3d_13 , zpronewd    => wrk_3d_14
      !
      INTEGER, INTENT(in) :: kt, jnt
      !
      INTEGER  ::   ji, jj, jk
      REAL(wp) ::   zsilfac, zfact
      REAL(wp) ::   zprdiachl, zprbiochl, zsilim, ztn, zadap, zadap2
      REAL(wp) ::   zlim, zsilfac2, zsiborn, zprod, zetot2, zmax, zproreg, zproreg2
      REAL(wp) ::   zmxltst, zmxlday, zlim1
      REAL(wp) ::   zpislopen  , zpislope2n
      REAL(wp) ::   zrum, zcodel, zargu, zval, zvol
#if defined key_diatrc
      REAL(wp) ::   zrfact2
#endif
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------

      IF( wrk_in_use(2, 1,2,3)                             .OR.  &
          wrk_in_use(3, 2,3,4,5,6,7,8,9,10,11,12,13,14)  ) THEN
          CALL ctl_stop('p4z_prod: requested workspace arrays unavailable')   ;   RETURN
      ENDIF

      zprorca (:,:,:) = 0._wp
      zprorcad(:,:,:) = 0._wp
      zprofed (:,:,:) = 0._wp
      zprofen (:,:,:) = 0._wp
      zprochln(:,:,:) = 0._wp
      zprochld(:,:,:) = 0._wp
      zpronew (:,:,:) = 0._wp
      zpronewd(:,:,:) = 0._wp
      zprdia  (:,:,:) = 0._wp
      zprbio  (:,:,:) = 0._wp
      zysopt  (:,:,:) = 0._wp

      ! Computation of the optimal production
# if defined key_degrad
      prmax(:,:,:) = rday1 * tgfunc(:,:,:) * facvol(:,:,:)
# else
      prmax(:,:,:) = rday1 * tgfunc(:,:,:)
# endif

      ! compute the day length depending on latitude and the day
      zrum = REAL( nday_year - 80, wp ) / REAL( nyear_len(1), wp )
      zcodel = ASIN(  SIN( zrum * rpi * 2._wp ) * SIN( rad * 23.5_wp )  )

      ! day length in hours
      zstrn(:,:) = 0._wp
      DO jj = 1, jpj
         DO ji = 1, jpi
            zargu = TAN( zcodel ) * TAN( gphit(ji,jj) * rad )
            zargu = MAX( -1., MIN(  1., zargu ) )
            zval  = MAX( 0.0, 24. - 2. * ACOS( zargu ) / rad / 15. )
            IF( zval < 1.e0 )   zval = 24.
            zstrn(ji,jj) = 24. / zval
         END DO
      END DO


!CDIR NOVERRCHK
      DO jk = 1, jpkm1
!CDIR NOVERRCHK
         DO jj = 1, jpj
!CDIR NOVERRCHK
            DO ji = 1, jpi

               ! Computation of the P-I slope for nanos and diatoms
               IF( etot(ji,jj,jk) > 1.E-3 ) THEN
                   ztn    = MAX( 0., tsn(ji,jj,jk,jp_tem) - 15. )
                   zadap  = 0.+ 1.* ztn / ( 2.+ ztn )
                   zadap2 = 0.e0

                   zfact  = EXP( -0.21 * emoy(ji,jj,jk) )

                   zpislopead (ji,jj,jk) = pislope  * ( 1.+ zadap  * zfact )
                   zpislopead2(ji,jj,jk) = pislope2 * ( 1.+ zadap2 * zfact )

                   zpislopen = zpislopead(ji,jj,jk) * trn(ji,jj,jk,jpnch)                 &
                     &         / ( trn(ji,jj,jk,jpphy) * 12.                   + rtrn )   &
                     &         / ( prmax(ji,jj,jk) * rday * xlimphy(ji,jj,jk) + rtrn )

                   zpislope2n = zpislopead2(ji,jj,jk) * trn(ji,jj,jk,jpdch)                &
                     &          / ( trn(ji,jj,jk,jpdia) * 12.                   + rtrn )   &
                     &          / ( prmax(ji,jj,jk) * rday * xlimdia(ji,jj,jk) + rtrn )

                   ! Computation of production function
                   zprbio(ji,jj,jk) = prmax(ji,jj,jk) * &
                     &                (  1.- EXP( -zpislopen * enano(ji,jj,jk) )  )
                   zprdia(ji,jj,jk) = prmax(ji,jj,jk) * &
                     &                (  1.- EXP( -zpislope2n * ediat(ji,jj,jk) )  )
               ENDIF
            END DO
         END DO
      END DO


      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi

                IF( etot(ji,jj,jk) > 1.E-3 ) THEN
                   !    Si/C of diatoms
                   !    ------------------------
                   !    Si/C increases with iron stress and silicate availability
                   !    Si/C is arbitrariliy increased for very high Si concentrations
                   !    to mimic the very high ratios observed in the Southern Ocean (silpot2)

                  zlim1  = trn(ji,jj,jk,jpsil) / ( trn(ji,jj,jk,jpsil) + xksi1 )
                  zlim   = xdiatno3(ji,jj,jk) + xdiatnh4(ji,jj,jk)

                  zsilim = MIN( zprdia(ji,jj,jk)    / ( rtrn + prmax(ji,jj,jk) ),                 &
                  &          trn(ji,jj,jk,jpfer) / ( concdfe(ji,jj,jk) + trn(ji,jj,jk,jpfer) ),   &
                  &          trn(ji,jj,jk,jppo4) / ( concdnh4 + trn(ji,jj,jk,jppo4) ),            &
                  &          zlim )
                  zsilfac = 5.4 * EXP( -4.23 * zsilim ) * MAX( 0.e0, MIN( 1., 2.2 * ( zlim1 - 0.5 ) )  ) + 1.e0
                  zsiborn = MAX( 0.e0, ( trn(ji,jj,jk,jpsil) - 15.e-6 ) )
                  zsilfac2 = 1.+ 3.* zsiborn / ( zsiborn + xksi2 )
                  zsilfac = MIN( 6.4,zsilfac * zsilfac2)
                  zysopt(ji,jj,jk) = grosip * zlim1 * zsilfac
              ENDIF
            END DO
         END DO
      END DO

      !  Computation of the limitation term due to
      !  A mixed layer deeper than the euphotic depth
      DO jj = 1, jpj
         DO ji = 1, jpi
            zmxltst = MAX( 0.e0, hmld(ji,jj) - heup(ji,jj) )
            zmxlday = zmxltst**2 / rday
            zmixnano(ji,jj) = 1.- zmxlday / ( 1.+ zmxlday )
            zmixdiat(ji,jj) = 1.- zmxlday / ( 3.+ zmxlday )
         END DO
      END DO
 
      !  Mixed-layer effect on production                                                                               
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( fsdepw(ji,jj,jk+1) <= hmld(ji,jj) ) THEN
                  zprbio(ji,jj,jk) = zprbio(ji,jj,jk) * zmixnano(ji,jj)
                  zprdia(ji,jj,jk) = zprdia(ji,jj,jk) * zmixdiat(ji,jj)
               ENDIF
            END DO
         END DO
      END DO


!CDIR NOVERRCHK
      DO jk = 1, jpkm1
!CDIR NOVERRCHK
         DO jj = 1, jpj
!CDIR NOVERRCHK
            DO ji = 1, jpi

               IF( etot(ji,jj,jk) > 1.E-3 ) THEN
                  !     Computation of the various production terms for nanophyto.
                  zetot2 = enano(ji,jj,jk) * zstrn(ji,jj)
                  zmax = MAX( 0.1, xlimphy(ji,jj,jk) )
                  zpislopen = zpislopead(ji,jj,jk)          &
                  &         * trn(ji,jj,jk,jpnch) / ( rtrn + trn(ji,jj,jk,jpphy) * 12.)         &
                  &         / ( prmax(ji,jj,jk) * rday * zmax + rtrn )

                  zprbiochl = prmax(ji,jj,jk) * (  1.- EXP( -zpislopen * zetot2 )  )

                  zprorca(ji,jj,jk) = zprbio(ji,jj,jk)  * xlimphy(ji,jj,jk) * trn(ji,jj,jk,jpphy) * rfact2

                  zpronew(ji,jj,jk) = zprorca(ji,jj,jk) * xnanono3(ji,jj,jk)    &
                  &             / ( xnanono3(ji,jj,jk) + xnanonh4(ji,jj,jk) + rtrn )
                  zprod = rday * zprorca(ji,jj,jk) * zprbiochl * trn(ji,jj,jk,jpphy) *zmax

                  zprofen(ji,jj,jk) = (fecnm)**2 * zprod / chlcnm            &
                  &              / (  zpislopead(ji,jj,jk) * zetot2 * trn(ji,jj,jk,jpnfe) + rtrn  )

                  zprochln(ji,jj,jk) = chlcnm * 144. * zprod                  &
                  &              / (  zpislopead(ji,jj,jk) * zetot2 * trn(ji,jj,jk,jpnch) + rtrn  )
               ENDIF
            END DO
         END DO
      END DO

!CDIR NOVERRCHK
      DO jk = 1, jpkm1
!CDIR NOVERRCHK
         DO jj = 1, jpj
!CDIR NOVERRCHK
            DO ji = 1, jpi
               IF( etot(ji,jj,jk) > 1.E-3 ) THEN
                  !  Computation of the various production terms for diatoms
                  zetot2 = ediat(ji,jj,jk) * zstrn(ji,jj)
                  zmax = MAX( 0.1, xlimdia(ji,jj,jk) )
                  zpislope2n = zpislopead2(ji,jj,jk) * trn(ji,jj,jk,jpdch)        &
                  &           / ( rtrn + trn(ji,jj,jk,jpdia) * 12.)        &
                  &           / ( prmax(ji,jj,jk) * rday * zmax + rtrn )

                  zprdiachl = prmax(ji,jj,jk) * (  1.- EXP( -zetot2 * zpislope2n )  )

                  zprorcad(ji,jj,jk) = zprdia(ji,jj,jk) * xlimdia(ji,jj,jk) * trn(ji,jj,jk,jpdia) * rfact2

                  zpronewd(ji,jj,jk) = zprorcad(ji,jj,jk) * xdiatno3(ji,jj,jk)     &
                  &              / ( xdiatno3(ji,jj,jk) + xdiatnh4(ji,jj,jk) + rtrn )

                  zprod = rday * zprorcad(ji,jj,jk) * zprdiachl * trn(ji,jj,jk,jpdia) * zmax

                  zprofed(ji,jj,jk) = (fecdm)**2 * zprod / chlcdm                   &
                  &              / ( zpislopead2(ji,jj,jk) * zetot2 * trn(ji,jj,jk,jpdfe) + rtrn )

                  zprochld(ji,jj,jk) = chlcdm * 144. * zprod       &
                  &              / ( zpislopead2(ji,jj,jk) * zetot2 * trn(ji,jj,jk,jpdch) + rtrn )

               ENDIF
            END DO
         END DO
      END DO
      !

      !   Update the arrays TRA which contain the biological sources and sinks
      DO jk = 1, jpkm1
         DO jj = 1, jpj
           DO ji =1 ,jpi
              zproreg  = zprorca(ji,jj,jk) - zpronew(ji,jj,jk)
              zproreg2 = zprorcad(ji,jj,jk) - zpronewd(ji,jj,jk)
              tra(ji,jj,jk,jppo4) = tra(ji,jj,jk,jppo4) - zprorca(ji,jj,jk) - zprorcad(ji,jj,jk)
              tra(ji,jj,jk,jpno3) = tra(ji,jj,jk,jpno3) - zpronew(ji,jj,jk) - zpronewd(ji,jj,jk)
              tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) - zproreg - zproreg2
              tra(ji,jj,jk,jpphy) = tra(ji,jj,jk,jpphy) + zprorca(ji,jj,jk) * texcret
              tra(ji,jj,jk,jpnch) = tra(ji,jj,jk,jpnch) + zprochln(ji,jj,jk) * texcret
              tra(ji,jj,jk,jpnfe) = tra(ji,jj,jk,jpnfe) + zprofen(ji,jj,jk) * texcret
              tra(ji,jj,jk,jpdia) = tra(ji,jj,jk,jpdia) + zprorcad(ji,jj,jk) * texcret2
              tra(ji,jj,jk,jpdch) = tra(ji,jj,jk,jpdch) + zprochld(ji,jj,jk) * texcret2
              tra(ji,jj,jk,jpdfe) = tra(ji,jj,jk,jpdfe) + zprofed(ji,jj,jk) * texcret2
              tra(ji,jj,jk,jpbsi) = tra(ji,jj,jk,jpbsi) + zprorcad(ji,jj,jk) * zysopt(ji,jj,jk) * texcret2
              tra(ji,jj,jk,jpdoc) = tra(ji,jj,jk,jpdoc) + &
              &                     excret2 * zprorcad(ji,jj,jk) + excret * zprorca(ji,jj,jk)
              tra(ji,jj,jk,jpoxy) = tra(ji,jj,jk,jpoxy) + o2ut * ( zproreg + zproreg2) &
              &                    + ( o2ut + o2nit ) * ( zpronew(ji,jj,jk) + zpronewd(ji,jj,jk) )
              tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) &
              &                     - texcret * zprofen(ji,jj,jk) - texcret2 * zprofed(ji,jj,jk)
              tra(ji,jj,jk,jpsil) = tra(ji,jj,jk,jpsil) &
              &                     - texcret2 * zprorcad(ji,jj,jk) * zysopt(ji,jj,jk)
              tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) - zprorca(ji,jj,jk) - zprorcad(ji,jj,jk)
              tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) &
              &                    + rno3 * ( zpronew(ji,jj,jk) + zpronewd(ji,jj,jk) )
          END DO
        END DO
     END DO

     ! Total primary production per year

#if defined key_degrad
     tpp = tpp + glob_sum( ( zprorca(:,:,:) + zprorcad(:,:,:) ) * cvol(:,:,:) * facvol(:,:,:) )
#else
     tpp = tpp + glob_sum( ( zprorca(:,:,:) + zprorcad(:,:,:) ) * cvol(:,:,:) )
#endif

     IF( kt == nitend .AND. jnt == nrdttrc .AND. lwp ) THEN
        WRITE(numout,*) 'Total PP (Gtc) :'
        WRITE(numout,*) '-------------------- : ',tpp * 12. / 1.E12
        WRITE(numout,*) 
      ENDIF

#if defined key_diatrc && ! defined key_iomput
      !   Supplementary diagnostics
      zrfact2 = 1.e3 * rfact2r
      trc3d(:,:,:,jp_pcs0_3d + 4)  = zprorca (:,:,:) * zrfact2 * tmask(:,:,:)
      trc3d(:,:,:,jp_pcs0_3d + 5)  = zprorcad(:,:,:) * zrfact2 * tmask(:,:,:)
      trc3d(:,:,:,jp_pcs0_3d + 6)  = zpronew (:,:,:) * zrfact2 * tmask(:,:,:)
      trc3d(:,:,:,jp_pcs0_3d + 7)  = zpronewd(:,:,:) * zrfact2 * tmask(:,:,:)
      trc3d(:,:,:,jp_pcs0_3d + 8)  = zprorcad(:,:,:) * zrfact2 * tmask(:,:,:) * zysopt(:,:,:)
      trc3d(:,:,:,jp_pcs0_3d + 9)  = zprofed (:,:,:) * zrfact2 * tmask(:,:,:)
#  if ! defined key_kriest
      trc3d(:,:,:,jp_pcs0_3d + 10) = zprofen (:,:,:) * zrfact2 * tmask(:,:,:)
#  endif
#endif

#if defined key_diatrc && defined key_iomput
      zrfact2 = 1.e3 * rfact2r
      IF ( jnt == nrdttrc ) then
         CALL iom_put( "PPPHY" , zprorca (:,:,:) * zrfact2 * tmask(:,:,:) )  ! primary production by nanophyto
         CALL iom_put( "PPPHY2", zprorcad(:,:,:) * zrfact2 * tmask(:,:,:) )  ! primary production by diatom
         CALL iom_put( "PPNEWN", zpronew (:,:,:) * zrfact2 * tmask(:,:,:) )  ! new primary production by nanophyto
         CALL iom_put( "PPNEWD", zpronewd(:,:,:) * zrfact2 * tmask(:,:,:) )  ! new primary production by diatom
         CALL iom_put( "PBSi"  , zprorcad(:,:,:) * zrfact2 * tmask(:,:,:) * zysopt(:,:,:) ) ! biogenic silica production
         CALL iom_put( "PFeD"  , zprofed (:,:,:) * zrfact2 * tmask(:,:,:) )  ! biogenic iron production by diatom
         CALL iom_put( "PFeN"  , zprofen (:,:,:) * zrfact2 * tmask(:,:,:) )  ! biogenic iron production by nanophyto
      ENDIF
#endif

      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('prod')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF

      IF(  wrk_not_released(2, 1,2,3)                          .OR.  &
           wrk_not_released(3, 2,3,4,5,6,7,8,9,10,11,12,13,14)   )   &
           CALL ctl_stop('p4z_prod: failed to release workspace arrays')
      !
   END SUBROUTINE p4z_prod


   SUBROUTINE p4z_prod_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_prod_init  ***
      !!
      !! ** Purpose :   Initialization of phytoplankton production parameters
      !!
      !! ** Method  :   Read the nampisprod namelist and check the parameters
      !!      called at the first timestep (nit000)
      !!
      !! ** input   :   Namelist nampisprod
      !!----------------------------------------------------------------------
      NAMELIST/nampisprod/ pislope, pislope2, excret, excret2, chlcnm, chlcdm,   &
         &              fecnm, fecdm, grosip
      !!----------------------------------------------------------------------

      REWIND( numnat )                     ! read numnat
      READ  ( numnat, nampisprod )

      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for phytoplankton growth, nampisprod'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    mean Si/C ratio                           grosip    =', grosip
         WRITE(numout,*) '    P-I slope                                 pislope   =', pislope
         WRITE(numout,*) '    excretion ratio of nanophytoplankton      excret    =', excret
         WRITE(numout,*) '    excretion ratio of diatoms                excret2   =', excret2
         WRITE(numout,*) '    P-I slope  for diatoms                    pislope2  =', pislope2
         WRITE(numout,*) '    Minimum Chl/C in nanophytoplankton        chlcnm    =', chlcnm
         WRITE(numout,*) '    Minimum Chl/C in diatoms                  chlcdm    =', chlcdm
         WRITE(numout,*) '    Maximum Fe/C in nanophytoplankton         fecnm     =', fecnm
         WRITE(numout,*) '    Minimum Fe/C in diatoms                   fecdm     =', fecdm
      ENDIF
      !
      rday1     = 0.6 / rday 
      texcret   = 1.0 - excret
      texcret2  = 1.0 - excret2
      tpp       = 0.
      !
   END SUBROUTINE p4z_prod_init


   INTEGER FUNCTION p4z_prod_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_prod_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( prmax(jpi,jpj,jpk), STAT=p4z_prod_alloc )
      !
      IF( p4z_prod_alloc /= 0 ) CALL ctl_warn('p4z_prod_alloc : failed to allocate arrays.')
      !
   END FUNCTION p4z_prod_alloc

#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_prod                    ! Empty routine
   END SUBROUTINE p4z_prod
#endif 

   !!======================================================================
END MODULE  p4zprod
