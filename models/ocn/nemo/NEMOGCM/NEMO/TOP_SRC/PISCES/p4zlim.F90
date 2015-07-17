MODULE p4zlim
   !!======================================================================
   !!                         ***  MODULE p4zlim  ***
   !! TOP :   PISCES 
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   'key_pisces'                                       PISCES bio-model
   !!----------------------------------------------------------------------
   !!   p4z_lim        :   Compute the nutrients limitation terms 
   !!   p4z_lim_init   :   Read the namelist 
   !!----------------------------------------------------------------------
   USE trc
   USE oce_trc         !
   USE trc         ! 
   USE sms_pisces      ! 

   IMPLICIT NONE
   PRIVATE

   PUBLIC p4z_lim    
   PUBLIC p4z_lim_init    

   !! * Shared module variables
   REAL(wp), PUBLIC ::   &
     conc0     = 2.e-6_wp      ,  &  !:
     conc1     = 10.e-6_wp     ,  &  !:
     conc2     = 2.e-11_wp     ,  &  !:
     conc2m    = 8.E-11_wp     ,  &  !:
     conc3     = 1.e-10_wp     ,  &  !:
     conc3m    = 4.e-10_wp     ,  &  !:
     concnnh4  = 1.e-7_wp      ,  &  !:
     concdnh4  = 5.e-7_wp      ,  &  !:
     xksi1     = 2.E-6_wp      ,  &  !:
     xksi2     = 3.33E-6_wp    ,  &  !:
     xkdoc     = 417.E-6_wp    ,  &  !:
     caco3r    = 0.3_wp              !:


   !!* Substitution
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: p4zlim.F90 2528 2010-12-27 17:33:53Z rblod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE p4z_lim( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_lim  ***
      !!
      !! ** Purpose :   Compute the co-limitations by the various nutrients
      !!              for the various phytoplankton species
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in)  :: kt
      INTEGER  ::   ji, jj, jk
      REAL(wp) ::   zlim1, zlim2, zlim3, zlim4, zno3, zferlim
      REAL(wp) ::   zconctemp, zconctemp2, zconctempn, zconctempn2
      REAL(wp) ::   ztemp, zdenom
      !!---------------------------------------------------------------------


      !  Tuning of the iron concentration to a minimum
      !  level that is set to the detection limit
      !  -------------------------------------

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               zno3=trn(ji,jj,jk,jpno3)
               zferlim = MAX( 1.5e-11*(zno3/40E-6)**2, 3e-12 )
               zferlim = MIN( zferlim, 1.5e-11 )
               trn(ji,jj,jk,jpfer) = MAX( trn(ji,jj,jk,jpfer), zferlim )
            END DO
         END DO
      END DO

      !  Computation of a variable Ks for iron on diatoms taking into account
      !  that increasing biomass is made of generally bigger cells
      !  ------------------------------------------------

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               zconctemp   = MAX( 0.e0 , trn(ji,jj,jk,jpdia)-5e-7 )
               zconctemp2  = trn(ji,jj,jk,jpdia) - zconctemp
               zconctempn  = MAX( 0.e0 , trn(ji,jj,jk,jpphy)-1e-6 )
               zconctempn2 = trn(ji,jj,jk,jpphy) - zconctempn
               concdfe(ji,jj,jk) = ( zconctemp2 * conc3 + conc3m * zconctemp)   &
                  &              / ( trn(ji,jj,jk,jpdia) + rtrn )
               concdfe(ji,jj,jk) = MAX( conc3, concdfe(ji,jj,jk) )
               concnfe(ji,jj,jk) = ( zconctempn2 * conc2 + conc2m * zconctempn)   &
                  &              / ( trn(ji,jj,jk,jpphy) + rtrn )
               concnfe(ji,jj,jk) = MAX( conc2, concnfe(ji,jj,jk) )
            END DO
         END DO
      END DO

     !  Michaelis-Menten Limitation term for nutrients Small flagellates
     !      -----------------------------------------------
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
              zdenom = 1. / &
                  & ( conc0 * concnnh4 + concnnh4 * trn(ji,jj,jk,jpno3) + conc0 * trn(ji,jj,jk,jpnh4) )
               xnanono3(ji,jj,jk) = trn(ji,jj,jk,jpno3) * concnnh4 * zdenom
               xnanonh4(ji,jj,jk) = trn(ji,jj,jk,jpnh4) * conc0    * zdenom

               zlim1 = xnanono3(ji,jj,jk) + xnanonh4(ji,jj,jk)
               zlim2 = trn(ji,jj,jk,jppo4) / ( trn(ji,jj,jk,jppo4) + concnnh4          ) 
               zlim3 = trn(ji,jj,jk,jpfer) / ( trn(ji,jj,jk,jpfer) + concnfe(ji,jj,jk) )
               xlimphy(ji,jj,jk) = MIN( zlim1, zlim2, zlim3 )
               zlim1 = trn(ji,jj,jk,jpnh4) / ( concnnh4 + trn(ji,jj,jk,jpnh4) )
               zlim3 = trn(ji,jj,jk,jpfer) / ( conc2    + trn(ji,jj,jk,jpfer) )
               zlim4 = trn(ji,jj,jk,jpdoc) / ( xkdoc   + trn(ji,jj,jk,jpdoc) )
               xlimbac(ji,jj,jk) = MIN( zlim1, zlim2, zlim3 ) * zlim4

            END DO
         END DO
      END DO

      !   Michaelis-Menten Limitation term for nutrients Diatoms
      !   ----------------------------------------------
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
              zdenom = 1. / &
                  & ( conc1 * concdnh4 + concdnh4 * trn(ji,jj,jk,jpno3) + conc1 * trn(ji,jj,jk,jpnh4) )

               xdiatno3(ji,jj,jk) = trn(ji,jj,jk,jpno3) * concdnh4 * zdenom
               xdiatnh4(ji,jj,jk) = trn(ji,jj,jk,jpnh4) * conc1    * zdenom 

               zlim1 = xdiatno3(ji,jj,jk) + xdiatnh4(ji,jj,jk)
               zlim2 = trn(ji,jj,jk,jppo4) / ( trn(ji,jj,jk,jppo4) + concdnh4          )
               zlim3 = trn(ji,jj,jk,jpsil) / ( trn(ji,jj,jk,jpsil) + xksi   (ji,jj)    )
               zlim4 = trn(ji,jj,jk,jpfer) / ( trn(ji,jj,jk,jpfer) + concdfe(ji,jj,jk) )
               xlimdia(ji,jj,jk) = MIN( zlim1, zlim2, zlim3, zlim4 )

            END DO
         END DO
      END DO


      ! Compute the fraction of nanophytoplankton that is made of calcifiers
      ! --------------------------------------------------------------------

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               ztemp = MAX( 0., tsn(ji,jj,jk,jp_tem) )
               xfracal(ji,jj,jk) = caco3r * xlimphy(ji,jj,jk)   &
                  &                       * MAX( 0.0001, ztemp / ( 2.+ ztemp ) )   &
                  &                       * MAX( 1., trn(ji,jj,jk,jpphy) * 1.e6 / 2. )
               xfracal(ji,jj,jk) = MIN( 0.8 , xfracal(ji,jj,jk) )
               xfracal(ji,jj,jk) = MAX( 0.01, xfracal(ji,jj,jk) )
            END DO
         END DO
      END DO
      !
   END SUBROUTINE p4z_lim

   SUBROUTINE p4z_lim_init

      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_lim_init  ***
      !!
      !! ** Purpose :   Initialization of nutrient limitation parameters
      !!
      !! ** Method  :   Read the nampislim namelist and check the parameters
      !!      called at the first timestep (nit000)
      !!
      !! ** input   :   Namelist nampislim
      !!
      !!----------------------------------------------------------------------

      NAMELIST/nampislim/ conc0, conc1, conc2, conc2m, conc3, conc3m,   &
         &             concnnh4, concdnh4, xksi1, xksi2, xkdoc, caco3r

      REWIND( numnat )                     ! read numnat
      READ  ( numnat, nampislim )

      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for nutrient limitations, nampislim'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    mean rainratio                            caco3r    =', caco3r
         WRITE(numout,*) '    NO3, PO4 half saturation                  conc0      =', conc0
         WRITE(numout,*) '    half saturation constant for Si uptake    xksi1     =', xksi1
         WRITE(numout,*) '    half saturation constant for Si/C         xksi2     =', xksi2
         WRITE(numout,*) '    2nd half-sat. of DOC remineralization     xkdoc    =', xkdoc
         WRITE(numout,*) '    Phosphate half saturation for diatoms     conc1     =', conc1
         WRITE(numout,*) '    Iron half saturation for phyto            conc2     =', conc2
         WRITE(numout,*) '    Max iron half saturation for phyto        conc2m    =', conc2m
         WRITE(numout,*) '    Iron half saturation for diatoms          conc3     =', conc3
         WRITE(numout,*) '    Maxi iron half saturation for diatoms     conc3m    =', conc3m
         WRITE(numout,*) '    NH4 half saturation for phyto             concnnh4  =', concnnh4
         WRITE(numout,*) '    NH4 half saturation for diatoms           concdnh4  =', concdnh4
      ENDIF

   END SUBROUTINE p4z_lim_init

#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_lim                   ! Empty routine
   END SUBROUTINE p4z_lim
#endif 

   !!======================================================================
END MODULE  p4zlim
