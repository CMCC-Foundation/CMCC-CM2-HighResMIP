MODULE flodom
   !!======================================================================
   !!                       ***  MODULE  flodom  ***
   !! Ocean floats :   domain
   !!======================================================================
   !! History :  OPA  ! 1998-07 (Y.Drillet, CLIPPER)  Original code
   !!----------------------------------------------------------------------
#if   defined key_floats   ||   defined key_esopa
   !!----------------------------------------------------------------------
   !!   'key_floats'                                     float trajectories
   !!----------------------------------------------------------------------
   !!   flo_dom        : initialization of floats
   !!   findmesh       : compute index of position 
   !!   dstnce         : compute distance between face mesh and floats 
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE flo_oce         ! ocean drifting floats
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! distribued memory computing library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   flo_dom    ! routine called by floats.F90

   !! * Substitutions
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: flodom.F90 2528 2010-12-27 17:33:53Z rblod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE flo_dom
      !! ---------------------------------------------------------------------
      !!                  ***  ROUTINE flo_dom  ***
      !!                 
      !!  ** Purpose :   Initialisation of floats
      !!
      !!  ** Method  :   We put the floats  in the domain with the latitude,
      !!               the longitude (degree) and the depth (m).
      !!----------------------------------------------------------------------      
      LOGICAL  ::   llinmesh
      INTEGER  ::   ji, jj, jk   ! DO loop index on 3 directions
      INTEGER  ::   jfl, jfl1    ! number of floats   
      INTEGER  ::   inum         ! logical unit for file read
      INTEGER, DIMENSION(jpnfl) ::   iimfl, ijmfl, ikmfl       ! index mesh of floats
      INTEGER, DIMENSION(jpnfl) ::   idomfl,  ivtest, ihtest   !   -             -
      REAL(wp) ::   zdxab, zdyad
      REAL(wp), DIMENSION(jpnnewflo+1)  :: zgifl, zgjfl,  zgkfl
      !!---------------------------------------------------------------------
      
      ! Initialisation with the geographical position or restart
      
      IF(lwp) WRITE(numout,*) 'flo_dom : compute initial position of floats'
      IF(lwp) WRITE(numout,*) '~~~~~~~~'
      IF(lwp) WRITE(numout,*) '           jpnfl = ',jpnfl
      
      IF(ln_rstflo) THEN
         IF(lwp) WRITE(numout,*) '        float restart file read'
         
         ! open the restart file 
         CALL ctl_opn( inum, 'restart_float', 'OLD', 'FORMATTED', 'SEQUENTIAL', -1, numout, lwp )

         ! read of the restart file
         READ(inum) ( tpifl  (jfl), jfl=1, jpnrstflo),   & 
                        ( tpjfl  (jfl), jfl=1, jpnrstflo),   &
                        ( tpkfl  (jfl), jfl=1, jpnrstflo),   &
                        ( nisobfl(jfl), jfl=1, jpnrstflo),   &
                        ( ngrpfl (jfl), jfl=1, jpnrstflo)    
         CLOSE(inum)

         ! if we want a  surface drift  ( like PROVOR floats )
         IF( ln_argo ) THEN
            DO jfl = 1, jpnrstflo
               nisobfl(jfl) = 0
            END DO
         ENDIF

         IF(lwp) WRITE(numout,*)' flo_dom: END of florstlec'
         
         ! It is possible to add new floats.          
         IF(lwp) WRITE(numout,*)' flo_dom:jpnfl jpnrstflo ',jpnfl,jpnrstflo
         IF( jpnfl > jpnrstflo ) THEN
            ! open the init file 
            CALL ctl_opn( inum, 'init_float', 'OLD', 'FORMATTED', 'SEQUENTIAL', -1, numout, lwp )
            DO jfl = jpnrstflo+1, jpnfl
               READ(inum,*) flxx(jfl),flyy(jfl),flzz(jfl), nisobfl(jfl),ngrpfl(jfl),jfl1
            END DO
            CLOSE(inum)
            IF(lwp) WRITE(numout,*)' flodom: END reading init_float file'
            
            ! Test to find the grid point coordonate with the geographical position            
            DO jfl = jpnrstflo+1, jpnfl
               ihtest(jfl) = 0
               ivtest(jfl) = 0
               ikmfl(jfl) = 0
# if   defined key_mpp_mpi
               DO ji = MAX(nldi,2), nlei
                  DO jj = MAX(nldj,2), nlej   ! NO vector opt.
# else
               DO ji = 2, jpi
                  DO jj = 2, jpj   ! NO vector opt.
# endif                     
                     ! For each float we find the indexes of the mesh                      
                     CALL findmesh(glamf(ji-1,jj-1),gphif(ji-1,jj-1),   &
                                   glamf(ji-1,jj  ),gphif(ji-1,jj  ),   &
                                   glamf(ji  ,jj  ),gphif(ji  ,jj  ),   &
                                   glamf(ji  ,jj-1),gphif(ji  ,jj-1),   &
                                   flxx(jfl)       ,flyy(jfl)       ,   &
                                   glamt(ji  ,jj  ),gphit(ji  ,jj  ), llinmesh)
                     IF(llinmesh) THEN
                        iimfl(jfl) = ji
                        ijmfl(jfl) = jj
                        ihtest(jfl) = ihtest(jfl)+1
                        DO jk = 1, jpk-1
                           IF( (fsdepw(ji,jj,jk) <= flzz(jfl)) .AND. (fsdepw(ji,jj,jk+1) > flzz(jfl)) ) THEN
                              ikmfl(jfl) = jk
                              ivtest(jfl) = ivtest(jfl) + 1
                           ENDIF
                        END DO
                     ENDIF
                  END DO
               END DO
               IF(lwp) WRITE(numout,*)'   flo_dom: END findmesh'
               
               ! If the float is in a mesh computed by an other processor we put iimfl=ijmfl=-1               
               IF( ihtest(jfl) ==  0 ) THEN
                  iimfl(jfl) = -1
                  ijmfl(jfl) = -1
               ENDIF
            END DO
            
            ! A zero in the sum of the arrays "ihtest" and "ivtest"             
# if   defined key_mpp_mpi
            CALL mpp_sum(ihtest,jpnfl)
            CALL mpp_sum(ivtest,jpnfl)
# endif 
            DO jfl = jpnrstflo+1, jpnfl
               IF( (ihtest(jfl) > 1 ) .OR. ( ivtest(jfl) > 1) ) THEN
                  IF(lwp) WRITE(numout,*) 'THE FLOAT',jfl,' IS NOT IN ONLY ONE MESH'
                  STOP
               ENDIF
               IF( (ihtest(jfl) == 0) ) THEN
                  IF(lwp) WRITE(numout,*)'THE FLOAT',jfl,' IS IN NO MESH'
                  STOP
               ENDIF
            END DO
            
            ! We compute the distance between the float and the face of the mesh            
            DO jfl = jpnrstflo+1, jpnfl               
               ! Made only if the float is in the domain of the processor              
               IF( (iimfl(jfl) >= 0) .AND. (ijmfl(jfl) >= 0) ) THEN
                  
                  ! TEST TO KNOW IF THE FLOAT IS NOT INITIALISED IN THE COAST
                  
                  idomfl(jfl) = 0
                  IF( tmask(iimfl(jfl),ijmfl(jfl),ikmfl(jfl)) == 0. ) idomfl(jfl) = 1
                                           
                  ! Computation of the distance between the float and the faces of the mesh
                  !            zdxab
                  !             .
                  !        B----.---------C
                  !        |    .         |
                  !        |<------>flo   |
                  !        |        ^     |
                  !        |        |.....|....zdyad
                  !        |        |     |
                  !        A--------|-----D
                  !
             
                  zdxab = dstnce( flxx(jfl), flyy(jfl), glamf(iimfl(jfl)-1,ijmfl(jfl)-1), flyy(jfl) )
                  zdyad = dstnce( flxx(jfl), flyy(jfl), flxx(jfl), gphif(iimfl(jfl)-1,ijmfl(jfl)-1) )
                  
                  ! Translation of this distances (in meter) in indexes
                  
                  zgifl(jfl-jpnrstflo)= (iimfl(jfl)-0.5) + zdxab/e1u(iimfl(jfl)-1,ijmfl(jfl)) + (mig(1)-jpizoom)
                  zgjfl(jfl-jpnrstflo)= (ijmfl(jfl)-0.5) + zdyad/e2v(iimfl(jfl),ijmfl(jfl)-1) + (mjg(1)-jpjzoom)
                  zgkfl(jfl-jpnrstflo) = (( fsdepw(iimfl(jfl),ijmfl(jfl),ikmfl(jfl)+1) - flzz(jfl) )* ikmfl(jfl))   &
                     &                 / (  fsdepw(iimfl(jfl),ijmfl(jfl),ikmfl(jfl)+1)                              &
                     &                    - fsdepw(iimfl(jfl),ijmfl(jfl),ikmfl(jfl) ) )                             &
                     &                 + (( flzz(jfl)-fsdepw(iimfl(jfl),ijmfl(jfl),ikmfl(jfl)) ) *(ikmfl(jfl)+1))   &
                     &                 / (  fsdepw(iimfl(jfl),ijmfl(jfl),ikmfl(jfl)+1)                              &
                     &                    - fsdepw(iimfl(jfl),ijmfl(jfl),ikmfl(jfl)) )
               ELSE
                  zgifl(jfl-jpnrstflo) = 0.e0
                  zgjfl(jfl-jpnrstflo) = 0.e0
                  zgkfl(jfl-jpnrstflo) = 0.e0
               ENDIF
            END DO
            
            ! The sum of all the arrays zgifl, zgjfl, zgkfl give 3 arrays with the positions of all the floats.
            IF( lk_mpp )   THEN
               CALL mpp_sum( zgjfl, jpnnewflo )   ! sums over the global domain
               CALL mpp_sum( zgkfl, jpnnewflo )
               IF(lwp) WRITE(numout,*) (zgifl(jfl),jfl=1,jpnnewflo)
               IF(lwp) WRITE(numout,*) (zgjfl(jfl),jfl=1,jpnnewflo)
               IF(lwp) WRITE(numout,*) (zgkfl(jfl),jfl=1,jpnnewflo) 
            ENDIF
           
            DO jfl = jpnrstflo+1, jpnfl
               tpifl(jfl) = zgifl(jfl-jpnrstflo)
               tpjfl(jfl) = zgjfl(jfl-jpnrstflo)
               tpkfl(jfl) = zgkfl(jfl-jpnrstflo)
            END DO
         ENDIF
      ELSE
         IF(lwp) WRITE(numout,*) '                     init_float read '
         
         ! First initialisation of floats
         ! the initials positions of floats are written in a file
         ! with a variable to know if it is a isobar float a number 
         ! to identified who want the trajectories of this float and 
         ! an index for the number of the float         
         ! open the init file 
         CALL ctl_opn( inum, 'init_float', 'OLD', 'FORMATTED', 'SEQUENTIAL', -1, numout, lwp )
         READ(inum) (flxx(jfl)   , jfl=1, jpnfl),   &
                    (flyy(jfl)   , jfl=1, jpnfl),   &
                    (flzz(jfl)   , jfl=1, jpnfl),   &
                    (nisobfl(jfl), jfl=1, jpnfl),   &
                    (ngrpfl(jfl) , jfl=1, jpnfl)
         CLOSE(inum)
            
         ! Test to find the grid point coordonate with the geographical position         
         DO jfl = 1, jpnfl
            ihtest(jfl) = 0
            ivtest(jfl) = 0
            ikmfl(jfl) = 0
# if   defined key_mpp_mpi
            DO ji = MAX(nldi,2), nlei
               DO jj = MAX(nldj,2), nlej   ! NO vector opt.
# else
            DO ji = 2, jpi
               DO jj = 2, jpj   ! NO vector opt.
# endif                  
                  ! for each float we find the indexes of the mesh 
                  
                  CALL findmesh(glamf(ji-1,jj-1),gphif(ji-1,jj-1),   &
                                glamf(ji-1,jj  ),gphif(ji-1,jj  ),   &
                                glamf(ji  ,jj  ),gphif(ji  ,jj  ),   &
                                glamf(ji  ,jj-1),gphif(ji  ,jj-1),   &
                                flxx(jfl)       ,flyy(jfl)       ,   &
                                glamt(ji  ,jj  ),gphit(ji  ,jj  ), llinmesh)
                  IF(llinmesh) THEN
                     iimfl(jfl)  = ji
                     ijmfl(jfl)  = jj
                     ihtest(jfl) = ihtest(jfl)+1
                     DO jk = 1, jpk-1
                        IF( (fsdepw(ji,jj,jk) <= flzz(jfl)) .AND. (fsdepw(ji,jj,jk+1) >  flzz(jfl)) ) THEN
                           ikmfl(jfl)  = jk
                           ivtest(jfl) = ivtest(jfl) + 1
                        ENDIF
                     END DO
                  ENDIF
               END DO
            END DO
            
            ! If the float is in a mesh computed by an other processor we put iimfl=ijmfl=-1            
            IF( ihtest(jfl) == 0 ) THEN
               iimfl(jfl) = -1
               ijmfl(jfl) = -1
            ENDIF
         END DO
         
         ! A zero in the sum of the arrays "ihtest" and "ivtest"          
         IF( lk_mpp )   CALL mpp_sum(ihtest,jpnfl)   ! sums over the global domain
         IF( lk_mpp )   CALL mpp_sum(ivtest,jpnfl)

         DO jfl = 1, jpnfl
            IF( (ihtest(jfl) > 1 ) .OR. ( ivtest(jfl) > 1 )) THEN
               IF(lwp) WRITE(numout,*) 'THE FLOAT',jfl,' IS NOT IN ONLY ONE MESH'
            ENDIF
            IF( ihtest(jfl) == 0 ) THEN 
               IF(lwp) WRITE(numout,*)'THE FLOAT',jfl,' IS IN NO MESH'
            ENDIF
         END DO
        
         ! We compute the distance between the float and the face of  the mesh         
         DO jfl = 1, jpnfl
            ! Made only if the float is in the domain of the processor
            IF( (iimfl(jfl) >= 0 ) .AND. ( ijmfl(jfl) >= 0 ) ) THEN
               
               ! TEST TO KNOW IF THE FLOAT IS NOT INITIALISED IN THE COAST
               
               idomfl(jfl) = 0
               IF( tmask(iimfl(jfl),ijmfl(jfl),ikmfl(jfl)) == 0. ) idomfl(jfl)=1
               
               ! Computation of the distance between the float
               ! and the faces of the mesh
               !            zdxab
               !             .
               !        B----.---------C
               !        |    .         |
               !        |<------>flo   |
               !        |        ^     |
               !        |        |.....|....zdyad
               !        |        |     |
               !        A--------|-----D
               
               zdxab = dstnce(flxx(jfl),flyy(jfl),glamf(iimfl(jfl)-1,ijmfl(jfl)-1),flyy(jfl))                
               zdyad = dstnce(flxx(jfl),flyy(jfl),flxx(jfl),gphif(iimfl(jfl)-1,ijmfl(jfl)-1))
               
               ! Translation of this distances (in meter) in indexes
               
               tpifl(jfl) = (iimfl(jfl)-0.5)+zdxab/ e1u(iimfl(jfl)-1,ijmfl(jfl))+(mig(1)-jpizoom)
               tpjfl(jfl) = (ijmfl(jfl)-0.5)+zdyad/ e2v(iimfl(jfl),ijmfl(jfl)-1)+(mjg(1)-jpjzoom)
               tpkfl(jfl) = (fsdepw(iimfl(jfl),ijmfl(jfl),ikmfl(jfl)+1) - flzz(jfl))*(ikmfl(jfl))                     &
                          / (fsdepw(iimfl(jfl),ijmfl(jfl),ikmfl(jfl)+1) - fsdepw(iimfl(jfl),ijmfl(jfl),ikmfl(jfl)))   &
                          + (flzz(jfl) - fsdepw(iimfl(jfl),ijmfl(jfl),ikmfl(jfl)))*(ikmfl(jfl)+1)                     &
                          / (fsdepw(iimfl(jfl),ijmfl(jfl),ikmfl(jfl)+1) - fsdepw(iimfl(jfl),ijmfl(jfl),ikmfl(jfl)))
            ELSE
               tpifl (jfl) = 0.e0
               tpjfl (jfl) = 0.e0
               tpkfl (jfl) = 0.e0
               idomfl(jfl) = 0
            ENDIF
         END DO
         
         ! The sum of all the arrays tpifl, tpjfl, tpkfl give 3 arrays with the positions of all the floats. 
         IF( lk_mpp )   CALL mpp_sum( tpifl , jpnfl )   ! sums over the global domain
         IF( lk_mpp )   CALL mpp_sum( tpjfl , jpnfl )
         IF( lk_mpp )   CALL mpp_sum( tpkfl , jpnfl )
         IF( lk_mpp )   CALL mpp_sum( idomfl, jpnfl )
      ENDIF
            
      ! Print the initial positions of the floats
      IF( .NOT. ln_rstflo ) THEN 
         ! WARNING : initial position not in the sea         
         DO jfl = 1, jpnfl
            IF( idomfl(jfl) == 1 ) THEN
               IF(lwp) WRITE(numout,*)'*****************************'
               IF(lwp) WRITE(numout,*)'!!!!!!!  WARNING   !!!!!!!!!!'
               IF(lwp) WRITE(numout,*)'*****************************'
               IF(lwp) WRITE(numout,*)'The float number',jfl,'is out of the sea.'
               IF(lwp) WRITE(numout,*)'geographical position',flxx(jfl),flyy(jfl),flzz(jfl)
               IF(lwp) WRITE(numout,*)'index position',tpifl(jfl),tpjfl(jfl),tpkfl(jfl)
            ENDIF
         END DO
      ENDIF

   END SUBROUTINE flo_dom


   SUBROUTINE findmesh( pax, pay, pbx, pby,   &
                        pcx, pcy, pdx, pdy,   &
                        px  ,py  ,ptx, pty, ldinmesh )
      !! -------------------------------------------------------------
      !!                ***  ROUTINE findmesh  ***
      !!      
      !! ** Purpose :   Find the index of mesh for the point spx spy
      !!
      !! ** Method  : 
      !!----------------------------------------------------------------------
      REAL(wp) ::   &
         pax, pay, pbx, pby,    &     ! ???
         pcx, pcy, pdx, pdy,    &     ! ???
         px, py,                &     ! longitude and latitude
         ptx, pty                     ! ???
      LOGICAL ::  ldinmesh            ! ???
      !!
      REAL(wp) ::   zabt, zbct, zcdt, zdat, zabpt, zbcpt, zcdpt, zdapt
      !!---------------------------------------------------------------------
      !! Statement function
      REAL(wp) ::   fsline
      REAL(wp) ::   psax, psay, psbx, psby, psx, psy
      fsline( psax, psay, psbx, psby, psx, psy ) = psy  * ( psbx - psax )   &
         &                                       - psx  * ( psby - psay )   &
         &                                       + psax *   psby - psay * psbx
      !!---------------------------------------------------------------------
      
      ! 4 semi plane defined by the 4 points and including the T point
      zabt = fsline(pax,pay,pbx,pby,ptx,pty)
      zbct = fsline(pbx,pby,pcx,pcy,ptx,pty)
      zcdt = fsline(pcx,pcy,pdx,pdy,ptx,pty)
      zdat = fsline(pdx,pdy,pax,pay,ptx,pty)
      
      ! 4 semi plane defined by the 4 points and including the extrememity
      zabpt = fsline(pax,pay,pbx,pby,px,py)
      zbcpt = fsline(pbx,pby,pcx,pcy,px,py)
      zcdpt = fsline(pcx,pcy,pdx,pdy,px,py)
      zdapt = fsline(pdx,pdy,pax,pay,px,py)
       
      ! We compare the semi plane T with the semi plane including the point
      ! to know if it is in this  mesh.
      ! For numerical reasons it is possible that for a point which is on
      ! the line we don't have exactly zero with fsline function. We want 
      ! that a point can't be in 2 mesh in the same time, so we put the 
      ! coefficient to zero if it is smaller than 1.E-12
      
      IF( ABS(zabpt) <= 1.E-12 ) zabpt = 0.
      IF( ABS(zbcpt) <= 1.E-12 ) zbcpt = 0.
      IF( ABS(zcdpt) <= 1.E-12 ) zcdpt = 0.
      IF( ABS(zdapt) <= 1.E-12 ) zdapt = 0.
      IF( (zabt*zabpt >  0.) .AND. (zbct*zbcpt >= 0. ) .AND. ( zcdt*zcdpt >= 0. ) .AND. ( zdat*zdapt > 0. )   &
         .AND. ( px <= MAX(pcx,pdx) ) .AND. ( px >= MIN(pax,pbx) )    &
         .AND. ( py <= MAX(pby,pcy) ) .AND. ( py >= MIN(pay,pdy) ) ) THEN
         ldinmesh=.TRUE.
      ELSE
         ldinmesh=.FALSE.
      ENDIF
      !
   END SUBROUTINE findmesh


   FUNCTION dstnce( pla1, phi1, pla2, phi2 )
      !! -------------------------------------------------------------
      !!                 ***  Function dstnce  ***
      !!          
      !! ** Purpose :   returns distance (in m) between two geographical
      !!                points
      !! ** Method  : 
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(in) ::   pla1, phi1, pla2, phi2   ! ???
      !!
      REAL(wp) ::   dly1, dly2, dlx1, dlx2, dlx, dls, dld, dpi
      REAL(wp) ::   dstnce
      !!---------------------------------------------------------------------
      !
      dpi  = 2.* ASIN(1.)
      dls  = dpi / 180.
      dly1 = phi1 * dls
      dly2 = phi2 * dls
      dlx1 = pla1 * dls
      dlx2 = pla2 * dls
      !
      dlx = SIN(dly1) * SIN(dly2) + COS(dly1) * COS(dly2) * COS(dlx2-dlx1)
      !
      IF( ABS(dlx) > 1.0 ) dlx = 1.0
      !
      dld = ATAN(DSQRT( ( 1-dlx )/( 1+dlx ) )) * 222.24 / dls
      dstnce = dld * 1000.
      !
   END FUNCTION dstnce

#  else
   !!----------------------------------------------------------------------
   !!   Default option                                         Empty module
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE flo_dom                 ! Empty routine
   END SUBROUTINE flo_dom
#endif

   !!======================================================================
END MODULE flodom
