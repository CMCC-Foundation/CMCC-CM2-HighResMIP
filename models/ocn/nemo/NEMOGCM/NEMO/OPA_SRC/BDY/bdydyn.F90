MODULE bdydyn
   !!======================================================================
   !!                       ***  MODULE  bdydyn  ***
   !! Unstructured Open Boundary Cond. :   Flow relaxation scheme on velocities
   !!======================================================================
   !! History :  1.0  !  2005-02  (J. Chanut, A. Sellar)  Original code
   !!             -   !  2007-07  (D. Storkey) Move Flather implementation to separate routine.
   !!            3.0  !  2008-04  (NEMO team)  add in the reference version
   !!            3.2  !  2008-04  (R. Benshila) consider velocity instead of transport 
   !!            3.3  !  2010-09  (E.O'Dea) modifications for Shelf configurations 
   !!            3.3  !  2010-09  (D.Storkey) add ice boundary conditions
   !!----------------------------------------------------------------------
#if defined key_bdy 
   !!----------------------------------------------------------------------
   !!   'key_bdy' :                    Unstructured Open Boundary Condition
   !!----------------------------------------------------------------------
   !!   bdy_dyn_frs    : relaxation of velocities on unstructured open boundary
   !!   bdy_dyn_fla    : Flather condition for barotropic solution
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers 
   USE dom_oce         ! ocean space and time domain
   USE bdy_oce         ! ocean open boundary conditions
   USE dynspg_oce      ! for barotropic variables
   USE phycst          ! physical constants
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE bdytides        ! for tidal harmonic forcing at boundary
   USE in_out_manager  !

   IMPLICIT NONE
   PRIVATE

   PUBLIC   bdy_dyn_frs   ! routine called in dynspg_flt (free surface case ONLY)
# if defined key_dynspg_exp || defined key_dynspg_ts
   PUBLIC   bdy_dyn_fla   ! routine called in dynspg_flt (free surface case ONLY)
# endif

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: bdydyn.F90 2528 2010-12-27 17:33:53Z rblod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE bdy_dyn_frs( kt )
      !!----------------------------------------------------------------------
      !!                  ***  SUBROUTINE bdy_dyn_frs  ***
      !!
      !! ** Purpose : - Apply the Flow Relaxation Scheme for dynamic in the  
      !!                case of unstructured open boundaries.
      !!
      !! References :- Engedahl H., 1995: Use of the flow relaxation scheme in 
      !!               a three-dimensional baroclinic ocean model with realistic
      !!               topography. Tellus, 365-382.
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt   ! Main time step counter
      !!
      INTEGER  ::   jb, jk         ! dummy loop indices
      INTEGER  ::   ii, ij, igrd   ! local integers
      REAL(wp) ::   zwgt           ! boundary weight
      !!----------------------------------------------------------------------
      !
      IF(ln_dyn_frs) THEN       ! If this is false, then this routine does nothing. 
         !
         IF( kt == nit000 ) THEN
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'bdy_dyn_frs : Flow Relaxation Scheme on momentum'
            IF(lwp) WRITE(numout,*) '~~~~~~~'
         ENDIF
         !
         igrd = 2                      ! Relaxation of zonal velocity
         DO jb = 1, nblen(igrd)
            DO jk = 1, jpkm1
               ii   = nbi(jb,igrd)
               ij   = nbj(jb,igrd)
               zwgt = nbw(jb,igrd)
               ua(ii,ij,jk) = ( ua(ii,ij,jk) * ( 1.- zwgt ) + ubdy(jb,jk) * zwgt ) * umask(ii,ij,jk)
            END DO
         END DO
         !
         igrd = 3                      ! Relaxation of meridional velocity
         DO jb = 1, nblen(igrd)
            DO jk = 1, jpkm1
               ii   = nbi(jb,igrd)
               ij   = nbj(jb,igrd)
               zwgt = nbw(jb,igrd)
               va(ii,ij,jk) = ( va(ii,ij,jk) * ( 1.- zwgt ) + vbdy(jb,jk) * zwgt ) * vmask(ii,ij,jk)
            END DO
         END DO 
         CALL lbc_lnk( ua, 'U', -1. )   ;   CALL lbc_lnk( va, 'V', -1. )   ! Boundary points should be updated
         !
      ENDIF ! ln_dyn_frs
      !
   END SUBROUTINE bdy_dyn_frs


# if defined   key_dynspg_exp   ||   defined key_dynspg_ts
   !!----------------------------------------------------------------------
   !!   'key_dynspg_exp'        OR              explicit sea surface height
   !!   'key_dynspg_ts '                  split-explicit sea surface height
   !!----------------------------------------------------------------------
   
!! Option to use Flather with dynspg_flt not coded yet...

   SUBROUTINE bdy_dyn_fla( pssh )
      !!----------------------------------------------------------------------
      !!                 ***  SUBROUTINE bdy_dyn_fla  ***
      !!             
      !!              - Apply Flather boundary conditions on normal barotropic velocities 
      !!                (ln_dyn_fla=.true. or ln_tides=.true.)
      !!
      !! ** WARNINGS about FLATHER implementation:
      !!1. According to Palma and Matano, 1998 "after ssh" is used. 
      !!   In ROMS and POM implementations, it is "now ssh". In the current 
      !!   implementation (tested only in the EEL-R5 conf.), both cases were unstable. 
      !!   So I use "before ssh" in the following.
      !!
      !!2. We assume that the normal ssh gradient at the bdy is zero. As a matter of 
      !!   fact, the model ssh just inside the dynamical boundary is used (the outside  
      !!   ssh in the code is not updated).
      !!
      !! References:  Flather, R. A., 1976: A tidal model of the northwest European
      !!              continental shelf. Mem. Soc. R. Sci. Liege, Ser. 6,10, 141-164.     
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   pssh

      INTEGER  ::   jb, igrd                         ! dummy loop indices
      INTEGER  ::   ii, ij, iim1, iip1, ijm1, ijp1   ! 2D addresses
      REAL(wp) ::   zcorr                            ! Flather correction
      REAL(wp) ::   zforc                            ! temporary scalar
      !!----------------------------------------------------------------------

      ! ---------------------------------!
      ! Flather boundary conditions     :!
      ! ---------------------------------! 
     
      IF(ln_dyn_fla .OR. ln_tides) THEN ! If these are both false, then this routine does nothing. 

         ! Fill temporary array with ssh data (here spgu):
         igrd = 4
         spgu(:,:) = 0.0
         DO jb = 1, nblenrim(igrd)
            ii = nbi(jb,igrd)
            ij = nbj(jb,igrd)
            IF( ln_dyn_fla ) spgu(ii, ij) = sshbdy(jb)
            IF( ln_tides )   spgu(ii, ij) = spgu(ii, ij) + sshtide(jb)
         END DO
         !
         igrd = 5      ! Flather bc on u-velocity; 
         !             ! remember that flagu=-1 if normal velocity direction is outward
         !             ! I think we should rather use after ssh ?
         DO jb = 1, nblenrim(igrd)
            ii  = nbi(jb,igrd)
            ij  = nbj(jb,igrd) 
            iim1 = ii + MAX( 0, INT( flagu(jb) ) )   ! T pts i-indice inside the boundary
            iip1 = ii - MIN( 0, INT( flagu(jb) ) )   ! T pts i-indice outside the boundary 
            !
            zcorr = - flagu(jb) * SQRT( grav * hur_e(ii, ij) ) * ( pssh(iim1, ij) - spgu(iip1,ij) )
            zforc = ubtbdy(jb) + utide(jb)
            ua_e(ii,ij) = zforc + zcorr * umask(ii,ij,1) 
         END DO
         !
         igrd = 6      ! Flather bc on v-velocity
         !             ! remember that flagv=-1 if normal velocity direction is outward
         DO jb = 1, nblenrim(igrd)
            ii  = nbi(jb,igrd)
            ij  = nbj(jb,igrd) 
            ijm1 = ij + MAX( 0, INT( flagv(jb) ) )   ! T pts j-indice inside the boundary
            ijp1 = ij - MIN( 0, INT( flagv(jb) ) )   ! T pts j-indice outside the boundary 
            !
            zcorr = - flagv(jb) * SQRT( grav * hvr_e(ii, ij) ) * ( pssh(ii, ijm1) - spgu(ii,ijp1) )
            zforc = vbtbdy(jb) + vtide(jb)
            va_e(ii,ij) = zforc + zcorr * vmask(ii,ij,1)
         END DO
         CALL lbc_lnk( ua_e, 'U', -1. )   ! Boundary points should be updated
         CALL lbc_lnk( va_e, 'V', -1. )   !
         !
      ENDIF ! ln_dyn_fla .or. ln_tides
      !
   END SUBROUTINE bdy_dyn_fla
#endif

#else
   !!----------------------------------------------------------------------
   !!   Dummy module                   NO Unstruct Open Boundary Conditions
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE bdy_dyn_frs( kt )      ! Empty routine
      WRITE(*,*) 'bdy_dyn_frs: You should not have seen this print! error?', kt
   END SUBROUTINE bdy_dyn_frs
   SUBROUTINE bdy_dyn_fla( pssh )    ! Empty routine
      REAL :: pssh(:,:)
      WRITE(*,*) 'bdy_dyn_fla: You should not have seen this print! error?', pssh(1,1)
   END SUBROUTINE bdy_dyn_fla
#endif

   !!======================================================================
END MODULE bdydyn
