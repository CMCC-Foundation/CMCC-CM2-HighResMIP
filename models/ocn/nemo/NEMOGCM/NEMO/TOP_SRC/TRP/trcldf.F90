MODULE trcldf
   !!======================================================================
   !!                       ***  MODULE  trcldf  ***
   !! Ocean Passive tracers : lateral diffusive trends 
   !!=====================================================================
   !! History :  9.0  ! 2005-11 (G. Madec)  Original code
   !!       NEMO 3.0  ! 2008-01  (C. Ethe, G. Madec)  merge TRC-TRA 
   !!----------------------------------------------------------------------
#if defined key_top
   !!----------------------------------------------------------------------
   !!   'key_top'                                                TOP models
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   trc_ldf     : update the tracer trend with the lateral diffusion
   !!       ldf_ctl : initialization, namelist read, and parameters control
   !!----------------------------------------------------------------------
   USE oce_trc         ! ocean dynamics and active tracers
   USE trc             ! ocean passive tracers variables
   USE trcnam_trp      ! passive tracers transport namelist variables
   USE ldftra_oce      ! lateral diffusion coefficient on tracers
   USE ldfslp          ! ???
   USE traldf_bilapg   ! lateral mixing            (tra_ldf_bilapg routine)
   USE traldf_bilap    ! lateral mixing            (tra_ldf_bilap routine)
   USE traldf_iso      ! lateral mixing            (tra_ldf_iso routine)
   USE traldf_lap      ! lateral mixing            (tra_ldf_lap routine)
   USE trdmod_oce
   USE trdtra
   USE prtctl_trc      ! Print control

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_ldf    ! called by step.F90 
   !                                                 !!: ** lateral mixing namelist (nam_trcldf) **
   INTEGER ::   nldf = 0   ! type of lateral diffusion used defined from ln_trcldf_... namlist logicals)
   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcldf.F90 2715 2011-03-30 15:58:35Z rblod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trc_ldf( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_ldf  ***
      !! 
      !! ** Purpose :   compute the lateral ocean tracer physics.
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step index
      !!
      INTEGER            :: jn
      CHARACTER (len=22) :: charout
      REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::   ztrtrd
      !!----------------------------------------------------------------------

      IF( kt == nit000 )   CALL ldf_ctl          ! initialisation & control of options

      IF( l_trdtrc )  THEN 
         ALLOCATE( ztrtrd(jpi,jpj,jpk,jptra) )  ! temporary save of trends
         ztrtrd(:,:,:,:)  = tra(:,:,:,:)
      ENDIF

      SELECT CASE ( nldf )                       ! compute lateral mixing trend and add it to the general trend
      CASE ( 0 )   ;   CALL tra_ldf_lap   ( kt, 'TRC', gtru, gtrv, trb, tra, jptra            )  ! iso-level laplacian
      CASE ( 1 )   ;   CALL tra_ldf_iso   ( kt, 'TRC', gtru, gtrv, trb, tra, jptra, rn_ahtb_0 )  ! rotated laplacian 
      CASE ( 2 )   ;   CALL tra_ldf_bilap ( kt, 'TRC', gtru, gtrv, trb, tra, jptra            )  ! iso-level bilaplacian
      CASE ( 3 )   ;   CALL tra_ldf_bilapg( kt, 'TRC',             trb, tra, jptra            )  ! s-coord. horizontal bilaplacian
         !
      CASE ( -1 )                                     ! esopa: test all possibility with control print
         CALL tra_ldf_lap   ( kt, 'TRC', gtru, gtrv, trb, tra, jptra            )
         WRITE(charout, FMT="('ldf0 ')") ;  CALL prt_ctl_trc_info(charout)
                                            CALL prt_ctl_trc( tab4d=tra, mask=tmask, clinfo=ctrcnm, clinfo2='trd' )
         CALL tra_ldf_iso   ( kt, 'TRC', gtru, gtrv, trb, tra, jptra, rn_ahtb_0 )
         WRITE(charout, FMT="('ldf1 ')") ;  CALL prt_ctl_trc_info(charout)
                                            CALL prt_ctl_trc( tab4d=tra, mask=tmask, clinfo=ctrcnm, clinfo2='trd' )
         CALL tra_ldf_bilap ( kt, 'TRC', gtru, gtrv, trb, tra, jptra            )
         WRITE(charout, FMT="('ldf2 ')") ;  CALL prt_ctl_trc_info(charout)
                                            CALL prt_ctl_trc( tab4d=tra, mask=tmask, clinfo=ctrcnm, clinfo2='trd' )
         CALL tra_ldf_bilapg( kt, 'TRC',             trb, tra, jptra            )
         WRITE(charout, FMT="('ldf3 ')") ;  CALL prt_ctl_trc_info(charout)
                                            CALL prt_ctl_trc( tab4d=tra, mask=tmask, clinfo=ctrcnm, clinfo2='trd' )
      END SELECT
      !
      IF( l_trdtrc )   THEN                      ! save the horizontal diffusive trends for further diagnostics
        DO jn = 1, jptra
           ztrtrd(:,:,:,jn) = tra(:,:,:,jn) - ztrtrd(:,:,:,jn)
           CALL trd_tra( kt, 'TRC', jn, jptra_trd_ldf, ztrtrd(:,:,:,jn) )
        END DO
        DEALLOCATE( ztrtrd ) 
      ENDIF
      !                                          ! print mean trends (used for debugging)
      IF( ln_ctl )   THEN
         WRITE(charout, FMT="('ldf ')") ;  CALL prt_ctl_trc_info(charout)
                                           CALL prt_ctl_trc( tab4d=tra, mask=tmask, clinfo=ctrcnm, clinfo2='trd' )
      ENDIF
      !
   END SUBROUTINE trc_ldf


   SUBROUTINE ldf_ctl
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ldf_ctl  ***
      !! 
      !! ** Purpose :   Choice of the operator for the lateral tracer diffusion
      !!
      !! ** Method  :   set nldf from the namtra_ldf logicals
      !!      nldf == -2   No lateral diffusion  
      !!      nldf == -1   ESOPA test: ALL operators are used
      !!      nldf ==  0   laplacian operator
      !!      nldf ==  1   Rotated laplacian operator
      !!      nldf ==  2   bilaplacian operator
      !!      nldf ==  3   Rotated bilaplacian
      !!----------------------------------------------------------------------
      INTEGER ::   ioptio, ierr         ! temporary integers 
      !!----------------------------------------------------------------------

      !  Define the lateral mixing oparator for tracers
      ! ===============================================
    
      !                               ! control the input
      ioptio = 0
      IF( ln_trcldf_lap   )   ioptio = ioptio + 1
      IF( ln_trcldf_bilap )   ioptio = ioptio + 1
      IF( ioptio >  1 )   CALL ctl_stop( '          use ONE or NONE of the 2 lap/bilap operator type on tracer' )
      IF( ioptio == 0 )   nldf = -2   ! No lateral diffusion
      ioptio = 0
      IF( ln_trcldf_level )   ioptio = ioptio + 1
      IF( ln_trcldf_hor   )   ioptio = ioptio + 1
      IF( ln_trcldf_iso   )   ioptio = ioptio + 1
      IF( ioptio /= 1 )   CALL ctl_stop( '          use only ONE direction (level/hor/iso)' )

      ! defined the type of lateral diffusion from ln_trcldf_... logicals
      ! CAUTION : nldf = 1 is used in trazdf_imp, change it carefully
      ierr = 0
      IF( ln_trcldf_lap ) THEN       ! laplacian operator
         IF ( ln_zco ) THEN                ! z-coordinate
            IF ( ln_trcldf_level )   nldf = 0      ! iso-level  (no rotation)
            IF ( ln_trcldf_hor   )   nldf = 0      ! horizontal (no rotation)
            IF ( ln_trcldf_iso   )   nldf = 1      ! isoneutral (   rotation)
         ENDIF
         IF ( ln_zps ) THEN             ! z-coordinate
            IF ( ln_trcldf_level )   ierr = 1      ! iso-level not allowed
            IF ( ln_trcldf_hor   )   nldf = 0      ! horizontal (no rotation)
            IF ( ln_trcldf_iso   )   nldf = 1      ! isoneutral (   rotation)
         ENDIF
         IF ( ln_sco ) THEN             ! z-coordinate
            IF ( ln_trcldf_level )   nldf = 0      ! iso-level  (no rotation)
            IF ( ln_trcldf_hor   )   nldf = 1      ! horizontal (   rotation)
            IF ( ln_trcldf_iso   )   nldf = 1      ! isoneutral (   rotation)
         ENDIF
      ENDIF

      IF( ln_trcldf_bilap ) THEN      ! bilaplacian operator
         IF ( ln_zco ) THEN                ! z-coordinate
            IF ( ln_trcldf_level )   nldf = 2      ! iso-level  (no rotation)
            IF ( ln_trcldf_hor   )   nldf = 2      ! horizontal (no rotation)
            IF ( ln_trcldf_iso   )   ierr = 2      ! isoneutral (   rotation)
         ENDIF
         IF ( ln_zps ) THEN             ! z-coordinate
            IF ( ln_trcldf_level )   ierr = 1      ! iso-level not allowed 
            IF ( ln_trcldf_hor   )   nldf = 2      ! horizontal (no rotation)
            IF ( ln_trcldf_iso   )   ierr = 2      ! isoneutral (   rotation)
         ENDIF
         IF ( ln_sco ) THEN             ! z-coordinate
            IF ( ln_trcldf_level )   nldf = 2      ! iso-level  (no rotation)
            IF ( ln_trcldf_hor   )   nldf = 3      ! horizontal (   rotation)
            IF ( ln_trcldf_iso   )   ierr = 2      ! isoneutral (   rotation)
         ENDIF
      ENDIF

      IF( ierr == 1 )   CALL ctl_stop( ' iso-level in z-coordinate - partial step, not allowed' )
      IF( ierr == 2 )   CALL ctl_stop( ' isoneutral bilaplacian operator does not exist' )
      IF( lk_traldf_eiv .AND. .NOT.ln_trcldf_iso )   &
           CALL ctl_stop( '          eddy induced velocity on tracers',   &
           &              ' the eddy induced velocity on tracers requires isopycnal laplacian diffusion' )
      IF( nldf == 1 .OR. nldf == 3 ) THEN      ! rotation
         IF( .NOT.lk_ldfslp )   CALL ctl_stop( '          the rotation of the diffusive tensor require key_ldfslp' )
#if defined key_offline
         l_traldf_rot = .TRUE.                 ! needed for trazdf_imp
#endif
      ENDIF

      IF( lk_esopa ) THEN
         IF(lwp) WRITE(numout,*) '          esopa control: use all lateral physics options'
         nldf = -1
      ENDIF

      IF( .NOT. ln_trcldf_diff ) THEN
         IF(lwp) WRITE(numout,*) '          No lateral diffusion on passive tracers'
         nldf = -2
      ENDIF

      IF(lwp) THEN
         WRITE(numout,*)
         IF( nldf == -2 )   WRITE(numout,*) '          NO lateral diffusion'
         IF( nldf == -1 )   WRITE(numout,*) '          ESOPA test All scheme used'
         IF( nldf ==  0 )   WRITE(numout,*) '          laplacian operator'
         IF( nldf ==  1 )   WRITE(numout,*) '          Rotated laplacian operator'
         IF( nldf ==  2 )   WRITE(numout,*) '          bilaplacian operator'
         IF( nldf ==  3 )   WRITE(numout,*) '          Rotated bilaplacian'
      ENDIF

      !
   END SUBROUTINE ldf_ctl
#else
   !!----------------------------------------------------------------------
   !!   Default option                                         Empty module
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_ldf( kt )
      INTEGER, INTENT(in) :: kt
      WRITE(*,*) 'trc_ldf: You should not have seen this print! error?', kt
   END SUBROUTINE trc_ldf
#endif
   !!======================================================================
END MODULE trcldf
