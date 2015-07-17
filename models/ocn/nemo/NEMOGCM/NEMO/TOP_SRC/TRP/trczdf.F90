MODULE trczdf
   !!==============================================================================
   !!                 ***  MODULE  trczdf  ***
   !! Ocean Passive tracers : vertical diffusive trends 
   !!=====================================================================
   !! History :  9.0  ! 2005-11 (G. Madec)  Original code
   !!       NEMO 3.0  ! 2008-01  (C. Ethe, G. Madec)  merge TRC-TRA 
   !!----------------------------------------------------------------------
#if defined key_top
   !!----------------------------------------------------------------------
   !!   'key_top'                                                TOP models
   !!----------------------------------------------------------------------
   !!   trc_ldf     : update the tracer trend with the lateral diffusion
   !!       ldf_ctl : initialization, namelist read, and parameters control
   !!----------------------------------------------------------------------
   USE oce_trc         ! ocean dynamics and active tracers
   USE trc             ! ocean passive tracers variables
   USE trcnam_trp      ! passive tracers transport namelist variables
   USE trazdf_exp      ! vertical diffusion: explicit (tra_zdf_exp     routine)
   USE trazdf_imp      ! vertical diffusion: implicit (tra_zdf_imp     routine)
   USE trdmod_oce
   USE trdtra
   USE prtctl_trc      ! Print control

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_zdf          ! called by step.F90 
   PUBLIC   trc_zdf_alloc    ! called by nemogcm.F90 

   INTEGER ::   nzdf = 0               ! type vertical diffusion algorithm used
      !                                ! defined from ln_zdf...  namlist logicals)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:) ::  r2dt   ! vertical profile time-step, = 2 rdttra
      !                                                 ! except at nit000 (=rdttra) if neuler=0

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "zdfddm_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trczdf.F90 2715 2011-03-30 15:58:35Z rblod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS
   
   INTEGER FUNCTION trc_zdf_alloc()
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trc_zdf_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( r2dt(jpk) , STAT=trc_zdf_alloc )
      !
      IF( trc_zdf_alloc /= 0 )   CALL ctl_warn('trc_zdf_alloc : failed to allocate array.')
      !
   END FUNCTION trc_zdf_alloc


   SUBROUTINE trc_zdf( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trc_zdf  ***
      !!
      !! ** Purpose :   compute the vertical ocean tracer physics.
      !!---------------------------------------------------------------------
      INTEGER, INTENT( in ) ::  kt      ! ocean time-step index
      !
      INTEGER               ::  jk, jn
      CHARACTER (len=22)    :: charout
      REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::   ztrtrd   ! 4D workspace
      !!---------------------------------------------------------------------

      IF( kt == nit000 )   CALL zdf_ctl          ! initialisation & control of options

#if ! defined key_pisces && ! defined key_bfm
      IF( neuler == 0 .AND. kt == nit000 ) THEN     ! at nit000
         r2dt(:) =  rdttrc(:)           ! = rdttrc (restarting with Euler time stepping)
      ELSEIF( kt <= nit000 + nn_dttrc ) THEN          ! at nit000 or nit000+1
         r2dt(:) = 2. * rdttrc(:)       ! = 2 rdttrc (leapfrog)
      ENDIF
#else
      r2dt(:) =  rdttrc(:)              ! = rdttrc (for PISCES and BFM use Euler time stepping)
#endif

      IF( l_trdtrc )  THEN
         ALLOCATE( ztrtrd(jpi,jpj,jpk,jptra) )   ! temporary save of trends
         ztrtrd(:,:,:,:)  = tra(:,:,:,:)
      ENDIF

      SELECT CASE ( nzdf )                       ! compute lateral mixing trend and add it to the general trend
      CASE ( -1 )                                       ! esopa: test all possibility with control print
         CALL tra_zdf_exp( kt, 'TRC', r2dt, nn_trczdf_exp, trb, tra, jptra ) 
         WRITE(charout, FMT="('zdf1 ')") ;  CALL prt_ctl_trc_info(charout)
                                            CALL prt_ctl_trc( tab4d=tra, mask=tmask, clinfo=ctrcnm, clinfo2='trd' )
         CALL tra_zdf_imp( kt, 'TRC', r2dt,                trb, tra, jptra ) 
         WRITE(charout, FMT="('zdf2 ')") ;  CALL prt_ctl_trc_info(charout)
                                            CALL prt_ctl_trc( tab4d=tra, mask=tmask, clinfo=ctrcnm, clinfo2='trd' )
      CASE ( 0 ) ;  CALL tra_zdf_exp( kt, 'TRC', r2dt, nn_trczdf_exp, trb, tra, jptra )    !   explicit scheme 
      CASE ( 1 ) ;  CALL tra_zdf_imp( kt, 'TRC', r2dt,                trb, tra, jptra )    !   implicit scheme          

      END SELECT

      IF( l_trdtra )   THEN                      ! save the vertical diffusive trends for further diagnostics
         DO jn = 1, jptra
            DO jk = 1, jpkm1
               ztrtrd(:,:,jk,jn) = ( ( tra(:,:,jk,jn) - trb(:,:,jk,jn) ) / r2dt(jk) ) - ztrtrd(:,:,jk,jn)
            END DO
            CALL trd_tra( kt, 'TRC', jn, jptra_trd_zdf, ztrtrd(:,:,:,jn) )
         END DO
         DEALLOCATE( ztrtrd )
      ENDIF

      !                                          ! print mean trends (used for debugging)
      IF( ln_ctl )   THEN
         WRITE(charout, FMT="('zdf ')") ;  CALL prt_ctl_trc_info(charout)
                                           CALL prt_ctl_trc( tab4d=tra, mask=tmask, clinfo=ctrcnm, clinfo2='trd' )
      END IF
      !
   END SUBROUTINE trc_zdf


   SUBROUTINE zdf_ctl
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE zdf_ctl  ***
      !!
      !! ** Purpose :   Choose the vertical mixing scheme
      !!
      !! ** Method  :   Set nzdf from ln_zdfexp
      !!      nzdf = 0   explicit (time-splitting) scheme (ln_trczdf_exp=T)
      !!           = 1   implicit (euler backward) scheme (ln_trczdf_exp=F)
      !!      NB: The implicit scheme is required when using : 
      !!             - rotated lateral mixing operator
      !!             - TKE, GLS or KPP vertical mixing scheme
      !!----------------------------------------------------------------------

      !  Define the vertical tracer physics scheme
      ! ==========================================

      ! Choice from ln_zdfexp already read in namelist in zdfini module
      IF( ln_trczdf_exp ) THEN           ! use explicit scheme
         nzdf = 0
      ELSE                               ! use implicit scheme
         nzdf = 1
      ENDIF

      ! Force implicit schemes
      IF( ln_trcldf_iso                               )   nzdf = 1      ! iso-neutral lateral physics
      IF( ln_trcldf_hor .AND. ln_sco                  )   nzdf = 1      ! horizontal lateral physics in s-coordinate
#if defined key_zdftke || defined key_zdfgls || defined key_zdfkpp
                                                          nzdf = 1      ! TKE, GLS or KPP physics       
#endif
      IF( ln_trczdf_exp .AND. nzdf == 1 )   THEN
         CALL ctl_stop( 'trc_zdf : If using the rotated lateral mixing operator or TKE, GLS or KPP vertical scheme ', &
            &           '          the implicit scheme is required, set ln_trczdf_exp = .false.' )
      ENDIF

      ! Test: esopa
      IF( lk_esopa )    nzdf = -1                      ! All schemes used

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'trc:zdf_ctl : vertical passive tracer physics scheme'
         WRITE(numout,*) '~~~~~~~~~~~'
         IF( nzdf == -1 )   WRITE(numout,*) '              ESOPA test All scheme used'
         IF( nzdf ==  0 )   WRITE(numout,*) '              Explicit time-splitting scheme'
         IF( nzdf ==  1 )   WRITE(numout,*) '              Implicit (euler backward) scheme'
      ENDIF

   END SUBROUTINE zdf_ctl
#else
   !!----------------------------------------------------------------------
   !!   Default option                                         Empty module
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_zdf( kt )
      INTEGER, INTENT(in) :: kt  
      WRITE(*,*) 'trc_zdf: You should not have seen this print! error?', kt
   END SUBROUTINE trc_zdf
#endif
   !!==============================================================================
END MODULE trczdf
