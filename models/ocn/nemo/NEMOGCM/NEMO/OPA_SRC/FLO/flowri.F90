MODULE flowri
   !!======================================================================
   !!                       ***  MODULE  flowri  ***
   !! lagrangian floats :   outputs
   !!======================================================================
   !! History :   OPA  ! 1999-09  (Y. Drillet)  Original code
   !!                  ! 2000-06  (J.-M. Molines)  Profiling floats for CLS 
   !!   NEMO      1.0  ! 2002-11  (G. Madec, A. Bozec)  F90: Free form and module
   !!----------------------------------------------------------------------
#if   defined key_floats   ||   defined key_esopa
   !!----------------------------------------------------------------------
   !!   'key_floats'                                     float trajectories
   !!----------------------------------------------------------------------
   !!    flowri     : write trajectories of floats in file 
   !!----------------------------------------------------------------------
   USE flo_oce         ! ocean drifting floats
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE lib_mpp         ! distribued memory computing library
   USE in_out_manager  ! I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   flo_wri         ! routine called by floats.F90
   PUBLIC   flo_wri_alloc   ! routine called by floats.F90

   INTEGER ::   jfl      ! number of floats
   INTEGER ::   numflo   ! logical unit for drifting floats

   ! Following are only workspace arrays but shape is not (jpi,jpj) and
   ! therefore make them module arrays rather than replacing with wrk_nemo
   ! member arrays.
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   ztemp, zsal   ! 2D workspace

   !! * Substitutions
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 4.0 , NEMO Consortium (2011)
   !! $Id: flowri.F90 2715 2011-03-30 15:58:35Z rblod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION flo_wri_alloc
      !!-------------------------------------------------------------------
      !!                ***  FUNCTION flo_wri_alloc  ***
      !!-------------------------------------------------------------------
      ALLOCATE( ztemp(jpk,jpnfl) , zsal(jpk,jpnfl) , STAT=flo_wri_alloc)
      !
      IF( lk_mpp             )   CALL mpp_sum ( flo_wri_alloc )
      IF( flo_wri_alloc /= 0 )   CALL ctl_warn('flo_wri_alloc: failed to allocate arrays.')
   END FUNCTION flo_wri_alloc


   SUBROUTINE flo_wri( kt )
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE flo_wri  ***
      !!             
      !! ** Purpose :   Write position of floats in "trajec_float" file
      !!      and the temperature and salinity at this position
      !!      
      !! ** Method  :   The frequency is nn_writefl
      !!----------------------------------------------------------------------
      INTEGER ::   kt   ! time step
      !!
      CHARACTER (len=21) ::  clname
      INTEGER ::   inum   ! temporary logical unit for restart file
      INTEGER ::   iafl, ibfl, icfl, ia1fl, ib1fl, ic1fl, jfl, irecflo
      INTEGER ::   iafloc, ibfloc, ia1floc, ib1floc, iafln, ibfln
      INTEGER  ::    ic, jc , jpn
      INTEGER, DIMENSION ( jpnij )  ::   iproc
      REAL(wp) ::   zafl, zbfl, zcfl, zdtj
      REAL(wp) ::   zxxu, zxxu_01,zxxu_10, zxxu_11
      !!---------------------------------------------------------------------
      
      IF( kt == nit000 .OR. MOD( kt,nn_writefl) == 0 ) THEN 

         ! header of output floats file
      
         IF(lwp) THEN
            WRITE(numout,*)
            WRITE(numout,*) 'flo_wri : write in trajec_float file '
            WRITE(numout,*) '~~~~~~~    '
         ENDIF

         ! open the file numflo 
         CALL ctl_opn( numflo, 'trajec_float', 'REPLACE', 'UNFORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )

         IF( kt == nit000 ) THEN
            irecflo = NINT( (nitend-nit000) / FLOAT(nn_writefl) )
            IF(lwp) WRITE(numflo)cexper,no,irecflo,jpnfl,nn_writefl
         ENDIF
         zdtj = rdt / 86400._wp

         ! translation of index position in geographical position

         IF( lk_mpp ) THEN
            DO jfl = 1, jpnfl
               iafl  = INT ( tpifl(jfl) )
               ibfl  = INT ( tpjfl(jfl) )
               icfl  = INT ( tpkfl(jfl) )
               iafln = NINT( tpifl(jfl) )
               ibfln = NINT( tpjfl(jfl) )
               ia1fl = iafl + 1
               ib1fl = ibfl + 1
               ic1fl = icfl + 1
               zafl  = tpifl(jfl) - FLOAT( iafl )
               zbfl  = tpjfl(jfl) - FLOAT( ibfl )
               zcfl  = tpkfl(jfl) - FLOAT( icfl )
               IF(   iafl >= mig(nldi)-jpizoom+1 .AND. iafl <= mig(nlei)-jpizoom+1 .AND.   &
                  &  ibfl >= mjg(nldj)-jpjzoom+1 .AND. ibfl <= mjg(nlej)-jpjzoom+1       ) THEN

                  ! local index

                  iafloc  = iafl -(mig(1)-jpizoom+1) + 1
                  ibfloc  = ibfl -(mjg(1)-jpjzoom+1) + 1
                  ia1floc = iafloc + 1
                  ib1floc = ibfloc + 1

                  flyy(jfl) = (1.-zafl)*(1.-zbfl)*gphit(iafloc ,ibfloc ) + (1.-zafl) * zbfl * gphit(iafloc ,ib1floc)   &
                     &      +     zafl *(1.-zbfl)*gphit(ia1floc,ibfloc ) +     zafl  * zbfl * gphit(ia1floc,ib1floc)
                  flxx(jfl) = (1.-zafl)*(1.-zbfl)*glamt(iafloc ,ibfloc ) + (1.-zafl) * zbfl * glamt(iafloc ,ib1floc)   &
                     &      +     zafl *(1.-zbfl)*glamt(ia1floc,ibfloc ) +     zafl  * zbfl * glamt(ia1floc,ib1floc)
                  flzz(jfl) = (1.-zcfl)*fsdepw(iafloc,ibfloc,icfl ) + zcfl * fsdepw(iafloc,ibfloc,ic1fl)

                  ! Change  by Alexandra Bozec et Jean-Philippe Boulanger
                  ! We save  the instantaneous profile of T and S of the column     
                  ! ztemp(jfl)=tn(iafloc,ibfloc,icfl)
                  ! zsal(jfl)=sn(iafloc,ibfloc,icfl)
                  ztemp(1:jpk,jfl) = tn(iafloc,ibfloc,1:jpk)
                  zsal (1:jpk,jfl) = sn(iafloc,ibfloc,1:jpk)            
               ELSE
                  flxx(jfl) = 0.
                  flyy(jfl) = 0.
                  flzz(jfl) = 0.
                  ztemp(1:jpk,jfl) = 0.
                  zsal (1:jpk,jfl) = 0.
               ENDIF
            END DO

            CALL mpp_sum( flxx, jpnfl )   ! sums over the global domain
            CALL mpp_sum( flyy, jpnfl )
            CALL mpp_sum( flzz, jpnfl )
            ! these 2 lines have accendentaly been removed from ATL6-V8 run hence
            ! giving 0 salinity and temperature on the float trajectory
!bug RB
!compilation failed in mpp
!            CALL mpp_sum( ztemp, jpk*jpnfl )
!            CALL mpp_sum( zsal , jpk*jpnfl )

         ELSE
            DO jfl = 1, jpnfl
               iafl  = INT (tpifl(jfl))
               ibfl  = INT (tpjfl(jfl))
               icfl  = INT (tpkfl(jfl))
               iafln = NINT(tpifl(jfl))
               ibfln = NINT(tpjfl(jfl))
               ia1fl = iafl+1
               ib1fl = ibfl+1
               ic1fl = icfl+1
               zafl  = tpifl(jfl) - FLOAT(iafl)
               zbfl  = tpjfl(jfl) - FLOAT(ibfl)
               zcfl  = tpkfl(jfl) - FLOAT(icfl)
               iafloc  = iafl
               ibfloc  = ibfl
               ia1floc = iafloc + 1
               ib1floc = ibfloc + 1
               !
               flyy(jfl) = (1.-zafl)*(1.-zbfl)*gphit(iafloc ,ibfloc ) + (1.-zafl) * zbfl * gphit(iafloc ,ib1floc)   &
                         +     zafl *(1.-zbfl)*gphit(ia1floc,ibfloc ) +     zafl  * zbfl * gphit(ia1floc,ib1floc)
               flxx(jfl) = (1.-zafl)*(1.-zbfl)*glamt(iafloc ,ibfloc ) + (1.-zafl) * zbfl * glamt(iafloc ,ib1floc)   &
                         +     zafl *(1.-zbfl)*glamt(ia1floc,ibfloc ) +     zafl  * zbfl * glamt(ia1floc,ib1floc)
               flzz(jfl) = (1.-zcfl)*fsdepw(iafloc,ibfloc,icfl ) + zcfl * fsdepw(iafloc,ibfloc,ic1fl)
               !ALEX
               ! Astuce pour ne pas avoir des flotteurs qui se baladent sur IDL
               zxxu_11 = glamt(iafloc ,ibfloc )
               zxxu_10 = glamt(iafloc ,ib1floc)
               zxxu_01 = glamt(ia1floc,ibfloc )
               zxxu    = glamt(ia1floc,ib1floc)

               IF( iafloc == 52 )  zxxu_10 = -181
               IF( iafloc == 52 )  zxxu_11 = -181
               flxx(jfl)=(1.-zafl)*(1.-zbfl)* zxxu_11 + (1.-zafl)*    zbfl * zxxu_10   &
                        +    zafl *(1.-zbfl)* zxxu_01 +     zafl *    zbfl * zxxu
               !ALEX         
               ! Change  by Alexandra Bozec et Jean-Philippe Boulanger
               ! We save  the instantaneous profile of T and S of the column     
               !     ztemp(jfl)=tn(iafloc,ibfloc,icfl)
               !     zsal(jfl)=sn(iafloc,ibfloc,icfl)
               ztemp(1:jpk,jfl) = tn(iafloc,ibfloc,1:jpk)
               zsal (1:jpk,jfl) = sn(iafloc,ibfloc,1:jpk)
            END DO
         ENDIF

         !
         WRITE(numflo) flxx,flyy,flzz,nisobfl,ngrpfl,ztemp,zsal, FLOAT(ndastp)
      !!
      !! case when profiles are dumped. In order to save memory, dumps are
      !! done level by level.
      !      IF (mod(kt,nflclean) == 0.) THEN
      !!     IF ( nwflo == nwprofil ) THEN
      !        DO jk = 1,jpk
      !         DO jfl=1,jpnfl
      !         iafl= INT(tpifl(jfl))
      !         ibfl=INT(tpjfl(jfl))
      !         iafln=NINT(tpifl(jfl))
      !         ibfln=NINT(tpjfl(jfl))
      !# if defined key_mpp_mpi   
      !        IF ( (iafl >= (mig(nldi)-jpizoom+1)) .AND.
      !     $       (iafl <= (mig(nlei)-jpizoom+1)) .AND.
      !     $       (ibfl >= (mjg(nldj)-jpjzoom+1)) .AND.
      !     $       (ibfl <= (mjg(nlej)-jpjzoom+1)) ) THEN
      !!
      !! local index
      !!
      !         iafloc=iafln-(mig(1)-jpizoom+1)+1
      !         ibfloc=ibfln-(mjg(1)-jpjzoom+1)+1
      !!         IF (jk == 1 ) THEN
      !!      PRINT *,'<<<>>> ',jfl,narea, iafloc ,ibfloc, iafln, ibfln,adatrj
      !!         ENDIF
      !# else
      !         iafloc=iafln
      !         ibfloc=ibfln
      !# endif
      !         ztemp(jfl)=tn(iafloc,ibfloc,jk)
      !         zsal(jfl)=sn(iaflo!,ibfloc,jk)
      !# if defined key_mpp_mpi   
      !        ELSE
      !         ztemp(jfl) = 0.
      !         zsal(jfl) = 0.
      !        ENDIF
      !# endif
      !! ... next float
      !        END DO
      !      IF( lk_mpp )   CALL mpp_sum( ztemp, jpnfl )
      !      IF( lk_mpp )   CALL mpp_sum( zsal , jpnfl )
      !
      !      IF (lwp) THEN 
      !         WRITE(numflo) ztemp, zsal
      !      ENDIF
      !! ... next level jk
      !      END DO
      !! ... reset nwflo to 0 for ALL processors, if profile has been written
      !!       nwflo = 0
      !      ENDIF
      !!
      !      CALL flush (numflo)
      !! ... time of dumping floats
      !!      END IF
      ENDIF
      
      IF( (MOD(kt,nn_stockfl) == 0) .OR. ( kt == nitend ) ) THEN 
         ! Writing the restart file 
         IF(lwp) THEN
            WRITE(numout,*)
            WRITE(numout,*) 'flo_wri : write in  restart_float file '
            WRITE(numout,*) '~~~~~~~    '
         ENDIF

         ! file is opened and closed every time it is used.

         clname = 'restart.float.'
         ic = 1
         DO jc = 1, 16
            IF( cexper(jc:jc) /= ' ' ) ic = jc
         END DO
         clname = clname(1:14)//cexper(1:ic)
         ic = 1
         DO jc = 1, 48
            IF( clname(jc:jc) /= ' ' ) ic = jc
         END DO

         CALL ctl_opn( inum, clname, 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )
         REWIND inum
         !
         DO jpn = 1, jpnij
            iproc(jpn) = 0
         END DO
         !
         IF(lwp) THEN
            REWIND(inum)
            WRITE (inum) tpifl,tpjfl,tpkfl,nisobfl,ngrpfl
            CLOSE (inum) 
         ENDIF
         !
         ! Compute the number of trajectories for each processor
         !
         IF( lk_mpp ) THEN
            DO jfl = 1, jpnfl
               IF( (INT(tpifl(jfl)) >= (mig(nldi)-jpizoom+1)) .AND.   &
                  &(INT(tpifl(jfl)) <= (mig(nlei)-jpizoom+1)) .AND.   &
                  &(INT(tpjfl(jfl)) >= (mjg(nldj)-jpjzoom+1)) .AND.   &
                  &(INT(tpjfl(jfl)) <= (mjg(nlej)-jpjzoom+1)) ) THEN
                  iproc(narea) = iproc(narea)+1
               ENDIF
            END DO
            CALL mpp_sum( iproc, jpnij )
            !
            IF(lwp) THEN 
               WRITE(numout,*) 'DATE',adatrj
               DO jpn = 1, jpnij
                  IF( iproc(jpn) /= 0 ) THEN
                     WRITE(numout,*)'PROCESSOR',jpn-1,'compute',iproc(jpn), 'trajectories.'
                  ENDIF
               END DO
            ENDIF
         ENDIF
      ENDIF 

      IF( kt == nitend )   CLOSE( numflo ) 
      !
   END SUBROUTINE flo_wri

#  else
   !!----------------------------------------------------------------------
   !!   Default option                                         Empty module
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE flo_wri                 ! Empty routine
   END SUBROUTINE flo_wri
#endif
   
   !!======================================================================
END MODULE flowri
