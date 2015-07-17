MODULE bdyini
   !!======================================================================
   !!                       ***  MODULE  bdyini  ***
   !! Unstructured open boundaries : initialisation
   !!======================================================================
   !! History :  1.0  !  2005-01  (J. Chanut, A. Sellar)  Original code
   !!             -   !  2007-01  (D. Storkey) Update to use IOM module
   !!             -   !  2007-01  (D. Storkey) Tidal forcing
   !!            3.0  !  2008-04  (NEMO team)  add in the reference version
   !!            3.3  !  2010-09  (E.O'Dea) updates for Shelf configurations
   !!            3.3  !  2010-09  (D.Storkey) add ice boundary conditions
   !!----------------------------------------------------------------------
#if defined key_bdy
   !!----------------------------------------------------------------------
   !!   'key_bdy'                     Unstructured Open Boundary Conditions
   !!----------------------------------------------------------------------
   !!   bdy_init       : Initialization of unstructured open boundaries
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain
   USE obc_par         ! ocean open boundary conditions
   USE bdy_oce         ! unstructured open boundary conditions
   USE bdydta, ONLY: bdy_dta_alloc ! open boundary data
   USE bdytides        ! tides at open boundaries initialization (tide_init routine)
   USE in_out_manager  ! I/O units
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp         ! for mpp_sum  
   USE iom             ! I/O

   IMPLICIT NONE
   PRIVATE

   PUBLIC   bdy_init   ! routine called by opa.F90

   !!----------------------------------------------------------------------
   !! NEMO/OPA 4.0 , NEMO Consortium (2011)
   !! $Id: bdyini.F90 2715 2011-03-30 15:58:35Z rblod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS
   
   SUBROUTINE bdy_init
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE bdy_init  ***
      !!         
      !! ** Purpose :   Initialization of the dynamics and tracer fields with 
      !!              unstructured open boundaries.
      !!
      !! ** Method  :   Read initialization arrays (mask, indices) to identify 
      !!              an unstructured open boundary
      !!
      !! ** Input   :  bdy_init.nc, input file for unstructured open boundaries
      !!----------------------------------------------------------------------      
      INTEGER  ::   ii, ij, ik, igrd, ib, ir   ! dummy loop indices
      INTEGER  ::   icount, icountr, ib_len, ibr_max   ! local integers
      INTEGER  ::   iw, ie, is, in, inum, id_dummy     !   -       -
      INTEGER  ::   igrd_start, igrd_end               !   -       -
      REAL(wp) ::   zefl, zwfl, znfl, zsfl              ! local scalars
      INTEGER, DIMENSION (2)             ::   kdimsz
      INTEGER, DIMENSION(jpbdta, jpbgrd) ::   nbidta, nbjdta   ! Index arrays: i and j indices of bdy dta
      INTEGER, DIMENSION(jpbdta, jpbgrd) ::   nbrdta           ! Discrete distance from rim points
      REAL(wp), DIMENSION(jpidta,jpjdta) ::   zmask            ! global domain mask
      REAL(wp), DIMENSION(jpbdta,1)      ::   zdta             ! temporary array 
      CHARACTER(LEN=80),DIMENSION(6)     ::   clfile
      !!
      NAMELIST/nambdy/cn_mask, cn_dta_frs_T, cn_dta_frs_U, cn_dta_frs_V,   &
         &            cn_dta_fla_T, cn_dta_fla_U, cn_dta_fla_V,            &
         &            ln_tides, ln_clim, ln_vol, ln_mask,                  &
         &            ln_dyn_fla, ln_dyn_frs, ln_tra_frs,ln_ice_frs,       &
         &            nn_dtactl, nn_rimwidth, nn_volctl
      !!----------------------------------------------------------------------

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'bdy_init : initialization of unstructured open boundaries'
      IF(lwp) WRITE(numout,*) '~~~~~~~~'
      !
      !                                      ! allocate bdy_oce arrays
      IF( bdy_oce_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'bdy_init : unable to allocate oce arrays' )
      IF( bdy_dta_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'bdy_init : unable to allocate dta arrays' )

      IF( jperio /= 0 )   CALL ctl_stop( 'Cyclic or symmetric,',   &
         &                               ' and unstructured open boundary condition are not compatible' )

      IF( lk_obc      )   CALL ctl_stop( 'Straight open boundaries,',   &
         &                               ' and unstructured open boundaries are not compatible' )

      ! ---------------------------
      REWIND( numnam )                    ! Read namelist parameters
      READ  ( numnam, nambdy )

      !                                   ! control prints
      IF(lwp) WRITE(numout,*) '         nambdy'

      !                                         ! check type of data used (nn_dtactl value)
      IF(lwp) WRITE(numout,*) 'nn_dtactl =', nn_dtactl      
      IF(lwp) WRITE(numout,*)
      SELECT CASE( nn_dtactl )                   ! 
      CASE( 0 )      ;   IF(lwp) WRITE(numout,*) '      initial state used for bdy data'        
      CASE( 1 )      ;   IF(lwp) WRITE(numout,*) '      boundary data taken from file'
      CASE DEFAULT   ;   CALL ctl_stop( 'nn_dtactl must be 0 or 1' )
      END SELECT

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'Boundary rim width for the FRS nn_rimwidth = ', nn_rimwidth

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) '      nn_volctl = ', nn_volctl

      IF( ln_vol ) THEN                     ! check volume conservation (nn_volctl value)
         SELECT CASE ( nn_volctl )
         CASE( 1 )      ;   IF(lwp) WRITE(numout,*) '      The total volume will be constant'
         CASE( 0 )      ;   IF(lwp) WRITE(numout,*) '      The total volume will vary according to the surface E-P flux'
         CASE DEFAULT   ;   CALL ctl_stop( 'nn_volctl must be 0 or 1' )
         END SELECT
         IF(lwp) WRITE(numout,*)
      ELSE
         IF(lwp) WRITE(numout,*) 'No volume correction with unstructured open boundaries'
         IF(lwp) WRITE(numout,*)
      ENDIF

      IF( ln_tides ) THEN
        IF(lwp) WRITE(numout,*) 'Tidal harmonic forcing at unstructured open boundaries'
        IF(lwp) WRITE(numout,*)
      ENDIF

      IF( ln_dyn_fla ) THEN
        IF(lwp) WRITE(numout,*) 'Flather condition on U, V at unstructured open boundaries'
        IF(lwp) WRITE(numout,*)
      ENDIF

      IF( ln_dyn_frs ) THEN
        IF(lwp) WRITE(numout,*) 'FRS condition on U and V at unstructured open boundaries'
        IF(lwp) WRITE(numout,*)
      ENDIF

      IF( ln_tra_frs ) THEN
        IF(lwp) WRITE(numout,*) 'FRS condition on T & S fields at unstructured open boundaries'
        IF(lwp) WRITE(numout,*)
      ENDIF

      IF( ln_ice_frs ) THEN
        IF(lwp) WRITE(numout,*) 'FRS condition on ice fields at unstructured open boundaries'
        IF(lwp) WRITE(numout,*)
      ENDIF

      IF( ln_tides )   CALL tide_init      ! Read tides namelist 


      ! Read arrays defining unstructured open boundaries
      ! -------------------------------------------------

      ! Read global 2D mask at T-points: bdytmask
      ! *****************************************
      ! bdytmask = 1  on the computational domain AND on open boundaries
      !          = 0  elsewhere   
 
      IF( cp_cfg == "eel" .AND. jp_cfg == 5 ) THEN          ! EEL configuration at 5km resolution
         zmask(         :                ,:) = 0.e0
         zmask(jpizoom+1:jpizoom+jpiglo-2,:) = 1.e0          
      ELSE IF( ln_mask ) THEN
         CALL iom_open( cn_mask, inum )
         CALL iom_get ( inum, jpdom_data, 'bdy_msk', zmask(:,:) )
         CALL iom_close( inum )
      ELSE
         zmask(:,:) = 1.e0
      ENDIF

      DO ij = 1, nlcj      ! Save mask over local domain      
         DO ii = 1, nlci
            bdytmask(ii,ij) = zmask( mig(ii), mjg(ij) )
         END DO
      END DO

      ! Derive mask on U and V grid from mask on T grid
      bdyumask(:,:) = 0.e0
      bdyvmask(:,:) = 0.e0
      DO ij=1, jpjm1
         DO ii=1, jpim1
            bdyumask(ii,ij)=bdytmask(ii,ij)*bdytmask(ii+1, ij )
            bdyvmask(ii,ij)=bdytmask(ii,ij)*bdytmask(ii  ,ij+1)  
         END DO
      END DO
      CALL lbc_lnk( bdyumask(:,:), 'U', 1. )   ;   CALL lbc_lnk( bdyvmask(:,:), 'V', 1. )      ! Lateral boundary cond.


      ! Read discrete distance and mapping indices
      ! ******************************************
      nbidta(:,:) = 0.e0
      nbjdta(:,:) = 0.e0
      nbrdta(:,:) = 0.e0

      IF( cp_cfg == "eel" .AND. jp_cfg == 5 ) THEN
         icount = 0
         DO ir = 1, nn_rimwidth                  ! Define west boundary (from ii=2 to ii=1+nn_rimwidth):
            DO ij = 3, jpjglo-2
               icount = icount + 1
               nbidta(icount,:) = ir + 1 + (jpizoom-1)
               nbjdta(icount,:) = ij     + (jpjzoom-1) 
               nbrdta(icount,:) = ir
            END DO
         END DO
         !
         DO ir = 1, nn_rimwidth                  ! Define east boundary (from ii=jpiglo-1 to ii=jpiglo-nn_rimwidth):
            DO ij=3,jpjglo-2
               icount = icount + 1
               nbidta(icount,:) = jpiglo-ir + (jpizoom-1)
               nbidta(icount,2) = jpiglo-ir-1 + (jpizoom-1) ! special case for u points
               nbjdta(icount,:) = ij + (jpjzoom-1)
               nbrdta(icount,:) = ir
            END DO
         END DO
         !       
      ELSE            ! Read indices and distances in unstructured boundary data files 
         !
         IF( ln_tides ) THEN             ! Read tides input files for preference in case there are no bdydata files
            clfile(4) = TRIM(filtide)//TRIM(tide_cpt(1))//'_grid_T.nc'
            clfile(5) = TRIM(filtide)//TRIM(tide_cpt(1))//'_grid_U.nc'
            clfile(6) = TRIM(filtide)//TRIM(tide_cpt(1))//'_grid_V.nc'
         ENDIF
         IF( ln_dyn_fla .AND. .NOT. ln_tides ) THEN 
            clfile(4) = cn_dta_fla_T
            clfile(5) = cn_dta_fla_U
            clfile(6) = cn_dta_fla_V
         ENDIF

         IF( ln_tra_frs ) THEN 
            clfile(1) = cn_dta_frs_T
            IF( .NOT. ln_dyn_frs ) THEN 
               clfile(2) = cn_dta_frs_T     ! Dummy read re read T file for sake of 6 files
               clfile(3) = cn_dta_frs_T     !
            ENDIF
         ENDIF          
         IF( ln_dyn_frs ) THEN 
            IF( .NOT. ln_tra_frs )   clfile(1) = cn_dta_frs_U      ! Dummy Read 
            clfile(2) = cn_dta_frs_U
            clfile(3) = cn_dta_frs_V 
         ENDIF

         !                                   ! how many files are we to read in?
         IF(ln_tides .OR. ln_dyn_fla)   igrd_start = 4
         !
         IF(ln_tra_frs    ) THEN   ;   igrd_start = 1
         ELSEIF(ln_dyn_frs) THEN   ;   igrd_start = 2
         ENDIF
         !
         IF( ln_tra_frs   )   igrd_end = 1
         !
         IF(ln_dyn_fla .OR. ln_tides) THEN   ;   igrd_end = 6
         ELSEIF( ln_dyn_frs             ) THEN   ;   igrd_end = 3
         ENDIF

         DO igrd = igrd_start, igrd_end
            CALL iom_open( clfile(igrd), inum )
            id_dummy = iom_varid( inum, 'nbidta', kdimsz=kdimsz )  
            IF(lwp) WRITE(numout,*) 'kdimsz : ',kdimsz
            ib_len = kdimsz(1)
            IF( ib_len > jpbdta)   CALL ctl_stop(  'Boundary data array in file too long.',                  &
                &                                  'File :', TRIM(clfile(igrd)),'increase parameter jpbdta.' )

            CALL iom_get( inum, jpdom_unknown, 'nbidta', zdta(1:ib_len,:) )
            DO ii = 1,ib_len
               nbidta(ii,igrd) = INT( zdta(ii,1) )
            END DO
            CALL iom_get( inum, jpdom_unknown, 'nbjdta', zdta(1:ib_len,:) )
            DO ii = 1,ib_len
               nbjdta(ii,igrd) = INT( zdta(ii,1) )
            END DO
            CALL iom_get( inum, jpdom_unknown, 'nbrdta', zdta(1:ib_len,:) )
            DO ii = 1,ib_len
               nbrdta(ii,igrd) = INT( zdta(ii,1) )
            END DO
            CALL iom_close( inum )

            IF( igrd < 4) THEN            ! Check that rimwidth in file is big enough for Frs case(barotropic is one):
               ibr_max = MAXVAL( nbrdta(:,igrd) )
               IF(lwp) WRITE(numout,*)
               IF(lwp) WRITE(numout,*) ' Maximum rimwidth in file is ', ibr_max
               IF(lwp) WRITE(numout,*) ' nn_rimwidth from namelist is ', nn_rimwidth
               IF (ibr_max < nn_rimwidth)   CALL ctl_stop( 'nn_rimwidth is larger than maximum rimwidth in file' )
            ENDIF !Check igrd < 4
            !
         END DO
         !
      ENDIF 

      ! Dispatch mapping indices and discrete distances on each processor
      ! *****************************************************************
     
      iw = mig(1) + 1            ! if monotasking and no zoom, iw=2
      ie = mig(1) + nlci-1 - 1   ! if monotasking and no zoom, ie=jpim1
      is = mjg(1) + 1            ! if monotasking and no zoom, is=2
      in = mjg(1) + nlcj-1 - 1   ! if monotasking and no zoom, in=jpjm1

      DO igrd = igrd_start, igrd_end
         icount  = 0
         icountr = 0
         nblen   (igrd) = 0
         nblenrim(igrd) = 0
         nblendta(igrd) = 0
         DO ir=1, nn_rimwidth
            DO ib = 1, jpbdta
               ! check if point is in local domain and equals ir
               IF(  nbidta(ib,igrd) >= iw .AND. nbidta(ib,igrd) <= ie .AND.   &
                  & nbjdta(ib,igrd) >= is .AND. nbjdta(ib,igrd) <= in .AND.   &
                  & nbrdta(ib,igrd) == ir  ) THEN
                  !
                  icount = icount  + 1
                  !
                  IF( ir == 1 )   icountr = icountr+1
                  IF (icount > jpbdim) THEN
                     IF(lwp) WRITE(numout,*) 'bdy_ini: jpbdim too small'
                     nstop = nstop + 1
                  ELSE
                     nbi(icount, igrd)  = nbidta(ib,igrd)- mig(1)+1
                     nbj(icount, igrd)  = nbjdta(ib,igrd)- mjg(1)+1
                     nbr(icount, igrd)  = nbrdta(ib,igrd)
                     nbmap(icount,igrd) = ib
                  ENDIF            
               ENDIF
            END DO
         END DO
         nblenrim(igrd) = icountr !: length of rim boundary data on each proc
         nblen   (igrd) = icount  !: length of boundary data on each proc        
      END DO 

      ! Compute rim weights
      ! -------------------
      DO igrd = igrd_start, igrd_end
         DO ib = 1, nblen(igrd)
            nbw(ib,igrd) = 1.- TANH( FLOAT( nbr(ib,igrd) - 1 ) *0.5 )                     ! tanh formulation
!           nbw(ib,igrd) = (FLOAT(nn_rimwidth+1-nbr(ib,igrd))/FLOAT(nn_rimwidth))**2      ! quadratic
!           nbw(ib,igrd) =  FLOAT(nn_rimwidth+1-nbr(ib,igrd))/FLOAT(nn_rimwidth)          ! linear
         END DO
      END DO 
   
      ! Mask corrections
      ! ----------------
      DO ik = 1, jpkm1
         DO ij = 1, jpj
            DO ii = 1, jpi
               tmask(ii,ij,ik) = tmask(ii,ij,ik) * bdytmask(ii,ij)
               umask(ii,ij,ik) = umask(ii,ij,ik) * bdyumask(ii,ij)
               vmask(ii,ij,ik) = vmask(ii,ij,ik) * bdyvmask(ii,ij)
               bmask(ii,ij)    = bmask(ii,ij)    * bdytmask(ii,ij)
            END DO      
         END DO
      END DO

      DO ik = 1, jpkm1
         DO ij = 2, jpjm1
            DO ii = 2, jpim1
               fmask(ii,ij,ik) = fmask(ii,ij,ik) * bdytmask(ii,ij  ) * bdytmask(ii+1,ij  )   &
                  &                              * bdytmask(ii,ij+1) * bdytmask(ii+1,ij+1)
            END DO      
         END DO
      END DO

      tmask_i (:,:) = tmask(:,:,1) * tmask_i(:,:)             
      bdytmask(:,:) = tmask(:,:,1)

      ! bdy masks and bmask are now set to zero on boundary points:
      igrd = 1       ! In the free surface case, bmask is at T-points
      DO ib = 1, nblenrim(igrd)     
        bmask(nbi(ib,igrd), nbj(ib,igrd)) = 0.e0
      END DO
      !
      igrd = 1
      DO ib = 1, nblenrim(igrd)      
        bdytmask(nbi(ib,igrd), nbj(ib,igrd)) = 0.e0
      END DO
      !
      igrd = 2
      DO ib = 1, nblenrim(igrd)
        bdyumask(nbi(ib,igrd), nbj(ib,igrd)) = 0.e0
      END DO
      !
      igrd = 3
      DO ib = 1, nblenrim(igrd)
        bdyvmask(nbi(ib,igrd), nbj(ib,igrd)) = 0.e0
      END DO

      ! Lateral boundary conditions
      CALL lbc_lnk( fmask        , 'F', 1. )   ;   CALL lbc_lnk( bdytmask(:,:), 'T', 1. )
      CALL lbc_lnk( bdyumask(:,:), 'U', 1. )   ;   CALL lbc_lnk( bdyvmask(:,:), 'V', 1. )

      IF( ln_vol .OR. ln_dyn_fla ) THEN      ! Indices and directions of rim velocity components
         !
         !flagu = -1 : u component is normal to the dynamical boundary but its direction is outward
         !flagu =  0 : u is tangential
         !flagu =  1 : u is normal to the boundary and is direction is inward
         icount = 0 
         flagu(:) = 0.e0
 
         igrd = 2      ! u-component 
         DO ib = 1, nblenrim(igrd)  
            zefl=bdytmask(nbi(ib,igrd)  , nbj(ib,igrd))
            zwfl=bdytmask(nbi(ib,igrd)+1, nbj(ib,igrd))
            IF( zefl + zwfl ==2 ) THEN
               icount = icount +1
            ELSE
               flagu(ib)=-zefl+zwfl
            ENDIF
         END DO

         !flagv = -1 : u component is normal to the dynamical boundary but its direction is outward
         !flagv =  0 : u is tangential
         !flagv =  1 : u is normal to the boundary and is direction is inward
         flagv(:) = 0.e0

         igrd = 3      ! v-component
         DO ib = 1, nblenrim(igrd)  
            znfl = bdytmask(nbi(ib,igrd), nbj(ib,igrd))
            zsfl = bdytmask(nbi(ib,igrd), nbj(ib,igrd)+1)
            IF( znfl + zsfl ==2 ) THEN
               icount = icount + 1
            ELSE
               flagv(ib) = -znfl + zsfl
            END IF
         END DO
 
         IF( icount /= 0 ) THEN
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) ' E R R O R : Some data velocity points,',   &
               ' are not boundary points. Check nbi, nbj, indices.'
            IF(lwp) WRITE(numout,*) ' ========== '
            IF(lwp) WRITE(numout,*)
            nstop = nstop + 1
         ENDIF 
    
      ENDIF

      ! Compute total lateral surface for volume correction:
      ! ----------------------------------------------------
      bdysurftot = 0.e0 
      IF( ln_vol ) THEN  
         igrd = 2      ! Lateral surface at U-points
         DO ib = 1, nblenrim(igrd)
            bdysurftot = bdysurftot + hu     (nbi(ib,igrd)  ,nbj(ib,igrd))                      &
               &                    * e2u    (nbi(ib,igrd)  ,nbj(ib,igrd)) * ABS( flagu(ib) )   &
               &                    * tmask_i(nbi(ib,igrd)  ,nbj(ib,igrd))                      &
               &                    * tmask_i(nbi(ib,igrd)+1,nbj(ib,igrd))                   
         END DO

         igrd=3 ! Add lateral surface at V-points
         DO ib = 1, nblenrim(igrd)
            bdysurftot = bdysurftot + hv     (nbi(ib,igrd),nbj(ib,igrd)  )                      &
               &                    * e1v    (nbi(ib,igrd),nbj(ib,igrd)  ) * ABS( flagv(ib) )   &
               &                    * tmask_i(nbi(ib,igrd),nbj(ib,igrd)  )                      &
               &                    * tmask_i(nbi(ib,igrd),nbj(ib,igrd)+1)
         END DO
         !
         IF( lk_mpp )   CALL mpp_sum( bdysurftot )      ! sum over the global domain
      END IF   

      ! Initialise bdy data arrays
      ! --------------------------
      tbdy(:,:) = 0.e0
      sbdy(:,:) = 0.e0
      ubdy(:,:) = 0.e0
      vbdy(:,:) = 0.e0
      sshbdy(:) = 0.e0
      ubtbdy(:) = 0.e0
      vbtbdy(:) = 0.e0
#if defined key_lim2
      frld_bdy(:) = 0.e0
      hicif_bdy(:) = 0.e0
      hsnif_bdy(:) = 0.e0
#endif

      ! Read in tidal constituents and adjust for model start time
      ! ----------------------------------------------------------
      IF( ln_tides )   CALL tide_data
      !
   END SUBROUTINE bdy_init

#else
   !!---------------------------------------------------------------------------------
   !!   Dummy module                                   NO unstructured open boundaries
   !!---------------------------------------------------------------------------------
CONTAINS
   SUBROUTINE bdy_init      ! Dummy routine
   END SUBROUTINE bdy_init
#endif

   !!=================================================================================
END MODULE bdyini
