PROGRAM mpp_domain_decomposition
 !!---------------------------------------------------------------------
 !!
 !!                       PROGRAM MPP_OPTIMIZ_NC
 !!                     ***********************
 !!
 !!  PURPOSE :
 !!  ---------
 !!              This program is build to optimize the domain beakdown into
 !!              subdomain for mpp computing.
 !!              Once the grid size, and the land/sea mask is known, it looks
 !!              for all the possibilities within a range of setting parameters
 !!              and determine the optimal.
 !!
 !!              Optimization is done with respect to the maximum number of
 !!              sea processors and to the maximum numbers of procs (jprocx)
 !!                     
 !! history:
 !! --------
 !!       original  : 95-12 (Imbard M) for OPA8.1, CLIPPER
 !!       f90       : 03-06 (Molines JM), namelist as input
 !!                 : 05-05 (Molines JM), bathy in ncdf
 !!                 : 2015  (Lovato T), revise structure of the tool
 !!                               and correct domain assignment
 !!----------------------------------------------------------------------
 !! * modules used
  USE netcdf

  IMPLICIT NONE

  INTEGER, PARAMETER :: dp = selected_real_kind(15,307) ! REAL*8
  INTEGER, PARAMETER :: sp = selected_real_kind(6,37) ! REAL*4
  INTEGER, PARAMETER :: wp = dp 
  INTEGER ::  jprocmin=1, jprocmax=1   !: min and max number of proc. to test along I and J
  !
  INTEGER      ::  &
       jpiglo ,    & !: I-size of the model (namelist)
       jpjglo ,    & !: J-size of the model (namelist)
       jpkglo        !: vertical levels (namelist)
  !
  INTEGER,PARAMETER :: jpreci=1, &         !: number of columns for overlap 
                &      jprecj=1            !: number of rows for overlap
  !
  CHARACTER(LEN=120) :: gridname,      &   !: Name of the input grid
                &      cbathy, clvar       !: File name of the netcdf bathymetry (namelist)
  LOGICAL      :: ln_mesh                  !: Logical flag if input is bathymetry or mesh mask.
  LOGICAL      :: ln_order                 !: Logical flag to output PEs combinantion by ocean domains
  !
  LOGICAL      :: ln_zcord                 !: use to mask minimum depth if z-coordinate is used in the configuration
  REAL(dp)     :: rn_hmin                  !: minimum depth threshold
  !
  ! OTHERS
  INTEGER :: jpnix ,jpnjx  
  INTEGER :: ji,jj,jn,jni,jnj,jni2,jnj2
  INTEGER :: ii,iim,ij,ijm,iwet,iost,iresti,irestj,isurf,iempty
  INTEGER :: iilb,ijlb,ireci,irecj,in
  INTEGER :: ipi,ipj
  INTEGER :: iii,iij
  !
  INTEGER :: totocepts, iocepts, iiocepts, itotocepts
  INTEGER :: jproc2, icnt, icnt2
  INTEGER,DIMENSION(:),ALLOCATABLE     :: jpnia, jpnja, jpsea, jpoce
  REAL(dp),DIMENSION(:),ALLOCATABLE    :: rto, rtoce
  !
  INTEGER,DIMENSION(:,:),ALLOCATABLE   ::  ibathy, ocepts
  INTEGER,DIMENSION(:,:),ALLOCATABLE   ::  ippdi, ippdj ,iidom, ijdom
  INTEGER,DIMENSION(:,:,:),ALLOCATABLE ::  tmask
  !
  REAL(dp)                             ::  rt, zper
  REAL(dp),DIMENSION(:,:),ALLOCATABLE  ::  zmask
  REAL(sp),DIMENSION(:,:),ALLOCATABLE  ::  zdta
  !
  INTEGER                              :: numnam, numout
  ! CDF stuff
  INTEGER            :: ncid, ivarid
  LOGICAL            :: llok
  CHARACTER(LEN=120) :: outfilename
  !
  ! Namelist
  NAMELIST /nam_mpprep/ jprocmin, jprocmax, jpiglo, jpjglo, jpkglo, &
           & gridname, cbathy, clvar, ln_mesh, ln_zcord, rn_hmin, ln_order
  !
  ! 0. Initialisation
  ! -----------------
  ! I/O channels
  numout = 1
  numnam = 4
  llok = .false.

  ! Read namelist content
  OPEN(numnam,FILE='namelist_mpp')
  REWIND(numnam)
  READ(numnam,nam_mpprep)
  !
  ! Check input values for testing decomposition  
  IF ( jprocmin <= 0 ) jprocmin = 1 
  jpnix = jprocmax ; jpnjx = jprocmax
  jproc2 = jprocmax * jprocmax

  ALLOCATE ( ibathy(jpiglo,jpjglo), zmask(jpiglo,jpjglo), zdta(jpiglo,jpjglo) )
  ALLOCATE ( ippdi(jpnix,jpnjx), ippdj(jpnix,jpnjx) )
  ALLOCATE ( iidom(jpnix,jpnjx), ijdom(jpnix,jpnjx) )
  ALLOCATE ( tmask(jpiglo,jpjglo,jpkglo) )
  ALLOCATE ( ocepts(jpiglo,jpjglo) )
  ALLOCATE ( jpnia(jproc2), jpnja(jproc2), jpsea(jproc2), jpoce(jproc2) )
  ALLOCATE ( rto(jproc2), rtoce(jproc2) )
  !
  ! Display analysis information
  WRITE(*,'(3a)')           'Begin processing of PEs decomposition for ',TRIM(gridname),' grid'
  WRITE(*,'(a,i5,a,i5)')    'Test processors in the range : ',jprocmin,' - ', jprocmax
  WRITE(*,'(4a)')           'Input field is ',TRIM(clvar),' from file ',TRIM(cbathy)
  WRITE(*,*)

  IF ( ln_zcord .AND. .NOT. ln_mesh ) THEN
    WRITE(*,'(a)')         'Vertical grid discretization is in z-coordinates (either zco or zps )'
    WRITE(*,'(a,f6.2,a)')  'Apply mininum value to bathymetry. rn_hmin = ',rn_hmin ,' meters'
  ELSEIF ( ln_mesh ) THEN
    WRITE(*,'(a)')         'Input field is the mesh_mask. Use tmask variable by default'
    WRITE(*,'(a)')         'Find also the PE with the largest value ocean points in 3D'
  ENDIF

  ! Set output filename
  WRITE(outfilename,'(a)') TRIM(gridname)//"_PEs_decomposition.layout" 
  !
  ! --------------------------------
  ! Read input data
  ! --------------------------------
  ! Force to read tmask if input meshmask is selected
  IF ( ln_mesh ) clvar = 'tmask' 
  !
  INQUIRE( FILE=cbathy, EXIST=llok )
  !
  IF( llok ) THEN
        call check ( NF90_OPEN(cbathy,NF90_NOWRITE,ncid) )
        call check ( NF90_INQ_VARID(ncid,clvar,ivarid) )
        IF ( ln_mesh ) THEN
          call check ( NF90_GET_VAR(ncid,ivarid,tmask) )
          ibathy(:,:) = tmask(:,:,1)
        ELSE
          call check ( NF90_GET_VAR(ncid,ivarid,zdta) )
          IF ( ln_zcord ) THEN
            WHERE (zdta <= 0.0 )   ; zdta = 0.0
            ELSEWHERE              ; zdta = MAX(  rn_hmin , zdta(:,:)  )
            END WHERE
          ENDIF
          ibathy(:,:) = int ( zdta(:,:) )
        ENDIF
        call check ( NF90_CLOSE(ncid) )
  ELSE
      PRINT *,' File missing : ', trim(cbathy)
      STOP
  ENDIF
  !
  WRITE(*,'(a)') 'Input data loaded ...'
  !
  ! --------------------------------
  ! Building the mask
  ! --------------------------------
  WHERE (ibathy <= 0 ) ; zmask = 0._wp
  ELSEWHERE            ; zmask = 1._wp
  END WHERE
  !
  ! Build colum wise total of ocean points
  IF (ln_mesh) THEN 
    totocepts = 0
    DO jj = 1 , jpjglo
       DO ji = 1 , jpiglo
          iocepts=0
          IF ( tmask(ji,jj,1) .GT. 0 ) iocepts = SUM( tmask(ji,jj,:) )
          ocepts(ji,jj) = float(iocepts)
          totocepts = totocepts + iocepts
       ENDDO
    ENDDO
  ENDIF
  !
  ! --------------------------------
  !  1. Main loop on PEs
  ! --------------------------------
  !
  WRITE(*,'(a)') 'Start loop on PEs decomposition'
  ! 
  icnt = 0
  DO jni = jprocmin , jprocmax
     DO jnj = jprocmin , jprocmax
        !
        ! Limitation on the maximum/minimum number of PE's
        IF(jni*jnj >  jprocmax) goto 1000
        IF(jni*jnj <  jprocmin) goto 1000
        !
        ! Partition
        ipi=(jpiglo-2*jpreci + (jni-1))/jni + 2*jpreci
        ipj=(jpjglo-2*jprecj + (jnj-1))/jnj + 2*jprecj
        !
        ! Ratio jpnij*domain/global domain
        zper=(jni*jnj*ipi*ipj)/float(jpiglo*jpjglo)
        !
        ! Bottom Left corner of each processor
        !
        iilb=1
        ijlb=1
        ireci=2*jpreci
        irecj=2*jprecj
        iresti = MOD ( jpiglo - ireci , jni )
        irestj = MOD ( jpjglo - irecj , jnj )
        !
        ! single sub-domain dimension along i-direction
        IF (iresti.EQ.0) iresti = jni
        DO jj=1,jnj
           DO ji=1,iresti
              ippdi(ji,jj) = ipi
           END DO
           DO ji=iresti+1,jni
              ippdi(ji,jj) = ipi -1
           END DO
        END DO
        ! single sub-domain dimension along j-direction
        IF (irestj.EQ.0) irestj = jnj
        DO ji=1,jni
           DO jj=1,irestj
              ippdj(ji,jj) = ipj
           END DO
           DO jj=irestj+1,jnj
              ippdj(ji,jj) = ipj -1
           END DO
        END DO
        DO jj=1,jnj
           DO ji=1,jni
              iidom(ji,jj)=iilb
              ijdom(ji,jj)=ijlb
           END DO
        END DO
        !
        ! --------------------------------
        !  2. Loop this decomposition PEs
        ! --------------------------------
        !
        iempty=0
        iwet=0
        iiocepts = 0
        itotocepts = 0
        icnt2 = 0
        !
        DO jnj2=1,jnj
          DO jni2=1,jni
              icnt2 = icnt2 + 1
              ! Global domain index along i-direction
              IF(jni.GT.1)THEN
                 DO jj=1,jnj
                    DO ji=2,jni
                       iidom(ji,jj)=iidom(ji-1,jj)+ippdi(ji-1,jj)-ireci
                    END DO
                 END DO
                 iilb=iidom(jni2,jnj2)
              ENDIF
              ! Global domain index along j-direction
              IF(jnj.GT.1)THEN
                 DO jj=2,jnj
                    DO ji=1,jni
                       ijdom(ji,jj)=ijdom(ji,jj-1)+ippdj(ji,jj-1)-irecj
                    END DO
                 END DO
                 ijlb=ijdom(jni2,jnj2)
              ENDIF

              ! Check wet points over the entire domain to preserve the MPI communication stencil
              isurf=0
              DO jj = 1, ippdj(jni2,jnj2)
                 DO  ji = 1, ippdi(jni2,jnj2)
                    IF ( zmask(ji+iilb-1,jj+ijlb-1) .EQ. 1. ) isurf = isurf + 1
                 END DO
              END DO

              IF(isurf.EQ.0) THEN
                 iempty=iempty+1
              ELSE
                 iwet=iwet+isurf
              ENDIF
              
              !
              ! Store 3D ocean points informations
              IF ( ln_mesh ) then
                iocepts=0
                IF (isurf.NE.0) THEN
                   DO jj=1+jprecj-1,ippdj(jni2,jnj2)-jprecj+1
                     DO  ji=1+jpreci-1,ippdi(jni2,jnj2)-jpreci+1
                       IF (ji.NE.0 .or. ji.NE.jpiglo+1 ) then
                         IF (jj.NE.0 .or. jj.NE.jpjglo+1 ) then
                           IF(ocepts(ji+iilb-1,jj+ijlb-1).GT.0.) iocepts=iocepts+int(ocepts(ji+iilb-1,jj+ijlb-1))
                         ENDIF
                       ENDIF
                     END DO
                   END DO
                   ! Store memory and ocepoints information
                   iiocepts = max (iiocepts , iocepts) 
                   itotocepts = itotocepts + iocepts
                ENDIF
              ENDIF

           END DO
        END DO
        !
        ! --------------------------------
        !  3. Store decomposition PEs
        ! --------------------------------
        icnt = icnt + 1
        jpnia(icnt) = jni
        jpnja(icnt) = jnj
        jpsea(icnt) = jni*jnj-iempty
        !rto(icnt) = real(jni*jnj-iempty)/real(jni*jnj)
        rto(icnt) = float((jni*jnj-iempty)) * float(ipi*ipj) / float(jpiglo*jpjglo)

        IF ( ln_mesh ) THEN
           jpoce(icnt) = float(iiocepts)
           rtoce(icnt) = float(itotocepts)/totocepts*100
        ENDIF
      
      1000 continue

     END DO
     WRITE(*,*) ' Processed row for procx = ',jni
  END DO
  !
  ! --------------------------------
  !  4. Reorder decomposition data
  ! --------------------------------
  !
  IF ( ln_order ) then
  !
  WRITE(*,*)
  WRITE(*,*) ' Reorder PEs decompositions by total ocean PEs and decomposed-to-global ratio'
  WRITE(*,*)
  !
  ! Rank the combinations by ocean domain only
  DO ji = 1, icnt-1
     DO jj = ji+1, icnt
        !if(abs(jpsea(ji)-jpnjx).gt.abs(jpsea(jj)-jpnjx))then
        IF ( jpsea(ji) .GT. jpsea(jj) ) THEN
           jn = jpnia(ji)
           jpnia(ji) = jpnia(jj)
           jpnia(jj) = jn
           jn = jpnja(ji)
           jpnja(ji) = jpnja(jj)
           jpnja(jj) = jn
           jn = jpsea(ji)
           jpsea(ji) = jpsea(jj)
           jpsea(jj) = jn
           rt = rto(ji)
           rto(ji) = rto(jj)
           rto(jj) = rt

           IF ( ln_mesh ) THEN
              jn = jpoce(ji)
              jpoce(ji) = jpoce(jj)
              jpoce(jj) = jn
              rt = rtoce(ji)
              rtoce(ji) = rtoce(jj)
              rtoce(jj) = rt      
           ENDIF
        ENDIF
     ENDDO
  ENDDO
  ! 
  ! Rank each PEs combination by 2D decomposed/original ratios
   DO ji = 1, icnt-1
      DO jj = ji+1, icnt
         IF( jpsea(ji) .EQ. jpsea(jj) .AND. rto(ji) .GT. rto(jj) ) THEN
            jn = jpnia(ji)
            jpnia(ji) = jpnia(jj)
            jpnia(jj) = jn
            jn = jpnja(ji)
            jpnja(ji) = jpnja(jj)
            jpnja(jj) = jn
            jn = jpsea(ji)
            jpsea(ji) = jpsea(jj)
            jpsea(jj) = jn
            rt = rto(ji)
            rto(ji) = rto(jj)
            rto(jj) = rt
            IF ( ln_mesh ) THEN
               jn = jpoce(ji)
               jpoce(ji) = jpoce(jj)
               jpoce(jj) = jn
               rt = rtoce(ji)
               rtoce(ji) = rtoce(jj)
               rtoce(jj) = rt
            ENDIF
         ENDIF
      ENDDO
   ENDDO
  !
  ENDIF
  !
  ! --------------------------------
  !  5. Save data to output file
  ! --------------------------------
  !
  OPEN(numout,FILE=TRIM(outfilename))
  WRITE(numout,*) 
  WRITE(numout,'(a)') TRIM(gridname)//' - NEMO Decompositon layout (listed by increasing number of ocean only PEs'
  WRITE(numout,*) 

  jn = 0
  DO ji = 1 , icnt
    IF ( ( jpnia(ji) .NE. 0 ) .AND. ( jpnja(ji) .NE. 0) ) THEN
       jn = jn + 1 
       IF ( ln_mesh ) THEN
         IF (jn .EQ. 1 ) WRITE(numout,'(a)') ' index | iproc | jproc | ixj(tot) | ixj(oce) | sea/total |  MAX # Ocepts | 3D Decomp/Glob '
         WRITE(numout,"(2(i6,2x),i6,1x,2(i10,1x),f10.2,'%',1x,i12,2x,f14.2,'%')") jn,jpnia(ji),jpnja(ji),jpnia(ji)*jpnja(ji),jpsea(ji),rto(ji) *100,jpoce(ji),rtoce(ji)
       ELSE
         IF (jn .EQ. 1 ) WRITE(numout,'(a)') ' index | iproc | jproc | ixj(tot) | ixj(oce) | sea/total '
         WRITE(numout,"(2(i6,2x),i6,1x,2(i10,1x),f10.2,'%')") jn,jpnia(ji),jpnja(ji),jpnia(ji)*jpnja(ji),jpsea(ji),rto(ji) *100
       ENDIF
    ENDIF
  ENDDO
  !
  WRITE(*,*) ' Decompositon analysis finished.'
  !
  STOP
  ! 
CONTAINS
  !
  ! Subroutine to handle NetCDF error message
  subroutine check(status)
    integer, intent ( in) :: status
  
    if(status /= nf90_noerr) then
      write(*,*) ' Input file error - ', trim(nf90_strerror(status))
      stop 2
    end if
  end subroutine check
  !
END PROGRAM mpp_domain_decomposition
