MODULE boxed_domain

  INTEGER,SAVE :: ni_glo
  INTEGER,SAVE :: nj_glo
  INTEGER,SAVE :: nb_proc
  
  TYPE local_domain
    INTEGER :: ni
    INTEGER :: nj
    INTEGER :: ibegin
    INTEGER :: iend
    INTEGER :: jbegin
    INTEGER :: jend
    REAL,POINTER :: lon(:,:)
    REAL,POINTER :: lat(:,:)
  END TYPE local_domain
  
  TYPE(local_domain),ALLOCATABLE,SAVE,TARGET :: domains(:)
  
  
  INTEGER,SAVE :: ni
  INTEGER,SAVE :: nj
  INTEGER,SAVE :: ibegin
  INTEGER,SAVE :: iend
  INTEGER,SAVE :: jbegin
  INTEGER,SAVE :: jend
  
  REAL,SAVE,ALLOCATABLE :: lon(:,:)
  REAL,SAVE,ALLOCATABLE :: lat(:,:)
  
  INTEGER,SAVE :: max_vert_size = 100
  REAL,SAVE,ALLOCATABLE :: Field_buffer(:,:,:)  
  
CONTAINS

  SUBROUTINE Init_global_domain(ni_glo_,nj_glo_,nb_proc_)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ni_glo_
  INTEGER, INTENT(IN) :: nj_glo_
  INTEGER, INTENT(IN) :: nb_proc_
  
    ni_glo=ni_glo_
    nj_glo=nj_glo_
    nb_proc=nb_proc_
    
    ALLOCATE(domains(nb_proc))
 
  END SUBROUTINE Init_global_domain


  SUBROUTINE Init_local_domain(rank,ni,nj,ibegin,jbegin,lon,lat)
  IMPLICIT NONE
    INTEGER, INTENT(IN) :: rank
    INTEGER, INTENT(IN) :: ni
    INTEGER, INTENT(IN) :: nj
    INTEGER, INTENT(IN) :: ibegin    
    INTEGER, INTENT(IN) :: jbegin
    REAL, INTENT(IN)    :: lon(ni,nj)
    REAL, INTENT(IN)    :: lat(ni,nj)
    
    domains(rank)%ni=ni
    domains(rank)%nj=nj
    domains(rank)%ibegin=ibegin
    domains(rank)%jbegin=jbegin   
    domains(rank)%iend=ibegin+ni-1
    domains(rank)%jend=jbegin+nj-1
    ALLOCATE(domains(rank)%lon(ni,nj))
    ALLOCATE(domains(rank)%lat(ni,nj))
    domains(rank)%lon(:,:)=lon(:,:)
    domains(rank)%lat(:,:)=lat(:,:)
  
  END SUBROUTINE Init_local_domain
  
  
  SUBROUTINE Init_domain
  IMPLICIT NONE
    INTEGER :: rank
    INTEGER :: ib,ie,jb,je
    ibegin=domains(1)%ibegin
    iend=domains(1)%iend
    jbegin=domains(1)%jbegin
    jend=domains(1)%jend
    
    DO rank=1,nb_proc
      IF (domains(rank)%ibegin<ibegin) ibegin=domains(rank)%ibegin
      IF (domains(rank)%iend>iend) iend=domains(rank)%iend
      IF (domains(rank)%jbegin<jbegin) jbegin=domains(rank)%jbegin
      IF (domains(rank)%jend>jend) jend=domains(rank)%jend
    ENDDO
    
    ni=iend-ibegin+1
    nj=jend-jbegin+1
    ALLOCATE(lon(ni,nj))
    ALLOCATE(lat(ni,nj))
    
    DO rank=1,nb_proc
      ib=domains(rank)%ibegin-ibegin+1
      ie=domains(rank)%iend-ibegin+1
      jb=domains(rank)%jbegin-jbegin+1
      je=domains(rank)%jend-jbegin+1

      domains(rank)%ibegin=ib
      domains(rank)%iend=ie
      domains(rank)%jbegin=jb
      domains(rank)%jend=je
      
      lon(ib:ie,jb:je)=domains(rank)%lon(:,:)
      lat(ib:ie,jb:je)=domains(rank)%lat(:,:)
    ENDDO
  
    ALLOCATE(field_buffer(ni,nj,max_vert_size))
    
  END SUBROUTINE 
  
  SUBROUTINE bufferize_field2d(field,rank)
  IMPLICIT NONE
    REAL,INTENT(IN) :: field(:,:)
    INTEGER,INTENT(IN) :: rank
    TYPE(local_domain), POINTER :: Pt

    Pt=>domains(rank)
    field_buffer(Pt%ibegin:Pt%iend,Pt%jbegin:Pt%jend,1)=field(:,:)
      
  END SUBROUTINE bufferize_field2d
  
  
  SUBROUTINE bufferize_field3d(field,rank)
  IMPLICIT NONE
    REAL,INTENT(IN) :: field(:,:,:)
    INTEGER,INTENT(IN) :: rank
    TYPE(local_domain), POINTER :: Pt
   
    Pt=>domains(rank)
    IF (size(field,3)>max_vert_size) THEN
      DEALLOCATE(field_buffer)
      max_vert_size=max_vert_size*2
      ALLOCATE(field_buffer(ni,nj,max_vert_size))
    ENDIF
    
    field_buffer(Pt%ibegin:Pt%iend,Pt%jbegin:Pt%jend,1:size(field,3))=field(:,:,:)
      
  END SUBROUTINE bufferize_field3d

    
END MODULE boxed_domain
