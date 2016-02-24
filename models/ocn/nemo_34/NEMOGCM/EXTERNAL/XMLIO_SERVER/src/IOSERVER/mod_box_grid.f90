MODULE box_grid
  
  INTEGER, PARAMETER :: nb_grid_max=100
  INTEGER            :: nb_grid=0

  TYPE domain
    INTEGER :: ni
    INTEGER :: nj
    INTEGER :: ibegin
    INTEGER :: iend
    INTEGER :: jbegin
    INTEGER :: jend
    REAL,POINTER :: lon(:,:)
    REAL,POINTER :: lat(:,:)
  END TYPE domain
  
  TYPE grid
    TYPE(domain)         :: local_domain
    TYPE(domain),POINTER :: task_domain(:)
    INTEGER                          :: ni
    INTEGER                          :: nj
    INTEGER                          :: nb_task
  END TYPE grid
  
  TYPE(grid),SAVE,TARGET :: grids(nb_grid_max)
  
  INTEGER,SAVE :: max_ni
  INTEGER,SAVE :: max_nj
  INTEGER,SAVE :: max_vert_size = 100
  REAL,SAVE,ALLOCATABLE :: Field_buffer(:,:,:) 
  
CONTAINS

  SUBROUTINE create_new_grid(ni_glo,nj_glo,nb_task,id)
  IMPLICIT NONE
    INTEGER,INTENT(IN)  :: ni_glo
    INTEGER,INTENT(IN)  :: nj_glo
    INTEGER,INTENT(IN)  :: nb_task
    INTEGER,INTENT(OUT) :: id
  
    nb_grid=nb_grid+1
    id=nb_grid
    
    grids(id)%ni=ni_glo
    grids(id)%nj=nj_glo
    grids(id)%nb_task=nb_task
   
    ALLOCATE(grids(id)%task_domain(nb_task))
 
  END SUBROUTINE create_new_grid
 
  SUBROUTINE Init_task_domain(id,task,ni,nj,ibegin,jbegin,lon,lat)
  IMPLICIT NONE
    INTEGER, INTENT(IN) :: id
    INTEGER, INTENT(IN) :: task
    INTEGER, INTENT(IN) :: ni
    INTEGER, INTENT(IN) :: nj
    INTEGER, INTENT(IN) :: ibegin
    INTEGER, INTENT(IN) :: jbegin
    REAL,    INTENT(IN) :: lon(ni,nj)
    REAL,    INTENT(IN) :: lat(ni,nj)
 
    TYPE(domain), POINTER :: Pt_dom
    
    Pt_dom=>grids(id)%task_domain(task)
    
    Pt_dom%ni=ni
    Pt_dom%nj=nj
    Pt_dom%ibegin=ibegin
    Pt_dom%jbegin=jbegin   
    Pt_dom%iend=ibegin+ni-1
    Pt_dom%jend=jbegin+nj-1
    ALLOCATE(Pt_dom%lon(ni,nj))
    ALLOCATE(Pt_dom%lat(ni,nj))
    Pt_dom%lon(:,:)=lon(:,:)
    Pt_dom%lat(:,:)=lat(:,:)
  
  END SUBROUTINE Init_task_domain

  
  SUBROUTINE Init_grids
  IMPLICIT NONE
    INTEGER :: rank
    INTEGER :: ib,ie,jb,je
    TYPE(grid),POINTER :: Pt_grid
    TYPE(domain),POINTER :: local_domain
    TYPE(domain),POINTER :: task_domain(:)
    INTEGER :: id
    
      max_ni=0
      max_nj=0
    
      DO id=1,nb_grid
        pt_grid=>grids(id)
        local_domain=>pt_grid%local_domain
        task_domain=>pt_grid%task_domain
        
        local_domain%ibegin=task_domain(1)%ibegin
        local_domain%iend=task_domain(1)%iend
        local_domain%jbegin=task_domain(1)%jbegin
        local_domain%jend=task_domain(1)%jend
      
      DO rank=1,pt_grid%nb_task
        IF (task_domain(rank)%ibegin < local_domain%ibegin)   local_domain%ibegin = task_domain(rank)%ibegin
        IF (task_domain(rank)%iend   > local_domain%iend  )   local_domain%iend   = task_domain(rank)%iend
        IF (task_domain(rank)%jbegin < local_domain%jbegin)   local_domain%jbegin = task_domain(rank)%jbegin
        IF (task_domain(rank)%jend   > local_domain%jend  )   local_domain%jend   = task_domain(rank)%jend
      ENDDO
    
      local_domain%ni = local_domain%iend-local_domain%ibegin+1
      local_domain%nj = local_domain%jend-local_domain%jbegin+1
      ALLOCATE(local_domain%lon(local_domain%ni,local_domain%nj))
      ALLOCATE(local_domain%lat(local_domain%ni,local_domain%nj))
      
      DO rank=1,pt_grid%nb_task
        ib=task_domain(rank)%ibegin-local_domain%ibegin+1
        ie=task_domain(rank)%iend-local_domain%ibegin+1
        jb=task_domain(rank)%jbegin-local_domain%jbegin+1
        je=task_domain(rank)%jend-local_domain%jbegin+1

        task_domain(rank)%ibegin=ib
        task_domain(rank)%iend=ie
        task_domain(rank)%jbegin=jb
        task_domain(rank)%jend=je
      
        local_domain%lon(ib:ie,jb:je)=task_domain(rank)%lon(:,:)
        local_domain%lat(ib:ie,jb:je)=task_domain(rank)%lat(:,:)
      ENDDO
  
      IF (local_domain%ni > max_ni) max_ni=local_domain%ni
      IF (local_domain%nj > max_nj) max_nj=local_domain%nj
    ENDDO
  
    ALLOCATE(field_buffer(max_ni,max_nj,max_vert_size))
    
  END SUBROUTINE Init_grids
  
  SUBROUTINE bufferize_field2d(field,id,task)
  IMPLICIT NONE
    REAL,INTENT(IN) :: field(:,:)
    INTEGER,INTENT(IN) :: id
    INTEGER,INTENT(IN) :: task
    TYPE(domain), POINTER :: Pt

    Pt=>grids(id)%task_domain(task)
    field_buffer(Pt%ibegin:Pt%iend,Pt%jbegin:Pt%jend,1)=field(:,:)
      
  END SUBROUTINE bufferize_field2d
  
  
  SUBROUTINE bufferize_field3d(field,id,task)
  IMPLICIT NONE
    REAL,INTENT(IN) :: field(:,:,:)
    INTEGER,INTENT(IN) :: id
    INTEGER,INTENT(IN) :: task
    TYPE(domain), POINTER :: Pt
   
    Pt=>grids(id)%task_domain(task)

    IF (size(field,3)>max_vert_size) THEN
      DEALLOCATE(field_buffer)
      max_vert_size=max_vert_size*2
      ALLOCATE(field_buffer(max_ni,max_nj,max_vert_size))
    ENDIF
    
    field_buffer(Pt%ibegin:Pt%iend,Pt%jbegin:Pt%jend,1:size(field,3))=field(:,:,:)
      
  END SUBROUTINE bufferize_field3d
    
END MODULE box_grid
