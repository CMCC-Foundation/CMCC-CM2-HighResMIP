MODULE field_bufferize
  
  INTEGER,SAVE      :: nx=0
  INTEGER,SAVE      :: ny=0
  INTEGER,SAVE      :: nz=0
    
  REAL,POINTER,SAVE :: field_buffer(:,:,:)
  

  
CONTAINS

  SUBROUTINE bufferize_field(ni,ibegin,nj,jbegin,nk,kbegin,nbp,field,i_index,j_index,mask)
  IMPLICIT NONE
    INTEGER,INTENT(IN) :: ni
    INTEGER,INTENT(IN) :: ibegin
    INTEGER,INTENT(IN) :: nj
    INTEGER,INTENT(IN) :: jbegin
    INTEGER,INTENT(IN) :: nk
    INTEGER,INTENT(IN) :: kbegin
    INTEGER,INTENT(IN) :: nbp
    REAL,INTENT(IN)    :: field(nbp,nk)
    INTEGER,INTENT(IN) :: i_index(nbp)
    INTEGER,INTENT(IN) :: j_index(nbp)
    LOGICAL,INTENT(IN),OPTIONAL :: mask(nbp)  
    
    INTEGER :: iend
    INTEGER :: jend
    INTEGER :: kend 

    LOGICAL      :: need_reallocate
    REAL,POINTER :: tmp_buffer(:,:,:)
    INTEGER      :: i,j,k,n
    
    iend=ibegin+ni-1
    jend=jbegin+nj-1
    kend=kbegin+nk-1
    
    IF (PRESENT(mask)) THEN
      DO k=kbegin,kend
        DO n=1,nbp
          i=i_index(n)+ibegin-1
          j=j_index(n)+jbegin-1
          IF (mask(n)) field_buffer(i,j,k)=field(n,k)
        ENDDO
      ENDDO
    ELSE
      DO k=kbegin,kend
        DO n=1,nbp
          i=i_index(n)+ibegin-1
          j=j_index(n)+jbegin-1
          field_buffer(i,j,k)=field(n,k)
        ENDDO
      ENDDO
    ENDIF
    
  END SUBROUTINE bufferize_field

  SUBROUTINE Init_field_Bufferize(ni,nj,nk)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: ni
  INTEGER,INTENT(IN) :: nj
  INTEGER,INTENT(IN) :: nk
     
    LOGICAL :: need_reallocate
    
    
    need_reallocate=.FALSE.
    
    IF ( ni > nx) THEN
      nx=ni
      need_reallocate=.TRUE.
    ENDIF
    
    IF ( nj > ny) THEN
      ny=nj
      need_reallocate=.TRUE.
    ENDIF

    IF ( nk > nz) THEN
      nz=nk
      need_reallocate=.TRUE.
    ENDIF
    
    IF (need_reallocate) THEN
      IF (ASSOCIATED(field_buffer)) DEALLOCATE(field_buffer)
      ALLOCATE(field_buffer(nx,ny,nz))
    ENDIF

    field_buffer(1:nx,1:ny,1:nz)=0.

  END SUBROUTINE  init_field_bufferize         

END MODULE field_bufferize
