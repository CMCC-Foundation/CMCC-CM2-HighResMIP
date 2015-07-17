MODULE mod_domain
   USE mod_xmlio_parameters
   INTEGER, PARAMETER :: box=1
   INTEGER, PARAMETER :: orange=2
   
  TYPE domain
    INTEGER :: rank
    INTEGER :: ni
    INTEGER :: nj
    INTEGER :: ibegin
    INTEGER :: iend
    INTEGER :: jbegin
    INTEGER :: jend
    REAL,POINTER :: lon(:,:)
    REAL,POINTER :: lat(:,:)
    
    INTEGER         :: type
    LOGICAL,POINTER :: mask(:)
    INTEGER,POINTER :: i_index(:)
    INTEGER,POINTER :: j_index(:)
    INTEGER         :: nbp
    LOGICAL         :: is_defined
  END TYPE domain
 
  INCLUDE 'vector_domain_def.inc'

CONTAINS

  INCLUDE 'vector_domain_contains.inc'

  SUBROUTINE domain__new(Pt_domain)
  IMPLICIT NONE
    TYPE(domain),POINTER  :: Pt_domain
    
    Pt_domain%is_defined=.FALSE.
    
  END SUBROUTINE domain__new
  
  SUBROUTINE domain__set(Pt_domain,rank,ni,nj,ibegin,jbegin,lon,lat)
  IMPLICIT NONE
    TYPE(domain),POINTER  :: Pt_domain
    INTEGER,INTENT(IN) :: rank
    INTEGER,INTENT(IN) :: ni
    INTEGER,INTENT(IN) :: nj
    INTEGER,INTENT(IN) :: ibegin
    INTEGER,INTENT(IN) :: jbegin
    REAL,INTENT(IN)    :: lon(ni,nj)
    REAL,INTENT(IN)    :: lat(ni,nj)

    Pt_domain%rank=rank
    Pt_domain%ni=ni
    Pt_domain%nj=nj
    Pt_domain%ibegin=ibegin
    Pt_domain%jbegin=jbegin   
    Pt_domain%iend=ibegin+ni-1
    Pt_domain%jend=jbegin+nj-1
    ALLOCATE(Pt_domain%lon(ni,nj))
    ALLOCATE(Pt_domain%lat(ni,nj))
    Pt_domain%lon(:,:)=lon(:,:)
    Pt_domain%lat(:,:)=lat(:,:)
    Pt_domain%is_defined=.TRUE.
  END SUBROUTINE domain__set


  SUBROUTINE domain__set_type_box(Pt_domain,mask)
  USE error_msg
  IMPLICIT NONE
    TYPE(domain),POINTER  :: Pt_domain
    LOGICAL,INTENT(IN),OPTIONAL :: mask(:,:)
    INTEGER                     :: i,j
  
    Pt_domain%type=box
    pt_domain%nbp=pt_domain%ni*pt_domain%nj
    ALLOCATE(Pt_domain%i_index(pt_domain%nbp))
    ALLOCATE(Pt_domain%j_index(pt_domain%nbp))
    ALLOCATE(pt_domain%mask(pt_domain%nbp))

    IF (PRESENT(mask)) THEN
      IF ( size(mask,1)==pt_domain%ni .AND. size(mask,2)==pt_domain%nj) THEN
        DO i=1,pt_domain%ni
          DO j=1,pt_domain%nj
            pt_domain%mask((j-1)*pt_domain%ni+i)=mask(i,j)
          ENDDO
        ENDDO
      ELSE
        WRITE (message,*) "mask dimensions are not compliant with domain dimensions :",   &
                           size(mask,1),",",size(mask,2),"/=",Pt_domain%ni,",",Pt_domain%nj
        CALL error("domain__set_type_box")
      ENDIF
    ELSE
      pt_domain%mask(:)=.TRUE.
    ENDIF
    
    DO i=1,pt_domain%ni
      DO j=1,pt_domain%nj
        pt_domain%i_index((j-1)*pt_domain%ni+i) = i
        pt_domain%j_index((j-1)*pt_domain%ni+i) = j
      ENDDO
    ENDDO

  END SUBROUTINE domain__set_type_box


  SUBROUTINE domain__set_type_orange(Pt_domain,nbp,offset,opt_index)
  IMPLICIT NONE
    TYPE(domain),POINTER        :: Pt_domain
    INTEGER,INTENT(IN)          :: nbp
    INTEGER,INTENT(IN)          :: offset
    INTEGER,INTENT(IN),OPTIONAL :: opt_index(:)
    
    INTEGER  :: index(nbp)
    INTEGER  :: i
    
    Pt_domain%type=orange
    
    IF (PRESENT(opt_index)) THEN
      index(:)=opt_index(:)       
    ELSE
      DO i=1,nbp
        index(i)=i
      ENDDO
    ENDIF

    pt_domain%nbp=nbp
    ALLOCATE(pt_domain%i_index(nbp))
    ALLOCATE(pt_domain%j_index(nbp))
    ALLOCATE(pt_domain%mask(nbp))
        
    DO i=1,nbp
!      Pt_domain%i_index(i)=(index(i)+offset)/pt_domain%ni+1+pt_domain%ibegin-1
!      Pt_domain%j_index(i)=MOD(index(i)+offset,pt_domain%ni)+1+pt_domain%jbegin-1
      Pt_domain%i_index(i)=MOD(index(i)+offset-1,pt_domain%ni)+1
      Pt_domain%j_index(i)=(index(i)+offset-1)/pt_domain%ni+1
      
    ENDDO

    Pt_domain%mask(:)=.TRUE.

  END SUBROUTINE domain__set_type_orange


    
  SUBROUTINE domain__print(pt_domain)
  IMPLICIT NONE
    TYPE(domain),POINTER   :: Pt_domain

    PRINT *,"---- DOMAIN ----"
    
    IF (pt_domain%is_defined) THEN
      PRINT *,"rank :",pt_domain%rank
      PRINT *,"ni :",pt_domain%ni
      PRINT *,"nj :",pt_domain%nj
      PRINT *,"ibegin",pt_domain%ibegin
      PRINT *,"iend",pt_domain%iend
    ELSE
      PRINT *,"  ---> domain undefined"
    ENDIF
    
    PRINT *,"-----------------"
  END SUBROUTINE domain__print
  
  SUBROUTINE domain__copy(Pt_in,pt_out)
  IMPLICIT NONE
    TYPE(domain),POINTER   :: Pt_in
    TYPE(domain),POINTER   :: Pt_out
  
    IF (pt_in%is_defined) THEN
      CALL domain__set(pt_out,pt_in%rank,pt_in%ni,pt_in%nj,pt_in%ibegin,pt_in%jbegin,pt_in%lon,pt_in%lat)
    ENDIF
  END SUBROUTINE domain__copy
  
END MODULE mod_domain
