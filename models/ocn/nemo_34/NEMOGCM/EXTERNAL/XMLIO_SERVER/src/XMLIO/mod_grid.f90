MODULE mod_grid
  USE mod_xmlio_parameters
  USE mod_sorted_list
  USE mod_domain
  USE mod_zoom

  IMPLICIT NONE

  TYPE, PUBLIC :: grid
    CHARACTER(len=str_len)      :: id
    LOGICAL                     :: has_id
    CHARACTER(len=str_len)      :: name
    LOGICAL                     :: has_name
    CHARACTER(len=str_len)      :: description
    LOGICAL                     :: has_description
    TYPE(domain),POINTER        :: domain
    TYPE(vector_domain),POINTER :: subdomain
    TYPE(sorted_list),POINTER   :: rank_ids
    INTEGER,POINTER             :: ranks(:)
    INTEGER                     :: ni
    INTEGER                     :: nj
    LOGICAL                     :: has_dimension
    TYPE(vector_zoom),POINTER   :: associated_zoom
    TYPE(zoom),POINTER          :: global_zoom
  END TYPE grid

  INCLUDE 'vector_grid_def.inc'
  
  TYPE(vector_grid),POINTER,SAVE             :: grid_Ids
  TYPE(sorted_list),POINTER,SAVE,PRIVATE     :: Ids 

  INTERFACE grid__set_attribut
    MODULE PROCEDURE grid__set_attribut_id,grid__set_attribut_pt
  END INTERFACE

CONTAINS
  INCLUDE 'vector_grid_contains.inc'

  SUBROUTINE grid__swap_context(saved_grid_Ids,saved_Ids)
  IMPLICIT NONE
    TYPE(vector_grid),POINTER          :: saved_grid_Ids
    TYPE(sorted_list),POINTER          :: saved_Ids 
    
    grid_ids=>saved_grid_ids
    ids=>saved_ids
  END SUBROUTINE grid__swap_context
  
  
  SUBROUTINE grid__init
  IMPLICIT NONE
    
    CALL vector_grid__new(grid_Ids)
    CALL sorted_list__new(Ids)
   
  END SUBROUTINE grid__init
  
  SUBROUTINE grid__get(Id,Pt_grid)
  USE string_function
  IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN)     :: Id
    TYPE(grid),POINTER              :: Pt_grid

    INTEGER                         :: Pos
    LOGICAL                         :: success
    
    CALL sorted_list__find(Ids,hash(Id),Pos,success)
    IF (success) THEN
      Pt_grid=>grid_ids%at(Pos)%Pt
    ELSE
      Pt_grid=>NULL()
    ENDIF
    
  END SUBROUTINE grid__get
  
  SUBROUTINE grid__new(pt_grid,Id)
  USE string_function
  IMPLICIT NONE
   TYPE(grid), POINTER           :: pt_grid
   CHARACTER(LEN=*),OPTIONAL     :: Id
   INTEGER                       :: Pos
   
   ALLOCATE(pt_grid%domain)
   ALLOCATE(pt_grid%subdomain)
   ALLOCATE(pt_grid%rank_ids)
   ALLOCATE(pt_grid%associated_zoom)
   
   CALL domain__new(pt_grid%domain)
   CALL vector_domain__new(pt_grid%subdomain)
   CALL sorted_list__new(pt_grid%rank_ids)
   CALL vector_zoom__new(pt_grid%associated_zoom)
   
   pt_grid%has_id          = .FALSE.
   pt_grid%has_name        = .FALSE.
   pt_grid%has_description = .FALSE.
   pt_grid%has_dimension   = .FALSE.
   
   IF (PRESENT(Id)) THEN
     Pt_grid%id=TRIM(ADJUSTL(Id))
     Pt_grid%has_id=.TRUE.
     CALL vector_grid__set_new(grid_Ids,Pt_grid,Pos)
     CALL sorted_list__Add(Ids,hash(id),Pos)
   ENDIF
   
   CALL grid__get_new_zoom(pt_grid,pt_grid%global_zoom,id)

 END SUBROUTINE grid__new

  SUBROUTINE grid__set(pt_grid, name, description)
  IMPLICIT NONE
    TYPE(grid), POINTER :: pt_grid
    CHARACTER(len=*)  ,OPTIONAL :: name
    CHARACTER(len=*)  ,OPTIONAL :: description

    IF (PRESENT(name)) THEN
        pt_grid%name=TRIM(ADJUSTL(name))
        pt_grid%has_name = .TRUE.
    ENDIF

    IF (PRESENT(description)) THEN
        pt_grid%description=TRIM(ADJUSTL(description))
        pt_grid%has_description = .TRUE.
    ENDIF

  END SUBROUTINE grid__set


  SUBROUTINE grid__set_attribut_id(id,attrib,ok)
  USE mod_attribut
  USE error_msg
  IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN)   :: id
    TYPE(attribut),INTENT(IN)     :: attrib
    LOGICAL,OPTIONAL,INTENT(OUT)  :: ok
    
    TYPE(grid),POINTER              :: Pt_grid
    INTEGER                         :: Pos
    LOGICAL                         :: success
    
    CALL sorted_list__find(Ids,hash(Id),Pos,success)
    IF (success) THEN
      Pt_grid=>grid_ids%at(Pos)%Pt
      CALL grid__set_attribut_pt(Pt_grid,attrib)
      IF (PRESENT(OK)) OK=.TRUE.
    ELSE
      IF (.NOT.PRESENT(OK)) THEN
        WRITE(message,*) 'grid id :',id,'is undefined'
        CALL error('mod_grid::grid__set_attribut')
      ELSE
        OK=.FALSE.
      ENDIF
    ENDIF  
  
  END SUBROUTINE grid__set_attribut_id
      
  SUBROUTINE grid__set_attribut_pt(Pt_grid,attrib)
  USE mod_attribut
  USE mod_grid_attribut
  USE error_msg
  IMPLICIT NONE
    TYPE(grid),POINTER        :: Pt_grid
    TYPE(attribut),INTENT(IN) :: attrib
    
    SELECT CASE(attrib%name)
      CASE (grid__name)
        IF (attrib%type==string0) CALL  grid__set(pt_grid,name=attrib%string0_ptr) ; RETURN
      CASE (grid__description)
        IF (attrib%type==string0) CALL  grid__set(pt_grid,description=attrib%string0_ptr) ; RETURN
     END SELECT

     WRITE(message,*) 'grid attribut ',attrib%name,' : type :',attrib%type,      &
                      ' : Attribute type is incompatible with the provided value'
     CALL error('mod_grid::grid__set_attribut')
    
  END SUBROUTINE grid__set_attribut_pt
  
  SUBROUTINE grid__set_dimension(pt_grid, ni, nj)
  IMPLICIT NONE
    TYPE(grid), POINTER   :: pt_grid
    INTEGER,INTENT(IN)    :: ni
    INTEGER,INTENT(IN)    :: nj
    
    pt_grid%ni=ni
    pt_grid%nj=nj
    pt_grid%has_dimension=.TRUE.
    
  END SUBROUTINE grid__set_dimension
    

  SUBROUTINE grid__get_new_subdomain(Pt_grid,rank,pt_domain)
  IMPLICIT NONE
    TYPE(grid), POINTER   :: pt_grid
    INTEGER,INTENT(IN)    :: rank
    TYPE(domain), POINTER :: Pt_domain
    
    INTEGER :: Pos
    
    CALL vector_domain__get_new(pt_grid%subdomain,pt_domain,Pos)
    CALL sorted_list__add(pt_grid%rank_ids,rank,Pos)
    CALL domain__new(pt_domain)
    
  END SUBROUTINE grid__get_new_subdomain    

  SUBROUTINE grid__get_subdomain(Pt_grid,rank,pt_domain)
  IMPLICIT NONE
    TYPE(grid), POINTER     :: pt_grid
    INTEGER,INTENT(IN)      :: rank
    TYPE(domain), POINTER   :: Pt_domain

    INTEGER    :: rank_id
    LOGICAL    :: success

    CALL sorted_list__find(pt_grid%rank_ids,rank,rank_id,success)
    IF (success) THEN
      pt_domain=>pt_grid%subdomain%at(rank_id)%pt
    ELSE
      !! message d'erreur 
   ENDIF

  END SUBROUTINE grid__get_subdomain
   
  SUBROUTINE grid__process_domain(Pt_grid)
  IMPLICIT NONE
    TYPE(grid), POINTER  :: pt_grid
    TYPE(domain),POINTER :: subdomain
    TYPE(zoom),POINTER :: pt_zoom
    
    REAL,ALLOCATABLE :: lon(:,:)
    REAL,ALLOCATABLE :: lat(:,:)
    INTEGER :: ib,ie,jb,je,ni,nj,ibegin,jbegin,iend,jend
    INTEGER :: i
    
      ALLOCATE(pt_grid%ranks(1:pt_grid%subdomain%size))

      DO i=1,pt_grid%subdomain%size
        subdomain=>pt_grid%subdomain%at(i)%pt
        IF (i==1) THEN
          ib=subdomain%ibegin
          ie=subdomain%iend   
          jb=subdomain%jbegin
          je=subdomain%jend
        ELSE
          IF (subdomain%ibegin<ib) ib=subdomain%ibegin
          IF (subdomain%iend>ie) ie=subdomain%iend
          IF (subdomain%jbegin<jb) jb=subdomain%jbegin
          IF (subdomain%jend>je) je=subdomain%jend
        ENDIF
        pt_grid%ranks(subdomain%rank)=i
      ENDDO
      
      ni=ie-ib+1
      nj=je-jb+1
      ibegin=ib
      jbegin=jb
      
      ALLOCATE(lon(ni,nj))
      ALLOCATE(lat(ni,nj))
      
      DO i=1,pt_grid%subdomain%size
        subdomain=>pt_grid%subdomain%at(i)%pt
        ib=subdomain%ibegin-ibegin+1
        ie=subdomain%iend-ibegin+1   
        jb=subdomain%jbegin-jbegin+1
        je=subdomain%jend-jbegin+1
        lon(ib:ie,jb:je)=subdomain%lon(:,:)
        lat(ib:ie,jb:je)=subdomain%lat(:,:)
      ENDDO
      
      CALL domain__set(pt_grid%domain,0,ni,nj,ibegin,jbegin,lon,lat)
      iend=ibegin+ni-1
      jend=jbegin+nj-1
      
      
      pt_grid%global_zoom%ni_glo=pt_grid%ni
      pt_grid%global_zoom%nj_glo=pt_grid%nj
      pt_grid%global_zoom%ibegin_glo=1      
      pt_grid%global_zoom%jbegin_glo=1
      
      DO i=1,pt_grid%associated_zoom%size
        pt_zoom=>pt_grid%associated_zoom%at(i)%pt
        
        ib=MAX(pt_zoom%ibegin_glo-ibegin+1,1)
        ie=MIN(pt_zoom%ibegin_glo+pt_zoom%ni_glo-ibegin,ni)
        pt_zoom%ni_loc=MAX(ie-ib+1,0)
        pt_zoom%ibegin_loc=ib

        jb=MAX(pt_zoom%jbegin_glo-jbegin+1,1)
        je=MIN(pt_zoom%jbegin_glo+pt_zoom%nj_glo-jbegin,nj)
        pt_zoom%nj_loc=MAX(je-jb+1,0)
        pt_zoom%jbegin_loc=jb
      ENDDO
                 
            
      DEALLOCATE(lon)
      DEALLOCATE(lat)
      
    END SUBROUTINE grid__process_domain
      
            
  SUBROUTINE grid__get_new_zoom(pt_grid,pt_zoom,zoom_id)
  USE string_function
  IMPLICIT NONE
    TYPE(grid), POINTER                      :: pt_grid
    TYPE(zoom),POINTER                       :: pt_zoom
    CHARACTER(LEN=*),INTENT(IN),OPTIONAL     :: zoom_id
    LOGICAL                                  :: success
    
     CALL vector_zoom__get_new(pt_grid%associated_zoom,Pt_zoom)
     CALL zoom__new(Pt_zoom,zoom_id)
    
   END SUBROUTINE grid__get_new_zoom
     
  SUBROUTINE grid__print(pt_grid)
  IMPLICIT NONE
    TYPE(grid), POINTER  :: pt_grid
    INTEGER              :: i
    
    PRINT *,"---- GRID ----"
    
    IF (pt_grid%has_id) THEN
      PRINT *,"id = ",TRIM(pt_grid%id)
    ELSE
      PRINT *,"id undefined"
    ENDIF
    
    IF (pt_grid%has_name) THEN
      PRINT *,"name = ",TRIM(pt_grid%name)
    ELSE
      PRINT *,"name undefined"
    ENDIF
    
    IF (pt_grid%has_description) THEN
      PRINT *,"description = ",TRIM(pt_grid%description)
    ELSE
      PRINT *,"description undefined"
    ENDIF
    
    IF (pt_grid%has_dimension) THEN
      PRINT *,"Global grid dimension :"
      PRINT *,"   ni =",pt_grid%ni
      PRINT *,"   nj =",pt_grid%nj
    ELSE
      PRINT *,"Global dimension ni, nj undefined"
    ENDIF
    
    PRINT *,"grid domain :"
    CALL domain__print(pt_grid%domain)
    
    PRINT *,"grid subdomain :",pt_grid%subdomain%size
    
    DO i=1,pt_grid%subdomain%size
      CALL domain__print(pt_grid%subdomain%at(i)%pt)
    ENDDO
    
    PRINT *,"--------------"
    
  END SUBROUTINE grid__print

  SUBROUTINE grid__apply_default(pt_grid_default, pt_grid_in, pt_grid_out)

    TYPE(grid), POINTER :: pt_grid_default, pt_grid_in, pt_grid_out

    IF (pt_grid_in%has_name) THEN
        pt_grid_out%name=pt_grid_in%name
        pt_grid_out%has_name=.TRUE.
    ELSE IF ( pt_grid_default%has_name) THEN
        pt_grid_out%name=pt_grid_default%name
        pt_grid_out%has_name=.TRUE.
    ELSE
        pt_grid_out%has_name=.FALSE.
    ENDIF
        
    IF (pt_grid_in%has_description) THEN
        pt_grid_out%description=pt_grid_in%description
        pt_grid_out%has_description=.TRUE.
    ELSE IF ( pt_grid_default%has_description ) THEN
        pt_grid_out%description=pt_grid_default%description
        pt_grid_out%has_description=.TRUE.
    ELSE
        pt_grid_out%has_description=.FALSE.
    ENDIF
  
  END SUBROUTINE grid__apply_default

END MODULE mod_grid
