MODULE mod_zoom
  USE mod_xmlio_parameters
  USE mod_sorted_list

  IMPLICIT NONE

  TYPE,PUBLIC :: zoom
    CHARACTER(len=str_len)      :: id
    LOGICAL                     :: has_id
    CHARACTER(len=str_len)      :: name
    LOGICAL                     :: has_name
    CHARACTER(len=str_len)      :: description
    LOGICAL                     :: has_description
    INTEGER                     :: ni_glo
    LOGICAL                     :: has_ni_glo
    INTEGER                     :: nj_glo
    LOGICAL                     :: has_nj_glo
    INTEGER                     :: ibegin_glo
    LOGICAL                     :: has_ibegin_glo
    INTEGER                     :: jbegin_glo
    LOGICAL                     :: has_jbegin_glo
    INTEGER                     :: ni_loc
    INTEGER                     :: nj_loc
    INTEGER                     :: ibegin_loc
    INTEGER                     :: jbegin_loc
  END TYPE zoom    

  INCLUDE 'vector_zoom_def.inc'

  TYPE(vector_zoom),POINTER,SAVE             :: zoom_Ids
  TYPE(sorted_list),POINTER,SAVE,PRIVATE     :: Ids 

  INTERFACE zoom__set_attribut
    MODULE PROCEDURE zoom__set_attribut_id,zoom__set_attribut_pt
  END INTERFACE
  
CONTAINS
  
  INCLUDE 'vector_zoom_contains.inc'

  SUBROUTINE zoom__swap_context(saved_zoom_Ids,saved_Ids)
  IMPLICIT NONE
    TYPE(vector_zoom),POINTER          :: saved_zoom_Ids
    TYPE(sorted_list),POINTER          :: saved_Ids 
    
    zoom_ids=>saved_zoom_ids
    ids=>saved_ids
  END SUBROUTINE zoom__swap_context


  SUBROUTINE zoom__init
  IMPLICIT NONE

    CALL vector_zoom__new(zoom_Ids)
    CALL sorted_list__new(Ids)
    
  END SUBROUTINE zoom__init
  

  SUBROUTINE zoom__new(pt_zoom,Id)
  USE string_function
  IMPLICIT NONE
    TYPE(zoom),POINTER :: pt_zoom
    CHARACTER(LEN=*),OPTIONAL     :: Id
    INTEGER   :: Pos 
      
    pt_zoom%has_id=.FALSE.
    pt_zoom%has_name=.FALSE.
    pt_zoom%has_description=.FALSE.
    pt_zoom%has_ni_glo=.FALSE.
    pt_zoom%has_nj_glo=.FALSE.
    pt_zoom%has_ibegin_glo=.FALSE.
    pt_zoom%has_jbegin_glo=.FALSE.
     
    IF (PRESENT(Id)) THEN
      Pt_zoom%id=TRIM(ADJUSTL(Id))
      Pt_zoom%has_id=.TRUE.
     CALL vector_zoom__set_new(zoom_Ids,Pt_zoom,Pos)
     CALL sorted_list__Add(Ids,hash(id),Pos)
    ENDIF
   
  END SUBROUTINE zoom__new
   
  SUBROUTINE zoom__set(pt_zoom,name,description,ni_glo,nj_glo,ibegin_glo,jbegin_glo)
  IMPLICIT NONE
    TYPE(zoom),POINTER                   :: pt_zoom
    CHARACTER(len=*),OPTIONAL            :: name
    CHARACTER(len=*),OPTIONAL            :: description
    INTEGER,OPTIONAL                     :: ni_glo
    INTEGER,OPTIONAL                     :: nj_glo
    INTEGER,OPTIONAL                     :: ibegin_glo
    INTEGER,OPTIONAL                     :: jbegin_glo
    
    IF (PRESENT(name)) THEN
      pt_zoom%name=TRIM(ADJUSTL(name))
      pt_zoom%has_name=.TRUE.
    ENDIF
    
    IF (PRESENT(description)) THEN
      pt_zoom%name=TRIM(ADJUSTL(description))
      pt_zoom%has_description=.TRUE.
    ENDIF
    
    IF (PRESENT(ni_glo)) THEN
      pt_zoom%ni_glo=ni_glo
      pt_zoom%has_ni_glo=.TRUE.
    ENDIF

    IF (PRESENT(nj_glo)) THEN
      pt_zoom%nj_glo=nj_glo
      pt_zoom%has_nj_glo=.TRUE.
    ENDIF

    IF (PRESENT(ibegin_glo)) THEN
      pt_zoom%ibegin_glo=ibegin_glo
      pt_zoom%has_ibegin_glo=.TRUE.
    ENDIF
    
    IF (PRESENT(jbegin_glo)) THEN
      pt_zoom%jbegin_glo=jbegin_glo
      pt_zoom%has_ibegin_glo=.TRUE.
    ENDIF
               
   END SUBROUTINE zoom__set   

  SUBROUTINE zoom__set_attribut_id(id,attrib,ok)
  USE mod_attribut
  USE error_msg
  IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN)   :: id
    TYPE(attribut),INTENT(IN)     :: attrib
    LOGICAL,OPTIONAL,INTENT(OUT)  :: ok
    
    TYPE(zoom),POINTER              :: Pt_zoom
    INTEGER                         :: Pos
    LOGICAL                         :: success
    
    CALL sorted_list__find(Ids,hash(Id),Pos,success)
    IF (success) THEN
      Pt_zoom=>zoom_ids%at(Pos)%Pt
      CALL zoom__set_attribut_pt(Pt_zoom,attrib)
      IF (PRESENT(OK)) OK=.TRUE.
    ELSE
      IF (.NOT.PRESENT(OK)) THEN
        WRITE(message,*) 'zoom id :',id,'is undefined'
        CALL error('mod_zoom::zoom__set_attribut')
      ELSE
        OK=.FALSE.
      ENDIF
    ENDIF  
  
  END SUBROUTINE zoom__set_attribut_id
      
  SUBROUTINE zoom__set_attribut_pt(Pt_zoom,attrib)
  USE mod_attribut
  USE mod_zoom_attribut
  USE error_msg
  IMPLICIT NONE
    TYPE(zoom),POINTER        :: Pt_zoom
    TYPE(attribut),INTENT(IN) :: attrib
    
    SELECT CASE(attrib%name)
      CASE (zoom__name)
        IF (attrib%type==string0) CALL  zoom__set(pt_zoom,name=attrib%string0_ptr) ; RETURN
      CASE (zoom__description)
        IF (attrib%type==string0) CALL  zoom__set(pt_zoom,description=attrib%string0_ptr) ; RETURN
      CASE (zoom__ni)
        IF (attrib%type==integer0) CALL  zoom__set(pt_zoom,ni_glo=attrib%integer0_ptr) ; RETURN
      CASE (zoom__nj)
        IF (attrib%type==integer0) CALL  zoom__set(pt_zoom,nj_glo=attrib%integer0_ptr) ; RETURN
      CASE (zoom__ibegin)
        IF (attrib%type==integer0) CALL  zoom__set(pt_zoom,ibegin_glo=attrib%integer0_ptr) ; RETURN
      CASE (zoom__jbegin)
        IF (attrib%type==integer0) CALL  zoom__set(pt_zoom,jbegin_glo=attrib%integer0_ptr) ; RETURN
     END SELECT

     WRITE(message,*) 'zoom attribut ',attrib%name,' : type :',attrib%type,   &
                      ' : Attribute type is incompatible with the provided value'
     CALL error('mod_zoom::zoom__set_attribut')
    
  END SUBROUTINE zoom__set_attribut_pt

  SUBROUTINE zoom__get(Id,pt_zoom)
  USE string_function
  IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN)     :: Id
    TYPE(zoom),POINTER              :: Pt_zoom

    INTEGER                         :: Pos
    LOGICAL                         :: success
    
    CALL sorted_list__find(Ids,hash(Id),Pos,success)
    IF (success) THEN
      Pt_zoom=>zoom_ids%at(Pos)%Pt
    ELSE
      Pt_zoom=>NULL()
    ENDIF
    
  END SUBROUTINE zoom__get

END MODULE mod_zoom
