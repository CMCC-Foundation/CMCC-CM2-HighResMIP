MODULE mod_axis_group
  USE mod_axis
  USE mod_xmlio_parameters

  IMPLICIT NONE

  TYPE axis_group
    CHARACTER(LEN=str_len)                    :: id
    LOGICAL                                   :: has_id
    TYPE(vector_axis_group), POINTER          :: groups
    TYPE(vector_axis),POINTER                 :: axis     
    TYPE(axis), POINTER                       :: default_attribut
  END TYPE axis_group

  INCLUDE "vector_axis_group_def.inc"  

  TYPE(vector_axis_group),POINTER       :: axis_group_Ids
  TYPE(sorted_list),POINTER,PRIVATE     :: Ids 

  INTERFACE axis_group__set_attribut
    MODULE PROCEDURE axis_group__set_attribut_id,axis_group__set_attribut_pt
  END INTERFACE

CONTAINS

  INCLUDE "vector_axis_group_contains.inc"

  SUBROUTINE axis_group__swap_context(saved_axis_group_Ids,saved_ids)
  IMPLICIT NONE
  TYPE(vector_axis_group),POINTER       :: saved_axis_group_Ids
  TYPE(sorted_list),POINTER             :: saved_Ids 
   
   axis_group_ids=>saved_axis_group_ids
   ids=>saved_ids
   
  END SUBROUTINE axis_group__swap_context

  SUBROUTINE axis_group__init
  IMPLICIT NONE
    
    CALL vector_axis_group__new(axis_group_Ids)
    CALL sorted_list__new(Ids)
   
  END SUBROUTINE axis_group__init

  SUBROUTINE axis_group__get(Id,Pt_ag)
  USE string_function
  IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN)     :: Id
    TYPE(axis_group),POINTER        :: Pt_ag

    INTEGER                         :: Pos
    LOGICAL                         :: success
    
    CALL sorted_list__find(Ids,hash(Id),Pos,success)
    IF (success) THEN
      Pt_ag=>axis_group_ids%at(Pos)%Pt
    ELSE
      Pt_ag=>NULL()
    ENDIF
    
  END SUBROUTINE axis_group__get


  SUBROUTINE axis_group__set_attribut_id(id,attrib,Ok)
  USE mod_attribut
  USE error_msg
  IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN) :: id
    TYPE(attribut),INTENT(IN) :: attrib
    LOGICAL,OPTIONAL,INTENT(out)  :: Ok
    
    TYPE(axis_group),POINTER             :: Pt_ag
    INTEGER                         :: Pos
    LOGICAL                         :: success
    
    CALL sorted_list__find(Ids,hash(Id),Pos,success)
    IF (success) THEN
      Pt_ag=>axis_group_ids%at(Pos)%Pt
      CALL axis_group__set_attribut(Pt_ag,attrib)
      IF (PRESENT(OK)) ok=.TRUE.
    ELSE
      IF (.NOT.PRESENT(OK)) THEN
        WRITE(message,*) 'axis group id :',id,'is undefined'
        CALL error('mod_axis_group::axis_group__set_attribut')
      ELSE
        OK=.FALSE.
      ENDIF
    ENDIF 
    
  END SUBROUTINE axis_group__set_attribut_id
  
  SUBROUTINE axis_group__set_attribut_pt(pt_ag,attrib)
  USE mod_attribut
  USE mod_object
  IMPLICIT NONE
    TYPE(axis_group),POINTER :: Pt_ag
    TYPE(attribut),INTENT(IN) :: attrib
     
    IF (attrib%object==axis_object) THEN
      CALL axis__set_attribut(pt_ag%default_attribut,attrib)
    ENDIF
    
  END SUBROUTINE axis_group__set_attribut_pt
      
  RECURSIVE SUBROUTINE axis_group__new(Pt_ag,Id)
  USE string_function
  IMPLICIT NONE
    TYPE(axis_group),POINTER                :: Pt_ag
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN)     :: Id
    
    INTEGER :: Pos
    
    ALLOCATE(Pt_ag%groups)
    ALLOCATE(Pt_ag%axis)
    ALLOCATE(Pt_ag%default_attribut)
    
    CALL vector_axis_group__new(Pt_ag%groups)
    CALL vector_axis__new(Pt_ag%axis)
    CALL axis__new(Pt_ag%default_attribut)
    Pt_ag%has_id=.FALSE.
      
    IF (PRESENT(Id)) THEN
      Pt_ag%id=TRIM(Id)
      Pt_ag%has_id=.TRUE.
      CALL vector_axis_group__set_new(axis_group_Ids,Pt_ag,Pos)
      CALL sorted_list__Add(Ids,hash(id),Pos)
    ENDIF

  END SUBROUTINE axis_group__new

      
  SUBROUTINE axis_group__get_new_group(Pt_ag,Pt_ag_out,Id)
  IMPLICIT NONE
    TYPE(axis_group),POINTER             :: Pt_ag
    TYPE(axis_group),POINTER             :: Pt_ag_out
    CHARACTER(LEN=*),OPTIONAL      :: Id
    
    CALL vector_axis_group__get_new(Pt_ag%groups,Pt_ag_out)
    CALL axis_group__new(Pt_ag_out)

    IF (PRESENT(id)) THEN
      CALL axis_group__new(Pt_ag_out,Id)
    ELSE
      CALL axis_group__new(Pt_ag_out)
    ENDIF
    
  END SUBROUTINE axis_group__get_new_group

  
  SUBROUTINE axis_group__get_new_axis(Pt_ag,Pt_a_out,Id)
  IMPLICIT NONE
    TYPE(axis_group),POINTER            :: Pt_ag
    TYPE(axis),POINTER                  :: Pt_a_out
    CHARACTER(LEN=*),OPTIONAL      :: Id
    
    CALL vector_axis__get_new(Pt_ag%axis,Pt_a_out)
    
    IF (PRESENT(id)) THEN
      CALL axis__new(Pt_a_out,Id)
    ELSE
      CALL axis__new(Pt_a_out)
    ENDIF
    
  END SUBROUTINE axis_group__get_new_axis
  
  
  SUBROUTINE axis_group__get_default_attrib(Pt_ag,Pt_a)
  IMPLICIT NONE
    TYPE(axis_group),POINTER  :: Pt_ag
    TYPE(axis),POINTER        :: Pt_a
    
    Pt_a=>Pt_ag%default_attribut
  END SUBROUTINE axis_group__get_default_attrib
  
  RECURSIVE SUBROUTINE axis_group__apply_default(Pt_ag,default)
  IMPLICIT NONE
    TYPE(axis_group),POINTER           :: Pt_ag
    TYPE(axis),POINTER,OPTIONAL        :: default
    INTEGER :: i
    
    IF (PRESENT(default)) THEN
      CALL axis__apply_default(default,Pt_ag%default_attribut,Pt_ag%default_attribut)
    ENDIF
      
    DO i=1,Pt_ag%groups%size
      CALL axis_group__apply_default(Pt_ag%groups%at(i)%pt,Pt_ag%default_attribut)
    ENDDO
    
    DO i=1,Pt_ag%axis%size
      CALL axis__apply_default(Pt_ag%default_attribut,Pt_ag%axis%at(i)%pt,Pt_ag%axis%at(i)%pt)
    ENDDO
  
  END SUBROUTINE axis_group__apply_default

  RECURSIVE SUBROUTINE axis_group__print(Pt_ag)
  IMPLICIT NONE
    TYPE(axis_group),POINTER  :: Pt_ag
   
    INTEGER :: i
    
    PRINT *,"--- AXIS GROUP ---"
    IF (pt_ag%has_id) THEN
      PRINT *,"id = ",TRIM(pt_ag%id)
    ELSE
      PRINT *,"id undefined"
    ENDIF
    
    PRINT *,"axis default attribut :"
    CALL axis__print(Pt_ag%default_attribut)    

    PRINT *,"owned axis groups :",Pt_ag%groups%size      
    DO i=1,Pt_ag%groups%size
      CALL axis_group__print(Pt_ag%groups%at(i)%pt)
    ENDDO

    PRINT *,"owned axis :",Pt_ag%axis%size      
    DO i=1,Pt_ag%axis%size
      CALL axis__print(Pt_ag%axis%at(i)%pt)
    ENDDO
    PRINT *,"------------"
    
  END SUBROUTINE axis_group__print

  RECURSIVE SUBROUTINE axis_group__Check(Pt_ag)
  IMPLICIT NONE
  
    TYPE(axis_group),POINTER  :: Pt_ag
    INTEGER :: i
   
    DO i=1,Pt_ag%groups%size
      CALL axis_group__check(pt_ag%groups%at(i)%pt)
    ENDDO

    DO i=1,Pt_ag%axis%size
      CALL axis__check(pt_ag%axis%at(i)%pt)
    ENDDO
  
  END SUBROUTINE axis_group__check     
   
END MODULE mod_axis_group

