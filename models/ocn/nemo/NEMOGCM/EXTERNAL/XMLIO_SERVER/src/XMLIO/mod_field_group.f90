MODULE mod_field_group

  USE mod_field
  USE mod_xmlio_parameters

  IMPLICIT NONE

  TYPE field_group
    CHARACTER(LEN=str_len)                    :: id
    LOGICAL                                   :: has_id
    TYPE(vector_field_group), POINTER         :: groups
    TYPE(vector_field),POINTER                :: fields     
    TYPE(field), POINTER                      :: default_attribut
  END TYPE field_group

  INCLUDE "vector_field_group_def.inc"  

  TYPE(vector_field_group),SAVE,POINTER      :: field_group_Ids
  TYPE(sorted_list),POINTER,SAVE,PRIVATE     :: Ids 

  INTERFACE field_group__set_attribut
    MODULE PROCEDURE field_group__set_attribut_id,field_group__set_attribut_pt
  END INTERFACE

CONTAINS

  INCLUDE "vector_field_group_contains.inc"


  SUBROUTINE field_group__swap_context(saved_field_group_ids, saved_ids)
  IMPLICIT NONE
    TYPE(vector_field_group),POINTER   :: saved_field_group_Ids
    TYPE(sorted_list),POINTER          :: saved_Ids
    
    field_group_ids=>saved_field_group_ids 
    ids=>saved_ids
    
  END SUBROUTINE field_group__swap_context

  SUBROUTINE field_group__init
  IMPLICIT NONE
    
    CALL vector_field_group__new(field_group_Ids)
    CALL sorted_list__new(Ids)
   
  END SUBROUTINE field_group__init

  SUBROUTINE field_group__get(Id,Pt_fg)
  USE string_function
  IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN)     :: Id
    TYPE(field_group),POINTER       :: Pt_fg

    INTEGER                         :: Pos
    LOGICAL                         :: success
    
    CALL sorted_list__find(Ids,hash(Id),Pos,success)
    IF (success) THEN
      Pt_fg=>field_group_ids%at(Pos)%Pt
    ELSE
      Pt_fg=>NULL()
    ENDIF
    
  END SUBROUTINE field_group__get

  SUBROUTINE field_group__set_attribut_id(id,attrib,Ok)
  USE mod_attribut
  USE error_msg
  IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN) :: id
    TYPE(attribut),INTENT(IN) :: attrib
    LOGICAL,OPTIONAL,INTENT(out)  :: Ok
    
    TYPE(field_group),POINTER             :: Pt_fg
    INTEGER                         :: Pos
    LOGICAL                         :: success
    
    CALL sorted_list__find(Ids,hash(Id),Pos,success)
    IF (success) THEN
      Pt_fg=>field_group_ids%at(Pos)%Pt
      CALL field_group__set_attribut(Pt_fg,attrib)
      IF (PRESENT(OK)) ok=.TRUE.
    ELSE
      IF (.NOT.PRESENT(OK)) THEN
        WRITE(message,*) 'Field group id :',id,'is undefined'
        CALL error('mod_field_group::field_group__set_attribut')
      ELSE
        OK=.FALSE.
      ENDIF
    ENDIF 
    
  END SUBROUTINE field_group__set_attribut_id
  
  SUBROUTINE field_group__set_attribut_pt(pt_fg,attrib)
  USE mod_attribut
  USE mod_object
  IMPLICIT NONE
    TYPE(field_group),POINTER :: Pt_fg
    TYPE(attribut),INTENT(IN) :: attrib
     
    IF (attrib%object==field_object) THEN
      CALL field__set_attribut(Pt_fg%default_attribut,attrib)
    ENDIF
    
  END SUBROUTINE field_group__set_attribut_pt
    
      
  RECURSIVE SUBROUTINE field_group__new(Pt_fg,Id)
  USE string_function
  IMPLICIT NONE
    TYPE(field_group),POINTER                :: Pt_fg
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN)     :: Id
    
    INTEGER :: Pos
    
    ALLOCATE(Pt_fg%groups)
    ALLOCATE(Pt_fg%fields)
    ALLOCATE(Pt_fg%default_attribut)
    
    CALL vector_field_group__new(Pt_fg%groups)
    CALL vector_field__new(Pt_fg%fields)
    CALL field__new(Pt_fg%default_attribut)
    Pt_fg%has_id=.FALSE.
      
    IF (PRESENT(Id)) THEN
      Pt_fg%id=TRIM(Id)
      Pt_fg%has_id=.TRUE.
      CALL vector_field_group__set_new(field_group_Ids,Pt_fg,Pos)
      CALL sorted_list__Add(Ids,hash(id),Pos)
    ENDIF

  END SUBROUTINE field_group__new

      
  SUBROUTINE field_group__get_new_group(Pt_fg,Pt_fg_out,Id)
  IMPLICIT NONE
    TYPE(field_group),POINTER            :: Pt_fg
    TYPE(field_group),POINTER            :: Pt_fg_out
    CHARACTER(LEN=*),OPTIONAL      :: Id
    
    CALL vector_field_group__get_new(Pt_fg%groups,Pt_fg_out)
 
    IF (PRESENT(id)) THEN
      CALL field_group__new(Pt_fg_out,Id)
    ELSE
      CALL field_group__new(Pt_fg_out)
    ENDIF
    
  END SUBROUTINE field_group__get_new_group

  
  SUBROUTINE field_group__get_new_field(Pt_fg,Pt_f_out,Id)
  IMPLICIT NONE
    TYPE(field_group),POINTER            :: Pt_fg
    TYPE(field),POINTER                  :: Pt_f_out
    CHARACTER(LEN=*),OPTIONAL      :: Id
    
    CALL vector_field__get_new(Pt_fg%fields,Pt_f_out)
    
    IF (PRESENT(id)) THEN
      CALL field__new(Pt_f_out,Id)
    ELSE
      CALL field__new(Pt_f_out)
    ENDIF
    
  END SUBROUTINE field_group__get_new_field
  
  
  SUBROUTINE field_group__get_default_attrib(Pt_fg,Pt_f)
  IMPLICIT NONE
    TYPE(field_group),POINTER  :: Pt_fg
    TYPE(field),POINTER        :: Pt_f
    
    Pt_f=>Pt_fg%default_attribut
  END SUBROUTINE field_group__get_default_attrib

  
  RECURSIVE SUBROUTINE field_group__apply_default(Pt_fg,default)
  IMPLICIT NONE
    TYPE(field_group),POINTER           :: Pt_fg
    TYPE(field),POINTER,OPTIONAL        :: default
    
    INTEGER :: i
    
    IF (PRESENT(default)) THEN
      CALL field__apply_default(default,Pt_fg%default_attribut,Pt_fg%default_attribut)
    ENDIF
      
    DO i=1,Pt_fg%groups%size
      CALL field_group__apply_default(Pt_fg%groups%at(i)%pt,Pt_fg%default_attribut)
    ENDDO
    
    DO i=1,Pt_fg%fields%size
      CALL field__apply_default(Pt_fg%default_attribut,Pt_fg%fields%at(i)%pt,Pt_fg%fields%at(i)%pt)
    ENDDO
  
  END SUBROUTINE field_group__apply_default
  
  SUBROUTINE field_group__solve_ref(pt_fg)
  IMPLICIT NONE
    TYPE(field_group),POINTER  :: Pt_fg

    CALL field_group__solve_field_ref(Pt_fg)
    CALL field_group__solve_axis_ref(Pt_fg)
    CALL field_group__solve_grid_ref(Pt_fg)
    CALL field_group__solve_zoom_ref(Pt_fg)

  END SUBROUTINE  field_group__solve_ref
     
  RECURSIVE SUBROUTINE field_group__solve_field_ref(Pt_fg)
  IMPLICIT NONE
    TYPE(field_group),POINTER  :: Pt_fg
    
    INTEGER :: i

    DO i=1,Pt_fg%groups%size
      CALL field_group__solve_field_ref(Pt_fg%groups%at(i)%pt)
    ENDDO
    
    DO i=1,Pt_fg%fields%size
      CALL field__solve_field_ref(Pt_fg%fields%at(i)%pt)
    ENDDO
  
  END SUBROUTINE field_group__solve_field_ref

  RECURSIVE SUBROUTINE field_group__solve_axis_ref(Pt_fg)
  IMPLICIT NONE
    TYPE(field_group),POINTER  :: Pt_fg
    
    INTEGER :: i

    DO i=1,Pt_fg%groups%size
      CALL field_group__solve_axis_ref(Pt_fg%groups%at(i)%pt)
    ENDDO
    
    DO i=1,Pt_fg%fields%size
      CALL field__solve_axis_ref(Pt_fg%fields%at(i)%pt)
    ENDDO
  
  END SUBROUTINE field_group__solve_axis_ref

  RECURSIVE SUBROUTINE field_group__solve_grid_ref(Pt_fg)
  IMPLICIT NONE
    TYPE(field_group),POINTER  :: Pt_fg
     
    INTEGER :: i
 
    DO i=1,Pt_fg%groups%size
      CALL field_group__solve_grid_ref(Pt_fg%groups%at(i)%pt)
    ENDDO
       
    DO i=1,Pt_fg%fields%size
      CALL field__solve_grid_ref(Pt_fg%fields%at(i)%pt)
    ENDDO
  
  END SUBROUTINE field_group__solve_grid_ref

  RECURSIVE SUBROUTINE field_group__solve_zoom_ref(Pt_fg)
  IMPLICIT NONE
    TYPE(field_group),POINTER  :: Pt_fg
     
    INTEGER :: i
 
    DO i=1,Pt_fg%groups%size
      CALL field_group__solve_zoom_ref(Pt_fg%groups%at(i)%pt)
    ENDDO
       
    DO i=1,Pt_fg%fields%size
      CALL field__solve_zoom_ref(Pt_fg%fields%at(i)%pt)
    ENDDO
  
  END SUBROUTINE field_group__solve_zoom_ref
  
  RECURSIVE SUBROUTINE field_group__print(Pt_fg)
  IMPLICIT NONE
    TYPE(field_group),POINTER  :: Pt_fg
   
    INTEGER :: i
    
    PRINT *,"--- FIELD GROUP ---"
    IF (Pt_fg%has_id) THEN
      PRINT *,"id :",TRIM(Pt_fg%id)
    ELSE
      PRINT *,"id undefined"
    ENDIF
    
    PRINT *,"field default attribut :"
    CALL field__print(Pt_fg%default_attribut)    

    PRINT *,"owned field groups :",Pt_fg%groups%size      
    DO i=1,Pt_fg%groups%size
      CALL field_group__print(Pt_fg%groups%at(i)%pt)
    ENDDO

    PRINT *,"owned field :",Pt_fg%fields%size      
    DO i=1,Pt_fg%fields%size
      CALL field__print(Pt_fg%fields%at(i)%pt)
    ENDDO
    
    PRINT *,"------------"
    
  END SUBROUTINE field_group__print     

END MODULE mod_field_group
