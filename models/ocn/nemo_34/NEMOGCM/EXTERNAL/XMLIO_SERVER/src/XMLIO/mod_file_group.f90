MODULE mod_file_group
  USE mod_file
  USE mod_xmlio_parameters

  IMPLICIT NONE

  TYPE file_group
    CHARACTER(LEN=str_len)                    :: id
    LOGICAL                                   :: has_id
    TYPE(vector_file_group), POINTER          :: groups
    TYPE(vector_file),POINTER                 :: files     
    TYPE(file), POINTER                       :: default_attribut
  END TYPE file_group

  INCLUDE "vector_file_group_def.inc"  

  TYPE(vector_file_group),POINTER,SAVE       :: file_group_Ids
  TYPE(sorted_list),POINTER,PRIVATE,SAVE     :: Ids 

  INTERFACE file_group__set_attribut
    MODULE PROCEDURE file_group__set_attribut_id,file_group__set_attribut_pt
  END INTERFACE

CONTAINS

  INCLUDE "vector_file_group_contains.inc"

  SUBROUTINE file_group__swap_context(saved_file_group_ids,saved_ids)
  IMPLICIT NONE
    TYPE(vector_file_group),POINTER      :: saved_file_group_Ids
    TYPE(sorted_list),POINTER            :: saved_Ids
    
    file_group_ids=>saved_file_group_ids
    ids=>saved_ids 
 
  END SUBROUTINE file_group__swap_context

  SUBROUTINE file_group__init
  IMPLICIT NONE
    
    CALL vector_file_group__new(file_group_Ids)
    CALL sorted_list__new(Ids)
   
  END SUBROUTINE file_group__init

  SUBROUTINE file_group__get(Id,Pt_fg)
  USE string_function
  IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN)     :: Id
    TYPE(file_group),POINTER        :: Pt_fg

    INTEGER                         :: Pos
    LOGICAL                         :: success
    
    CALL sorted_list__find(Ids,hash(Id),Pos,success)
    IF (success) THEN
      Pt_fg=>file_group_ids%at(Pos)%Pt
    ELSE
      Pt_fg=>NULL()
    ENDIF
    
  END SUBROUTINE file_group__get

  SUBROUTINE file_group__set_attribut_id(id,attrib,Ok)
  USE mod_attribut
  USE error_msg
  IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN) :: id
    TYPE(attribut),INTENT(IN) :: attrib
    LOGICAL,OPTIONAL,INTENT(out)  :: Ok
    
    TYPE(file_group),POINTER             :: Pt_fg
    INTEGER                         :: Pos
    LOGICAL                         :: success
    
    CALL sorted_list__find(Ids,hash(Id),Pos,success)
    IF (success) THEN
      Pt_fg=>file_group_ids%at(Pos)%Pt
      CALL file_group__set_attribut(Pt_fg,attrib)
      IF (PRESENT(OK)) ok=.TRUE.
    ELSE
      IF (.NOT.PRESENT(OK)) THEN
        WRITE(message,*) 'file group id :',id,' is undefined'
        CALL error('mod_file_group::file_group__set_attribut')
      ELSE
        OK=.FALSE.
      ENDIF
    ENDIF 
    
  END SUBROUTINE file_group__set_attribut_id
  
  SUBROUTINE file_group__set_attribut_pt(pt_fg,attrib)
  USE mod_attribut
  USE mod_object
  IMPLICIT NONE
    TYPE(file_group),POINTER :: Pt_fg
    TYPE(attribut),INTENT(IN) :: attrib
     
    IF (attrib%object==file_object) THEN
      CALL file__set_attribut(Pt_fg%default_attribut,attrib)
    ENDIF
    
  END SUBROUTINE file_group__set_attribut_pt
  
      
  RECURSIVE SUBROUTINE file_group__new(Pt_fg,Id)
  USE string_function
  IMPLICIT NONE
    TYPE(file_group),POINTER                :: Pt_fg
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN)     :: Id
    
    INTEGER :: Pos
    
    ALLOCATE(Pt_fg%groups)
    ALLOCATE(Pt_fg%files)
    ALLOCATE(Pt_fg%default_attribut)
    
    CALL vector_file_group__new(Pt_fg%groups)
    CALL vector_file__new(Pt_fg%files)
    CALL file__new(Pt_fg%default_attribut)
    Pt_fg%has_id=.FALSE.
      
    IF (PRESENT(Id)) THEN
      Pt_fg%id=TRIM(Id)
      Pt_fg%has_id=.TRUE.
      CALL vector_file_group__set_new(file_group_Ids,Pt_fg,Pos)
      CALL sorted_list__Add(Ids,hash(id),Pos)
    ENDIF

  END SUBROUTINE file_group__new

      
  SUBROUTINE file_group__get_new_group(Pt_fg,Pt_fg_out,Id)
  IMPLICIT NONE
    TYPE(file_group),POINTER             :: Pt_fg
    TYPE(file_group),POINTER             :: Pt_fg_out
    CHARACTER(LEN=*),OPTIONAL      :: Id
    
    CALL vector_file_group__get_new(Pt_fg%groups,Pt_fg_out)
    CALL file_group__new(Pt_fg_out)

    IF (PRESENT(id)) THEN
      CALL file_group__new(Pt_fg_out,Id)
    ELSE
      CALL file_group__new(Pt_fg_out)
    ENDIF
    
  END SUBROUTINE file_group__get_new_group

  
  SUBROUTINE file_group__get_new_file(Pt_fg,Pt_f_out,Id)
  IMPLICIT NONE
    TYPE(file_group),POINTER            :: Pt_fg
    TYPE(file),POINTER                  :: Pt_f_out
    CHARACTER(LEN=*),OPTIONAL      :: Id
    
    CALL vector_file__get_new(Pt_fg%files,Pt_f_out)
    
    IF (PRESENT(id)) THEN
      CALL file__new(Pt_f_out,Id)
    ELSE
      CALL file__new(Pt_f_out)
    ENDIF
    
  END SUBROUTINE file_group__get_new_file
  
  
  SUBROUTINE file_group__get_default_attrib(Pt_fg,Pt_f)
  IMPLICIT NONE
    TYPE(file_group),POINTER  :: Pt_fg
    TYPE(file),POINTER        :: Pt_f
    
    Pt_f=>Pt_fg%default_attribut
  END SUBROUTINE file_group__get_default_attrib
  
  RECURSIVE SUBROUTINE file_group__apply_default(Pt_fg,default)
  IMPLICIT NONE
    TYPE(file_group),POINTER  :: Pt_fg
    TYPE(file),POINTER,OPTIONAL        :: default
    
    INTEGER :: i
    
    IF (PRESENT(default)) THEN
      CALL file__apply_default(default,Pt_fg%default_attribut,Pt_fg%default_attribut)
    ENDIF
      
    DO i=1,Pt_fg%groups%size
      CALL file_group__apply_default(Pt_fg%groups%at(i)%pt,Pt_fg%default_attribut)
    ENDDO
    
    DO i=1,Pt_fg%files%size
      CALL file__apply_default(Pt_fg%default_attribut,Pt_fg%files%at(i)%pt,Pt_fg%files%at(i)%pt)
    ENDDO
  
  END SUBROUTINE file_group__apply_default

  RECURSIVE SUBROUTINE file_group__solve_field_ref(pt_fg)
  IMPLICIT NONE
    TYPE(file_group),POINTER  :: Pt_fg

    INTEGER :: i

    DO i=1,Pt_fg%groups%size
      CALL file_group__solve_field_ref(Pt_fg%groups%at(i)%pt)
    ENDDO
    
    DO i=1,Pt_fg%files%size
      CALL file__solve_field_ref(Pt_fg%files%at(i)%pt)
    ENDDO
  
  END SUBROUTINE file_group__solve_field_ref
  
  
  RECURSIVE SUBROUTINE file_group__print(Pt_fg)
  IMPLICIT NONE
    TYPE(file_group),POINTER  :: Pt_fg
   
    INTEGER :: i
    
    PRINT *,"--- FILE GROUP ---"
    IF (Pt_fg%has_id) THEN
      PRINT *,"id :",TRIM(Pt_fg%id)
    ELSE
      PRINT *,"id undefined"
    ENDIF
    
    PRINT *,"file default attribut :"
    CALL file__print(Pt_fg%default_attribut)    

    PRINT *,"owned file groups :"      
    DO i=1,Pt_fg%groups%size
      CALL file_group__print(Pt_fg%groups%at(i)%pt)
    ENDDO

    PRINT *,"owned file :"      
    DO i=1,Pt_fg%files%size
      CALL file__print(Pt_fg%files%at(i)%pt)
    ENDDO
  
  END SUBROUTINE file_group__print     

  RECURSIVE SUBROUTINE file_group__Check(Pt_fg)
  IMPLICIT NONE
  
    TYPE(file_group),POINTER  :: Pt_fg
    INTEGER :: i
   
    DO i=1,Pt_fg%groups%size
      CALL file_group__check(pt_fg%groups%at(i)%pt)
    ENDDO

    DO i=1,Pt_fg%files%size
      CALL file__check(pt_fg%files%at(i)%pt)
    ENDDO
  
  END SUBROUTINE file_group__check     
          
END MODULE mod_file_group

