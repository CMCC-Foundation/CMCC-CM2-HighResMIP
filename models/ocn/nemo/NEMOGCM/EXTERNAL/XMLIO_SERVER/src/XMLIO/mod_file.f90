MODULE mod_file

  USE mod_xmlio_parameters
  USE mod_field_group
  USE mod_sorted_list

  TYPE, PUBLIC :: file
    CHARACTER(len=str_len)           :: id
    LOGICAL                          :: has_id
    CHARACTER(len=str_len)           :: name
    LOGICAL                          :: has_name
    CHARACTER(len=str_len)           :: name_suffix
    LOGICAL                          :: has_name_suffix
    CHARACTER(len=str_len)           :: description
    LOGICAL                          :: has_description
    INTEGER                          :: output_freq
    LOGICAL                          :: has_output_freq
    INTEGER                          :: output_level
    LOGICAL                          :: has_output_level
    LOGICAL                          :: enabled
    LOGICAL                          :: has_enabled
    INTEGER                          :: internal(internal_file)
    TYPE(field_group),POINTER        :: field_list 
  END TYPE file

  INCLUDE 'vector_file_def.inc'
  
  TYPE(vector_file),POINTER,SAVE             :: file_Ids
  TYPE(sorted_list),POINTER,SAVE,PRIVATE     :: Ids 

  INTERFACE file__set_attribut
    MODULE PROCEDURE file__set_attribut_id,file__set_attribut_pt
  END INTERFACE

CONTAINS
  INCLUDE 'vector_file_contains.inc'

  SUBROUTINE file__swap_context(saved_file_ids,saved_ids)
  IMPLICIT NONE
    TYPE(vector_file),POINTER      :: saved_file_Ids
    TYPE(sorted_list),POINTER      :: saved_Ids
    
    file_ids=>saved_file_ids
    ids=>saved_ids 
 
  END SUBROUTINE file__swap_context


  SUBROUTINE file__init
  IMPLICIT NONE
    
    CALL vector_file__new(file_Ids)
    CALL sorted_list__new(Ids)
   
  END SUBROUTINE file__init
  
  SUBROUTINE file__get(Id,Pt_file)
  USE string_function
  IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN)     :: Id
    TYPE(file),POINTER              :: Pt_file

    INTEGER                         :: Pos
    LOGICAL                         :: success
    
    CALL sorted_list__find(Ids,hash(Id),Pos,success)
    IF (success) THEN
      Pt_file=>file_ids%at(Pos)%Pt
    ELSE
      Pt_file=>NULL()
    ENDIF
    
  END SUBROUTINE file__get
  
  SUBROUTINE file__new(pt_file,Id)
  USE string_function
  IMPLICIT NONE
   TYPE(file), POINTER           :: pt_file
   CHARACTER(LEN=*),OPTIONAL     :: Id
   INTEGER                       :: Pos

   ALLOCATE(pt_file%field_list)
   CALL field_group__new(pt_file%field_list)
      
   pt_file%has_id           = .FALSE.
   pt_file%has_name         = .FALSE.
   pt_file%has_name_suffix  = .FALSE.
   pt_file%has_description  = .FALSE.
   pt_file%has_output_freq  = .FALSE.
   pt_file%has_output_level = .FALSE.
   pt_file%has_enabled = .FALSE.
   
   IF (PRESENT(Id)) THEN
     Pt_file%id=TRIM(ADJUSTL(Id))
     Pt_file%has_id=.TRUE.
     CALL vector_file__set_new(file_Ids,Pt_file,Pos)
     CALL sorted_list__Add(Ids,hash(id),Pos)
   ENDIF

  END SUBROUTINE file__new

  SUBROUTINE file__set(pt_file, name, name_suffix, description, output_freq, output_level,enabled)
  IMPLICIT NONE
    TYPE(file), POINTER         :: pt_file
    CHARACTER(len=*)  ,OPTIONAL :: name
    CHARACTER(len=*)  ,OPTIONAL :: name_suffix
    CHARACTER(len=*)  ,OPTIONAL :: description
    INTEGER           ,OPTIONAL :: output_freq
    INTEGER           ,OPTIONAL :: output_level
    LOGICAL           ,OPTIONAL :: enabled

    IF (PRESENT(name)) THEN
        pt_file%name=TRIM(ADJUSTL(name))
        pt_file%has_name = .TRUE.
    ENDIF

    IF (PRESENT(name_suffix)) THEN
        pt_file%name_suffix=TRIM(ADJUSTL(name_suffix))
        pt_file%has_name_suffix = .TRUE.
    ENDIF

    IF (PRESENT(description)) THEN
        pt_file%description=TRIM(ADJUSTL(description))
        pt_file%has_description = .TRUE.
    ENDIF
 
    IF (PRESENT(output_freq)) then
        pt_file%output_freq=output_freq
        pt_file%has_output_freq = .TRUE.
    ENDIF

    IF (PRESENT(output_level)) then
        pt_file%output_level = output_level
        pt_file%has_output_level = .TRUE.
    ENDIF

    IF (PRESENT(enabled)) then
        pt_file%enabled = enabled
        pt_file%has_enabled = .TRUE.
    ENDIF
    
  END SUBROUTINE file__set

  SUBROUTINE file__set_attribut_id(id,attrib,ok)
  USE mod_attribut
  USE error_msg
  IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN)   :: id
    TYPE(attribut),INTENT(IN)     :: attrib
    LOGICAL,OPTIONAL,INTENT(OUT)  :: ok
    
    TYPE(file),POINTER              :: Pt_file
    INTEGER                         :: Pos
    LOGICAL                         :: success
    
    CALL sorted_list__find(Ids,hash(Id),Pos,success)
    IF (success) THEN
      Pt_file=>file_ids%at(Pos)%Pt
      CALL file__set_attribut_pt(Pt_file,attrib)
      IF (PRESENT(OK)) OK=.TRUE.
    ELSE
      IF (.NOT.PRESENT(OK)) THEN
        WRITE(message,*) 'File id : ',id,' is undefined'
        CALL error('mod_file::file__set_attribut')
      ELSE
        OK=.FALSE.
      ENDIF
    ENDIF  
  
  END SUBROUTINE file__set_attribut_id
      
  SUBROUTINE file__set_attribut_pt(Pt_file,attrib)
  USE mod_attribut
  USE mod_file_attribut
  USE error_msg
  IMPLICIT NONE
    TYPE(file),POINTER        :: Pt_file
    TYPE(attribut),INTENT(IN) :: attrib

    
    SELECT CASE(attrib%name)
      CASE (file__name)
        IF (attrib%type==string0) CALL  file__set(pt_file,name=attrib%string0_ptr(1:attrib%string_len)) ; RETURN
      CASE (file__name_suffix)
        IF (attrib%type==string0) CALL  file__set(pt_file,name_suffix=attrib%string0_ptr(1:attrib%string_len)) ; RETURN
      CASE (file__description)
        IF (attrib%type==string0) CALL  file__set(pt_file,description=attrib%string0_ptr(1:attrib%string_len)) ; RETURN
      CASE (file__output_freq)
        IF (attrib%type==integer0) CALL  file__set(pt_file,output_freq=attrib%integer0_ptr) ; RETURN
      CASE (file__output_level)
        IF (attrib%type==integer0) CALL  file__set(pt_file,output_level=attrib%integer0_ptr) ; RETURN
      CASE (file__enabled)
        IF (attrib%type==logical0) CALL  file__set(pt_file,enabled=attrib%logical0_ptr) ; RETURN
     END SELECT

     WRITE(message,*) 'file attribut ',attrib%name,' : type : ',attrib%type,       &
                      ' : Attribute type is incompatible with the provided value'
     CALL error('mod_file::file__set_attribut')
    
  END SUBROUTINE file__set_attribut_pt
  
     
  SUBROUTINE file__get_field_list(pt_file,pt_field_list)
  IMPLICIT NONE
    TYPE(file),POINTER         :: pt_file
    TYPE(field_group),POINTER  :: pt_field_list
    
      pt_field_list=>pt_file%field_list
  
  END SUBROUTINE file__get_field_list
    
  SUBROUTINE file__print(pt_file)
  IMPLICIT NONE
    TYPE(file), POINTER         :: pt_file

    PRINT *,"---- FILE ----"
    IF (pt_file%has_id) THEN
      PRINT *,"id = ",TRIM(pt_file%id)
    ELSE
      PRINT *,"id undefined"
    ENDIF
    
    IF (pt_file%has_name) THEN
      PRINT *,"name = ",TRIM(pt_file%name)
    ELSE
      PRINT *,"name undefined"
    ENDIF

    IF (pt_file%has_name_suffix) THEN
      PRINT *,"name_suffix = ",TRIM(pt_file%name_suffix)
    ELSE
      PRINT *,"name_suffix undefined"
    ENDIF
    
    IF (pt_file%has_description) THEN
      PRINT *,"description = ",TRIM(pt_file%description)
    ELSE
      PRINT *,"description undefined"
    ENDIF
 
    IF (pt_file%has_output_freq) THEN
      PRINT *,"output_freq = ",pt_file%output_freq
    ELSE
      PRINT *,"output_freq undefined"
    ENDIF

    IF (pt_file%has_output_level) THEN
      PRINT *,"output_level = ",pt_file%output_level
    ELSE
      PRINT *,"output_level undefined"
    ENDIF

    IF (pt_file%has_enabled) THEN
      PRINT *,"enabled = ",pt_file%enabled
    ELSE
      PRINT *,"enabled undefined"
    ENDIF

    PRINT *,"field_list :"
    CALL field_group__print(pt_file%field_list)

    PRINT *,"--------------"

  END SUBROUTINE file__print


  SUBROUTINE file__apply_default(pt_file_default, pt_file_in, pt_file_out)

    TYPE(file), POINTER :: pt_file_default, pt_file_in, pt_file_out

    IF (pt_file_in%has_name) THEN
        pt_file_out%name=pt_file_in%name
        pt_file_out%has_name=.TRUE.
    ELSE IF ( pt_file_default%has_name) THEN
        pt_file_out%name=pt_file_default%name
        pt_file_out%has_name=.TRUE.
    ELSE
        pt_file_out%has_name=.FALSE.
    ENDIF

    IF (pt_file_in%has_name_suffix) THEN
        pt_file_out%name_suffix=pt_file_in%name_suffix
        pt_file_out%has_name_suffix=.TRUE.
    ELSE IF ( pt_file_default%has_name_suffix) THEN
        pt_file_out%name_suffix=pt_file_default%name_suffix
        pt_file_out%has_name_suffix=.TRUE.
    ELSE
        pt_file_out%has_name_suffix=.FALSE.
    ENDIF
        
    IF (pt_file_in%has_description) THEN
        pt_file_out%description=pt_file_in%description
        pt_file_out%has_description=.TRUE.
    ELSE IF ( pt_file_default%has_description ) THEN
        pt_file_out%description=pt_file_default%description
        pt_file_out%has_description=.TRUE.
    ELSE
        pt_file_out%has_description=.FALSE.
    ENDIF

    IF (pt_file_in%has_output_freq) THEN
        pt_file_out%output_freq=pt_file_in%output_freq
        pt_file_out%has_output_freq=.TRUE.
    ELSE IF ( pt_file_default%has_output_freq ) THEN
        pt_file_out%output_freq=pt_file_default%output_freq
        pt_file_out%has_output_freq=.TRUE.
    ELSE
        pt_file_out%has_output_freq=.FALSE.
    ENDIF

    IF (pt_file_in%has_output_level) THEN
        pt_file_out%output_level=pt_file_in%output_level
        pt_file_out%has_output_level=.TRUE.
    ELSE IF ( pt_file_default%has_output_level ) THEN
        pt_file_out%output_level=pt_file_default%output_level
        pt_file_out%has_output_level=.TRUE.
    ELSE
        pt_file_out%has_output_level=.FALSE.
    ENDIF

    IF (pt_file_in%has_enabled) THEN
        pt_file_out%enabled=pt_file_in%enabled
        pt_file_out%has_enabled=.TRUE.
    ELSE IF ( pt_file_default%has_enabled ) THEN
        pt_file_out%enabled=pt_file_default%enabled
        pt_file_out%has_enabled=.TRUE.
    ELSE
        pt_file_out%has_enabled=.FALSE.
    ENDIF
    
    CALL field_group__apply_default(pt_file_out%field_list)


  END SUBROUTINE file__apply_default

   
  SUBROUTINE file__solve_field_ref(pt_file)
  IMPLICIT NONE
    TYPE(file), POINTER :: pt_file
    
    CALL field_group__solve_ref(pt_file%field_list)
  
  END SUBROUTINE file__solve_field_ref
  
  
  SUBROUTINE file__Check(pt_file)
  USE error_msg
  IMPLICIT NONE
    TYPE(file), POINTER :: pt_file
      
    IF (.NOT. pt_file%has_name) THEN
      IF (pt_file%has_id) THEN
        pt_file%name=TRIM(pt_file%id)
      ELSE
        WRITE(message,*) "File has no name and no id" 
        CALL error("mod_file::file__check")
      ENDIF 
    ENDIF
 
 END SUBROUTINE file__Check
    
END MODULE mod_file
