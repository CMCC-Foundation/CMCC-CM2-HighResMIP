MODULE mod_attribut
  USE mod_xmlio_parameters
  USE mod_attribut_list
  USE mod_stdtype
  USE string_function  
  
  TYPE, PUBLIC :: attribut
    INTEGER :: object
    INTEGER :: name
    INTEGER :: type
    INTEGER :: dim(7)
    INTEGER :: ndim
    INTEGER :: integer0_ptr
    INTEGER,POINTER :: integer1_ptr(:)
    INTEGER,POINTER :: integer2_ptr(:,:)
    REAL    :: real0_ptr
    REAL,POINTER    :: real1_ptr(:)
    REAL,POINTER    :: real2_ptr(:,:)
    LOGICAL :: logical0_ptr
    LOGICAL,POINTER :: logical1_ptr(:)
    LOGICAL,POINTER :: logical2_ptr(:,:)
    
    CHARACTER(LEN=str_len) :: string0_ptr
    CHARACTER(LEN=str_len),POINTER :: string1_ptr(:)
    CHARACTER(LEN=str_len),POINTER :: string2_ptr(:,:)
    INTEGER                        :: string_len
    
  END TYPE attribut
  
  INTERFACE attr
    MODULE PROCEDURE attr_int0,attr_int1,attr_int2,              &
                     attr_real0,attr_real1,attr_real2,           &
                     attr_logical0,attr_logical1,attr_logical2,  &
                     attr_string0,attr_string1,attr_string2
  END INTERFACE
                     
CONTAINS
  
  FUNCTION attr_get_object(attr_name)
  USE mod_attribut_list
  USE mod_object
  USE error_msg
  IMPLICIT NONE
    INTEGER, INTENT(IN) :: attr_name
    INTEGER             :: attr_get_object
    
    IF (attr_name > field__begin .AND. attr_name < field__end) THEN
      attr_get_object=field_object
      RETURN
    ENDIF
    
    IF (attr_name > file__begin .AND. attr_name < file__end) THEN
      attr_get_object=file_object
      RETURN
    ENDIF

    IF (attr_name > grid__begin .AND. attr_name < grid__end) THEN
      attr_get_object=grid_object
      RETURN
    ENDIF

    IF (attr_name > axis__begin .AND. attr_name < axis__end) THEN
      attr_get_object=axis_object
      RETURN
    ENDIF
    
    IF (attr_name > zoom__begin .AND. attr_name < zoom__end) THEN
      attr_get_object=zoom_object
      RETURN
    ENDIF

    WRITE (message,*) 'Attribut name value :',attr_name,'is undefined'
    CALL error("mod_attributd::attr_get_object") 
   
  END FUNCTION attr_get_object 

  SUBROUTINE attr_deallocate(attrib)
  IMPLICIT NONE
  TYPE(attribut) :: attrib

    SELECT CASE(attrib%type)
       CASE (integer0)
!!$         DEALLOCATE(attrib%integer0_ptr)
       CASE (integer1)
         DEALLOCATE(attrib%integer1_ptr)
       CASE (integer2)
         DEALLOCATE(attrib%integer2_ptr)
       CASE (real0)
!!$         DEALLOCATE(attrib%real0_ptr)
       CASE (real1)
         DEALLOCATE(attrib%real1_ptr)
       CASE (real2)
         DEALLOCATE(attrib%real2_ptr)
       CASE (logical0)
!!$         DEALLOCATE(attrib%logical0_ptr)
       CASE (logical1)
         DEALLOCATE(attrib%logical1_ptr)
       CASE (logical2)
         DEALLOCATE(attrib%logical2_ptr)
       CASE (string0)
!!$         DEALLOCATE(attrib%string0_ptr)
       CASE (string1)
         DEALLOCATE(attrib%string1_ptr)
       CASE (string2)
         DEALLOCATE(attrib%string2_ptr)
    END SELECT         
  
  END SUBROUTINE attr_deallocate

  FUNCTION attr_int0(attr_name,value) 
  USE mod_stdtype
    IMPLICIT NONE
    INTEGER :: attr_name
    INTEGER :: value
    
    TYPE(attribut) :: attr_int0
    attr_int0%object=attr_get_object(attr_name) 
    attr_int0%name=attr_name
    attr_int0%type=integer0
    attr_int0%integer0_ptr=value
    
  END FUNCTION attr_int0

  FUNCTION attr_int1(attr_name,value) 
  USE mod_stdtype
    IMPLICIT NONE
    INTEGER :: attr_name
    INTEGER,TARGET, INTENT(INOUT) :: value(:)
    
    TYPE(attribut) :: attr_int1
    attr_int1%object=attr_get_object(attr_name) 
    attr_int1%name=attr_name
    attr_int1%type=integer1
    attr_int1%integer1_ptr=>value
    
  END FUNCTION attr_int1    

  FUNCTION attr_int2(attr_name,value) 
  USE mod_stdtype
    IMPLICIT NONE
    INTEGER :: attr_name
    INTEGER,TARGET, INTENT(INOUT) :: value(:,:)
    
    TYPE(attribut) :: attr_int2
    attr_int2%object=attr_get_object(attr_name) 
    attr_int2%name=attr_name
    attr_int2%type=integer2
    attr_int2%integer2_ptr=>value
    
  END FUNCTION attr_int2
          

  FUNCTION attr_real0(attr_name,value) 
  USE mod_stdtype
    IMPLICIT NONE
    INTEGER :: attr_name
    REAL :: value
    
    TYPE(attribut) :: attr_real0
    attr_real0%object=attr_get_object(attr_name) 
    attr_real0%name=attr_name
    attr_real0%type=real0
    attr_real0%real0_ptr=value
    
  END FUNCTION attr_real0

  FUNCTION attr_real1(attr_name,value) 
  USE mod_stdtype
    IMPLICIT NONE
    INTEGER :: attr_name
    REAL,TARGET, INTENT(INOUT) :: value(:)
    
    TYPE(attribut) :: attr_real1
    attr_real1%object=attr_get_object(attr_name) 
    attr_real1%name=attr_name
    attr_real1%type=real1
    attr_real1%real1_ptr=>value
    
  END FUNCTION attr_real1    

  FUNCTION attr_real2(attr_name,value) 
  USE mod_stdtype
    IMPLICIT NONE
    INTEGER :: attr_name
    REAL,TARGET, INTENT(INOUT) :: value(:,:)
    
    TYPE(attribut) :: attr_real2
    attr_real2%object=attr_get_object(attr_name) 
    attr_real2%name=attr_name
    attr_real2%type=REAL2
    attr_real2%real2_ptr=>value
    
  END FUNCTION attr_real2
  
    FUNCTION attr_logical0(attr_name,value) 
  USE mod_stdtype
    IMPLICIT NONE
    INTEGER :: attr_name
    LOGICAL :: value
    
    TYPE(attribut) :: attr_logical0
    attr_logical0%object=attr_get_object(attr_name) 
    attr_logical0%name=attr_name
    attr_logical0%type=logical0
    attr_logical0%logical0_ptr=value
    
  END FUNCTION attr_logical0

  FUNCTION attr_logical1(attr_name,value) 
  USE mod_stdtype
    IMPLICIT NONE
    INTEGER :: attr_name
    LOGICAL,TARGET, INTENT(INOUT) :: value(:)
    
    TYPE(attribut) :: attr_logical1
    attr_logical1%object=attr_get_object(attr_name) 
    attr_logical1%name=attr_name
    attr_logical1%type=logical1
    attr_logical1%logical1_ptr=>value
    
  END FUNCTION attr_logical1    

  FUNCTION attr_logical2(attr_name,value) 
  USE mod_stdtype
    IMPLICIT NONE
    INTEGER :: attr_name
    LOGICAL,TARGET, INTENT(INOUT) :: value(:,:)
    
    TYPE(attribut) :: attr_logical2
    attr_logical2%object=attr_get_object(attr_name) 
    attr_logical2%name=attr_name
    attr_logical2%type=logical2
    attr_logical2%logical2_ptr=>value
    
  END FUNCTION attr_logical2
  

  FUNCTION attr_string0(attr_name,value) 
  USE mod_stdtype
    IMPLICIT NONE
    INTEGER :: attr_name
    CHARACTER(LEN=*) :: value
    TYPE(attribut) :: attr_string0

      attr_string0%object=attr_get_object(attr_name) 
      attr_string0%name=attr_name
      attr_string0%type=string0
      attr_string0%string_len=LEN_TRIM(value)
      attr_string0%string0_ptr=value
    
  END FUNCTION attr_string0

  FUNCTION attr_string1(attr_name,value) 
  USE mod_stdtype
    IMPLICIT NONE
    INTEGER :: attr_name
    CHARACTER(LEN=*),TARGET, INTENT(INOUT) :: value(:)
    
    TYPE(attribut) :: attr_string1
    attr_string1%object=attr_get_object(attr_name) 
    attr_string1%name=attr_name
    attr_string1%type=string1
    attr_string1%string_len=LEN(value)
    attr_string1%string1_ptr=>value
    
  END FUNCTION attr_string1    

  FUNCTION attr_string2(attr_name,value) 
  USE mod_stdtype
    IMPLICIT NONE
    INTEGER :: attr_name
    CHARACTER(LEN=*),TARGET, INTENT(INOUT) :: value(:,:)
    
    TYPE(attribut) :: attr_string2
    attr_string2%object=attr_get_object(attr_name) 
    attr_string2%name=attr_name
    attr_string2%type=string2
    attr_string2%string_len=LEN(value)
    attr_string2%string2_ptr=>value
    
  END FUNCTION attr_string2
 
END MODULE mod_attribut
