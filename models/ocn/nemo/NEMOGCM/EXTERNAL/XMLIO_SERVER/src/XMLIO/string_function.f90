MODULE string_function
 



CONTAINS

  FUNCTION stdstr(string)
  USE mod_xmlio_parameters
  IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN) :: string
    CHARACTER(LEN=str_len) :: stdstr
    
    stdstr=string
  END FUNCTION stdstr
    
  
  FUNCTION Hash(Str)
  IMPLICIT NONE  
    CHARACTER(LEN=*),INTENT(IN) :: Str
    CHARACTER(LEN=LEN(str)) :: new_str
    INTEGER :: Hash
    INTEGER :: c
    INTEGER :: i
    
    new_str=ADJUSTL(ADJUSTR(Str))
    
    Hash=0
    DO i=1,LEN_TRIM(new_Str)
      c = IACHAR(new_Str(i:i))
      Hash=c+ISHFT(Hash,6)+ISHFT(Hash,16)-Hash
    ENDDO
    
  END FUNCTION Hash
  
  FUNCTION String_to_integer(str,succed)
  USE error_msg
  IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN) :: Str
    LOGICAL,OPTIONAL            :: succed
    
    INTEGER :: String_to_integer
    INTEGER :: ierr
    
    READ(Str,FMT=*,iostat=ierr) String_to_integer
    
    IF (PRESENT(succed)) succed=.TRUE.
      
    IF (ierr/=0) THEN
      IF (PRESENT(succed)) THEN
        succed=.FALSE.
      ELSE 
        WRITE (message,*) 'Error when attempting to convert string :<<',TRIM(str),'>> to integer'
        CALL Error("string_function:string_to_integer")
      ENDIF
    ENDIF
  END FUNCTION String_to_integer
  
  FUNCTION String_to_real(str,succed)
  USE error_msg
  IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN) :: Str
    LOGICAL,OPTIONAL            :: succed
    
    REAL    :: String_to_real
    INTEGER :: ierr
    
    READ(Str,FMT=*,iostat=ierr) String_to_real
    
    IF (PRESENT(succed)) succed=.TRUE.
      
    IF (ierr/=0) THEN
      IF (PRESENT(succed)) THEN
        succed=.FALSE.
      ELSE 
        WRITE (message,*) 'Error when attempting to convert string :<<',TRIM(str),'>> to real'
        CALL Error("string_function:string_to_real")
      ENDIF
    ENDIF
    
  END FUNCTION String_to_real
  
  FUNCTION String_to_logical(str,succed)
  USE error_msg
  IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN) :: Str
    LOGICAL,OPTIONAL            :: succed
    
    LOGICAL :: String_to_logical
    INTEGER :: ierr
    
    READ(Str,FMT=*,iostat=ierr) String_to_logical
    
    IF (PRESENT(succed)) succed=.TRUE.
      
    IF (ierr/=0) THEN
      IF (PRESENT(succed)) THEN
        succed=.FALSE.
      ELSE 
        WRITE (message,*) 'Error when attempting to convert string :<<',TRIM(str),'>> to logical'
        CALL Error("string_function:string_to_logical")
      ENDIF
    ENDIF
    
  END FUNCTION String_to_logical
    
END MODULE string_function
