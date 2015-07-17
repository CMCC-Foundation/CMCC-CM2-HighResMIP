MODULE error_msg
  USE mod_xmlio_parameters,ONLY : str_msg_len

  CHARACTER(LEN=str_msg_len),SAVE :: message

CONTAINS

  SUBROUTINE Warning(sub_name)
  IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: sub_name
    
    IF (PRESENT(sub_name)) THEN
      PRINT *,"--* WARNING *-- : in subroutine "//sub_name
    ELSE
      PRINT *,"--* WARNING *-- "
    ENDIF
    
    PRINT *,"----> ",TRIM(message)
    message=REPEAT(' ',str_msg_len)
    
  END SUBROUTINE warning
  
  SUBROUTINE error(sub_name)
  IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: sub_name
    
    IF (PRESENT(sub_name)) THEN
      PRINT *,"--* WARNING *-- : in subroutine "//sub_name
    ELSE
      PRINT *,"--* WARNING *-- "
    ENDIF

    PRINT *,"----> ",TRIM(message)
    message=REPEAT(' ',str_msg_len)
    
    STOP 
  
  END SUBROUTINE error    

END MODULE error_msg
