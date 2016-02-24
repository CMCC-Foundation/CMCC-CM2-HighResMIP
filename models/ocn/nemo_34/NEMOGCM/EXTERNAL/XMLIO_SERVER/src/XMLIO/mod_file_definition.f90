MODULE mod_file_definition
  USE mod_file_group
  
  TYPE(file_group),SAVE,POINTER :: file_definition
  
CONTAINS

  SUBROUTINE file_definition__swap_context(saved_file_definition)
  IMPLICIT NONE
    TYPE(file_group),POINTER      :: saved_file_definition
     
    file_definition=>saved_file_definition
 
  END SUBROUTINE file_definition__swap_context

  SUBROUTINE file_definition__Init
  USE mod_file_group
  IMPLICIT NONE
  
    CALL file_group__new(file_definition,"file_definition")
    
  END SUBROUTINE file_definition__Init
        
END MODULE mod_file_definition
