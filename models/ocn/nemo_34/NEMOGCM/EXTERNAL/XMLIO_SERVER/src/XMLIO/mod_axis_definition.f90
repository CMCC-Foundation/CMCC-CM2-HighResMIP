MODULE mod_axis_definition
  USE mod_axis_group
  
  TYPE(axis_group),SAVE,POINTER :: axis_definition
  
CONTAINS

  SUBROUTINE axis_definition__swap_context(saved_axis_definition)      
  IMPLICIT NONE
    TYPE(axis_group),POINTER :: saved_axis_definition
    
    axis_definition=>saved_axis_definition
    
  END SUBROUTINE axis_definition__swap_context 

  SUBROUTINE axis_definition__Init
  USE mod_axis_group
  IMPLICIT NONE
  
    CALL axis_group__new(axis_definition,"axis_definition")
    
  END SUBROUTINE axis_definition__Init
        
END MODULE mod_axis_definition
