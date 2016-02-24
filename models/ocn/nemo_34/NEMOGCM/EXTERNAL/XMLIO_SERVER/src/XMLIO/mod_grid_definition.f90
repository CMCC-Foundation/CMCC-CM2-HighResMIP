MODULE mod_grid_definition
  USE mod_grid_group
  
  TYPE(grid_group),SAVE,POINTER :: grid_definition
  
CONTAINS

  SUBROUTINE grid_definition__swap_context(saved_grid_definition)      
  IMPLICIT NONE
    TYPE(grid_group),POINTER :: saved_grid_definition
    
    grid_definition=>saved_grid_definition
    
  END SUBROUTINE grid_definition__swap_context 

  SUBROUTINE grid_definition__Init
  USE mod_grid_group
  IMPLICIT NONE
  
    CALL grid_group__new(grid_definition,"grid_definition")
    
  END SUBROUTINE grid_definition__Init
 
END MODULE mod_grid_definition
