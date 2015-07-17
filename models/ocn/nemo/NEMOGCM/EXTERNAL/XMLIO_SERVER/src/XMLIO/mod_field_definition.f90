MODULE mod_field_definition
  USE mod_field_group
  
  TYPE(field_group),SAVE,POINTER :: field_definition
 
CONTAINS

  SUBROUTINE field_definition__swap_context(saved_field_definition)
  IMPLICIT NONE
    TYPE(field_group),POINTER :: saved_field_definition
  
    field_definition=>saved_field_definition

  END SUBROUTINE field_definition__swap_context


  SUBROUTINE field_definition__Init
  IMPLICIT NONE
  
   CALL field_group__new(field_definition,"field_definition")

  END SUBROUTINE field_definition__Init
  
        
END MODULE mod_field_definition

