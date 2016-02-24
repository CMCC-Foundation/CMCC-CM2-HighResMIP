MODULE xmlio

  USE mod_field
  USE mod_field_group
  USE mod_field_definition
  USE mod_grid
  USE mod_grid_group
  USE mod_grid_definition
  USE mod_axis
  USE mod_axis_group
  USE mod_axis_definition
  USE mod_file
  USE mod_file_group
  USE mod_file_definition
  USE mod_dependency
  USE string_function
  USE error_msg
  USE mod_context
  USE mod_time_parameters
  
CONTAINS

  SUBROUTINE xmlio__Init(xml_file)
  USE parsing_xml
  IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN) :: xml_file
    
    CALL context__init
    CALL parsing_xml_file(xml_file)
    
  END SUBROUTINE xmlio__init
  
  
  SUBROUTINE xmlio__close_definition
  IMPLICIT NONE
  
    CALL axis_group__apply_default(axis_definition)
    CALL grid_group__apply_default(grid_definition)
    CALL field_group__apply_default(field_definition)
    CALL file_group__apply_default(file_definition)
  
    CALL field_group__solve_ref(field_definition)
    CALL file_group__solve_field_ref(file_definition)

    CALL file_group__check(file_definition)
    CALL axis_group__check(axis_definition)
    CALL grid_group__process_domain(grid_definition)
     
    CALL set_dependency

  END SUBROUTINE xmlio__close_definition

END MODULE xmlio
