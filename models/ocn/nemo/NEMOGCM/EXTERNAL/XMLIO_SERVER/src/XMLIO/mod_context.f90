MODULE mod_context
  USE mod_xmlio_parameters
  USE mod_sorted_list
  USE mod_field
  USE mod_field_group
  USE mod_field_definition
  USE mod_file
  USE mod_file_group
  USE mod_file_definition
  USE mod_grid
  USE mod_zoom
  USE mod_grid_group
  USE mod_grid_definition
  USE mod_axis
  USE mod_axis_group
  USE mod_axis_definition
  USE mod_dependency
  USE mod_time_parameters


  INTERFACE context__swap
    MODULE PROCEDURE context__swap_id,context__swap_pt
  END INTERFACE  
  
  TYPE, PUBLIC :: context
    CHARACTER(len=str_len)         :: id
    INTEGER                        :: number
  
    TYPE(vector_field),POINTER        :: field__field_Ids
    TYPE(sorted_list),POINTER         :: field__Ids
    TYPE(vector_field_group),POINTER  :: field_group__field_group_ids
    TYPE(sorted_list),POINTER         :: field_group__Ids
    TYPE(field_group),POINTER         :: field_definition__field_def
    TYPE(vector_file),POINTER         :: file__file_Ids
    TYPE(sorted_list),POINTER         :: file__Ids
    TYPE(vector_file_group),POINTER   :: file_group__file_group_Ids
    TYPE(sorted_list),POINTER         :: file_group__Ids 
    TYPE(file_group),POINTER          :: file_definition__file_def 
    TYPE(vector_grid),POINTER         :: grid__grid_Ids
    TYPE(sorted_list),POINTER         :: grid__Ids 
    TYPE(vector_zoom),POINTER         :: zoom__zoom_Ids
    TYPE(sorted_list),POINTER         :: zoom__Ids 
    TYPE(vector_grid_group),POINTER   :: grid_group__grid_group_Ids
    TYPE(sorted_list),POINTER         :: grid_group__Ids 
    TYPE(grid_group),POINTER          :: grid_definition__grid_def
    TYPE(vector_axis),POINTER         :: axis__axis_Ids
    TYPE(sorted_list),POINTER         :: axis__Ids 
    TYPE(vector_axis_group),POINTER   :: axis_group__axis_group_Ids
    TYPE(sorted_list),POINTER         :: axis_group__Ids 
    TYPE(axis_group),POINTER          :: axis_definition__axis_def
    TYPE(vector_file_dep),POINTER     :: dependency__file_enabled
    TYPE(vector_field_out),POINTER    :: dependency__field_enabled 
    TYPE(vector_field_dep),POINTER    :: dependency__field_id
    TYPE(sorted_list),POINTER         :: dependency__sorted_id 
    INTEGER,POINTER                   :: time_param__initial_timestep
    REAL,POINTER                      :: time_param__initial_date
    REAL,POINTER                      :: time_param__timestep_value 
    INTEGER,POINTER                   :: time_param__timestep_number

  END TYPE context
  
  INCLUDE 'vector_context_def.inc'
  TYPE(vector_context),SAVE,POINTER  :: context_ids
  TYPE(sorted_list),POINTER,SAVE,PRIVATE :: Ids

 
  CONTAINS
    INCLUDE 'vector_context_contains.inc'
    
    SUBROUTINE context__init
    IMPLICIT NONE

      ALLOCATE(context_Ids)
      ALLOCATE(Ids)

      CALL vector_context__new(context_Ids)
      CALL sorted_list__new(ids)

    END SUBROUTINE context__init

    SUBROUTINE context__get(Id,Pt_context)
      USE string_function
      IMPLICIT NONE
      CHARACTER(LEN=*),INTENT(IN)     :: Id
      TYPE(context),POINTER              :: Pt_context

      INTEGER                         :: Pos
      LOGICAL                         :: success
    
      CALL sorted_list__find(Ids,hash(Id),Pos,success)
      IF (success) THEN
        Pt_context=>context_ids%at(Pos)%Pt
      ELSE
        Pt_context=>NULL()
      ENDIF
    
    END SUBROUTINE context__get

    SUBROUTINE context__get_new(Id,pt_context)
    USE string_function
    IMPLICIT NONE
      CHARACTER(LEN=*),INTENT(IN)     :: Id
      TYPE(context),POINTER              :: Pt_context

      INTEGER                         :: Pos
      LOGICAL                         :: success
      
      CALL sorted_list__find(Ids,hash(Id),Pos,success)
      IF (success) THEN
        Pt_context=>context_ids%at(Pos)%Pt
      ELSE
        CALL vector_context__get_new(context_ids,pt_context,Pos)
        pt_context%number=Pos
        Pt_context%id=Id
        CALL sorted_list__Add(Ids,hash(Id),Pos)
        
        ALLOCATE(pt_context%field__field_Ids )
        ALLOCATE(pt_context%field__Ids )
        ALLOCATE(pt_context%field_group__field_group_ids )
        ALLOCATE(pt_context%field_group__field_group_ids )
        ALLOCATE(pt_context%field_group__Ids )
        ALLOCATE(pt_context%field_definition__field_def )
        ALLOCATE(pt_context%file__file_Ids)
        ALLOCATE(pt_context%file__Ids)
        ALLOCATE(pt_context%file_group__file_group_Ids)
        ALLOCATE(pt_context%file_group__Ids) 
        ALLOCATE(pt_context%file_definition__file_def) 
        ALLOCATE(pt_context%grid__grid_Ids)
        ALLOCATE(pt_context%grid__Ids) 
        ALLOCATE(pt_context%zoom__zoom_Ids)
        ALLOCATE(pt_context%zoom__Ids) 
        ALLOCATE(pt_context%grid_group__grid_group_Ids)
        ALLOCATE(pt_context%grid_group__Ids) 
        ALLOCATE(pt_context%grid_definition__grid_def)
        ALLOCATE(pt_context%axis__axis_Ids)
        ALLOCATE(pt_context%axis__Ids) 
        ALLOCATE(pt_context%axis_group__axis_group_Ids)
        ALLOCATE(pt_context%axis_group__Ids) 
        ALLOCATE(pt_context%axis_definition__axis_def)
        ALLOCATE(pt_context%dependency__file_enabled)
        ALLOCATE(pt_context%dependency__field_enabled) 
        ALLOCATE(pt_context%dependency__field_id)
        ALLOCATE(pt_context%dependency__sorted_id)  
        ALLOCATE(pt_context%time_param__initial_timestep)
        ALLOCATE(pt_context%time_param__initial_date)
        ALLOCATE(pt_context%time_param__timestep_value) 
        ALLOCATE(pt_context%time_param__timestep_number)

        CALL context__swap(pt_context)
        
        CALL field__init
        CALL field_group__Init
        CALL field_definition__Init

        CALL axis__init
        CALL axis_group__Init
        CALL axis_definition__Init

        CALL grid__init
        CALL grid_group__Init
        CALL grid_definition__Init

        CALL zoom__init

        CALL file__init
        CALL file_group__Init
        CALL file_definition__Init

      ENDIF

    END SUBROUTINE context__get_new


    SUBROUTINE context__create(Id)
    IMPLICIT NONE
      CHARACTER(LEN=*),INTENT(IN)     :: Id
      
      TYPE(context),POINTER              :: Pt_context

      CALL context__get(Id,Pt_context)
      IF (.NOT. ASSOCIATED(Pt_context)) CALL context__get_new(Id,Pt_context)
      
!      CALL field__init
!      CALL field_group__init
!      CALL field_definition__init
      
     END SUBROUTINE context__create

     
    SUBROUTINE context__swap_id(Id)
    USE mod_field
      IMPLICIT NONE
      CHARACTER(LEN=*),INTENT(IN)     :: Id
      TYPE(context),POINTER       :: Pt_context

      INTEGER :: number
      
      CALL context__get(Id,Pt_context)
      IF (.NOT. ASSOCIATED(Pt_context)) THEN
!!      error message
      ENDIF
      
      CALL context__swap(pt_context)  

    END SUBROUTINE context__swap_id

    SUBROUTINE context__swap_pt(Pt_context)
    USE mod_field
      IMPLICIT NONE
      TYPE(context),POINTER       :: Pt_context

      CALL field__swap_context(Pt_context%field__field_Ids ,Pt_context%field__Ids)
      CALL field_group__swap_context(Pt_context%field_group__field_group_ids ,Pt_context%field_group__Ids)
      CALL field_definition__swap_context(Pt_context%field_definition__field_def)
      CALL file__swap_context(Pt_context%file__file_Ids,Pt_context%file__Ids)
      CALL file_group__swap_context(Pt_context%file_group__file_group_Ids,Pt_context%file_group__Ids)
      CALL file_definition__swap_context(Pt_context%file_definition__file_def) 
      CALL grid__swap_context(pt_context%grid__grid_Ids,pt_context%grid__Ids) 
      CALL zoom__swap_context(pt_context%zoom__zoom_Ids,pt_context%zoom__Ids) 
      CALL grid_group__swap_context(pt_context%grid_group__grid_group_Ids,pt_context%grid_group__Ids) 
      CALL grid_definition__swap_context(pt_context%grid_definition__grid_def)
      CALL axis__swap_context(pt_context%axis__axis_Ids,pt_context%axis__Ids) 
      CALL axis_group__swap_context(pt_context%axis_group__axis_group_Ids,pt_context%axis_group__Ids) 
      CALL axis_definition__swap_context(pt_context%axis_definition__axis_def)
      CALL dependency__swap_context(pt_context%dependency__file_enabled,pt_context%dependency__field_enabled,  &
                                    pt_context%dependency__field_id,pt_context%dependency__sorted_id)
      CALL time_parameters__swap_context(pt_context%time_param__initial_timestep,pt_context%time_param__initial_date,  &
                                         pt_context%time_param__timestep_value,pt_context%time_param__timestep_number)
                                    
    END SUBROUTINE context__swap_pt


END MODULE mod_context
