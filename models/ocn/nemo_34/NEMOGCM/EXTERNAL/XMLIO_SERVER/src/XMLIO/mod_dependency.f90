MODULE mod_dependency
USE mod_field
USE mod_file
USE mod_grid
USE mod_axis
USE mod_sorted_list
USE mod_zoom

TYPE file_dep
  TYPE(file),POINTER           :: file
  INTEGER                      :: hash
  TYPE(vector_field),POINTER   :: fields
  TYPE(vector_grid),POINTER    :: grids
  TYPE(vector_zoom),POINTER    :: zooms
  TYPE(vector_axis),POINTER    :: axis
END TYPE file_dep

TYPE field_dep
  TYPE(field),POINTER :: field
  INTEGER             :: hash
  
  TYPE(vector_field_out),POINTER   :: field_out
END TYPE field_dep

TYPE field_out
  TYPE(Field), POINTER  :: field
  TYPE(file),  POINTER  :: file
  TYPE(axis),  POINTER  :: axis
  TYPE(grid),  POINTER  :: grid
  TYPE(zoom),  POINTER  :: zoom
END TYPE field_out

INCLUDE 'vector_field_dep_def.inc'
INCLUDE 'vector_file_dep_def.inc'
INCLUDE 'vector_field_out_def.inc'

TYPE(vector_file_dep),POINTER,SAVE      :: file_enabled
TYPE(vector_field_out),POINTER,SAVE     :: field_enabled 
TYPE(vector_field_dep),POINTER,SAVE     :: field_id

 
TYPE(sorted_list),POINTER,SAVE :: sorted_id

CONTAINS

INCLUDE 'vector_field_dep_contains.inc'
INCLUDE 'vector_file_dep_contains.inc'
INCLUDE 'vector_field_out_contains.inc'

  
  SUBROUTINE dependency__swap_context(saved_file_enabled,saved_field_enabled,save_field_id,saved_sorted_id)  
  IMPLICIT NONE
    TYPE(vector_file_dep),POINTER      :: saved_file_enabled
    TYPE(vector_field_out),POINTER     :: saved_field_enabled 
    TYPE(vector_field_dep),POINTER     :: save_field_id
    TYPE(sorted_list),POINTER          :: saved_sorted_id  
    
    file_enabled=>saved_file_enabled   
    field_enabled=>saved_field_enabled
    field_id=>save_field_id
    sorted_id=>saved_sorted_id
      
  END SUBROUTINE dependency__swap_context
  
  SUBROUTINE set_dependency
  IMPLICIT NONE
  
   CALL set_file_dependency
   CALL set_field_enabled
   CALL set_field_dependency

  END SUBROUTINE set_dependency
  
  
  RECURSIVE SUBROUTINE set_file_dependency(Pt_file_group)
  USE mod_file_definition
  USE mod_file_group
  USE string_function
  USE mod_sorted_list
  IMPLICIT NONE
    TYPE (file_group),POINTER,OPTIONAL :: pt_file_group 
    TYPE (file_group),POINTER          :: Pt_fg
    TYPE (file)     ,POINTER           :: Pt_file
    TYPE (file_dep),POINTER            :: Pt_file_dep
    TYPE (sorted_list),POINTER         :: sorted_axis
    TYPE (sorted_list),POINTER         :: sorted_grid
    TYPE (sorted_list),POINTER         :: sorted_zoom
    INTEGER                            :: i
    INTEGER                            :: j
    
    ALLOCATE(sorted_axis)
    ALLOCATE(sorted_grid)
    ALLOCATE(sorted_zoom)
    
    IF (PRESENT(Pt_file_group)) THEN
      Pt_fg=>Pt_file_group
    ELSE
      CALL vector_file_dep__new(file_enabled)
      Pt_fg=>file_definition
    ENDIF
    
    DO i=1,Pt_fg%groups%size
      CALL set_file_dependency(Pt_fg%groups%at(i)%pt)
    ENDDO
    
    DO i=1,Pt_fg%files%size
      Pt_file=>pt_fg%files%at(i)%pt
      IF (Pt_file%enabled) THEN
        CALL vector_file_dep__get_new(file_enabled,Pt_file_dep)
        
        ALLOCATE(Pt_file_dep%fields)
        ALLOCATE(Pt_file_dep%grids)
        ALLOCATE(Pt_file_dep%zooms)
        ALLOCATE(Pt_file_dep%axis)
        pt_file_dep%file=>pt_file
        pt_file_dep%hash=hash(pt_file%id)
        CALL vector_field__new(Pt_file_dep%fields)
        CALL vector_grid__new(Pt_file_dep%grids)
        CALL vector_zoom__new(Pt_file_dep%zooms)
        CALL vector_axis__new(Pt_file_dep%axis)
        CALL sorted_list__new(sorted_axis)
        CALL sorted_list__new(sorted_grid)
        CALL sorted_list__new(sorted_zoom)
        
        CALL Treat_field_group(pt_file%field_list)
        
        CALL sorted_list__delete(sorted_axis)
        CALL sorted_list__delete(sorted_grid)
        CALL sorted_list__delete(sorted_zoom)
      ENDIF
    ENDDO
    
    CONTAINS 
      RECURSIVE SUBROUTINE treat_field_group(pt_fg)
      IMPLICIT NONE  
        TYPE(field_group),POINTER :: Pt_fg
        INTEGER                   :: i
        
        DO i=1,Pt_fg%groups%size
          CALL treat_field_group(Pt_fg%groups%at(i)%pt)
        ENDDO
        
        DO i=1,Pt_fg%fields%size
          CALL treat_field(Pt_fg%fields%at(i)%pt)
        ENDDO
      END SUBROUTINE treat_field_group
      
      
      SUBROUTINE treat_field(pt_field)
      IMPLICIT NONE  
        TYPE(field),POINTER :: Pt_field
        LOGICAL             :: found
        INTEGER             :: Pos

        IF (Pt_field%enabled .AND. Pt_field%level <= Pt_file%output_level) THEN
          CALL vector_field__set_new(Pt_file_dep%fields,Pt_field)
  
          IF (Pt_field%has_grid) THEN
            CALL sorted_list__find(sorted_grid,hash(Pt_field%grid%id),pos,found)
            IF (.NOT. found) THEN 
              CALL vector_grid__set_new(pt_file_dep%grids,Pt_field%grid,pos)
              CALL sorted_list__add(sorted_grid,hash(Pt_field%grid%id),pos)
            ENDIF
          ENDIF

          IF (Pt_field%has_zoom) THEN
            CALL sorted_list__find(sorted_zoom,hash(Pt_field%zoom%id),pos,found)
            IF (.NOT. found) THEN 
              CALL vector_zoom__set_new(pt_file_dep%zooms,Pt_field%zoom,pos)
              CALL sorted_list__add(sorted_zoom,hash(Pt_field%zoom%id),pos)
            ENDIF
          ENDIF
  
          IF (Pt_field%has_axis) THEN
            CALL sorted_list__find(sorted_axis,hash(Pt_field%axis%id),Pos,found)
            IF (.NOT. found) THEN 
              CALL vector_axis__set_new(Pt_file_dep%axis,Pt_field%axis,pos)
              CALL sorted_list__add(sorted_axis,hash(Pt_field%axis%id),pos)
            ENDIF
          ENDIF
        ENDIF
      
      END SUBROUTINE treat_field
      
  END SUBROUTINE set_file_dependency
  
  SUBROUTINE set_field_enabled
  IMPLICIT NONE
  TYPE(file_dep),POINTER   :: pt_file_dep
  TYPE(field_out),POINTER  :: pt_field_out
  INTEGER                  :: i
  INTEGER                  :: j
  
    CALL vector_field_out__new(field_enabled)
    
    DO i=1,file_enabled%size
      pt_file_dep=>file_enabled%at(i)%pt
      DO j=1,pt_file_dep%fields%size
        CALL vector_field_out__get_new(field_enabled,pt_field_out)
        pt_field_out%field=>pt_file_dep%fields%at(j)%pt
        pt_field_out%file=>pt_file_dep%file
        pt_field_out%axis=>pt_field_out%field%axis
        pt_field_out%grid=>pt_field_out%field%grid
        pt_field_out%zoom=>pt_field_out%field%zoom
      ENDDO
    ENDDO
   
      
  END SUBROUTINE set_field_enabled
      
  SUBROUTINE set_field_dependency
  USE string_function
  IMPLICIT NONE
  TYPE(field_out),POINTER :: pt_field_out
  TYPE(field_dep),POINTER :: pt_field_dep
  TYPE(field),POINTER     :: pt_field
  TYPE(field),POINTER     :: pt_field_base
  INTEGER :: pos
  LOGICAL :: found
  INTEGER :: i
  
    CALL vector_field_dep__new(field_id)
    CALL sorted_list__new(sorted_id)
    
    DO i=1,field_enabled%size
      pt_field_out=>field_enabled%at(i)%pt
      pt_field=>pt_field_out%field
      pt_field_base=>pt_field%field_base
      CALL sorted_list__find(sorted_id,hash(pt_field_base%id),pos,found)
      IF (.NOT. found) THEN
        CALL vector_field_dep__get_new(field_id,pt_field_dep,pos)
        ALLOCATE(pt_field_dep%field_out)
        CALL vector_field_out__new(pt_field_dep%field_out)
        pt_field_dep%field=>pt_field_base       
        CALL sorted_list__add(sorted_id,hash(pt_field_base%id),pos)
      ELSE
        pt_field_dep=>field_id%at(pos)%pt
      ENDIF
      
      CALL vector_field_out__set_new(pt_field_dep%field_out,pt_field_out) 
    ENDDO  
      
  END SUBROUTINE set_field_dependency     
      
    
    
END MODULE mod_dependency
