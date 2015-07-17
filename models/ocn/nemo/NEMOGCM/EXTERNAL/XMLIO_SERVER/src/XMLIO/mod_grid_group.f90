MODULE mod_grid_group
  USE mod_grid
  USE mod_xmlio_parameters

  IMPLICIT NONE

  TYPE grid_group
    CHARACTER(LEN=str_len)                    :: id
    LOGICAL                                   :: has_id
    TYPE(vector_grid_group), POINTER          :: groups
    TYPE(vector_grid),POINTER                 :: grids     
    TYPE(grid), POINTER                       :: default_attribut
  END TYPE grid_group

  INCLUDE "vector_grid_group_def.inc"  

  TYPE(vector_grid_group),POINTER,SAVE       :: grid_group_Ids
  TYPE(sorted_list),POINTER,SAVE,PRIVATE     :: Ids 

  INTERFACE grid_group__set_attribut
    MODULE PROCEDURE grid_group__set_attribut_id,grid_group__set_attribut_pt
  END INTERFACE

CONTAINS

  INCLUDE "vector_grid_group_contains.inc"

  SUBROUTINE grid_group__swap_context(saved_grid_group_Ids,saved_ids)
  IMPLICIT NONE
  TYPE(vector_grid_group),POINTER       :: saved_grid_group_Ids
  TYPE(sorted_list),POINTER             :: saved_Ids 
   
   grid_group_ids=>saved_grid_group_ids
   ids=>saved_ids
   
  END SUBROUTINE grid_group__swap_context
   
  SUBROUTINE grid_group__init
  IMPLICIT NONE
    
    CALL vector_grid_group__new(grid_group_Ids)
    CALL sorted_list__new(Ids)
   
  END SUBROUTINE grid_group__init

  SUBROUTINE grid_group__get(Id,Pt_gg)
  USE string_function
  IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN)     :: Id
    TYPE(grid_group),POINTER        :: Pt_gg

    INTEGER                         :: Pos
    LOGICAL                         :: success
    
    CALL sorted_list__find(Ids,hash(Id),Pos,success)
    IF (success) THEN
      Pt_gg=>grid_group_ids%at(Pos)%Pt
    ELSE
      Pt_gg=>NULL()
    ENDIF
    
  END SUBROUTINE grid_group__get


  SUBROUTINE grid_group__set_attribut_id(id,attrib,Ok)
  USE mod_attribut
  USE error_msg
  IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN) :: id
    TYPE(attribut),INTENT(IN) :: attrib
    LOGICAL,OPTIONAL,INTENT(out)  :: Ok
    
    TYPE(grid_group),POINTER             :: Pt_gg
    INTEGER                         :: Pos
    LOGICAL                         :: success
    
    CALL sorted_list__find(Ids,hash(Id),Pos,success)
    IF (success) THEN
      Pt_gg=>grid_group_ids%at(Pos)%Pt
      CALL grid_group__set_attribut(Pt_gg,attrib)
      IF (PRESENT(OK)) ok=.TRUE.
    ELSE
      IF (.NOT.PRESENT(OK)) THEN
        WRITE(message,*) 'grid group id :',id,'is undefined'
        CALL error('mod_grid_group::grid_group__set_attribut')
      ELSE
        OK=.FALSE.
      ENDIF
    ENDIF 
    
  END SUBROUTINE grid_group__set_attribut_id
  
  SUBROUTINE grid_group__set_attribut_pt(pt_gg,attrib)
  USE mod_attribut
  USE mod_object
  IMPLICIT NONE
    TYPE(grid_group),POINTER :: Pt_gg
    TYPE(attribut),INTENT(IN) :: attrib
     
    IF (attrib%object==grid_object) THEN
      CALL grid__set_attribut(pt_gg%default_attribut,attrib)
    ENDIF
    
  END SUBROUTINE grid_group__set_attribut_pt    



  RECURSIVE SUBROUTINE grid_group__new(Pt_gg,Id)
  USE string_function
  IMPLICIT NONE
    TYPE(grid_group),POINTER                :: Pt_gg
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN)     :: Id
    
    INTEGER :: Pos
    
    ALLOCATE(Pt_gg%groups)
    ALLOCATE(Pt_gg%grids)
    ALLOCATE(Pt_gg%default_attribut)
    
    CALL vector_grid_group__new(Pt_gg%groups)
    CALL vector_grid__new(Pt_gg%grids)
    CALL grid__new(Pt_gg%default_attribut)
    Pt_gg%has_id=.FALSE.
      
    IF (PRESENT(Id)) THEN
      Pt_gg%id=TRIM(Id)
      Pt_gg%has_id=.TRUE.
      CALL vector_grid_group__set_new(grid_group_Ids,Pt_gg,Pos)
      CALL sorted_list__Add(Ids,hash(id),Pos)
    ENDIF

  END SUBROUTINE grid_group__new

      
  SUBROUTINE grid_group__get_new_group(Pt_gg,Pt_gg_out,Id)
  IMPLICIT NONE
    TYPE(grid_group),POINTER             :: Pt_gg
    TYPE(grid_group),POINTER             :: Pt_gg_out
    CHARACTER(LEN=*),OPTIONAL            :: Id
    
    CALL vector_grid_group__get_new(Pt_gg%groups,Pt_gg_out)
    CALL grid_group__new(Pt_gg_out)

    IF (PRESENT(id)) THEN
      CALL grid_group__new(Pt_gg_out,Id)
    ELSE
      CALL grid_group__new(Pt_gg_out)
    ENDIF
    
  END SUBROUTINE grid_group__get_new_group

  
  SUBROUTINE grid_group__get_new_grid(Pt_gg,Pt_g_out,Id)
  IMPLICIT NONE
    TYPE(grid_group),POINTER            :: Pt_gg
    TYPE(grid),POINTER                  :: Pt_g_out
    CHARACTER(LEN=*),OPTIONAL      :: Id
    
    CALL vector_grid__get_new(Pt_gg%grids,Pt_g_out)
    
    IF (PRESENT(id)) THEN
      CALL grid__new(Pt_g_out,Id)
    ELSE
      CALL grid__new(Pt_g_out)
    ENDIF
    
  END SUBROUTINE grid_group__get_new_grid
  
  
  SUBROUTINE grid_group__get_default_attrib(Pt_gg,Pt_g)
  IMPLICIT NONE
    TYPE(grid_group),POINTER  :: Pt_gg
    TYPE(grid),POINTER        :: Pt_g
    
    Pt_g=>Pt_gg%default_attribut
  END SUBROUTINE grid_group__get_default_attrib
  
  RECURSIVE SUBROUTINE grid_group__apply_default(Pt_gg,default)
  IMPLICIT NONE
    TYPE(grid_group),POINTER      :: Pt_gg
    TYPE(grid),POINTER,OPTIONAL   :: default
    
    INTEGER :: i
    
    IF (PRESENT(default)) THEN
      CALL grid__apply_default(default,Pt_gg%default_attribut,Pt_gg%default_attribut)
    ENDIF
      
    DO i=1,Pt_gg%groups%size
      CALL grid_group__apply_default(Pt_gg%groups%at(i)%pt,Pt_gg%default_attribut)
    ENDDO
    
    DO i=1,Pt_gg%grids%size
      CALL grid__apply_default(Pt_gg%default_attribut,Pt_gg%grids%at(i)%pt,Pt_gg%grids%at(i)%pt)
    ENDDO
  
  END SUBROUTINE grid_group__apply_default

  RECURSIVE SUBROUTINE grid_group__Process_domain(Pt_gg)
  IMPLICIT NONE
  TYPE(grid_group),POINTER  :: Pt_gg
    INTEGER :: i
    
    DO i=1,pt_gg%groups%size
      CALL grid_group__process_domain(pt_gg%groups%at(i)%pt)
    ENDDO
    
    DO i=1,pt_gg%grids%size
      CALL grid__process_domain(pt_gg%grids%at(i)%pt)
    ENDDO
    
  END SUBROUTINE grid_group__Process_domain
      

  RECURSIVE SUBROUTINE grid_group__print(Pt_gg)
  IMPLICIT NONE
    TYPE(grid_group),POINTER  :: Pt_gg
   
    INTEGER :: i
    
    PRINT *,"--- GRID GROUP ---"
    IF (pt_gg%has_id) THEN
      PRINT *,"id  ",TRIM(pt_gg%id)
    ELSE
      PRINT *,"id undefined"
    ENDIF
    
    PRINT *,"grid default attribut :"
    CALL grid__print(Pt_gg%default_attribut)    

    PRINT *,"owned grid groups : ",Pt_gg%groups%size      
    DO i=1,Pt_gg%groups%size
      CALL grid_group__print(Pt_gg%groups%at(i)%pt)
    ENDDO

    PRINT *,"owned grid : ",Pt_gg%grids%size     
    DO i=1,Pt_gg%grids%size
      CALL grid__print(Pt_gg%grids%at(i)%pt)
    ENDDO
    PRINT *,"------------"
  
  END SUBROUTINE grid_group__print
     
END MODULE mod_grid_group
