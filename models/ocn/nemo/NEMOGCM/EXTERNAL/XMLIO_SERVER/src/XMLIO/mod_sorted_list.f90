MODULE mod_sorted_list

  PRIVATE
  
  TYPE node
   INTEGER :: key
   INTEGER :: value
   LOGICAL :: has_value
   TYPE(node),POINTER :: child1
   TYPE(node),POINTER :: child2
  END TYPE node

TYPE, PUBLIC :: sorted_list
  PRIVATE
  TYPE(vector_node),POINTER :: tree => NULL()
  TYPE(node),POINTER        :: first => NULL()
  INTEGER :: size = 0
  LOGICAL :: is_init=.FALSE.
END TYPE sorted_list 

  PUBLIC :: sorted_list__new, sorted_list__add, sorted_list__find, sorted_list__delete 
INCLUDE 'vector_node_def.inc'

CONTAINS

INCLUDE 'vector_node_contains.inc'
  
  SUBROUTINE node__new(Pt_n)
  IMPLICIT NONE
    TYPE(node), POINTER :: Pt_n
    
    Pt_n%has_value=.FALSE.
    Pt_n%child1=>NULL()
    pt_n%child2=>NULL()
  END SUBROUTINE node__new
  
  
  SUBROUTINE sorted_list__new(Pt_sl)
  IMPLICIT NONE
    TYPE(sorted_list), POINTER :: Pt_sl
    TYPE(node), POINTER :: Pt_n
    
    IF (Pt_sl%is_init) CALL sorted_list__delete(Pt_sl)
    
    ALLOCATE(Pt_sl%tree)
    CALL vector_node__new(Pt_sl%tree)
    CALL vector_node__get_new(Pt_sl%tree,Pt_n)
    CALL node__new(Pt_n)  
    Pt_sl%first=>Pt_n
    Pt_sl%Size=0
    Pt_sl%is_init=.TRUE.
    
  END SUBROUTINE  sorted_list__new
  
  SUBROUTINE sorted_list__delete(Pt_sl)
  IMPLICIT NONE
    TYPE(sorted_list), POINTER :: Pt_sl
    TYPE(node), POINTER :: Pt_n

    IF (Pt_sl%is_init) THEN
      CALL vector_node__delete(Pt_sl%tree)
      DEALLOCATE(Pt_sl%tree)
      Pt_sl%is_init=.FALSE.
    ENDIF
   
  END SUBROUTINE  sorted_list__delete

   
  SUBROUTINE sorted_list__Add(Pt_sl,key,value)
  IMPLICIT NONE
    TYPE(sorted_list), POINTER :: Pt_sl
    INTEGER,INTENT(IN)  :: key
    INTEGER,INTENT(IN)  :: value

    LOGICAL               :: Out
    TYPE(node), POINTER   :: current
    TYPE(node), POINTER   :: new 
    
    
    CALL vector_node__get_new(Pt_sl%tree,new)
    CALL node__new(new)
    current=>Pt_sl%first
    
      
    Out=.FALSE.
    
    DO WHILE (.NOT. out)
      IF (current%has_value) THEN
        IF (key > current%key) THEN
          IF (ASSOCIATED(current%child2)) THEN
            current=>current%child2
          ELSE
            current%child2=>new
            out=.FALSE.
          ENDIF
        ELSE
          IF (ASSOCIATED(current%child1)) THEN
            current=>current%child1
          ELSE
            current%child1=>new
            out=.FALSE.
          ENDIF
        ENDIF 
      ELSE
        current%has_value=.TRUE.
        current%key=key
        current%value=value
        out=.TRUE.
      ENDIF
    ENDDO
    
    Pt_sl%Size=Pt_sl%Size+1
       
  END SUBROUTINE sorted_list__Add

  
  SUBROUTINE sorted_list__find(Pt_sl,key,value,success)
  IMPLICIT NONE
    TYPE(sorted_list),POINTER   :: Pt_sl
    INTEGER,INTENT(IN)          :: key
    INTEGER,INTENT(OUT)         :: value
    LOGICAL,INTENT(OUT)         :: success
    
    LOGICAL               :: Out
    TYPE(node), POINTER   :: current
    
    current=>Pt_sl%first
    Out=.FALSE.
    Success=.FALSE.
    
    DO WHILE (.NOT. out)
      IF (current%has_value) THEN
        IF (key == current%key) THEN
          value=current%value
          success=.TRUE.
          out=.TRUE.
        ELSE IF (key > current%key) THEN
          IF (ASSOCIATED(current%child2)) THEN
            current=>current%child2
          ELSE
            out=.TRUE.
          ENDIF
        ELSE
          IF (ASSOCIATED(current%child1)) THEN
            current=>current%child1
          ELSE
            out=.TRUE.
          ENDIF
        ENDIF
      ELSE
        out=.TRUE. 
      ENDIF
    ENDDO

  END SUBROUTINE sorted_list__find

  
  SUBROUTINE sorted_list__get_first(Pt_sl,Pt_n)
  IMPLICIT NONE
    TYPE(sorted_list),POINTER   :: Pt_sl
    TYPE(node),POINTER          :: Pt_n

    Pt_n=>Pt_sl%tree%at(1)%Pt
  END SUBROUTINE sorted_list__get_first  


END MODULE mod_sorted_list
