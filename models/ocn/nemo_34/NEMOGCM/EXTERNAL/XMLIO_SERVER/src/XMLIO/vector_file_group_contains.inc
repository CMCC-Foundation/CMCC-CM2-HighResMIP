
  SUBROUTINE vector_file_group__new(vect,vect_size)
  use mod_xmlio_parameters
  IMPLICIT NONE
    TYPE(vector_file_group), INTENT(INOUT) :: vect
    INTEGER, OPTIONAL, INTENT(IN) :: vect_size
     
    IF (ASSOCIATED(vect%at)) CALL vector_file_group__delete(vect)
    
    IF (PRESENT(vect_size)) THEN 
      vect%size_max=MAX(1,vect_size)
    ELSE
      vect%size_max=default_vector_size
    ENDIF
    
    ALLOCATE(vect%at(vect%size_max))
    
    vect%size=0
    vect%grow_factor=default_vector_grow_factor
    
  END SUBROUTINE vector_file_group__new
  
  SUBROUTINE vector_file_group__delete(vect)
  IMPLICIT NONE
    TYPE(vector_file_group), INTENT(IN OUT) :: vect
    
    INTEGER                             :: i
    
    DO i=1,vect%size
      IF (vect%at(i)%owned) DEALLOCATE(vect%at(i)%Pt)
! doit-on avoir une subroutine du genre : CALL file_group__delete() ???
    ENDDO 
    
    IF (ASSOCIATED(vect%at)) DEALLOCATE(vect%at)
    vect%size_max=0
    
  END SUBROUTINE vector_file_group__delete
  
  SUBROUTINE vector_file_group__get_new(vect,pt_value,pos)
  IMPLICIT NONE
  
    TYPE(vector_file_group), INTENT(INOUT)   :: vect
    TYPE(file_group),POINTER                 :: pt_value
    INTEGER,OPTIONAL,INTENT(OUT)         :: Pos
    
    IF (vect%size==vect%size_max) CALL vector_file_group__increase(vect)
    
    vect%size=vect%size+1
    ALLOCATE(vect%at(vect%size)%pt)
    vect%at(vect%size)%owned=.TRUE.
    pt_value=>vect%at(vect%size)%pt
    IF (PRESENT(pos)) pos=vect%size
    
  END SUBROUTINE vector_file_group__get_new

  SUBROUTINE vector_file_group__set_new(vect,pt_value,pos)
  IMPLICIT NONE
  
    TYPE(vector_file_group), INTENT(INOUT)   :: vect
    TYPE(file_group),POINTER                 :: pt_value
    INTEGER,OPTIONAL,INTENT(OUT)         :: Pos
    
    IF (vect%size==vect%size_max) CALL vector_file_group__increase(vect)
    
    vect%size=vect%size+1
    vect%at(vect%size)%pt=>pt_value
    vect%at(vect%size)%owned=.FALSE.
    IF (PRESENT(pos)) pos=vect%size
    
  END SUBROUTINE vector_file_group__set_new
  
  SUBROUTINE vector_file_group__increase(vect)
    TYPE(vector_file_group), INTENT(INOUT) :: vect
    TYPE (Pt_file_group__), POINTER, DIMENSION(:) :: Pt_tmp
    INTEGER :: new_size
    INTEGER :: i
    
    vect%size_max=MAX(vect%size_max+1,INT(vect%size_max*vect%grow_factor))
    ALLOCATE(Pt_tmp(vect%size_max))
    
    DO i=1,vect%size
     Pt_tmp(i)%pt=>vect%at(i)%pt
    ENDDO
    
    DEALLOCATE(vect%at)
    vect%at=>Pt_tmp
    
  END SUBROUTINE vector_file_group__increase 
    
