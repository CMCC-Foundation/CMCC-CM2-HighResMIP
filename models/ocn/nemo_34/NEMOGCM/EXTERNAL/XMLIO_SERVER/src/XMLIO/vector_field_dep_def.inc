
  TYPE Pt_field_dep__ 
    TYPE (field_dep),POINTER :: Pt
    LOGICAL               :: owned
  END TYPE Pt_field_dep__

  TYPE vector_field_dep
  
    TYPE (Pt_field_dep__), POINTER, DIMENSION(:) :: at => NULL()
    INTEGER :: size
    
    INTEGER :: size_max
    REAL    :: grow_factor
    INTEGER :: chunk
  
  END TYPE vector_field_dep
  
