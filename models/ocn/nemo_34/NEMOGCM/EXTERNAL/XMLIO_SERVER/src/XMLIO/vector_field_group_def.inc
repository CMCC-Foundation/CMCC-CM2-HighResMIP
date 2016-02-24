
  TYPE Pt_field_group__ 
    TYPE (field_group),POINTER :: Pt
    LOGICAL               :: owned
  END TYPE Pt_field_group__

  TYPE vector_field_group
  
    TYPE (Pt_field_group__), POINTER, DIMENSION(:) :: at => NULL()
    INTEGER :: size
    
    INTEGER :: size_max
    REAL    :: grow_factor
    INTEGER :: chunk
  
  END TYPE vector_field_group
  
