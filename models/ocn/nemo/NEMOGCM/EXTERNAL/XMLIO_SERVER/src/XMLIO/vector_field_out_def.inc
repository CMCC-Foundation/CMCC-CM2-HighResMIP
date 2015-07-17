
  TYPE Pt_field_out__ 
    TYPE (field_out),POINTER :: Pt
    LOGICAL               :: owned
  END TYPE Pt_field_out__

  TYPE vector_field_out
  
    TYPE (Pt_field_out__), POINTER, DIMENSION(:) :: at => NULL()
    INTEGER :: size
    
    INTEGER :: size_max
    REAL    :: grow_factor
    INTEGER :: chunk
  
  END TYPE vector_field_out
  
