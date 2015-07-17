
  TYPE Pt_field__ 
    TYPE (field),POINTER :: Pt
    LOGICAL               :: owned
  END TYPE Pt_field__

  TYPE vector_field
  
    TYPE (Pt_field__), POINTER, DIMENSION(:) :: at => NULL()
    INTEGER :: size
    
    INTEGER :: size_max
    REAL    :: grow_factor
    INTEGER :: chunk
  
  END TYPE vector_field
  
