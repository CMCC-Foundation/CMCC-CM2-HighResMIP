
  TYPE Pt_context__ 
    TYPE (context),POINTER :: Pt
    LOGICAL               :: owned
  END TYPE Pt_context__

  TYPE vector_context
  
    TYPE (Pt_context__), POINTER, DIMENSION(:) :: at => NULL()
    INTEGER :: size
    
    INTEGER :: size_max
    REAL    :: grow_factor
    INTEGER :: chunk
  
  END TYPE vector_context
  
