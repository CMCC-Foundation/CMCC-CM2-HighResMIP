
  TYPE Pt_file__ 
    TYPE (file),POINTER :: Pt
    LOGICAL               :: owned
  END TYPE Pt_file__

  TYPE vector_file
  
    TYPE (Pt_file__), POINTER, DIMENSION(:) :: at => NULL()
    INTEGER :: size
    
    INTEGER :: size_max
    REAL    :: grow_factor
    INTEGER :: chunk
  
  END TYPE vector_file
  
