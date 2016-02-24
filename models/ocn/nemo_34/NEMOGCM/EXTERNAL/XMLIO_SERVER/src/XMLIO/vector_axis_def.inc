
  TYPE Pt_axis__ 
    TYPE (axis),POINTER :: Pt
    LOGICAL               :: owned
  END TYPE Pt_axis__

  TYPE vector_axis
  
    TYPE (Pt_axis__), POINTER, DIMENSION(:) :: at => NULL()
    INTEGER :: size
    
    INTEGER :: size_max
    REAL    :: grow_factor
    INTEGER :: chunk
  
  END TYPE vector_axis
  
