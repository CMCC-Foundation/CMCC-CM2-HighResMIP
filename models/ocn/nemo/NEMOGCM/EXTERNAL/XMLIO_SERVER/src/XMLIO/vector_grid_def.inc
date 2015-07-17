
  TYPE Pt_grid__ 
    TYPE (grid),POINTER :: Pt
    LOGICAL               :: owned
  END TYPE Pt_grid__

  TYPE vector_grid
  
    TYPE (Pt_grid__), POINTER, DIMENSION(:) :: at => NULL()
    INTEGER :: size
    
    INTEGER :: size_max
    REAL    :: grow_factor
    INTEGER :: chunk
  
  END TYPE vector_grid
  
