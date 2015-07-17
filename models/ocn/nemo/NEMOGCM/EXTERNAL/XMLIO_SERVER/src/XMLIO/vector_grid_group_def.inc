
  TYPE Pt_grid_group__ 
    TYPE (grid_group),POINTER :: Pt
    LOGICAL               :: owned
  END TYPE Pt_grid_group__

  TYPE vector_grid_group
  
    TYPE (Pt_grid_group__), POINTER, DIMENSION(:) :: at => NULL()
    INTEGER :: size
    
    INTEGER :: size_max
    REAL    :: grow_factor
    INTEGER :: chunk
  
  END TYPE vector_grid_group
  
