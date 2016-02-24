
  TYPE Pt_axis_group__ 
    TYPE (axis_group),POINTER :: Pt
    LOGICAL               :: owned
  END TYPE Pt_axis_group__

  TYPE vector_axis_group
  
    TYPE (Pt_axis_group__), POINTER, DIMENSION(:) :: at => NULL()
    INTEGER :: size
    
    INTEGER :: size_max
    REAL    :: grow_factor
    INTEGER :: chunk
  
  END TYPE vector_axis_group
  
