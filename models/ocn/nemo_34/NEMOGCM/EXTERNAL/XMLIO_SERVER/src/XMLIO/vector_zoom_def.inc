
  TYPE Pt_zoom__ 
    TYPE (zoom),POINTER :: Pt
    LOGICAL               :: owned
  END TYPE Pt_zoom__

  TYPE vector_zoom
  
    TYPE (Pt_zoom__), POINTER, DIMENSION(:) :: at => NULL()
    INTEGER :: size
    
    INTEGER :: size_max
    REAL    :: grow_factor
    INTEGER :: chunk
  
  END TYPE vector_zoom
  
