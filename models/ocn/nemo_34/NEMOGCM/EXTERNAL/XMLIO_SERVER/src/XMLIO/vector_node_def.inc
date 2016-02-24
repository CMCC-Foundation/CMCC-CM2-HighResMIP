
  TYPE Pt_node__ 
    TYPE (node),POINTER :: Pt
    LOGICAL               :: owned
  END TYPE Pt_node__

  TYPE vector_node
  
    TYPE (Pt_node__), POINTER, DIMENSION(:) :: at => NULL()
    INTEGER :: size
    
    INTEGER :: size_max
    REAL    :: grow_factor
    INTEGER :: chunk
  
  END TYPE vector_node
  
